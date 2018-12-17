/*
    FITS - Flexible Inference from Time-Series data
    (c) 2016-2018 by Tal Zinger
    tal.zinger@outlook.com

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/


#include "ActualDataFile.hpp"


ActualDataPositionData::ActualDataPositionData() :
_actual_data(),
_position(-1),
_wt_index(-1),
_num_alleles(-1),
_actual_generations(0),
_actual_frequencies(0.0f),
_init_frequencies(0.0f)
{}

ActualDataPositionData::ActualDataPositionData( const ActualDataPositionData& other) :
_actual_data(other._actual_data),
_position(other._position),
_wt_index(other._wt_index),
_num_alleles(other._num_alleles),
_actual_generations(other._actual_generations),
_actual_frequencies(other._actual_frequencies),
_init_frequencies(other._init_frequencies)
{}

void ActualDataPositionData::Clear()
{
    _actual_data.clear();
    _position = -1;
    _wt_index = -1;
    _num_alleles = -1;
    _actual_generations.clear();
    _actual_frequencies.clear();
    _init_frequencies.clear();
}

ActualDataFile::ActualDataFile() :
_multi_positions(false),
_position_data(0),
_is_initialized(false)
{}


ActualDataFile::ActualDataFile( const ActualDataFile& other ) :
_multi_positions(other._multi_positions),
_is_initialized(other._is_initialized)
{}

void ActualDataFile::ValidateMultiPosition(int position)
{
    if ( position >= 0 && !_multi_positions ) {
        throw "Only single position, but asked for position " + std::to_string(position);
    }
    
    if (_multi_positions && position < 0) {
        throw "Multiple positions detected but non specified";
    }
}

void ActualDataFile::LoadActualData( std::string filename )
{
    _position_data.clear();
    
    std::ifstream infile(filename);
    
    if (!infile.is_open()) {
        std::cerr << "error while opening data file: " << filename << std::endl;
        throw "error while opening data file";
    }
    
    std::string tmp_line;
    bool is_first_line = true;
    int current_position = -1;
    std::size_t expected_num_columns = 0;
    std::size_t current_line_num = 0;
    
    ActualDataPositionData tmp_position_data;
    
    // todo: check if newline needs to be normalized
    while (std::getline(infile, tmp_line)) {
        
        ++current_line_num;
        
        if (is_first_line) {
            
            // keep the number of columns - make sure rest of file contains same number of columns
            std::vector<std::string> line_fields;
            try {
                boost::split(line_fields, tmp_line, boost::is_any_of( fits_constants::FILE_FIELD_DELIMITER ));
                
                expected_num_columns = line_fields.size();
            }
            catch (...) {
                throw "Error while parsing data file to columns. Header line:\n" + tmp_line + "\n";
            }
            
            if ( line_fields.size() < ACTUAL_DATA_MINIMAL_COLS ) {
                throw "Not enough columns in data file.";
            }
            
            
            is_first_line = false;
            continue;
        }
        
        std::vector<std::string> line_fields;
        try {
            //boost::split(line_fields, tmp_line, boost::is_any_of("\t"));
            boost::split(line_fields, tmp_line, boost::is_any_of( fits_constants::FILE_FIELD_DELIMITER ));
        }
        catch (...) {
            throw "Error while parsing data file to columns. Line from file:\n" + tmp_line + "\n";
        }
        
        
        if ( line_fields.size() < ACTUAL_DATA_MINIMAL_COLS || line_fields.size() != expected_num_columns ) {
            std::cerr << "\nError while reading data file (line " << std::to_string(current_line_num)  << "): line contains " << line_fields.size() << " columns, expected " << expected_num_columns << "." << std::endl;
            throw "Line in data file contains inconsistent number of columns.";
        }
        
        ActualDataEntry tmp_data_entry;
        
        try {
            boost::trim(line_fields[ACTUAL_DATA_COLUMN_GENERATION]);
            boost::trim(line_fields[ACTUAL_DATA_COLUMN_ALLELE]);
            boost::trim(line_fields[ACTUAL_DATA_COLUMN_FREQ]);
        }
        catch (...) {
            std::cerr << "\nError while parsing data file (trim)." << std::endl;
            throw "Error while pasring data file (trim).";
        }
        
        try {
            tmp_data_entry.gen = boost::lexical_cast<int>(line_fields[ACTUAL_DATA_COLUMN_GENERATION]);
            tmp_data_entry.allele = boost::lexical_cast<int>(line_fields[ACTUAL_DATA_COLUMN_ALLELE]);
            tmp_data_entry.freq = boost::lexical_cast<FLOAT_TYPE>(line_fields[ACTUAL_DATA_COLUMN_FREQ]);
            tmp_data_entry.ref = -1;
            tmp_data_entry.read_count = -1;
        }
        catch (...) {
            std::cerr << "\nError while parsing data file (cast):" << std::endl << tmp_line << std::endl;
            throw "Error while pasring file (cast).";
        }
        
        
        // assuming position column exists
        int tmp_pos = -1;
        if (line_fields.size() > ACTUAL_DATA_MINIMAL_COLS) {
            try {
                tmp_pos = boost::lexical_cast<int>(line_fields[ACTUAL_DATA_COLUMN_POSITION]);
                if (tmp_pos < 0) {
                    throw "Generation number cannot be negative: " + std::to_string(tmp_pos);
                }
                
                
            }
            catch ( const char* str ) {
                throw "Error while loading data: " + std::string(str);
            }
            catch (...) {
                std::cerr << "Error while parsing data file - position column (cast). Data is: " << line_fields[ACTUAL_DATA_COLUMN_POSITION] << std::endl;
                throw "Error while pasring data file - position column (cast). Data is:" + line_fields[ACTUAL_DATA_COLUMN_POSITION];
            }
        }
        
        try {
            // have we started a new position?
            if ( (current_position != tmp_pos) && (current_position>0) ) {
                
                if ( !tmp_position_data._actual_data.empty() ) {
                    _is_initialized = true;
                }
                std::sort(tmp_position_data._actual_data.begin(), tmp_position_data._actual_data.end());
                tmp_position_data._position = current_position;
                _position_data.push_back(tmp_position_data);
                _multi_positions = true;
                
                tmp_position_data.Clear();
            }
            current_position = tmp_pos;
            tmp_position_data._actual_data.push_back( tmp_data_entry );
        }
        catch (...) {
            std::cerr << "Error while adding data entry" << std::endl;
            throw "Error while adding data entry";
        }
    }
    
    try {
        if ( !tmp_position_data._actual_data.empty() ) {
            _is_initialized = true;
        }
        
        std::sort(tmp_position_data._actual_data.begin(), tmp_position_data._actual_data.end());
        tmp_position_data._position = current_position;
        _position_data.push_back(tmp_position_data);
    }
    catch (...) {
        std::cerr << "Error while sorting data." << std::endl;
        throw "Error while sorting data.";
    }
    
    infile.close();
    
    // finished going through file - have we read any data?
    if ( !_is_initialized ) {
        throw "Data file appears to be empty.";
    }
    
    // do sorting and some consistency checks for all positions
    if ( _multi_positions ) {
        
        std::sort( _position_data.begin(), _position_data.end() );
        
        int global_num_alleles = -1;
        for ( auto current_position_data : _position_data ) {
            
            if ( global_num_alleles < 0 ) {
                global_num_alleles = current_position_data.GetNumberOfAlleles();
                continue;
            }
            
            auto current_num_alleles = current_position_data.GetNumberOfAlleles();
            
            if ( global_num_alleles != current_num_alleles ) {
                std::cerr << "Inconsistent number of alleles found in data file (" << global_num_alleles << " vs. " << current_num_alleles << ")" << std::endl;
                throw "Inconsistent number of alleles found in data file";
            }
        }
    }
    
}

