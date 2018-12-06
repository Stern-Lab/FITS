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
        std::cerr << "error while opening actual data file: " << filename << std::endl;
        throw "error while opening actual data file";
    }
    
    std::string tmp_line;
    bool is_first_line = true;
    int current_position = -1;
    
    ActualDataPositionData tmp_position_data;
    
    // todo: check if newline needs to be normalized
    while (std::getline(infile, tmp_line)) {
        
        if (is_first_line) {
            is_first_line = false;
            continue;
        }
        
        std::vector<std::string> line_fields;
        try {
            boost::split(line_fields, tmp_line, boost::is_any_of("\t"));
            
        }
        catch (std::exception& e) {
            std::cerr << "Error while parsing actual data file to columns. Line: " << std::endl << tmp_line << std::endl;
            throw e;
        }
        catch (...) {
            std::cerr << "Unknown error while parsing actual data file to columns. Line: " << std::endl << tmp_line << std::endl;
            std::cerr << "NOT TERMINATING" << std::endl;
        }
        
        
        if (line_fields.size() <=1) {
            std::cerr << "skipping line" << std::endl;
            continue;
        }
        
        if (line_fields.size() < ACTUAL_DATA_MINIMAL_COLS) {
            std::cerr << "incompatible actual data file with " << line_fields.size() << " columns, expected " << ACTUAL_DATA_MINIMAL_COLS << std::endl;
            throw "incompatible actual data file.";
        }
        
        
        ActualDataEntry tmp_data_entry;
        
        try {
            boost::trim(line_fields[ACTUAL_DATA_COLUMN_GENERATION]);
            boost::trim(line_fields[ACTUAL_DATA_COLUMN_ALLELE]);
            boost::trim(line_fields[ACTUAL_DATA_COLUMN_FREQ]);
        }
        catch (...) {
            std::cerr << "Error while parsing actual data file (trim)." << std::endl;
            throw "Error while pasring actual data file (trim).";
        }
        
        try {
            tmp_data_entry.gen = boost::lexical_cast<int>(line_fields[ACTUAL_DATA_COLUMN_GENERATION]);
            tmp_data_entry.allele = boost::lexical_cast<int>(line_fields[ACTUAL_DATA_COLUMN_ALLELE]);
            tmp_data_entry.freq = boost::lexical_cast<FLOAT_TYPE>(line_fields[ACTUAL_DATA_COLUMN_FREQ]);
            tmp_data_entry.ref = -1;
            tmp_data_entry.read_count = -1;
        }
        catch (...) {
            std::cerr << "Error while parsing actual data file (cast):" << std::endl << tmp_line << std::endl;
            throw "Error while pasring actual data file (cast).";
        }
        
        
        // assuming position column exists
        int tmp_pos = -1;
        if (line_fields.size() > ACTUAL_DATA_MINIMAL_COLS) {
            try {
                tmp_pos = boost::lexical_cast<int>(line_fields[ACTUAL_DATA_COLUMN_POSITION]);
                if (tmp_pos < 0) {
                    throw "Generation number cannot be negative: " + std::to_string(tmp_pos);
                }
                
                _multi_positions = true;
            }
            catch ( const char* str ) {
                throw "Error while loading data: " + std::string(str);
            }
            catch (...) {
                std::cerr << "Error while parsing actual data file - position column (cast). Data is: " << line_fields[ACTUAL_DATA_COLUMN_POSITION] << std::endl;
                throw "Error while pasring actual data file - position column (cast). Data is:" + line_fields[ACTUAL_DATA_COLUMN_POSITION];
            }
        }
        
        try {
            // have we started a new position?
            if ( (current_position != tmp_pos) && (current_position>0) ) {
                
                std::sort(tmp_position_data._actual_data.begin(), tmp_position_data._actual_data.end());
                tmp_position_data._position = current_position;
                _position_data.push_back(tmp_position_data);
                
                tmp_position_data.Clear();
            }
            current_position = tmp_pos;
            tmp_position_data._actual_data.push_back( tmp_data_entry );
        }
        catch (...) {
            std::cerr << "Error while adding actual data entry" << std::endl;
            throw "Error while adding actual data entry";
        }
    }
    
    try {
        std::sort(tmp_position_data._actual_data.begin(), tmp_position_data._actual_data.end());
        tmp_position_data._position = current_position;
        _position_data.push_back(tmp_position_data);
    }
    catch (...) {
        std::cerr << "Error while sorting actual data." << std::endl;
        throw "Error while sorting actual data.";
    }
    
    infile.close();
    
    _is_initialized = true;
}

