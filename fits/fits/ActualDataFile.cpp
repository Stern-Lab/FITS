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
{
}

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
_is_initialized(other._is_initialized),
_position_data(other._position_data)
{}

void ActualDataFile::ValidateMultiPosition(int position)
{
    if ( position >= 0 && !_multi_positions ) {
        std::string tmp_str = "Only single position, but asked for position " + std::to_string(position);
        throw tmp_str;
    }
    
    if (_multi_positions && position < 0) {
        std::string tmp_str = "Multiple positions detected but non specified";
        throw tmp_str;
    }
}



std::vector<ActualDataEntry> ActualDataFile::DataFileToEntries( std::string filename )
{
    std::ifstream infile(filename);
    
    if (!infile.is_open()) {
        std::string tmp_str = "error while opening data file: " + filename;
        throw tmp_str;
    }
    
    std::vector<ActualDataEntry> raw_entries_vec;
    
    std::string tmp_line;
    bool is_first_line = true;
    std::size_t expected_num_columns = 0;
    std::size_t current_line_num = 0;
    
    
    while ( std::getline(infile, tmp_line) ) {
        
        ++current_line_num;
        
        if (is_first_line) {
            
            // keep the number of columns - make sure rest of file contains same number of columns
            std::vector<std::string> line_fields;
            try {
                boost::split(line_fields, tmp_line, boost::is_any_of( fits_constants::FILE_FIELD_DELIMITER ));
                
                expected_num_columns = line_fields.size();
            }
            catch (...) {
                std::string tmp_str = "Error while parsing data file to columns. Header line:\n" + tmp_line;
                throw tmp_str;
            }
            
            if ( line_fields.size() < ACTUAL_DATA_MINIMAL_COLS ) {
                infile.close();
                std::string tmp_str = "Not enough columns in data file (" + std::to_string(line_fields.size()) + " instead of " + std::to_string(ACTUAL_DATA_MINIMAL_COLS) + ")";
                throw tmp_str;
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
            infile.close();
            std::string tmp_str = "Error while parsing data file to columns. Line from file:\n" + tmp_line;
            throw tmp_str;
        }
        
        
        if ( line_fields.size() < ACTUAL_DATA_MINIMAL_COLS || line_fields.size() != expected_num_columns ) {
            
            infile.close();
            
            std::string tmp_str = "Error while reading data file (line " + std::to_string(current_line_num) +
            "): line contains " + std::to_string(line_fields.size())  + " columns, expected " + std::to_string(expected_num_columns) + ".\n";
            
            throw tmp_str;
        }
        
        ActualDataEntry tmp_data_entry;
        try {
            boost::trim(line_fields[ACTUAL_DATA_COLUMN_GENERATION]);
            boost::trim(line_fields[ACTUAL_DATA_COLUMN_ALLELE]);
            boost::trim(line_fields[ACTUAL_DATA_COLUMN_FREQ]);
        }
        catch (...) {
            infile.close();
            std::string tmp_str = "Error while parsing data file while trimming (line " + std::to_string(current_line_num) + ")";
            throw tmp_str;
        }
        
        
        tmp_data_entry.ref = -1;
        tmp_data_entry.read_count = -1;
        try {
            tmp_data_entry.gen = boost::lexical_cast<int>(line_fields[ACTUAL_DATA_COLUMN_GENERATION]);
        }
        catch (...) {
            infile.close();
            std::string tmp_str = "Error while loading data: Generation column ("  + std::to_string(ACTUAL_DATA_COLUMN_GENERATION) + ") in line " + std::to_string(current_line_num) + " does not contain a valid value.";
            throw tmp_str;
        }
        
        try {
            tmp_data_entry.allele = boost::lexical_cast<int>(line_fields[ACTUAL_DATA_COLUMN_ALLELE]);
        }
        catch (...) {
            infile.close();
            std::string tmp_str = "Error while loading data: Allele column ("  + std::to_string(ACTUAL_DATA_COLUMN_GENERATION) + ") in line " + std::to_string(current_line_num) + " does not contain a valid value.";
            throw tmp_str;
        }
        
        try {
            tmp_data_entry.freq = boost::lexical_cast<FLOAT_TYPE>(line_fields[ACTUAL_DATA_COLUMN_FREQ]);
        }
        catch (...) {
            infile.close();
            std::string tmp_str = "Error while loading data: Frequency column ("  + std::to_string(ACTUAL_DATA_COLUMN_FREQ) + ") in line " + std::to_string(current_line_num) + " does not contain a valid value.";
            throw tmp_str;
        }
        
        
        
        // assuming position column exists
        int tmp_pos = -1;
        if (line_fields.size() > ACTUAL_DATA_MINIMAL_COLS) {
            
            try {
                tmp_pos = boost::lexical_cast<int>(line_fields[ACTUAL_DATA_COLUMN_POSITION]);
                if (tmp_pos < 0) {
                    std::string tmp_str = "Position number cannot be negative (line " + std::to_string(current_line_num) + "): " + std::to_string(tmp_pos);
                    throw tmp_str;
                }
                tmp_data_entry.pos = tmp_pos;
            }
            catch (...) {
                infile.close();
                std::string tmp_str = "Error while loading data: Position column ("  + std::to_string(ACTUAL_DATA_COLUMN_POSITION) + ") in line " + std::to_string(current_line_num) + " does not contain a valid value.";
                throw tmp_str;
            }
        }
        
        raw_entries_vec.push_back( tmp_data_entry );
    }
    
    infile.close();
    
    return raw_entries_vec;
}


// Load all data file
// Because data may not be properly sorted by the user, we first need to load all data
// and then divide it into positions (if applicable)
void ActualDataFile::LoadActualData( std::string filename )
{
    _position_data.clear();
    
    // not catching any exceptions from here, should be transparent that I use helper functions
    auto all_data_entries_vec = DataFileToEntries( filename );
    
    std::sort( all_data_entries_vec.begin(), all_data_entries_vec.end() );
    
    
    int current_position = -1;
    ActualDataPositionData tmp_position_data;
    for ( auto current_entry : all_data_entries_vec ) {
        
        // new position
        if ( current_entry.pos != current_position && current_position > 0 ) {
            
            // I like to make sure it's all sorted
            std::sort( tmp_position_data._actual_data.begin(), tmp_position_data._actual_data.end() );
            
            _position_data.push_back( tmp_position_data );
            
            // after copying into the position vector, clear the tmp object
            tmp_position_data.Clear();
            
            // if we identified additional position, we have multi-position data
            _multi_positions = true;
        }
        
        // store data in current position (may be the only one)
        tmp_position_data._position = current_entry.pos;
        current_position = current_entry.pos;
        tmp_position_data._actual_data.push_back( current_entry );
        
        // at least one data entry is stored
        _is_initialized = true;
    
    } // iterating all entries
    
    // add the final (or only) position
    std::sort( tmp_position_data._actual_data.begin(), tmp_position_data._actual_data.end() );
    _position_data.push_back( tmp_position_data );
    
    // finished going through file - have we read any data?
    if ( !_is_initialized ) {
        std::string tmp_str = "Data file appears to be empty.";
        throw tmp_str;
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
                std::string tmp_str = "Inconsistent number of alleles found in data file (" + std::to_string(global_num_alleles) + " vs. " + std::to_string(current_num_alleles) + ")";
                throw tmp_str;
            }
        } // iterating through positions
    }
    
    ValidateDataFile();
}

void ActualDataFile::ValidateDataFile()
{
    for ( auto current_position : _position_data ) {
 
        auto position_wt_idx = current_position.GetWTIndex();
        
        auto generation_vec = current_position.GetActualGenerations();
        
        // verify non-negative frequency and sum up mutant alleles
        for ( auto current_generation : generation_vec ) {
            
            FLOAT_TYPE current_freq_sum = 0;
            
            for ( auto current_entry : current_position._actual_data ) {
                
                if ( current_entry.gen == current_generation ) {
                    
                    if ( current_entry.freq < 0 ) {
                        std::string tmp_str = "negative frequency value (" + std::to_string(current_entry.freq)
                        + ") in generation " + std::to_string(current_generation)
                        + " at position " + std::to_string(current_position._position);
                        
                        throw tmp_str;
                    }
                    
                    if ( current_entry.allele != position_wt_idx ) {
                        current_freq_sum += current_entry.freq;
                    }
                }
            }
            
            // update wt allele
            for ( auto &current_entry : current_position._actual_data ) {
                
                if ( current_entry.gen == current_generation ) {
                    
                    if ( current_entry.allele == position_wt_idx ) {
                        // std::cout << "before = " << current_entry.freq;
                        current_entry.freq = 1 - current_freq_sum;
                        // std::cout << " ; after = " << current_entry.freq << std::endl;
                    }
                }
                
            }
            
            
            /*
            if ( std::fabs(1.0 - current_freq_sum) > 2.0f * std::numeric_limits<float>::epsilon() ) {
                
                std::string tmp_str = "frequency values do not sum up to 1 (" + std::to_string(current_freq_sum)
                + ") in generation " + std::to_string(current_generation)
                + " at position " + std::to_string(current_position._position);
                
                throw tmp_str;
            }
            */
        }
    }
}
