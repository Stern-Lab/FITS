/*
 FITS - Flexible Inference from Time-Series data
 (c) 2016-2019 by Tal Zinger
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

std::string ActualDataFile::GetAsRawText()
{
    std::string tmp_str = "";
    
    tmp_str += "gen\t";
    tmp_str += "allele\t";
    tmp_str += "freq\t";
    tmp_str += "pos\n";
    
    for ( auto current_position : _position_data ) {
        
        for ( auto current_data_entry : current_position._actual_data ) {
            
            auto current_generation = current_data_entry.gen;
            auto current_allele = current_data_entry.allele;
            auto current_frequency = current_data_entry.freq;
            auto current_position = current_data_entry.pos;
            
            tmp_str += std::to_string(current_generation) + "\t";
            tmp_str += std::to_string(current_allele) + "\t";
            tmp_str += std::to_string(current_frequency) + "\t";
            
            if ( current_position < 0 ) {
                tmp_str += "N/A";
            }
            else {
                tmp_str += std::to_string(current_position);
            }
            tmp_str += "\n";
        }
    }
    
    return tmp_str;
}

std::size_t ActualDataFile::GetNumberOfPositions()
{
    return _position_data.size();
}

std::vector<int> ActualDataFile::GetPositionNumbers()
{
    std::vector<int> tmp_vec;
    
    for ( auto current_pos : _position_data ) {
        tmp_vec.push_back( current_pos._position );
    }
    
    return tmp_vec;
}

ActualDataPositionData ActualDataFile::GetFirstPosition()
{
    return _position_data[0];
}

ActualDataPositionData ActualDataFile::GetPosition( int position )
{
    std::vector<ActualDataPositionData>::iterator tmp_iterator;
    // get to first entry of desired position
    for ( tmp_iterator = _position_data.begin(); tmp_iterator < _position_data.end(); ++tmp_iterator ) {
        if ( tmp_iterator->_position == position ) {
            break;
        }
    }
    
    // we actually got to the end of the vector, without finding the desired position
    if ( tmp_iterator == _position_data.end() ) {
        std::string tmp_str = "GetInitFreqs - position not found: " + std::to_string(position);
        throw tmp_str;
    }
    
    return *tmp_iterator;
}


MATRIX_TYPE ActualDataFile::GetActualFreqsAsMatrix(int position )
{
    ValidateMultiPosition(position);
    
    if ( !_multi_positions) {
        return _position_data[0].GetActualFreqsAsMatrix();
    }
    
    auto position_object = GetPosition(position);
    auto tmp_matrix = position_object.GetActualFreqsAsMatrix();
    
    return tmp_matrix;
}


std::vector<int> ActualDataFile::GetActualGenerations( bool only_unique, int position )
{
    ValidateMultiPosition(position);
    
    if ( !_multi_positions) {
        return _position_data[0].GetActualGenerations();
    }
    
    auto position_object = GetPosition(position);
    auto tmp_vec = position_object.GetActualGenerations(only_unique);
    
    return tmp_vec;
}


std::vector<FLOAT_TYPE> ActualDataFile::GetActualFrequencies(int position)
{
    
    ValidateMultiPosition(position);
    
    std::cout << "positiony " << position << std::endl;
    
    if ( !_multi_positions) {
        return _position_data[0].GetActualFrequencies();
    }
    
    auto position_object = GetPosition(position);
    auto tmp_vec = position_object.GetActualFrequencies();
    
    return tmp_vec;
}


int ActualDataFile::GetFirstGeneration( int position )
{
    ValidateMultiPosition(position);
    
    if ( !_multi_positions) {
        return _position_data[0].GetFirstGeneration();
    }
    
    auto position_object = GetPosition(position);
    auto tmp_val = position_object.GetFirstGeneration();
    
    return tmp_val;

}


int ActualDataFile::GetLastGeneration( int position )
{
    ValidateMultiPosition(position);
    
    if ( !_multi_positions) {
        return _position_data[0].GetLastGeneration();
    }
    
    auto position_object = GetPosition(position);
    auto tmp_val = position_object.GetLastGeneration();
    
    return tmp_val;
}


int ActualDataFile::GetNumberOfAlleles()
{
    return _position_data[0].GetNumberOfAlleles();
    
    // cannot be different between positions
    
    /*
    ValidateMultiPosition(position);
     
    if ( !_multi_positions) {
        return _position_data[0].GetNumberOfAlleles();
    }
    
    auto position_object = GetPosition(position);
    auto tmp_val = position_object.GetNumberOfAlleles();
    
    return tmp_val;
     */
}


int ActualDataFile::GetWTIndex(int position )
{
    ValidateMultiPosition(position);
    
    if ( !_multi_positions) {
        return _position_data[0].GetWTIndex();
    }
    
    auto position_object = GetPosition(position);
    auto tmp_val = position_object.GetWTIndex();
    
    return tmp_val;
}


std::vector<FLOAT_TYPE> ActualDataFile::GetInitFreqs(int position)
{
    ValidateMultiPosition(position);
    
    if ( !_multi_positions) {
        return _position_data[0].GetInitFreqs();
    }
    
    auto position_object = GetPosition(position);
    auto tmp_vec = position_object.GetInitFreqs();
    
    return tmp_vec;
}
