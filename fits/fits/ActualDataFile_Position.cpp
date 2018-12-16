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


bool ActualDataPositionData::operator<( const ActualDataPositionData& other ) const
{
    return _position < other._position;
}


std::vector<int> ActualDataPositionData::GetActualGenerations( bool only_unique )
{
    // returned cached result
    if ( _actual_generations.size() > 0 ) {
        return _actual_generations;
    }
    
    
    std::vector<int> tmp_gen_vec;
    
    for (auto current_entry : _actual_data) {
        tmp_gen_vec.push_back(current_entry.gen);
    }
    
    if (only_unique) {
        auto last = std::unique(tmp_gen_vec.begin(), tmp_gen_vec.end());
        tmp_gen_vec.erase(last, tmp_gen_vec.end());
    }
    
    _actual_generations = tmp_gen_vec;
    return tmp_gen_vec;
}



MATRIX_TYPE ActualDataPositionData::GetActualFreqsAsMatrix()
{
    auto num_alleles = GetNumberOfAlleles();
    auto actual_generations = GetActualGenerations();
    auto num_generations = actual_generations.size();
    
    MATRIX_TYPE freq_matrix( num_generations, num_alleles );
    
    auto actual_freqs = GetActualFrequencies();
    
    for ( auto i=0; i<actual_freqs.size(); ++i ) {
        
        auto current_row = i / num_alleles;
        auto current_col = i %  num_alleles;
        
        freq_matrix(current_row, current_col) = actual_freqs[i];
    }
    
    return freq_matrix;
}

int ActualDataPositionData::GetFirstGeneration()
{
    auto ret_gen = _actual_data[0].gen;
    return ret_gen;
}

int ActualDataPositionData::GetLastGeneration()
{
    auto ret_gen = _actual_data[ _actual_data.size() - 1 ].gen;
    return ret_gen;
}

int ActualDataPositionData::GetNumberOfAlleles()
{
    if ( _actual_data.size() <= 1 ) {
        throw "Get number of alleles - actual data vector contains only 1 entry.";
    }
    
    if ( _num_alleles > 0 ) {
        return _num_alleles;
    }
    
    auto first_gen = _actual_data[0].gen;
    int i=0;
    
    while ( _actual_data[i].gen == first_gen ) {
        ++i;
    }
    
    //std::cout << "GetNumberOfAlleles: " << i << std::endl;
    _num_alleles = i;
    
    return _num_alleles;
}


int ActualDataPositionData::GetWTIndex()
{
    
    if ( _wt_index > -1 ) {
        return _wt_index;
    }
    
    auto first_generation = _actual_data[0].gen;
    
    int current_idx = 0;
    auto current_wt_idx = -1;
    FLOAT_TYPE current_wt_freq = -1.0f;
    
    while (_actual_data[current_idx].gen == first_generation ) {
        if ( _actual_data[current_idx].freq > current_wt_freq ) {
            current_wt_freq = _actual_data[current_idx].freq;
            current_wt_idx = current_idx;
        }
        ++current_idx;
    }
    
    _wt_index = current_wt_idx;
    return current_wt_idx;
    
    /*
    for (auto current_entry : _actual_data) {
        if ( current_entry.gen == first_generation ) {
            if (current_entry.freq > current_wt_freq) {
                current_wt_idx = current_entry.allele;
                current_wt_freq = current_entry.freq;
            }
        }
    }
     */
}


std::vector<FLOAT_TYPE> ActualDataPositionData::GetActualFrequencies()
{
    
    if ( _actual_frequencies.size() > 0 ) {
        return _actual_frequencies;
    }
    std::vector<FLOAT_TYPE> tmp_freq_vec;
    
    
    for (auto current_entry : _actual_data) {
        tmp_freq_vec.push_back(current_entry.freq);
    }
    
    _actual_frequencies = tmp_freq_vec;
    return tmp_freq_vec;
}

std::vector<FLOAT_TYPE> ActualDataPositionData::GetInitFreqs()
{
    if ( _init_frequencies.size() > 0 ) {
        return _init_frequencies;
    }
    
    
    auto first_generation = _actual_data[0].gen;
    std::vector<FLOAT_TYPE> tmp_freqs(0);
    
    for (auto current_entry : _actual_data) {
        if (current_entry.gen == first_generation) {
            tmp_freqs.push_back(current_entry.freq);
        }
    }
    
    _init_frequencies = tmp_freqs;
    return tmp_freqs;
}
