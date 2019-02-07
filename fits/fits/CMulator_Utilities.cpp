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

#include "CMulator.h"

void CMulator::Reset_Soft()
{
    if ( _allele_init_freqs.size() != _num_alleles ) {
        throw " Reset soft: number of alleles changed unexpectedly, consider hard reset.";
    }
    
    if ( _allele_init_freqs.size() != _num_alleles ) {
        throw " Reset soft: number of alleles changed unexpectedly, consider hard reset.";
    }
    
    // resize the data matrix
    _all_simulated_data.resize( _num_generations+1, _num_alleles );
    _observed_simulated_data.resize( _num_generations+1, _num_alleles );
    
    // for consistency, update initial frequencies
    for (auto j=0; j<_num_alleles; ++j ) {
        _all_simulated_data(0,j) = _allele_init_freqs[j];
        _observed_simulated_data(0,j) = _allele_init_freqs[j];
    }
    
    // now zero all others
    for (auto i=1; i<_num_generations; ++i ) {
        for (auto j=0; j<_num_alleles; ++j ) {
            _all_simulated_data(i,j) = 0.0f;
            _observed_simulated_data(i,j) = 0.0f;
        }
    }
    
    _current_generation = 0;
}

/*
void CMulator::Reset_Soft( std::string new_sim_uid )
{
    Reset_Soft();
    SetSimUID( new_sim_uid );
}
*/
