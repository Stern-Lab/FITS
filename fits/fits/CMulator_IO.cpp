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


#include <iostream>
#include <fstream>

#include "CMulator.h"


MATRIX_TYPE CMulator::GetAllOutputAsMatrix() const
{
    MATRIX_TYPE tmp_matrix(_current_generation, _num_alleles);
    
    for ( auto gen=0; gen<_current_generation; ++gen ) {
        
        for ( auto allele=0; allele<_num_alleles; ++allele ) {
            
            if (_use_observed_data) {
                tmp_matrix(gen,allele) = _observed_simulated_data(gen,allele);
            }
            else {
                tmp_matrix(gen,allele) = _all_simulated_data(gen,allele);
            }
        }
    }
    
    return tmp_matrix;
}


MATRIX_TYPE CMulator::GetAllOutputAsMatrix( std::vector<int> actual_generations ) const
{
    auto rows_to_return = actual_generations.size();
    
    MATRIX_TYPE tmp_matrix(rows_to_return, _num_alleles);
    
    for ( auto gen_idx=0; gen_idx<actual_generations.size(); ++gen_idx ) {
        
        auto shifted_gen = actual_generations[gen_idx] - _generation_shift;
        
        for ( auto allele=0; allele<_num_alleles; ++allele ) {
            
            if (_use_observed_data) {
                tmp_matrix(gen_idx,allele) = _observed_simulated_data( shifted_gen, allele );
            }
            else {
                tmp_matrix(gen_idx,allele) = _all_simulated_data( shifted_gen, allele );
            }
        }
    }
    
    return tmp_matrix;
}


std::string CMulator::GetAllOutputAsText(bool header, std::string delimiter) const
{
	if (!_initialized_with_parameters) {
		throw " GetAllOutputAsText: object not initialized with parameters.";
	}
	
    std::string tmp_str = "";
	if (_current_generation==0) {
		return tmp_str;
	}
	
	// header
    if (header) {
        for (auto myAllele=0; myAllele<_num_alleles; myAllele++) {
            tmp_str +=  "allele" + std::to_string(myAllele) + delimiter;
        }
        tmp_str += "\n";
    }
	
	if (_current_generation>1) {
		
		for (auto myGeneration=0; myGeneration<_current_generation; ++myGeneration) {
			
			for (auto myAllele=0; myAllele<_num_alleles; ++myAllele) {
            
                if (_use_observed_data) {
                    tmp_str += std::to_string( _observed_simulated_data( myGeneration, myAllele)) + fits_constants::FILE_FIELD_DELIMITER;
                }
                else {
                    tmp_str += std::to_string( _all_simulated_data( myGeneration, myAllele)) + fits_constants::FILE_FIELD_DELIMITER;
                }
                
			}
			
			tmp_str += "\n";
		}
	}
	
	return tmp_str;
}


std::string CMulator::GetAllOutputAsTextForR( bool header ) const
{
    if (!_initialized_with_parameters) {
        throw " GetAllOutputAsText: object not initialized with parameters.";
    }
    
    std::string tmp_str = "";
    if (_current_generation==0) {
        return tmp_str;
    }
    
    // header
    if (header) {
        tmp_str += "gen" + fits_constants::FILE_FIELD_DELIMITER + "allele" + fits_constants::FILE_FIELD_DELIMITER + "freq\n";
    }
    
    if (_current_generation>1) {
        
        for (auto myGeneration=0; myGeneration<_current_generation; ++myGeneration) {
            
            for (auto myAllele=0; myAllele<_num_alleles; ++myAllele) {
                
                tmp_str += std::to_string(myGeneration) + fits_constants::FILE_FIELD_DELIMITER;
                tmp_str += std::to_string(myAllele) + fits_constants::FILE_FIELD_DELIMITER;
                
                tmp_str += std::to_string( _all_simulated_data(myGeneration,myAllele) ) + "\n";
            }
        }
    }
    
    return tmp_str;
}
