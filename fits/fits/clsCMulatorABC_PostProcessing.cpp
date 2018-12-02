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

#include "clsCMulatorABC.h"

std::vector<int> clsCMulatorABC::GetUniqueIndexSet( int num_items )
{
    std::cout << "Generating " << num_items << " unique items..." << std::endl;
    
    std::vector<bool> chosen_array(num_items, false);
    
    auto chosen_counter = 0;;
    boost::random::uniform_int_distribution<> local_int_distrib(0, num_items);
    
    while (chosen_counter<num_items) {
        
        auto chosen_index = local_int_distrib(_boost_gen);
        
        if (!chosen_array[chosen_index]) {
            
            chosen_array[chosen_index] = true;
            ++chosen_counter;
            
            //std::cout << chosen_index << std::endl;
        }
    }
    
    std::vector<int> return_vec(num_items, 0);
    
    for ( auto i=0; i<chosen_array.size(); ++i ) {
        
        if ( chosen_array[i] ) {
            return_vec.push_back(i);
        }
    }
    
    return return_vec;
    
}


// Calculate Median Absolute Deviation (MAD)
// This is used for scaling frequencies
std::vector<FLOAT_TYPE> clsCMulatorABC::GetMADPerAllele( std::size_t start_idx, std::size_t end_idx )
{
    
    // gather all frequencies per allele
    std::vector< std::vector<FLOAT_TYPE> > allele_freq_storage( _num_alleles );
    
    for ( auto current_idx=start_idx; current_idx<end_idx; ++current_idx ) {
        
        for ( auto current_allele=0; current_allele<_num_alleles; ++current_allele ) {
            
            boost::numeric::ublas::matrix_column<MATRIX_TYPE> current_col( _simulation_result_vector[current_idx].sim_data_matrix, current_allele );
            
            for ( auto current_freq : current_col ) {
                allele_freq_storage[current_allele].push_back(current_freq);
            }
        }
    }
    
    
    // calculate median per allele
    std::vector<FLOAT_TYPE> allele_median_vec( _num_alleles, 0.0f);
    
    for ( auto current_allele=0; current_allele<allele_freq_storage.size(); ++current_allele ) {
        allele_median_vec[current_allele] = GetMedian( allele_freq_storage[current_allele] );
    }
    
    
    // gather all distances per allele
    std::vector< std::vector<FLOAT_TYPE> > allele_distance_vec( _num_alleles );
    
    for ( auto current_allele=0; current_allele<allele_distance_vec.size(); ++current_allele ) {
        for ( auto val : allele_freq_storage[current_allele] ) {
            allele_distance_vec[current_allele].push_back( std::fabs(val - allele_median_vec[current_allele] ) );
        }
    }
    
    
    // calculate median distance
    std::vector<FLOAT_TYPE> allele_mad_vec( _num_alleles, 0.0f);
    
    for ( auto current_allele=0; current_allele<allele_freq_storage.size(); ++current_allele ) {
        allele_mad_vec[current_allele] = GetMedian( allele_distance_vec[current_allele] );
    }
    
    
    return allele_mad_vec;
}


std::vector<FLOAT_TYPE> clsCMulatorABC::GetSDPerAllele( std::size_t start_idx, std::size_t end_idx )
{
    
    std::vector<FLOAT_TYPE> allele_sd_vec( _num_alleles, 0.0f);
    
    std::vector< boost::accumulators::accumulator_set<
    FLOAT_TYPE,
    boost::accumulators::stats<
    boost::accumulators::tag::variance,
    boost::accumulators::tag::mean,
    boost::accumulators::tag::min,
    boost::accumulators::tag::max> > > allele_accumulator_vec;
    
    allele_accumulator_vec.resize(_num_alleles);
    
    // for each of the unchosen-ones [ ones that were not pseudo data
    for ( auto current_idx=start_idx; current_idx<end_idx; ++current_idx ) {
        
        for ( auto current_allele=0; current_allele<_num_alleles; ++current_allele ) {
            
            boost::numeric::ublas::matrix_column<MATRIX_TYPE> current_col( _simulation_result_vector[current_idx].sim_data_matrix, current_allele );
            
            for ( auto current_freq : current_col ) {
                allele_accumulator_vec[current_allele]( current_freq );
            }
        }
    }
    
    for ( auto current_allele=0; current_allele<_num_alleles; ++current_allele ) {
        allele_sd_vec[current_allele] = std::sqrt( boost::accumulators::variance(allele_accumulator_vec[current_allele]) );
    }
    
    return allele_sd_vec;
}


// for a given dataset, sum the number of results in which the inferred parameter value is smaller than the actual
// parameter values used to generate the dataset
// NOTE: assumes that pseudo-data sets are contiguous (e.g. first 100) and the reset of sims is contiguous
std::vector<unsigned int>
clsCMulatorABC::CoverageSingleDatasetFitness( std::size_t dataset_idx, std::size_t start_idx, std::size_t end_idx ) const
{
    std::vector<unsigned int> local_fitness_underestimate_count_vec(_num_alleles, 0);
    
    // process
    for ( auto current_scaled_sim=start_idx; current_scaled_sim<end_idx; ++current_scaled_sim ) {
        
        for ( auto current_allele=0; current_allele<local_fitness_underestimate_count_vec.size(); ++current_allele ) {
            if ( _simulation_result_vector[current_scaled_sim].fitness_values[current_allele] <
                _simulation_result_vector[dataset_idx].fitness_values[current_allele] ) {
                ++local_fitness_underestimate_count_vec[current_allele];
            }
        }
    }
    
    return local_fitness_underestimate_count_vec;
}
