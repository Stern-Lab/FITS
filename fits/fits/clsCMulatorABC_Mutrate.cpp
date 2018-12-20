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


//std::vector<SimulationResult> clsCMulatorABC::RunMutationInferenceBatch( const PRIOR_DISTRIB &prior_distrib )
std::vector<SimulationResult> clsCMulatorABC::RunMutationInferenceBatch( const PRIOR_DISTRIB &prior_distrib, std::size_t start_idx, std::size_t end_idx )
{
    VerifyIndece( prior_distrib, start_idx, end_idx );
    
    // initialize sim object
    CMulator local_sim_object(_zparams);
    
    /*
    if ( local_sim_object.IsObservedDataUsed()) {
        std::cout << "Using observed (sampled) data for comparison." << std::endl;
    }
    else {
        std::cout << "Using all data for comparison." << std::endl;
    }
    */
    
    if ( !local_sim_object.IsAbleToInferMutationRate() ) {
        std::string tmp_str = "Not enough parameters to infer mutation rate";
        throw tmp_str;
    }
    
    // set initial frequencies from actual data
    auto init_freq_vec = _actual_data_position.GetInitFreqs();
    for (auto i = 0; i < init_freq_vec.size(); ++i) {
        local_sim_object.SetAlleleInitFreq(i, init_freq_vec[i]);
    }
    
    auto first_generation = _actual_data_position.GetFirstGeneration();
    auto last_generation = _actual_data_position.GetLastGeneration();
    auto num_generations = last_generation - first_generation + 1;
    local_sim_object.SetGenerationShift(first_generation);
    local_sim_object.SetNumOfGeneration(num_generations);
    
    // identify wt
    auto wt_allele_it = std::max_element(init_freq_vec.begin(), init_freq_vec.end());
    auto wt_allele_idx = static_cast<unsigned int>(std::distance(init_freq_vec.begin(), wt_allele_it));
    local_sim_object.SetWTAllele(wt_allele_idx, 0.0f, 2.0f);
    
    //auto min_matrix = local_sim_object.GetMinMutationRateMatrix();
    //auto max_matrix = local_sim_object.GetMaxMutationRateMatrix();
    
    //PriorSampler<FLOAT_TYPE> sampler( min_matrix, max_matrix, PriorDistributionType::UNIFORM);
    
    
    //std::vector< std::vector<FLOAT_TYPE> >  mutrate_vector_list;
    
    
    std::vector<SimulationResult> tmp_res_vector;
    
    if ( _zparams.GetInt( "Debug", 0 ) > 0 ) {
        std::cout << "Prior distrib size : " << _global_prior.size() << std::endl;
    }
    
    // simulation for each set of parameters
    //for (auto current_mutrate_idx=0; current_mutrate_idx<prior_distrib.size(); ++current_mutrate_idx ) {
    for (auto current_prior_sample_idx=start_idx; current_prior_sample_idx<end_idx; ++current_prior_sample_idx ) {
        
        auto current_mutrate_vector = prior_distrib[current_prior_sample_idx];
        
        local_sim_object.Reset_Soft();
        
        MATRIX_TYPE tmp_mutrates_matrix( local_sim_object.GetAlleleNumber(), local_sim_object.GetAlleleNumber() );
        
        if ( _zparams.GetInt( "Debug", 0 ) > 0 ) {
            std::cout << "Putting sampled mutation rates in matrix." << std::endl;
        }
        
        std::vector<FLOAT_TYPE> tmp_line_sum( tmp_mutrates_matrix.size1(), 0.0f );
        
        // if it were 0, then it would be overridden
        auto current_single_mutation_rate = current_mutrate_vector[1];
        
        if ( _zparams.GetInt( "Debug", 0 ) > 0 ) {
            std::cout << "Single mutation rate:" << current_single_mutation_rate << std::endl;
        }

        for ( auto i=0; i<current_mutrate_vector.size(); ++i ) {
            
            auto row = i / local_sim_object.GetAlleleNumber();
            auto col = i % local_sim_object.GetAlleleNumber();
            
            if ( local_sim_object.IsSingleMutrateUsed() ) {
                if ( _zparams.GetInt( "Debug", 0 ) > 0 ) {
                    std::cout << "Using single mutation rate: " << current_single_mutation_rate << std::endl;
                }
                
                tmp_mutrates_matrix(row,col) = std::pow( 10, current_single_mutation_rate );
            }
            else {
                tmp_mutrates_matrix(row,col) = std::pow( 10, current_mutrate_vector[i] );
            }
            
         
            // to normalize such that row sums up to 1.0
            if ( row != col ) {
                tmp_line_sum[row] += tmp_mutrates_matrix(row, col);
            }
        }
        
        // normalize
        for ( auto row=0; row<tmp_mutrates_matrix.size1(); ++row ) {
            tmp_mutrates_matrix(row,row) = 1.0f - tmp_line_sum[row];
        }
        
        if ( _zparams.GetInt( "Debug", 0 ) > 0 ) {
            std::cout << "Mutation matrix processing finished:" << std::endl;
            std::cout << tmp_mutrates_matrix << std::endl;
        }
        
        local_sim_object.SetMutationRateMatrix(tmp_mutrates_matrix);
        
        local_sim_object.EvolveAllGenerations();
        
        //_float_prior_archive.push_back( current_mutrate_vector);
        
        
        // SimulationResult sim_result(local_sim_object);
        
        // keep only the generations we need to conserve memory
        auto tmp_actual_generations = _actual_data_position.GetActualGenerations();
        SimulationResult sim_result(local_sim_object, tmp_actual_generations);
        
        sim_result.prior_sample_index = current_prior_sample_idx;
        
        //tmp_res_vector.push_back(std::move(sim_result));
        tmp_res_vector.push_back(sim_result);
        
    } // mutation rate
    
    
    return tmp_res_vector;
    
} // RunMutrateInferenceBatch


