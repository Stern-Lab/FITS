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

#include "clsCMulatorABC.h"

std::vector<SimulationResult> clsCMulatorABC::RunPopulationSizeInferenceBatch( const PRIOR_DISTRIB_VECTOR &prior_distrib, std::size_t start_idx, std::size_t end_idx )
{
    VerifyIndece( prior_distrib, start_idx, end_idx );
    
    // initialization
    CMulator local_sim_object(_zparams);
    
    /*
    if ( local_sim_object.IsObservedDataUsed() ) {
        std::cout << "Using observed (sampled) data for comparison." << std::endl;
    }
    else {
        std::cout << "Using all data for comparison." << std::endl;
    }
    */
    
    local_sim_object.AssertAbleToInferPopulationSize();
    
    // set initial frequencies from actual data
    auto init_freq_vec = _actual_data_position.GetInitFreqs();
    for (auto i = 0; i < init_freq_vec.size(); ++i) {
        local_sim_object.SetAlleleInitFreq(i, init_freq_vec[i]);
    }
    
    //auto first_generation = _actual_data_position.GetFirstGeneration();
    //auto last_generation = _actual_data_position.GetLastGeneration();
    // auto num_generations = last_generation - first_generation + 1;
    // local_sim_object.SetGenerationShift(first_generation);
    // local_sim_object.SetNumOfGeneration(num_generations);
    
    // identify wt
    auto wt_allele_it = std::max_element(init_freq_vec.begin(), init_freq_vec.end());
    auto wt_allele_idx = static_cast<unsigned int>(std::distance(init_freq_vec.begin(), wt_allele_it));
    
    local_sim_object.SetWTAllele(wt_allele_idx, 0.0f, 2.0f);
    
    //std::vector<FLOAT_TYPE> minN {_zparams.GetDouble(fits_constants::PARAM_MIN_LOG_POPSIZE)};
    //std::vector<FLOAT_TYPE> maxN {_zparams.GetDouble(fits_constants::PARAM_MAX_LOG_POPSIZE)};
    
    //PriorSampler<FLOAT_TYPE> sampler( minN, maxN, PriorDistributionType::UNIFORM );
    
    //auto popsize_vector_list = sampler.SamplePrior(num_simulations);
    
    std::vector<SimulationResult> tmp_res_vector;
    
    // simulation for each set of parameters
    //for (auto current_popsize : prior_distrib) {
    for ( auto current_prior_sample_idx=start_idx; current_prior_sample_idx<end_idx; ++current_prior_sample_idx ) {
        
        if ( _expected_prior_sample_size != prior_distrib[current_prior_sample_idx].size() ) {
            std::string tmp_str = "Expected sample size is " + std::to_string(_expected_prior_sample_size)
            + " but current sample size is " + std::to_string(prior_distrib[current_prior_sample_idx].size()) +
            ". Sample: ";
            
            for ( auto val : prior_distrib[current_prior_sample_idx] ) {
                tmp_str += std::to_string(val) + " ";
            }
            
            throw tmp_str;
        }
        
        local_sim_object.Reset_Soft();
        //local_sim_object.SetPopulationSize( std::pow( 10, current_popsize[0]) );
        
        auto current_popsize_log10 = prior_distrib[current_prior_sample_idx][0];
        
        local_sim_object.SetPopulationSize( std::pow( 10, current_popsize_log10 ) );
        local_sim_object.EvolveAllGenerations();
        
        // 2018-12-03
        //std::vector<FLOAT_TYPE> dummy_popsize_vector(_num_alleles, 0);
        //std::vector<FLOAT_TYPE> dummy_popsize_vector(1, 0);
        //dummy_popsize_vector[0] = static_cast<FLOAT_TYPE>(current_popsize[0]);
        //_float_prior_archive.push_back( dummy_popsize_vector );
        
        // keep only the generations we need to conserve memory
        auto tmp_actual_generations = _actual_data_position.GetActualGenerations();
        SimulationResult sim_result(local_sim_object, tmp_actual_generations);
        
        sim_result.prior_sample_index = current_prior_sample_idx;
        
        //tmp_res_vector.push_back(std::move(sim_result));
        tmp_res_vector.push_back(sim_result);
        
        
    } // fitness
    
    
    return tmp_res_vector;
    
} // RunPopsizeInferenceBatch


