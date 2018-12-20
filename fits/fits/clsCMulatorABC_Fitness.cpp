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

std::vector<SimulationResult> clsCMulatorABC::RunFitnessInferenceBatch( const PRIOR_DISTRIB &prior_distrib, std::size_t start_idx, std::size_t end_idx  )
{
    
    VerifyIndece( prior_distrib, start_idx, end_idx );
    
    // initialization
    CMulator local_sim_object(_zparams);
    
    if ( !local_sim_object.IsAbleToInferFitness() ) {
        std::string tmp_str = "Not enough parameters to infer fitness";
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
    
    if ( _zparams.GetInt( "Debug", 0 ) > 0 ) {
        std::cout << "first generation (shift) = " << first_generation << std::endl;
        std::cout << "last generation = " << last_generation << std::endl;
        std::cout << "number of generations = " << num_generations << std::endl;
    }
    local_sim_object.SetWTAllele( _actual_data_position.GetWTIndex() );
    
    
    if ( _zparams.GetInt( "Debug", 0 ) > 0 ) {
        std::cout << "BEGIN Debug: Prior distribution" << std::endl;
        std::cout << "================================" << std::endl;
        
        for ( auto current_fitness_vector : prior_distrib ) {
            
            for ( auto current_fitness_val : current_fitness_vector ) {
                
                std::cout << current_fitness_val << fits_constants::FILE_FIELD_DELIMITER;
            }
            
            std::cout << std::endl;
        }
        
        std::cout << "END Debug: Prior distribution" << std::endl;
        std::cout << "==============================" << std::endl;
    }
    
    std::vector<SimulationResult> tmp_res_vector;
    
    
    // simulation for each set of parameters
    //for (auto current_fitness_vector : prior_distrib) {
    for ( auto current_prior_sample_idx=start_idx; current_prior_sample_idx<end_idx; ++current_prior_sample_idx ) {

        local_sim_object.Reset_Soft();
        //local_sim_object.SetFitnessValues( current_fitness_vector );
        local_sim_object.SetFitnessValues( prior_distrib[current_prior_sample_idx] );
        local_sim_object.EvolveAllGenerations();
        
        //_float_prior_archive.push_back( current_fitness_vector );
                
        // keep only the generations we need, to conserve memory 2017-04-02
        auto tmp_actual_generations = _actual_data_position.GetActualGenerations();
        SimulationResult sim_result(local_sim_object, tmp_actual_generations);
        
        sim_result.prior_sample_index = current_prior_sample_idx;
        
        tmp_res_vector.push_back(sim_result);
        //tmp_res_vector.push_back(std::move(sim_result));
        
    } // fitness
    
    
    return tmp_res_vector;
    
} // RunFitnessInferenceBatch


