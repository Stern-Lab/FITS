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

std::vector<SimulationResult> clsCMulatorABC::RunFitnessInferenceBatch( std::size_t num_simulations )
{
    // initialization
    CMulator local_sim_object(_zparams);
    
    if ( !local_sim_object.IsAbleToInferFitness() ) {
        std::cerr << "Not enough parameters to infer fitness" << std::endl;
        throw "Not enough parameters to infer fitness";
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
    
    
    // composite is the default for fitness
    auto tmp_prior = _zparams.GetString( fits_constants::PARAM_PRIOR_DISTRIB,
                                         fits_constants::PARAM_PRIOR_DISTRIB_DEFAULT );
    
    if ( tmp_prior.compare( fits_constants::PARAM_PRIOR_DISTRIB_UNIFORM ) == 0 ) {
        _prior_type = UNIFORM;
    }
    else if ( tmp_prior.compare( fits_constants::PARAM_PRIOR_DISTRIB_COMPOSITE ) == 0 ) {
        _prior_type = FITNESS_COMPOSITE;
    }
    else if ( tmp_prior.compare( fits_constants::PARAM_PRIOR_DISTRIB_SMOOTHED_COMPOSITE ) == 0 ) {
        _prior_type = SMOOTHED_COMPOSITE;
    }
    else if ( tmp_prior.compare( fits_constants::PARAM_PRIOR_DISTRIB_LOGNORMAL ) == 0 ) {
        _prior_type = FITNESS_LOGNORMAL;
    }
    else {
        std::cerr << "Unkown prior distribution: " << tmp_prior << ". Setting to uniform as default." << std::endl;
    }

    auto min_fitness_vec = local_sim_object.GetAlleleMinFitnessValues();
    auto max_fitness_vec = local_sim_object.GetAlleleMaxFitnessValues();
    
    PriorSampler<FLOAT_TYPE> sampler(min_fitness_vec, max_fitness_vec, _prior_type);

    auto fitness_vector_list = sampler.SamplePrior(num_simulations);


    
    if ( _zparams.GetInt( "Debug", 0 ) > 0 ) {
        std::cout << "BEGIN Debug: Prior distribution" << std::endl;
        std::cout << "================================" << std::endl;
        
        for ( auto current_fitness_vector : fitness_vector_list ) {
            
            for ( auto current_fitness_val : current_fitness_vector ) {
                
                std::cout << current_fitness_val << "\t";
            }
            
            std::cout << std::endl;
        }
        
        std::cout << "END Debug: Prior distribution" << std::endl;
        std::cout << "==============================" << std::endl;
    }
    
    std::vector<SimulationResult> tmp_res_vector;
    
    
    // simulation for each set of parameters
    for (auto current_fitness_vector : fitness_vector_list) {

        local_sim_object.Reset_Soft();
        local_sim_object.SetFitnessValues(current_fitness_vector);
        local_sim_object.EvolveAllGenerations();
        
        _float_prior_archive.push_back( current_fitness_vector );
                
        // keep only the generations we need, to conserve memory 2017-04-02
        auto tmp_actual_generations = _actual_data_position.GetActualGenerations();
        SimulationResult sim_result(local_sim_object, tmp_actual_generations);
        
        tmp_res_vector.push_back(sim_result);
        //tmp_res_vector.push_back(std::move(sim_result));
        
    } // fitness
    
    
    return tmp_res_vector;
    
} // RunFitnessInferenceBatch


