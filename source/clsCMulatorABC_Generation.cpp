
#include "clsCMulatorABC.h"



// old version - generate trajectory for each iteration
// new version - generate one trajectory for each batch, then compare using different sets of generations


// NEW
std::vector<SimulationResult> clsCMulatorABC::RunGenerationInferenceBatch( std::size_t num_simulations )
{
    // initialization
    CMulator local_sim_object(_zparams);
    
    if ( !local_sim_object.IsAbleToInferGeneration() ) {
        std::cerr << "Not enough parameters to infer generations" << std::endl;
        throw "Not enough parameters to infer generations";
    }
    
    // set initial frequencies from actual data
    auto init_freq_vec = _actual_data_file.GetInitFreqs();
    for (auto i = 0; i < init_freq_vec.size(); ++i) {
        local_sim_object.SetAlleleInitFreq(i, init_freq_vec[i]);
    }
    
    std::vector<int> min_vec;
    std::vector<int> max_vec;
    
    auto actual_generations = _actual_data_file.GetActualGenerations(true);
    
    if ( actual_generations.size() > 2 ) {
        std::cout << "Warning: More than two time points detected. Make sure the interval between them is uniform." << std::endl;
    }
    
    auto min_generation = _zparams.GetInt(fits_constants::PARAM_MIN_FIXED_INT_GENERATION);
    auto max_generation = _zparams.GetInt(fits_constants::PARAM_MAX_FIXED_INT_GENERATION);
    min_vec.push_back(min_generation);
    max_vec.push_back(max_generation);
    
    PriorSampler<int> sampler(min_vec, max_vec, UNIFORM);
    
    auto fitness_vector_list = sampler.SamplePrior(num_simulations);
    
    std::vector<SimulationResult> tmp_res_vector;
    
    auto first_generation = _actual_data_file.GetFirstGeneration();
    
    //auto num_given_generations = actual_generations.size();
    local_sim_object.SetWTAllele( _actual_data_file.GetWTIndex() );
    local_sim_object.SetGenerationShift(first_generation);
    
    
    // single simulation for this batch
    local_sim_object.Reset_Soft();
    local_sim_object.SetNumOfGeneration(max_generation);
    local_sim_object.EvolveAllGenerations();
    
    for (auto current_interval_vector : fitness_vector_list) {
        
        auto current_interval = current_interval_vector[0];
        auto current_generations_vec(actual_generations);
        std::vector<int> interval_vector;
        interval_vector.push_back(current_interval);
        
        for ( auto idx=1; idx<current_generations_vec.size(); ++idx ) {
            current_generations_vec[idx] = current_generations_vec[idx-1] + current_interval;
        }
        
        auto last_generation = *(current_generations_vec.cend()-1);
        auto num_generations = last_generation - first_generation + 1;
        
        /*
         std::cout << "Generations: ";
         for ( auto gen : current_generations_vec ) {
         std::cout << gen << ", ";
         }
         std::cout << std::endl;
         */
        
        _int_prior_archive.push_back( interval_vector );
        
        // keep only the generations we need, to conserve memory 2017-04-02
        SimulationResult sim_result(local_sim_object, current_generations_vec);
        sim_result.generation_interval = current_interval;
        sim_result.actual_generations = current_generations_vec;
        tmp_res_vector.push_back(sim_result);
        
    } // generations
    
    
    return tmp_res_vector;
    
} // RunFitnessInferenceBatch


/*
 Je Olde
 
std::vector<SimulationResult> clsCMulatorABC::RunGenerationInferenceBatch( std::size_t num_simulations )
{
    // initialization
    CMulator local_sim_object(_zparams);
    
    if ( !local_sim_object.IsAbleToInferGeneration() ) {
        std::cerr << "Not enough parameters to infer generations" << std::endl;
        throw "Not enough parameters to infer generations";
    }
    
    // set initial frequencies from actual data
    auto init_freq_vec = _actual_data_file.GetInitFreqs();
    for (auto i = 0; i < init_freq_vec.size(); ++i) {
        local_sim_object.SetAlleleInitFreq(i, init_freq_vec[i]);
    }
    
    std::vector<int> min_vec;
    std::vector<int> max_vec;
    
    auto actual_generations = _actual_data_file.GetActualGenerations(true);
    
    if ( actual_generations.size() > 2 ) {
        std::cout << "Warning: More than two time points detected. Make sure the interval between them is uniform." << std::endl;
    }
    
    min_vec.push_back(_zparams.GetInt(fits_constants::PARAM_MIN_FIXED_INT_GENERATION));
    max_vec.push_back(_zparams.GetInt(fits_constants::PARAM_MAX_FIXED_INT_GENERATION));
    
    
    PriorSampler<int> sampler(min_vec, max_vec, UNIFORM);
    
    auto fitness_vector_list = sampler.SamplePrior(num_simulations);
    
    std::vector<SimulationResult> tmp_res_vector;
    
    auto first_generation = _actual_data_file.GetFirstGeneration();
    //auto first_generation = 1;
    
    //auto num_given_generations = actual_generations.size();
    local_sim_object.SetWTAllele( _actual_data_file.GetWTIndex() );
    local_sim_object.SetGenerationShift(first_generation);
    
    // should be set when initializing cmulator from zparams
    //local_sim_object.SetFitnessValues(current_interval_vector);
    
    // simulation for each set of parameters
    for (auto current_interval_vector : fitness_vector_list) {
        
        local_sim_object.Reset_Soft();
        
        auto current_interval = current_interval_vector[0];
        auto current_generations_vec(actual_generations);
        std::vector<int> interval_vector;
        interval_vector.push_back(current_interval);
        
        for ( auto idx=1; idx<current_generations_vec.size(); ++idx ) {
            current_generations_vec[idx] = current_generations_vec[idx-1] + current_interval;
        }
        
        
        auto last_generation = *(current_generations_vec.cend()-1);
        auto num_generations = last_generation - first_generation + 1;
        //std::cout << "setting number of generations to " << num_generations << std::endl;
        local_sim_object.SetNumOfGeneration(num_generations+1);
        
 
        
        local_sim_object.EvolveAllGenerations();
        
        
        _int_prior_archive.push_back( interval_vector );
        
        // keep only the generations we need, to conserve memory 2017-04-02
        
        SimulationResult sim_result(local_sim_object, current_generations_vec);
        sim_result.generation_interval = current_interval;
        
        sim_result.actual_generations = current_generations_vec;
        //tmp_res_vector.push_back(std::move(sim_result));
        tmp_res_vector.push_back(sim_result);
        
    } // generations
    
    
    return tmp_res_vector;
    
} // RunFitnessInferenceBatch


*/
