
#include "clsCMulatorABC.h"

std::vector<SimulationResult> clsCMulatorABC::RunPopulationSizeInferenceBatch( std::size_t num_simulations )
{
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
    

    if ( !local_sim_object.IsAbleToInferPopulationSize() ) {
        std::cerr << "Not enough parameters to infer population size" << std::endl;
        throw "Not enough parameters to infer population size";
    }
    
    // set initial frequencies from actual data
    auto init_freq_vec = _actual_data_file.GetInitFreqs();
    for (auto i = 0; i < init_freq_vec.size(); ++i) {
        local_sim_object.SetAlleleInitFreq(i, init_freq_vec[i]);
    }
    
    auto first_generation = _actual_data_file.GetFirstGeneration();
    auto last_generation = _actual_data_file.GetLastGeneration();
    auto num_generations = last_generation - first_generation + 1;
    local_sim_object.SetGenerationShift(first_generation);
    local_sim_object.SetNumOfGeneration(num_generations);
    
    // identify wt
    auto wt_allele_it = std::max_element(init_freq_vec.begin(), init_freq_vec.end());
    auto wt_allele_idx = static_cast<unsigned int>(std::distance(init_freq_vec.begin(), wt_allele_it));
    
    local_sim_object.SetWTAllele(wt_allele_idx, 0.0f, 2.0f);
    //std::cout << " WT allele set to " << wt_allele_idx << std::endl;
    
    
    // std::vector<int> minN {_zparams.GetInt("_Nlog_min")};
    // std::vector<int> maxN {_zparams.GetInt("_Nlog_max")};
    
    std::vector<int> minN {_zparams.GetInt(fits_constants::PARAM_MIN_LOG_POPSIZE)};
    std::vector<int> maxN {_zparams.GetInt(fits_constants::PARAM_MAX_LOG_POPSIZE)};
    
    
    PriorSampler<int> sampler( minN, maxN, PriorDistributionType::UNIFORM );
    
    //auto fitness_vector_list = _fitness_range.GetRandomCombinations(num_simulations, false);
    auto popsize_vector_list = sampler.SamplePrior(num_simulations);
    
    
    std::vector<SimulationResult> tmp_res_vector;
    
    // simulation for each set of parameters
    for (auto current_popsize : popsize_vector_list) {
        
        local_sim_object.Reset_Soft();
        local_sim_object.SetPopulationSize( std::pow( 10, current_popsize[0]) );
        local_sim_object.EvolveAllGenerations();
        
        // patchy, I know, TODO make it better 2017-02-23
        std::vector<int> dummy_popsize_vector(_num_alleles, 0);
        dummy_popsize_vector[0] = static_cast<int>(current_popsize[0]);
        _int_prior_archive.push_back( dummy_popsize_vector );
        
        
        // SimulationResult sim_result(local_sim_object);
        
        // right here keep only the generations we need
        // to conserve memory 2017-04-02
        auto tmp_actual_generations = _actual_data_file.GetActualGenerations();
        SimulationResult sim_result(local_sim_object, tmp_actual_generations);
        
        //tmp_res_vector.push_back(std::move(sim_result));
        tmp_res_vector.push_back(sim_result);
        
        
    } // fitness
    
    
    return tmp_res_vector;
    
} // RunPopsizeInferenceBatch


