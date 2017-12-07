
#include "clsCMulatorABC.h"

std::vector<SimulationResult> clsCMulatorABC::RunMutationInferenceBatch( std::size_t num_simulations )
{
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
        std::cerr << "Not enough parameters to infer mutation rate" << std::endl;
        throw "Not enough parameters to infer mutation rate";
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
    
    auto min_matrix = local_sim_object.GetMinMutationRateMatrix();
    auto max_matrix = local_sim_object.GetMaxMutationRateMatrix();
    
    // 2017-03-20 - change sampling from float to int
    PriorSampler<int> sampler( min_matrix, max_matrix, PriorDistributionType::UNIFORM);
    
    /*
     // convert matrices to vectors in order to use same sampling mechanism for all inferred parameters
     std::vector<FLOAT_TYPE> min_mutrates;
     std::vector<FLOAT_TYPE> max_mutrates;
    for ( auto row=0; row<local_sim_object.GetAlleleNumber(); ++row ) {
        for ( auto col=0; col<local_sim_object.GetAlleleNumber(); ++col ) {
            
            min_mutrates.push_back( min_matrix(row,col) );
            max_mutrates.push_back( max_matrix(row,col) );
        }
    }
    
    PriorSampler<FLOAT_TYPE> sampler( min_mutrates, max_mutrates, PriorDistributionType::UNIFORM);
    */
    auto mutrate_vector_list = sampler.SamplePrior( num_simulations );
    
    
    std::vector<SimulationResult> tmp_res_vector;
    
    // std::cout << "Prior archive size : " << _prior_archive.size() << std::endl;
    // simulation for each set of parameters
    // 2017-02-07 changed to regular for loop to be able to record index
    //for (auto current_mutrate_vector : mutrate_vector_list) {
    for (auto current_mutrate_idx=0; current_mutrate_idx<mutrate_vector_list.size(); ++current_mutrate_idx ) {
        
        auto current_mutrate_vector = mutrate_vector_list[current_mutrate_idx];
        
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
            
            // power of the mutation rate
            
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
        
        _int_prior_archive.push_back( current_mutrate_vector);
        
        
        // SimulationResult sim_result(local_sim_object);
        
        // right here keep only the generations we need
        // to conserve memory 2017-04-02
        auto tmp_actual_generations = _actual_data_file.GetActualGenerations();
        SimulationResult sim_result(local_sim_object, tmp_actual_generations);
        
        //tmp_res_vector.push_back(std::move(sim_result));
        tmp_res_vector.push_back(sim_result);
        
    } // mutation rate
    
    
    return tmp_res_vector;
    
} // RunMutrateInferenceBatch


