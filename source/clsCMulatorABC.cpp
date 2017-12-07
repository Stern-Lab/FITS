/**********************
 CMulator Core Methods
 **********************/

// #define DEBUG_VERBOSE

#include "clsCMulatorABC.h"
#include "fits_constants.h"

clsCMulatorABC::clsCMulatorABC( ZParams sim_params, ActualDataFile actual_data_file ) :
_zparams(sim_params),
_total_running_time(0.0f),
_prior_type(PriorDistributionType::UNIFORM),
_actual_data_file(actual_data_file),
_simulation_result_vector(),
_float_prior_archive(),
_int_prior_archive(),
_use_rejection_threshold(true)
{
    ResetRejectionThreshold();
    
    _repeats = _zparams.GetInt( "_num_repeats" );
    
    _num_alleles = _zparams.GetUnsignedInt( "_num_alleles", 0 );
    
    _sims_to_keep = _zparams.GetFloat( fits_constants::PARAM_ACCEPTANCE_RATE,
                                       fits_constants::ACCEPTANCE_RATE_DEFAULT ) * _repeats;
    
    _rejection_threshold = _zparams.GetFloat( fits_constants::PARAM_REJECTION_THRESHOLD,
                                             -1.0f );
    if ( _sims_to_keep == 0.0f ) {
        _sims_to_keep = _zparams.GetFloat( fits_constants::PARAM_ACCEPTANCE_LIMIT,
                                        fits_constants::ACCEPTANCE_LIMIT_DEFAULT ) * _repeats;
    }
    
    // TODO: check what rejection method we use, verify that at least one is given
    
    _selected_actual_generations = actual_data_file.GetActualGenerations();
    _actual_data_raw_freqs = actual_data_file.GetActualFrequencies();
    
    auto _time_for_seeding = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    _boost_gen.seed(_time_for_seeding);
}


void clsCMulatorABC::RunABCInference( FactorToInfer factor, std::size_t number_of_batches )
{
    _simulation_result_vector.clear();
    _float_prior_archive.clear();
    _int_prior_archive.clear();
    
    
    
    // performace tracking
    boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::mean>> rate_stats;
    
    auto remaining_repeats = _repeats;
    auto repeats_in_batch = _repeats / number_of_batches;
    
    _simulation_result_vector.reserve(_repeats);
    
    // std::size_t threshold_update_counter = 0;
    
    
    auto start_global = std::chrono::high_resolution_clock::now();
    
    std::string completion_eta_str = "";
    
    while (remaining_repeats > 0) {
        
        if (repeats_in_batch > remaining_repeats) {
            repeats_in_batch = remaining_repeats;
        }
        remaining_repeats -= repeats_in_batch;
        
        std::vector<SimulationResult> tmp_result_vector;
        
        auto start = std::chrono::high_resolution_clock::now();
        
        _factor_to_infer = factor;
        
        // TODO: TRY TO USE ASYNC - doesn't work. inconsistent - one time works the other one not.
        switch (_factor_to_infer) {
            case Fitness:
                tmp_result_vector = RunFitnessInferenceBatch(repeats_in_batch);
                break;
                
            case PopulationSize:
                tmp_result_vector = RunPopulationSizeInferenceBatch(repeats_in_batch);
                break;
                
            case MutationRate:
                tmp_result_vector = RunMutationInferenceBatch(repeats_in_batch);
                break;
                
            case Generations:
                tmp_result_vector = RunGenerationInferenceBatch(repeats_in_batch);
                break;
                
            case None:
                throw "Factor to infer is not defined (none)";
        }
        auto end = std::chrono::high_resolution_clock::now();
        
        
        // move results from chunk to the main result vector
        // std::cout << "moving.... ";
        //for ( auto &&tmp_result : tmp_result_vector ) {
        //    _simulation_result_vector.push_back(std::move(tmp_result));
        //}
        for ( auto &tmp_result : tmp_result_vector ) {
            _simulation_result_vector.push_back(tmp_result);
        }
        // std::cout << " Done." << std::endl;
        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        
        
        auto sim_count = repeats_in_batch;
        
        auto calculated_speed = static_cast<double>(sim_count) / static_cast<double>(elapsed_ms.count());
        calculated_speed *= 1000; // 1/milisecond to 1/second
        rate_stats(calculated_speed);
        
        // estimate time for completion in seconds
        if (completion_eta_str.empty()) {
            
            auto seconds_remaining = remaining_repeats / static_cast<int>(std::round(calculated_speed));
            auto duration_remaining = std::chrono::seconds(seconds_remaining);
            auto current_time = std::chrono::system_clock::now();
            auto completion_ETA = current_time + duration_remaining;
            auto completion_ETA_timet = std::chrono::system_clock::to_time_t(completion_ETA);
            auto completion_ETA_tm = *std::localtime(&completion_ETA_timet);
        
            std::stringstream ss;
            ss << std::put_time(&completion_ETA_tm, "%c");
            completion_eta_str = ss.str();
        }
        
        
        std::cout << "Remaining " << remaining_repeats
        << " repeats (rate of "
        << std::round(calculated_speed)
        << "sim/sec) - ETA is "
        //<< std::put_time(&completion_ETA_tm, "%c")
        << completion_eta_str
        << std::endl;
        
    } // while (remaining_repeats > 0)
    
    
    /*
     // NOT RELEVANT ANYMORE
     
    std::cout << "Freeing up memory... ";
    std::cout << _simulation_result_vector.capacity() << " -> ";
    _simulation_result_vector.shrink_to_fit();
    std::cout << _simulation_result_vector.capacity();
    std::cout << "Done." << std::endl;
     */
    
    /*
    if (_simulation_result_vector.size() > _sims_to_keep) {
        
        std::cout << "Total results count is " << _simulation_result_vector.size()
        << ", capacity is " << _simulation_result_vector.capacity() << std::endl;
        
        // added 2016-09-04 maybe this will help to prevent page fault previouly caused by large vectors
        _simulation_result_vector.shrink_to_fit();
        
        std::cout << "Final results count is " << _simulation_result_vector.size() << ", capacity is " << _simulation_result_vector.capacity() << std::endl;
    }
    */
    
    std::cout << "Calculating distance... ";
    
    // scaling factor; default actualling means no scaling (divide by 1)
    std::vector<FLOAT_TYPE> scaling_vector(0);
    
    // std::cout << "scaling vector" << scaling_vector.size() << std::endl;
    
    auto scaling_option_str = _zparams.GetString( fits_constants::PARAM_SCALING,
                                             fits_constants::PARAM_SCALING_DEFAULT );
    
    if ( scaling_option_str.compare(fits_constants::PARAM_SCALING_SD) == 0 ) {
        scaling_vector = GetSDPerAllele(0, _simulation_result_vector.size());
    }
    
    if ( scaling_option_str.compare(fits_constants::PARAM_SCALING_MAD) == 0) {
        scaling_vector = GetMADPerAllele(0, _simulation_result_vector.size());
    }
    
    //std::cout << "scaling vector" << scaling_vector.size() << std::endl;
    //auto scaling_vector = GetMADPerAllele(0, _simulation_result_vector.size());
    
    
    // get actual data matrix
    auto actual_generations = _actual_data_file.GetActualGenerations();
    auto actual_matrix = _actual_data_file.GetActualFreqsAsMatrix();
    
    
    // from each simulation result, get only the relevant generations
    for ( auto current_idx=0; current_idx<_simulation_result_vector.size(); ++current_idx ) {
    //for ( auto& current_result : _simulation_result_vector )
        _simulation_result_vector[current_idx].distance_from_actual =
            GetDistanceSimActual(actual_matrix, _simulation_result_vector[current_idx].sim_data_matrix, scaling_vector );
        // std::cout << " distance in current result = " << current_result.distance_from_actual << std::endl;
    }

    std::cout << "Done." << std::endl;
    
    std::cout << "Sorting results... ";
    std::sort(_simulation_result_vector.begin(), _simulation_result_vector.end());
    std::cout << "Done." << std::endl;
    
    auto end_global = std::chrono::high_resolution_clock::now();
    auto global_elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_global - start_global);
    _total_running_time = static_cast<double>(global_elapsed_ms.count()) / 1000.0;
    
    std::cout << std::endl;
    std::cout << "Total running time: " << _total_running_time << " seconds" << std::endl;
    std::cout << "Average simulation rate: " << boost::accumulators::mean(rate_stats) << " simulations/second" << std::endl;
}


FLOAT_TYPE clsCMulatorABC::ResetRejectionThreshold()
{
    return SetRejectionThreshold(THRESHOLD_RESET_VALUE);
}


FLOAT_TYPE clsCMulatorABC::GetRejectionThreshold()
{
    return _rejection_threshold;
}


FLOAT_TYPE clsCMulatorABC::SetRejectionThreshold(FLOAT_TYPE new_threshold)
{
    FLOAT_TYPE previous_threshold = _rejection_threshold;
    
    _rejection_threshold = new_threshold;
    
    return previous_threshold;
}


std::size_t clsCMulatorABC::GetRepeats()
{
    return _repeats;
}

std::size_t clsCMulatorABC::GetNumberOfKeptResults()
{
    return _sims_to_keep;
}

std::vector<SimulationResult> clsCMulatorABC::GetResultsVector(bool only_accepted_results)
{
    if (only_accepted_results) {
        
        if ( _rejection_threshold > 0.0f ) {
            std::cout << "Using rejection threshold" << _rejection_threshold << std::endl;
            
            if (!std::is_sorted(_simulation_result_vector.cbegin(),
                                _simulation_result_vector.cend() ) ) {
                std::cerr << "GetResultsVector: results vector is not sorted!" << std::endl;
                throw "GetResultsVector: results vector is not sorted!";
            }
            
            std::cout << "test" << std::endl;
            int current_idx=0;
            while ( current_idx < _simulation_result_vector.size() &&
                   _simulation_result_vector[current_idx].distance_from_actual < _rejection_threshold ) {
                std::cout << "distance=" << _simulation_result_vector[current_idx].distance_from_actual << std::endl;
                ++current_idx;
            }
            _sims_to_keep = current_idx - 1;
            
            if ( _sims_to_keep > _simulation_result_vector.size() ) {
                _sims_to_keep = _simulation_result_vector.size();
            }
        }
        
        // find the first item to hold distance>0
        std::size_t first_nonzero_distance_idx = 0;
        
        /*
        while ( first_nonzero_distance_idx<_simulation_result_vector.size() &&
               _simulation_result_vector[first_nonzero_distance_idx].distance_from_actual==0.0f ) {
            ++first_nonzero_distance_idx;
        }
        */
        
        std::nth_element(_simulation_result_vector.begin() + first_nonzero_distance_idx,
                         _simulation_result_vector.begin() + first_nonzero_distance_idx + _sims_to_keep,
                         _simulation_result_vector.end());
        
        std::vector<SimulationResult> tmp_vec( _simulation_result_vector.begin() + first_nonzero_distance_idx,
                                               _simulation_result_vector.begin() + first_nonzero_distance_idx + _sims_to_keep );
        
        // std::cout << "sims " << _sims_to_keep << std::endl;
        // std::cout << "tmpvec " << tmp_vec.size() << std::endl;
        
        return tmp_vec;
    }
    
    return _simulation_result_vector;
}

void clsCMulatorABC::SetImmediateRejection(bool new_val)
{
    _use_rejection_threshold = new_val;
}


FLOAT_TYPE clsCMulatorABC::GetDistanceSimActual( const MATRIX_TYPE &actual_data, const MATRIX_TYPE &sim_data, const std::vector<FLOAT_TYPE> &scaling_vector )
{
    
    if ( actual_data.size1() != sim_data.size1() ) {
        std::cerr << "DistanceSimActual - matrices don't match in size1: actual " <<
        actual_data.size1() << "vs sim " << sim_data.size1() << std::endl;
        throw "DistanceSimActual - matrices don't match in size1";
    }
    
    if ( actual_data.size2() != sim_data.size2() ) {
        std::cerr << "DistanceSimActual - matrices don't match in size2: actual " <<
        actual_data.size2() << "vs sim " << sim_data.size2() << std::endl;
        throw "DistanceSimActual - matrices don't match in size2";
    }
    
    MATRIX_TYPE diff_matrix = actual_data - sim_data;
    
    //std::cout << "actual matrix: " << actual_data << std::endl;
    //std::cout << "sim matrix: " << sim_data << std::endl;
    
    diff_matrix = boost::numeric::ublas::abs(diff_matrix);
    
    // scale before summing up data
    if (scaling_vector.size() > 0 ) {
        for ( auto col=0; col<diff_matrix.size2(); ++col ) {
            boost::numeric::ublas::matrix_column<MATRIX_TYPE> current_col(diff_matrix, col);
            current_col = current_col / scaling_vector[col];
            //std::cout << "scaling by " << scaling_vector[col] << std::endl;
            
            // current_col = boost::numeric::ublas::element_prod( current_col, current_col );
        }
    }
    
    // sum each allele
    FLOAT_TYPE sum_diff = 0.0f;
    for ( auto col=0; col<diff_matrix.size2(); ++col ) {
        
        boost::numeric::ublas::matrix_column<MATRIX_TYPE> current_col(diff_matrix, col);
        auto current_allele_sum = boost::numeric::ublas::sum(current_col);
        
        sum_diff += current_allele_sum;
    }
    // trick - sum can be for vectors, not matrices. multiply matrix by unit vector to get a vector
    //auto sum_diff = sum(prod(boost::numeric::ublas::scalar_vector<float>(diff_matrix.size1()), diff_matrix));
    
    //sum_diff = std::sqrt(sum_diff);
    if ( sum_diff == 0.0f ) {
        std::cerr << "Warning: actual data identical to simulated data (distance=0)." << std::endl;
        std::cerr << "matrix data: " << actual_data << std::endl;
    }
    return sum_diff;
}

std::vector< std::vector<FLOAT_TYPE>> clsCMulatorABC::GetPriorFloat()
{
    return _float_prior_archive;
}

std::vector< std::vector<int>> clsCMulatorABC::GetPriorInt()
{
    return _int_prior_archive;
}

std::string clsCMulatorABC::GetPriorFloatAsString()
{
    std::string tmp_str = "";
    for ( auto current_vec : _float_prior_archive ) {
        
        for ( auto current_val : current_vec ) {
            
            tmp_str += std::to_string(current_val) + "\t";
        }
        tmp_str = tmp_str + "\n";
    }
    return tmp_str;
}

std::string clsCMulatorABC::GetPriorIntAsString()
{
    std::string tmp_str = "";
    for ( auto current_vec : _int_prior_archive ) {
        
        for ( auto current_val : current_vec ) {
            
            tmp_str += std::to_string(current_val) + "\t";
            std::cout << tmp_str << "\t";
        }
        tmp_str = tmp_str + "\n";
    }
    return tmp_str;
}
