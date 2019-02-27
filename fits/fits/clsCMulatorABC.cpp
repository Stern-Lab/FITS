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
#include "fits_constants.h"

clsCMulatorABC::clsCMulatorABC()
{}

clsCMulatorABC::clsCMulatorABC( ZParams sim_params, ActualDataPositionData actual_data_position, FactorToInfer factor_to_infer, const PRIOR_DISTRIB_VECTOR& prior_distribution_vec, const PRIOR_DISTRIB_MATRIX& prior_distribution_matrix  ) :
_zparams(sim_params),
_total_running_time_sec(0),
_prior_type(PriorDistributionType::UNIFORM),
_actual_data_position(actual_data_position),
_simulation_result_vector(),
//_float_prior_archive(),
_global_prior(),
_use_rejection_threshold(true),
_factor_to_infer(factor_to_infer),
_simulation_speed(0.0f),
_verbose_output(false)
{
    ResetRejectionThreshold();
    
    _wt_allele_idx = actual_data_position.GetWTIndex();
    
    _verbose_output = ( sim_params.GetInt( fits_constants::PARAM_VERBOSE_SWITCH ) == fits_constants::PARAM_VERBOSE_SWITCH_ON );
    
    
    // _repeats = _zparams.GetInt( fits_constants::PARAM_SIM_REPEATS, fits_constants::PARAM_SIM_REPEATS_DEFAULT );
    _repeats = _zparams.GetSize_t( fits_constants::PARAM_SIM_REPEATS, fits_constants::PARAM_SIM_REPEATS_DEFAULT );
    if ( _repeats <= 0 ) {
        std::string tmp_str = "Error: Number of repoeats must be positive.";
        throw tmp_str;
    }
    
    _num_alleles = _zparams.GetUnsignedInt( fits_constants::PARAM_NUM_ALLELES, -1 );
    
    _sims_to_keep = _zparams.GetFloat( fits_constants::PARAM_ACCEPTANCE_RATE,
                                      fits_constants::ACCEPTANCE_RATE_DEFAULT ) * _repeats;
    
    _rejection_threshold = _zparams.GetFloat( fits_constants::PARAM_REJECTION_THRESHOLD,
                                             -1.0f );
    if ( _sims_to_keep == 0.0f ) {
        _sims_to_keep = _zparams.GetFloat( fits_constants::PARAM_ACCEPTANCE_LIMIT,
                                          fits_constants::ACCEPTANCE_LIMIT_DEFAULT ) * _repeats;
    }
    
    
    _selected_actual_generations = _actual_data_position.GetActualGenerations();
    _actual_data_raw_freqs = _actual_data_position.GetActualFrequencies();
    
    auto _time_for_seeding = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    _boost_gen.seed(_time_for_seeding);
    
    
    // Create prior distribution
    CMulator local_sim_object(_zparams);
    local_sim_object.SetWTAllele( _actual_data_position.GetWTIndex() );
    
    _expected_prior_sample_size = 0;
    
    switch (_factor_to_infer) {
        case Fitness: {
            auto tmp_prior = _zparams.GetString( fits_constants::PARAM_FITNESS_PRIOR_DISTRIB,
                                                fits_constants::PARAM_PRIOR_DISTRIB_FITNESS_DEFAULT );
            
            _expected_prior_sample_size = _zparams.GetSize_t( fits_constants::PARAM_NUM_ALLELES );
            
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
                //std::cerr << "Unkown prior distribution: " << tmp_prior << ". Setting to uniform as default." << std::endl;
                std::string tmp_str = "Error: unknown prior selected (" + tmp_prior + ")";
                throw tmp_str;
            }
            
            
            auto min_fitness_vec = local_sim_object.GetAlleleMinFitnessValues();
            auto max_fitness_vec = local_sim_object.GetAlleleMaxFitnessValues();
            
            if ( prior_distribution_vec.empty() ) {
                
                std::cout << "Sampling from prior... " << std::flush;
                
                PriorSampler<FLOAT_TYPE> sampler(min_fitness_vec, max_fitness_vec, _prior_type);
                _global_prior = sampler.SamplePrior(_repeats);
                
                std::cout << "Done." << std::endl;
            }
            else {
                _global_prior = prior_distribution_vec;
                
                if (_verbose_output) {
                    std::cout << "Reusing prior." << std::endl;
                }
            }
            break;
        }
        case PopulationSize: {
            std::vector<FLOAT_TYPE> minN {_zparams.GetDouble(fits_constants::PARAM_MIN_LOG_POPSIZE)};
            std::vector<FLOAT_TYPE> maxN {_zparams.GetDouble(fits_constants::PARAM_MAX_LOG_POPSIZE)};
            
            _expected_prior_sample_size = 1;
            
            if ( prior_distribution_vec.empty() ) {
                if (_verbose_output) {
                    std::cout << "Sampling from prior... " << std::flush;
                }
                
                PriorSampler<FLOAT_TYPE> sampler( minN, maxN, PriorDistributionType::UNIFORM );
                _global_prior = sampler.SamplePrior(_repeats);
                
                if (_verbose_output) {
                    std::cout << "Done." << std::endl;
                }
            }
            else {
                _global_prior = prior_distribution_vec;
                
                if (_verbose_output) {
                    std::cout << "Reusing prior." << std::endl;
                }
            }
            break;
        }
        case MutationRate: {
            auto min_matrix = local_sim_object.GetMinMutationRateMatrix();
            auto max_matrix = local_sim_object.GetMaxMutationRateMatrix();
            
            _expected_prior_sample_size = _zparams.GetSize_t( fits_constants::PARAM_NUM_ALLELES ) * _zparams.GetSize_t( fits_constants::PARAM_NUM_ALLELES );
            
            if ( prior_distribution_vec.empty() ) {
                
                if (_verbose_output) {
                    std::cout << "Sampling from prior... " << std::flush;
                }
                
                PriorSampler<FLOAT_TYPE> sampler( min_matrix, max_matrix, PriorDistributionType::UNIFORM );
                _global_prior = sampler.SamplePrior(_repeats);
                NormalizePrior();
                
                if (_verbose_output) {
                    std::cout << "Done." << std::endl;
                }
            }
            else {
                _global_prior = prior_distribution_vec;
                
                if (_verbose_output) {
                    std::cout << "Reusing prior." << std::endl;
                }
            }
            
            break;
        }
    }
}


void clsCMulatorABC::NormalizePrior()
{
    // auto tmp_prior = _global_prior;
    
    /*
     for ( auto vec : _global_prior ) {
     for ( auto val : vec ) {
     std::cout << val << ", ";
     }
     std::cout << std::endl;
     }
     */
    
    for (auto current_prior_sample_idx=0; current_prior_sample_idx<_global_prior.size(); ++current_prior_sample_idx ) {
        
        auto current_mutrate_vector = _global_prior[current_prior_sample_idx];
        
        //std::coxut << "normalizing prior: alleles " << _num_alleles << std::endl;
        MATRIX_TYPE tmp_mutrates_matrix( _num_alleles, _num_alleles );
        
        std::vector<FLOAT_TYPE> tmp_line_sum( tmp_mutrates_matrix.size1(), 0.0f );
        
        // if it were 0, then it would be overridden
        auto current_single_mutation_rate = current_mutrate_vector[1];
        
        
        // fill the matrix with values
        for ( auto i=0; i<current_mutrate_vector.size(); ++i ) {
            
            auto row = i / _num_alleles;
            auto col = i % _num_alleles;
            
            if ( _zparams.GetInt( fits_constants::PARAM_SINGLE_MUTATION_RATE, -1 ) > 0 ) {
                
                std::cout << "single mutation rate: " << std::pow( 10, current_single_mutation_rate ) << std::endl;
                tmp_mutrates_matrix(row,col) = std::pow( 10, current_single_mutation_rate );
            }
            else {
                tmp_mutrates_matrix(row,col) = std::pow( 10, current_mutrate_vector[i] );
            }
        }
        
        // normalize
        for ( auto current_row=0; current_row<tmp_mutrates_matrix.size1(); ++current_row ) {
            
            boost::numeric::ublas::matrix_row<MATRIX_TYPE> current_mutrate_row( tmp_mutrates_matrix, current_row );
            
            auto sum_mutation_rates = sum(current_mutrate_row) - current_mutrate_row[current_row];
            
            if ( sum_mutation_rates > 1.0f ) {
                std::string tmp_str = "Sum of mutation rates from allele " + std::to_string(current_row) + " is greater than 1 (" + std::to_string(sum_mutation_rates) + ")";
                throw tmp_str;
            }
            
            // if the samples mutation rate sum up to a value smaller than 1
            // then we need to make sure non-mutation event probability would
            // complement to 1.0
            if ( sum_mutation_rates < 1.0f ) {
                current_mutrate_row[current_row] = 1.0f - sum_mutation_rates;
            }
            
            // fill the vector back with the normalized values
            for ( auto i=0; i<current_mutrate_vector.size(); ++i ) {
                
                auto row = i / _num_alleles;
                auto col = i % _num_alleles;
                
                _global_prior[current_prior_sample_idx][i] = std::log10( tmp_mutrates_matrix(row, col) );
            }
        }
        
    } // mutation rate
    
    /*
     std::cout << "normalized: " << std::endl;
     for ( auto vec : _global_prior ) {
     for ( auto val : vec ) {
     std::cout << val << ", ";
     }
     std::cout << std::endl;
     }
     */
}


void clsCMulatorABC::RunABCInference( FactorToInfer factor, std::size_t number_of_batches )
{

    
    
    _simulation_result_vector.clear();
    // _float_prior_archive.clear();
    // _int_prior_archive.clear();
    
    
    // performace tracking
    boost::accumulators::accumulator_set<std::size_t, boost::accumulators::stats<boost::accumulators::tag::mean>> rate_stats;
    
    auto remaining_repeats = _repeats;
    std::size_t current_repeat = 0;
    auto repeats_in_batch = _repeats / number_of_batches;
    
    _simulation_result_vector.reserve(_repeats);
    
    _factor_to_infer = factor;
    
    _total_running_time_sec = 0;
    
    auto start_global = std::chrono::high_resolution_clock::now();
    
    std::string completion_eta_str = "";
    
    
    if (_verbose_output) {
        std::cout << "Remaining " << remaining_repeats
        << " repeats (calculating speed...)"
        << std::endl;
    }
    else {
        std::cout << "." << std::flush;
    }
    
    while (remaining_repeats > 0) {
        
        if (repeats_in_batch > remaining_repeats) {
            repeats_in_batch = remaining_repeats;
        }
        remaining_repeats -= repeats_in_batch;
        
        std::vector<SimulationResult> tmp_result_vector;
        
        
        //PRIOR_DISTRIB prior_subset;
        //prior_subset.resize(repeats_in_batch);
        
        if ( current_repeat > _global_prior.size() ) {
            std::string tmp_str = "Error: current repeat (" + std::to_string(current_repeat) +
            ") is larger than prior size (" + std::to_string(_global_prior.size()) + ")";
            throw tmp_str;
        }
        
        
        auto start = std::chrono::high_resolution_clock::now();
        
        auto start_idx = current_repeat;
        auto end_idx = current_repeat + repeats_in_batch;
        
        if ( end_idx > _global_prior.size() ) {
            end_idx = _global_prior.size();
        }
        
        
        //std::cout << std::endl << "\tPrior start_index =" << start_idx << ", end_index=" << end_idx << std::endl;
        
        
        try {
            switch (_factor_to_infer) {
                case Fitness: {
                    // tmp_result_vector = RunFitnessInferenceBatch(prior_subset, current_repeat, current_repeat+repeats_in_batch);
                    tmp_result_vector = RunFitnessInferenceBatch( _global_prior, start_idx, end_idx );
                    break;
                }
                case PopulationSize: {
                    //tmp_result_vector = RunPopulationSizeInferenceBatch(prior_subset);
                    tmp_result_vector = RunPopulationSizeInferenceBatch( _global_prior, start_idx, end_idx );
                    break;
                }
                case MutationRate: {
                    //tmp_result_vector = RunMutationInferenceBatch(prior_subset);
                    tmp_result_vector = RunMutationInferenceBatch( _global_prior, start_idx, end_idx );
                    break;
                }
            }
        }
        catch (std::string str) {
            std::string tmp_str = "Exception while running inference batch: " + str;
            throw tmp_str;
        }
        catch (...) {
            std::cerr << "Unknown exception while running inference batch. ";
        }
        auto end = std::chrono::high_resolution_clock::now();
        
        
        // move results from chunk to the main result vector
        // std::cout << "moving.... ";
        //for ( auto &&tmp_result : tmp_result_vector ) {
        //    _simulation_result_vector.push_back(std::move(tmp_result));
        //}
        for ( auto tmp_result : tmp_result_vector ) {
            _simulation_result_vector.push_back(tmp_result);
        }
        // std::cout << " Done." << std::endl;
        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        
        
        auto sim_count = repeats_in_batch;
        
        auto calculated_speed = static_cast<double>(sim_count) / static_cast<double>(elapsed_ms.count());
        calculated_speed *= 1000; // 1/milisecond to 1/second
        rate_stats(calculated_speed);
        _simulation_speed = calculated_speed;
        
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
        
        if (_verbose_output) {
            std::cout << "Remaining " << remaining_repeats
            << " repeats (rate of "
            << std::round(calculated_speed)
            << "sim/sec)"
            // << "sim/sec) - estimated time for completion of position: "
            // std::put_time(&completion_ETA_tm, "%c")
            // << completion_eta_str
            << std::endl;
        }
        else {
            std::cout << "." << std::flush;
        }
        
        current_repeat += repeats_in_batch;
        
        
    } // while (remaining_repeats > 0)
    
    
    std::cout << std::endl << "Calculating distance... " << std::flush;
    
    // scaling factor; default actualling means no scaling (divide by 1)
    std::vector<FLOAT_TYPE> scaling_vector(0);
    
    // std::cout << "scaling vector" << scaling_vector.size() << std::endl;
    
    auto scaling_option_str = _zparams.GetString( fits_constants::PARAM_SCALING,
                                                 fits_constants::PARAM_SCALING_DEFAULT );
    
    
    if ( scaling_option_str.compare(fits_constants::PARAM_SCALING_SD) == 0 ) {
        //std::string tmp_str = "Not scaling SD";
        //throw tmp_str;
        
        scaling_vector = GetSDPerAllele(0, _simulation_result_vector.size());
    }
    
    if ( scaling_option_str.compare(fits_constants::PARAM_SCALING_MAD) == 0) {
        //std::string tmp_str = "Not scaling MAD";
        //throw tmp_str;
        
        scaling_vector = GetMADPerAllele(0, _simulation_result_vector.size());
    }
    
    
    // get actual data matrix
    auto actual_generations = _actual_data_position.GetActualGenerations();
    auto actual_matrix = _actual_data_position.GetActualFreqsAsMatrix();
    
    
    // from each simulation result, get only the relevant generations
    for ( auto current_idx=0; current_idx<_simulation_result_vector.size(); ++current_idx ) {
        //for ( auto& current_result : _simulation_result_vector )
        _simulation_result_vector[current_idx].distance_from_actual =
        GetDistanceSimActual( actual_matrix, _simulation_result_vector[current_idx].sim_data_matrix, scaling_vector );
        
        _simulation_result_vector[current_idx].distance_metric = _zparams.GetString( fits_constants::PARAM_DISTANCE, fits_constants::PARAM_DISTANCE_L1 );
        // std::cout << " distance in current result = " << current_result.distance_from_actual << std::endl;
    }
    
    std::cout << "Done." << std::endl;
    
    std::cout << "Sorting results... ";
    std::sort(_simulation_result_vector.begin(), _simulation_result_vector.end());
    std::cout << "Done." << std::endl;
    
    auto end_global = std::chrono::high_resolution_clock::now();
    auto global_elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_global - start_global);
    _total_running_time_sec = static_cast<std::size_t>(global_elapsed_ms.count()) / 1000;
    
    std::cout << "Running time for position: " << _total_running_time_sec << " seconds" << std::endl;
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
                std::string tmp_str = "GetResultsVector: results vector is not sorted!";
                throw tmp_str;
            }
            
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


// Assumes input valididty has been checked
FLOAT_TYPE clsCMulatorABC::DistanceL1( const MATRIX_TYPE &actual_data, const MATRIX_TYPE &sim_data, const std::vector<FLOAT_TYPE> &scaling_vector )
{
    
    if ( _zparams.GetInt( "Debug", 0 ) > 0 ) {
        std::cout << "Calculating L1 distance:" << std::endl;
    }
    
    MATRIX_TYPE diff_matrix = actual_data - sim_data;
    
    diff_matrix = boost::numeric::ublas::abs(diff_matrix);
    
    // scale before summing up data
    if (scaling_vector.size() > 0 ) {
        for ( auto col=0; col<diff_matrix.size2(); ++col ) {
            boost::numeric::ublas::matrix_column<MATRIX_TYPE> current_col(diff_matrix, col);
            current_col = current_col / scaling_vector[col];
            
            if ( _zparams.GetInt( "Debug", 0 ) > 0 ) {
                std::cout << "scaling by " << scaling_vector[col] << std::endl;
            }
        }
    }
    
    // sum each allele
    FLOAT_TYPE sum_diff = 0.0f;
    for ( auto col=0; col<diff_matrix.size2(); ++col ) {
        
        /*
        if ( col == _wt_allele_idx ) {
             std::cout << " distance skipping allele " << col << std::endl;
            continue;
        }
         */
        
        boost::numeric::ublas::matrix_column<MATRIX_TYPE> current_col(diff_matrix, col);
        auto current_allele_sum = boost::numeric::ublas::sum(current_col);
        
        sum_diff += current_allele_sum;
    }
    
    if ( _zparams.GetInt( "Debug", 0 ) > 0 ) {
        if ( sum_diff == 0.0f) {
            std::cerr << "Warning: data file identical to simulated data (distance=0)." << std::endl;
        }
    }
    
    return sum_diff;
}


// Assumes input valididty has been checked
FLOAT_TYPE clsCMulatorABC::DistanceL2( const MATRIX_TYPE &actual_data, const MATRIX_TYPE &sim_data, const std::vector<FLOAT_TYPE> &scaling_vector )
{
    if ( _zparams.GetInt( "Debug", 0 ) > 0 ) {
        std::cout << "Calculating L2 distance:" << std::endl;
    }
    
    MATRIX_TYPE diff_matrix = actual_data - sim_data;
    
    // instead of ABS as in L1, square the values
    MATRIX_TYPE squared_diff_matrix = boost::numeric::ublas::element_prod(diff_matrix, diff_matrix);
    
    // scale before summing up data
    if (scaling_vector.size() > 0 ) {
        std::cerr << "Warning: L2 distance used, scaling for this measure has not been tested and thus not used." << std::endl;
    }
    
    // sum all squared differences
    FLOAT_TYPE sum_diff = 0.0f;
    for ( auto col=0; col<squared_diff_matrix.size2(); ++col ) {
        
        /*
        if ( col == _wt_allele_idx ) {
            continue;
        }
        */
        
        boost::numeric::ublas::matrix_column<MATRIX_TYPE> current_col(squared_diff_matrix, col);
        auto current_allele_sum = boost::numeric::ublas::sum(current_col);
        
        sum_diff += current_allele_sum;
    }
    
    return std::sqrt(sum_diff);
}


FLOAT_TYPE clsCMulatorABC::GetDistanceSimActual( const MATRIX_TYPE &actual_data, const MATRIX_TYPE &sim_data, const std::vector<FLOAT_TYPE> &scaling_vector )
{
    
    FLOAT_TYPE calculated_distance = -1.0f;
    
    if ( actual_data.size1() != sim_data.size1() ) {
        std::cerr << "DistanceSimActual - matrices don't match in size1: actual " <<
        actual_data.size1() << "vs sim " << sim_data.size1() << std::endl;
        
        std::string tmp_str = "DistanceSimActual - matrices don't match in size1";
        throw tmp_str;
    }
    
    if ( actual_data.size2() != sim_data.size2() ) {
        std::string tmp_str = "DistanceSimActual - matrices don't match in size2: actual " +
        std::to_string(actual_data.size2()) + "vs sim " + std::to_string(sim_data.size2());
        throw tmp_str;
    }
    
    
    auto tmp_distance_metric = _zparams.GetString( fits_constants::PARAM_DISTANCE, fits_constants::PARAM_DISTANCE_L1 );
    
    if ( tmp_distance_metric.compare(fits_constants::PARAM_DISTANCE_L1) == 0 ) {
        
        calculated_distance = DistanceL1( actual_data, sim_data, scaling_vector );
    }
    else if ( tmp_distance_metric.compare(fits_constants::PARAM_DISTANCE_L2) == 0 ) {
        
        calculated_distance = DistanceL2( actual_data, sim_data, scaling_vector );
    }
    
    return calculated_distance;
}

PRIOR_DISTRIB_VECTOR clsCMulatorABC::GetPriorFloat()
{
    //return _float_prior_archive;
    return _global_prior;
}

std::string clsCMulatorABC::GetPriorFloatAsString()
{
    std::string tmp_str = "";
    for ( auto current_vec : _global_prior ) {
        
        for ( auto current_val : current_vec ) {
            
            tmp_str += std::to_string(current_val) + fits_constants::FILE_FIELD_DELIMITER;
        }
        tmp_str = tmp_str + "\n";
    }
    return tmp_str;
}

void clsCMulatorABC::SetPriorFloat( PRIOR_DISTRIB_VECTOR given_prior )
{
    _global_prior.clear();
    
    for ( auto current_vec : given_prior ) {
        _global_prior.push_back(current_vec);
    }
    
    _use_stored_prior = true;
}


void clsCMulatorABC::VerifyIndece( const PRIOR_DISTRIB_VECTOR &prior_distrib, std::size_t start_idx, std::size_t end_idx )
{
    if ( start_idx > end_idx ) {
        std::string tmp_str = "clsCMulatorABC: start index of prior is bigger than the end (" + std::to_string(start_idx) + "," + std::to_string(end_idx) + ")";
        throw tmp_str;
    }
    if ( start_idx > prior_distrib.size() ) {
        std::string tmp_str = "clsCMulatorABC: start index of prior is bigger prior size (" + std::to_string(start_idx) + "," + std::to_string(prior_distrib.size()) + ")";
        throw tmp_str;
    }
    if ( end_idx > prior_distrib.size() ) {
        std::string tmp_str = "clsCMulatorABC: end index of prior is bigger than prior size (" + std::to_string(end_idx) + "," + std::to_string(prior_distrib.size()) + ")";
        throw tmp_str;
    }
}
