//
//  CMulator_Initialization.cpp
//  fits
//
//  Created by Tal Zinger on 15/11/2016.
//  Copyright © 2016 Stern Lab. All rights reserved.
//

#include "CMulator.h"
using namespace fits_constants;

/*****************************************************************
 Initialization of simulator variables
 
 Each function throws an exception if its parameters are missing.
 *****************************************************************/

void CMulator::InitMutationRates( ZParams zparams, bool available_mutrate_range )
{
    // first, check if we need to use a single value
    auto single_mutation_rate = zparams.GetFloat( fits_constants::PARAM_SINGLE_MUTATION_RATE, 0.0f );
    
    if ( single_mutation_rate != 0.0f ) {
        
        _use_single_mutation_rate = true;
        
        if ( _debug_mode ) {
            std::cout << "Initialized with a single mutation rate u=" << single_mutation_rate << std::endl;
        }
        
        for (auto from_allele=0; from_allele<_num_alleles; ++from_allele) {
            
            FLOAT_TYPE complement_mutation_rate = 1.0f - static_cast<FLOAT_TYPE>(_num_alleles-1) * single_mutation_rate  ;
            
            for ( auto to_allele=0; to_allele<_num_alleles; ++to_allele ) {
                
                if ( from_allele == to_allele ) {
                    _mutation_rates_matrix(from_allele, to_allele) = complement_mutation_rate;
                }
                else {
                    _mutation_rates_matrix(from_allele, to_allele) = single_mutation_rate;
                }
            }
        }
        
        return;
    }
    
    // individual mutation rate
    for (auto from_allele=0; from_allele<_num_alleles; ++from_allele) {
        
        FLOAT_TYPE sanity_mutation_rate_sum = 0.0;
        
        if ( _debug_mode ) {
            std::cout << "Initialized mutation rates:" << std::endl;
        }
        
        for ( auto to_allele=0; to_allele<_num_alleles; ++to_allele ) {
            
            std::string current_mutation_rate_str = PARAM_MUTATION_RATE + std::to_string(from_allele) + "_" + std::to_string(to_allele);
            
            auto tmp_mutation_rate = zparams.GetFloat( current_mutation_rate_str );
            
            //_mutation_rates[from_allele][to_allele] = tmp_mutation_rate;
            
            // side by side while moving to ublast
            _mutation_rates_matrix.at_element(from_allele, to_allele) = tmp_mutation_rate;
            
            if ( _debug_mode ) {
                std::cout << "from:" << from_allele << " to: " << to_allele << " u=" << tmp_mutation_rate << std::endl;
            }
            

            
            sanity_mutation_rate_sum += tmp_mutation_rate;
        }
        
        
        if ( std::fabs(1.0 - sanity_mutation_rate_sum ) > _epsilon_float_compare ) {
            std::cerr << "Warning: Mutation rates don't sum up to 1.0: " << sanity_mutation_rate_sum
            << ". Delta is " << 1.0 - sanity_mutation_rate_sum << " and epsilon is " << _epsilon_float_compare << std::endl;
            //throw( _name_of_run + ": Mutation rates don't sum up to 1.0: " + std::to_string(sanity_mutation_rate_sum) );
        }
    }
}


// 2016-12-22 - Note that now this is log-uniform - i.e. the data read is the power x (as in 10^x)
// 20170116 why should we mention mutation 0_0 in this format? it turns out quite tricky
void CMulator::InitMutationInferenceVariables( ZParams zparams )
{
    // first check if we infer a single mutation rate
    auto infer_single_mutrate = zparams.GetInt(fits_constants::PARAM_MIN_LOG_SINGLE_MUTATION_RATE, 0);
    
    if ( infer_single_mutrate != 0 ) {
        
        if ( _debug_mode ) {
            std::cout << "Inferring a single mutation rate." << std::endl;
        }
        _use_single_mutation_rate = true;
        
        for (auto from_allele=0; from_allele<_num_alleles; from_allele++) {
            
            for (auto to_allele=0; to_allele<_num_alleles; to_allele++) {
                
                auto tmp_min = zparams.GetInt( PARAM_MIN_LOG_SINGLE_MUTATION_RATE );
                auto tmp_max = zparams.GetInt( PARAM_MAX_LOG_SINGLE_MUTATION_RATE );
                
                if (_debug_mode) {
                    std::cout << "single mutrate min=" << tmp_min << "; max=" << tmp_max << std::endl;
                }
                
                _min_log_mutation_rate_matrix(from_allele, to_allele) = tmp_min;
                _max_log_mutation_rate_matrix(from_allele, to_allele) = tmp_max;
            }
        }
        
        return;
    }
    
    // individual mutation rate
    for (auto from_allele=0; from_allele<_num_alleles; from_allele++) {
        
        for (auto to_allele=0; to_allele<_num_alleles; to_allele++) {
            
            std::string min_log_mutation_rate_str = PARAM_MIN_LOG_MUTATION_RATE +
                std::to_string(from_allele) + "_" + std::to_string(to_allele);
            
            std::string max_log_mutation_rate_str = PARAM_MAX_LOG_MUTATION_RATE +
            std::to_string(from_allele) + "_" + std::to_string(to_allele);
            
            _min_log_mutation_rate_matrix(from_allele, to_allele) = zparams.GetInt( min_log_mutation_rate_str );
            _max_log_mutation_rate_matrix(from_allele, to_allele) = zparams.GetInt( max_log_mutation_rate_str );
        }
    }
}


void CMulator::InitBasicVariables( ZParams zparams )
{
    /* Mandatory Parameters - No Defaults */
    
    try {
        if ( zparams.GetInt( "Debug", 0 ) > 0 ) {
            _debug_mode = true;
        }
        
        _N = zparams.GetInt(PARAM_POPULATION_SIZE);
        
        _sample_size = zparams.GetInt( PARAM_SAMPLE_SIZE, PARAM_SAMPLE_SIZE_DEFAULT );
        
        if ( _sample_size > 0 ) {
            _use_observed_data = true;
        }
        
        if ( !IsValid_PopulationSize(_N)  ) {
            throw "Parameters: Missing or invalid N (must be positive).";
        }
        
        _alt_N = zparams.GetInt( PARAM_ALT_POPULATION_SIZE, 0 );
        _alt_generation = zparams.GetInt( PARAM_ALT_GENERATION, -1 );
        
        _repeats = zparams.GetInt( PARAM_SIM_REPEATS, 1 );
        
        _N0 = _N;
        
        // _num_alleles = Parameters::getInt( PARAM_NUM_ALLELES, PARAM_DEFAULT_VAL_INT );
        _num_alleles = zparams.GetInt(PARAM_NUM_ALLELES);
        
        if (  !IsValid_NumAlleles(_num_alleles) ) {
            throw "Parameters: Missing or invalid number of alleles (must be >=2)";
        }
        
        // not required for inference but we do need it for simulation
        _num_generations = zparams.GetInt(PARAM_NUM_GENERATIONS, 0);
        
        
        //_sample_size = Parameters::getInt( PARAM_SAMPLE_SIZE, PARAM_DEFAULT_VAL_INT );
        //_sample_size = zparams.GetInt( PARAM_SAMPLE_SIZE, PARAM_DEFAULT_VAL_INT );
        
        _bottleneck_interval = zparams.GetInt( PARAM_BOTTLENECK_INTERVAL, PARAM_DEFAULT_VAL_INT );
        
        _bottleneck_size = zparams.GetInt( PARAM_BOTTLENECK_SIZE, PARAM_DEFAULT_VAL_INT );
    }
    catch (const char* exp_txt) {
        std::cerr << "Exception while reading parameters: " << exp_txt << std::endl;
        throw "Can't continue initializtion.";
    }
    catch (std::exception& exp) {
        std::cerr << "Exception while reading parameters: " << exp.what() << std::endl;
        throw "Can't continue initializtion.";
    }
    catch (...) {
        std::cerr << "Unknown exception while reading essential parameters." << std::endl;
        throw "Can't continue initializtion.";
    }
    
    /* Not Mandatory, General Parameters */
    // I fear that it would cause lots of errors, we'll see
    auto tmp_default_epsilon = 2.0f * std::numeric_limits<FLOAT_TYPE>::epsilon();
    //_epsilon_float_compare = zparams.GetFloat( PARAM_EPSILON_FLOAT_COMPARE, tmp_default_epsilon );
    _epsilon_float_compare = zparams.GetDouble( PARAM_EPSILON_FLOAT_COMPARE, tmp_default_epsilon );
    
    _name_of_run = zparams.GetString( PARAM_SIM_ID, std::string("SimX") );
    
    _logistic_growth = ( zparams.GetInt( PARAM_LOGISTIC_GROWTH, PARAM_DEFAULT_VAL_INT ) > 0 );
    _logistic_growth_K = zparams.GetInt( PARAM_LOGISTIC_GROWTH_K, PARAM_DEFAULT_VAL_INT );
    _logistic_growth_r = zparams.GetFloat( PARAM_LOGISTIC_GROWTH_r, PARAM_DEFAULT_VAL_FLOAT );
    _logistic_growth_t = 0;
    
    //_skip_stochastic_step = ( zparams.GetInt(PARAM_SKIP_STOCHASTIC_STEP, PARAM_SKIP_STOCHASTIC_STEP_DEFAULT) > 0);
    
    // seed already initialized, replace only if manual seed is given
    auto original_seed = _time_for_seeding;
    _time_for_seeding = zparams.GetUnsignedInt(PARAM_MANUAL_SEED, original_seed);
    
    
}



void CMulator::InitGeneralInferenceVariables( ZParams zparams )
{
    //_generation_shift = zparams.GetInt( PARAM_GENERATION_SHIFT, PARAM_DEFAULT_VAL_INT );
    
    _generation_shift = 0;
    //_fitness_increment = zparams.GetFloat( PARAM_ALLELE_FITNESS_INCRENEMT, PARAM_ALLELE_FITNESS_INCRENEMT_DEFAULT);
    
    _acceptance_rate = zparams.GetFloat( PARAM_ACCEPTANCE_RATE, 0.0 );
    //_top_percent_to_keep = zparams.GetFloat( PARAM_ACCEPTANCE_RATE, 0.0 ) * 100.0f;
    
    // 2017-01-15 can't allow this to happen that we don't have the percent.
    /* TODO: return to this thought when the CMulator object knows what is its purpose - simulate or inference.
     till then - this is passed to inference section */
    /*try {
        
    }
    catch (...) {
        std::cerr << "Error: cannot get top percent to keep in parameter file." << std::endl;
        throw "Error: cannot get top percent to keep in parameter file.";
    }
    */
    
    
    
    // TODO: make this obsolete, convert to actual file to take this data from
    // auto tmp_no_init_freqs_as_parameters = Parameters::getInt( PARAM_ALLELE_NO_INIT_FREQS, 0 );
    //auto tmp_no_init_freqs_as_parameters = zparams.GetInt( PARAM_ALLELE_NO_INIT_FREQS, 1 );
    //_no_init_freqs_as_parameters = ( tmp_no_init_freqs_as_parameters == 1 );
}


void CMulator::InitInitialAlleleFreqs( ZParams zparams )
{
    FLOAT_TYPE sanity_allele_init_freq_sum = 0.0;
    
    //FLOAT_TYPE highest_freq = -1.0;
    
    for (auto current_allele_num=0; current_allele_num<_num_alleles; current_allele_num++) {
        
        std::string current_allele_freq_str = PARAM_ALLELE_INIT_FREQ + std::to_string(current_allele_num);
        
        _allele_init_freqs[current_allele_num] = zparams.GetFloat(current_allele_freq_str, -1.0f);
        
        if ( _allele_init_freqs[current_allele_num] < 0.0f ) {
            
            throw "Missing initial frequency for allele " + std::to_string(current_allele_num);
        }
        
        sanity_allele_init_freq_sum += _allele_init_freqs[current_allele_num];
        
        //_sim_data[0][current_allele_num] = _allele_init_freqs[current_allele_num];
        _all_simulated_data(0, current_allele_num) = _allele_init_freqs[current_allele_num];
        _observed_simulated_data(0, current_allele_num) = _allele_init_freqs[current_allele_num];
        
    } // end fitness and init freq
    
    /*
    if ( std::fabs( 1.0 - sanity_allele_init_freq_sum ) > _epsilon_float_compare ) {
        std::cerr << _name_of_run << ": Allele initial frequencies don't sum up to 1.0: " << sanity_allele_init_freq_sum << std::endl;
        throw( _name_of_run + ": Allele initial frequencies don't sum up to 1.0: " + std::to_string(sanity_allele_init_freq_sum) );
    }
     */
}


void CMulator::InitFitnessValues( ZParams zparams )
{
    
    for (auto current_allele_num=0; current_allele_num<_num_alleles; current_allele_num++) {
        
        std::string current_allele_fitness_str = PARAM_ALLELE_FITNESS + std::to_string(current_allele_num);
        std::string current_allele_min_fitness = PARAM_ALLELE_MIN_FITNESS + std::to_string(current_allele_num);
        std::string current_allele_max_fitness = PARAM_ALLELE_MAX_FITNESS + std::to_string(current_allele_num);
        
        _allele_fitness[current_allele_num] = zparams.GetFloat(current_allele_fitness_str, -1.0);
        
        if ( _allele_fitness[current_allele_num] < 0 ) {
            //std::cerr << "Missing fitness value for allele " << current_allele_num << std::endl;
            throw "Missing fitness value for allele " + std::to_string(current_allele_num);
        }
    }
    
}


void CMulator::InitFitnessInference( ZParams zparams )
{

    for (auto current_allele_num=0; current_allele_num<_num_alleles; current_allele_num++) {
        
        std::string current_allele_min_fitness = PARAM_ALLELE_MIN_FITNESS + std::to_string(current_allele_num);
        std::string current_allele_max_fitness = PARAM_ALLELE_MAX_FITNESS + std::to_string(current_allele_num);
        
        _allele_min_fitness[current_allele_num] = zparams.GetFloat(current_allele_min_fitness, -1.0);
        if ( _allele_min_fitness[current_allele_num] < 0.0 ) {
            //std::cerr << "Missing minimum fitness value for allele " << current_allele_num << std::endl;
            throw "Missing minimum fitness value for allele " + std::to_string(current_allele_num);
        }
        
        
        _allele_max_fitness[current_allele_num] = zparams.GetFloat(current_allele_max_fitness, -1.0);
        if ( _allele_max_fitness[current_allele_num] < 0.0 ) {
            //std::cerr << "Missing maximum fitness value for allele " << current_allele_num << std::endl;
            throw "Missing maximum fitness value for allele " + std::to_string(current_allele_num);
        }
        
    }
}


void CMulator::InitPopulationSizeInference( ZParams zparams )
{
    _Nlog_min = zparams.GetInt("_Nlog_min");
    _Nlog_max = zparams.GetInt("_Nlog_max");

}

void CMulator::InitGenerationInference( ZParams zparams )
{
    // either
    try {
        zparams.GetInt( fits_constants::PARAM_MIN_FIXED_INT_GENERATION);
        zparams.GetInt( fits_constants::PARAM_MAX_FIXED_INT_GENERATION);
    }
    catch (...) {
        auto tmp_min_str = fits_constants::PARAM_MIN_INT_GENERATION + "0";
        auto tmp_max_str = fits_constants::PARAM_MAX_INT_GENERATION + "0";
        
        zparams.GetInt(tmp_min_str);
        zparams.GetInt(tmp_max_str);
    }
}

void CMulator::WarnAgainstDeprecatedParameters(ZParams zparams)
{
    std::string msg =  "Warning - deprecated parameter detected - ";
    
    if ( !zparams.GetString( "_infer_mutation_rate", "" ).empty()  ) {
        std::cout << msg << "_infer_mutation_rate" << std::endl;
    }
    
    if ( !zparams.GetString( "_generations_shift", "" ).empty() ) {
        std::cout << msg << "_generations_shift" << std::endl;
    }
    
    if ( !zparams.GetString( "_skip_stochastic_step", "" ).empty() ) {
        std::cout << msg << "_skip_stochastic_step" << std::endl;
    }
    
    if ( !zparams.GetString( "_num_generations", "" ).empty() ) {
        std::cout << "NOTE: _num_generations parameter not used for inference" << std::endl;
    }
}



void CMulator::InitMemberVariables( ZParams zparams )
{
    _available_mutrate = true;
    _available_fitness = true;
    _available_popsize = true;
    _available_essential = true;
    _available_inference = true;
    _available_fitness_range = true;
    _available_mutrate_range = true;
    _available_popsize_range = true;
    _available_init_freqs = true;
    _available_generation_range = true;
    
    // exception here will crash the simulator, this is sort of fine with me
    
    InitBasicVariables(zparams);
    
    // I want to make sure these are initialized
    // only possible after basic data was read
    //_mutation_rates.resize(boost::extents[_num_alleles][_num_alleles]);
    _mutation_rates_matrix.resize(_num_alleles, _num_alleles );
    _wt_allele_index = -1;
    
    // generations + 1 because 0 is given and is not the first processed generation.
    //_sim_data.resize(boost::extents[_num_generations+1][_num_alleles]);
    
    _all_simulated_data.resize( _num_generations+1, _num_alleles );
    _observed_simulated_data.resize( _num_generations+1, _num_alleles );
    
    _allele_init_freqs.resize(_num_alleles);
    _allele_fitness.resize(_num_alleles);
    //_allele_fitness_bar.resize(_num_alleles);
    //_allele_fitness_adjusted.resize(_num_alleles);
    
    _min_log_mutation_rate_matrix.resize(_num_alleles, _num_alleles);
    _max_log_mutation_rate_matrix.resize(_num_alleles, _num_alleles);
    
    _allele_fitness.resize(_num_alleles);
    _allele_max_fitness.resize(_num_alleles);
    _allele_min_fitness.resize(_num_alleles);
    
    
    try {
        InitInitialAlleleFreqs(zparams);
    }
    catch (const char* txt) {
        std::cerr << txt << std::endl;
    }
    catch (...) {
        _available_init_freqs = false;
    }
    
    try {
        InitGeneralInferenceVariables(zparams);
    }
    catch (...) {
        _available_essential = false;
        std::cerr << "Error: essential inference parameter missing." << std::endl;
        throw "Error: essential simulation parameter missing.";
    }
    
    
    try {
        InitMutationInferenceVariables(zparams);
    }
    catch (...) {
        _available_mutrate_range = false;
    }
    
    try {
        InitMutationRates(zparams, _available_mutrate_range);
    }
    catch (...) {
        _available_mutrate = false;
    }
    
    
    try {
        InitPopulationSizeInference(zparams);
    }
    catch (...) {
        _available_popsize_range = false;
    }
    
    try {
        InitGenerationInference(zparams);
    }
    catch (...) {
        _available_generation_range = false;
    }
    
    
    try {
        InitFitnessValues(zparams);
    }
    catch (...) {
        _available_fitness = false;
    }
    
    try {
        InitFitnessInference(zparams);
    }
    catch (...) {
        _available_fitness_range = false;
    }
    
    
    /* Internal State */
    _current_generation = 0;
    
    _initialized_with_parameters = true;
}


