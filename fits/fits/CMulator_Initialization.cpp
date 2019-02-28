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

#include "CMulator.h"
//using namespace fits_constants;


//Destroy and rewrite everything


/*****************************************************************
 Initialization of simulator variables
 
 Each function throws an exception if its parameters are missing.
 *****************************************************************/

bool CMulator::InitMutationRates( ZParams zparams )
{
    // first, check if we need to use a single value
    auto single_mutation_rate = zparams.GetDouble( fits_constants::PARAM_SINGLE_MUTATION_RATE, 0.0f );
    
    if ( single_mutation_rate > 0.0f ) {
        
        _use_single_mutation_rate = true;
        
        if ( _debug_mode ) {
            std::cout << "Initialized with a single mutation rate u=" << single_mutation_rate << std::endl;
        }
        
        //FLOAT_TYPE sanity_mutation_rate_sum = 0.0f;
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
        
        for ( auto current_row = 0; current_row < _mutation_rates_matrix.size2(); ++current_row ) {
            boost::numeric::ublas::matrix_row<MATRIX_TYPE> current_row_data( _mutation_rates_matrix, current_row );
            
            if ( sum(current_row_data) != 1 ) {
                std::string tmp_str = "Error: single mutation rate doesn't sum up to 1 for allele (row) " + std::to_string(current_row);
                throw tmp_str;
            }
            
        }
        
        return true;
    }
    
    // individual mutation rate
    _use_single_mutation_rate = false;
    for (auto from_allele=0; from_allele<_num_alleles; ++from_allele) {
        
        //FLOAT_TYPE sanity_mutation_rate_sum = 0.0f;
        
        if ( _debug_mode ) {
            std::cout << "Initialized mutation rates:" << std::endl;
        }
        

        for ( auto to_allele=0; to_allele<_num_alleles; ++to_allele ) {
            
            std::string current_mutation_rate_str = fits_constants::PARAM_MUTATION_RATE + std::to_string(from_allele) + "_" + std::to_string(to_allele);
            
            // we would later complete this to 1, don't actually need to read this value
            if ( from_allele == to_allele ) {
                _mutation_rates_matrix.at_element(from_allele, to_allele) = 0.0f;
                continue;
            }
            
            if ( !zparams.IsParameterFound( current_mutation_rate_str ) ) {
                // if we reached here, single is not found AND specific mutation rate is not found
                return false;
            }
            
            //std::cout << "mutatio rate str: " <<current_mutation_rate_str << std::endl;
            auto tmp_mutation_rate = zparams.GetDouble( current_mutation_rate_str );
            
            _mutation_rates_matrix.at_element(from_allele, to_allele) = tmp_mutation_rate;
            
            //sanity_mutation_rate_sum += tmp_mutation_rate;
        }
        
    }// finished loading individual mutation rates
    
    // instead of making sure this == 1 calculate the diagonal such that it does
    for ( auto current_row = 0; current_row < _mutation_rates_matrix.size2(); ++current_row ) {
        boost::numeric::ublas::matrix_row<MATRIX_TYPE> current_row_data( _mutation_rates_matrix, current_row );
        
        auto row_sum = sum(current_row_data) - current_row_data[current_row];
        auto diagonal_value = 1.0f - row_sum;
        
        current_row_data[current_row] = diagonal_value;
    }
    
    return true;
}


//void CMulator::InitMutationInferenceVariables( ZParams zparams )
bool CMulator::InitMutatioRateRange( ZParams zparams )
{
    // first, check if we need to use a single value
    if ( zparams.IsParameterFound(fits_constants::PARAM_MIN_LOG_SINGLE_MUTATION_RATE) ) {
        _use_single_mutation_rate = true;
        
        for (auto from_allele=0; from_allele<_num_alleles; from_allele++) {
            
            for (auto to_allele=0; to_allele<_num_alleles; to_allele++) {
                
                try {
                    auto tmp_min = zparams.GetDouble( fits_constants::PARAM_MIN_LOG_SINGLE_MUTATION_RATE );
                    auto tmp_max = zparams.GetDouble( fits_constants::PARAM_MAX_LOG_SINGLE_MUTATION_RATE );
                    
                    _min_log_mutation_rate_matrix(from_allele, to_allele) = tmp_min;
                    _max_log_mutation_rate_matrix(from_allele, to_allele) = tmp_max;
                }
                catch (...) {
                    return false;
                }
            }
        }
        
        return true;
    }
    
    
    // individual mutation rate
    for (auto from_allele=0; from_allele<_num_alleles; from_allele++) {
        
        for (auto to_allele=0; to_allele<_num_alleles; to_allele++) {
            
            std::string min_log_mutation_rate_str = fits_constants::PARAM_MIN_LOG_MUTATION_RATE +
                std::to_string(from_allele) + "_" + std::to_string(to_allele);
            
            std::string max_log_mutation_rate_str = fits_constants::PARAM_MAX_LOG_MUTATION_RATE +
            std::to_string(from_allele) + "_" + std::to_string(to_allele);
            
            try {
                _min_log_mutation_rate_matrix(from_allele, to_allele) = zparams.GetDouble( min_log_mutation_rate_str );
                _max_log_mutation_rate_matrix(from_allele, to_allele) = zparams.GetDouble( max_log_mutation_rate_str );
                
                if ( _min_log_mutation_rate_matrix(from_allele, to_allele) >
                    _max_log_mutation_rate_matrix(from_allele, to_allele) ) {
                    std::string tmp_str = "Error: Minimum mutation rate bigger than maximum between alleles " + std::to_string(from_allele) +
                    " and " + std::to_string(to_allele);
                    throw tmp_str;
                }
                
                if ( _min_log_mutation_rate_matrix(from_allele, to_allele) > 0 || _max_log_mutation_rate_matrix(from_allele, to_allele) > 0 ) {
                    std::string tmp_str = "Error: Mutation rate limits between alleles " + std::to_string(from_allele) +
                    " and " + std::to_string(to_allele) + " are bigger than 1.";
                    throw tmp_str;
                }
            }
            catch (...) {
                return false;
            }
        }
    }
    
    return true;
}

/*
void CMulator::InitBasicVariables( ZParams zparams )
{
    
    try {
        if ( zparams.GetInt( "Debug", 0 ) > 0 ) {
            _debug_mode = true;
        }
        
        _N = zparams.GetInt(fits_constants::PARAM_POPULATION_SIZE, -1);
        
        _sample_size = zparams.GetInt( fits_constants::PARAM_SAMPLE_SIZE, fits_constants::PARAM_SAMPLE_SIZE_DEFAULT );
        if ( _sample_size > 0 ) {
            _use_observed_data = true;
        }
        
        //_alt_N = zparams.GetInt( PARAM_ALT_POPULATION_SIZE, 0 );
        //_alt_generation = zparams.GetInt( PARAM_ALT_GENERATION, -1 );
        
        _repeats = zparams.GetInt( fits_constants::PARAM_SIM_REPEATS, fits_constants::PARAM_SIM_REPEATS_DEFAULT );
        if ( _repeats <= 0 ) {
            throw "Error: Number of repoeats must be positive.";
        }
        
        _N0 = _N;
        
        _num_alleles = zparams.GetInt(fits_constants::PARAM_NUM_ALLELES);
        
        // not required for inference but we do need it for simulation
        _num_generations = zparams.GetInt(fits_constants::PARAM_NUM_GENERATIONS, 0);
        
        _bottleneck_interval = zparams.GetInt( fits_constants::PARAM_BOTTLENECK_INTERVAL, fits_constants::PARAM_DEFAULT_VAL_INT );
        
        _bottleneck_size = zparams.GetInt( fits_constants::PARAM_BOTTLENECK_SIZE, fits_constants::PARAM_DEFAULT_VAL_INT );
    }
    catch (const char* exp_txt) {
        std::cerr << "Exception while reading parameters: " << exp_txt << std::endl;
        throw "Can't continue initializtion.";
    }
    catch (const std::string exp_txt) {
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
    
    /// Not Mandatory, General Parameters
    // I fear that it would cause lots of errors, we'll see
    auto tmp_default_epsilon = 2.0f * std::numeric_limits<float>::epsilon();
    //_epsilon_float_compare = zparams.GetFloat( PARAM_EPSILON_FLOAT_COMPARE, tmp_default_epsilon );
    _epsilon_float_compare = zparams.GetDouble( fits_constants::PARAM_EPSILON_FLOAT_COMPARE, tmp_default_epsilon );
    
    
    _logistic_growth = ( zparams.GetInt( fits_constants::PARAM_LOGISTIC_GROWTH, fits_constants::PARAM_DEFAULT_VAL_INT ) > 0 );
    _logistic_growth_K = zparams.GetInt( fits_constants::PARAM_LOGISTIC_GROWTH_K, fits_constants::PARAM_DEFAULT_VAL_INT );
    _logistic_growth_r = zparams.GetFloat( fits_constants::PARAM_LOGISTIC_GROWTH_r, fits_constants::PARAM_DEFAULT_VAL_FLOAT );
    _logistic_growth_t = 0;
    
    // seed already initialized, replace only if manual seed is given
    auto original_seed = _time_for_seeding;
    _time_for_seeding = zparams.GetUnsignedInt(fits_constants::PARAM_MANUAL_SEED, original_seed);
    
    
}
*/

/*
void CMulator::InitGeneralInferenceVariables( ZParams zparams )
{
    //_generation_shift = zparams.GetInt( PARAM_GENERATION_SHIFT, PARAM_DEFAULT_VAL_INT );
    
    _generation_shift = 0;
    //_fitness_increment = zparams.GetFloat( PARAM_ALLELE_FITNESS_INCRENEMT, PARAM_ALLELE_FITNESS_INCRENEMT_DEFAULT);
    
    _acceptance_rate = 0.0f;
    try {
        _acceptance_rate = zparams.GetFloat( fits_constants::PARAM_ACCEPTANCE_RATE );
    }
    catch (...) {
        throw "Cannot find a valid acceptance (>0) rate in parameters file";
    }
    
    if ( _acceptance_rate <= 0.0f ) {
        throw "Cannot find a valid acceptance rate (>0) in parameters file";
    }
    
    //_top_percent_to_keep = zparams.GetFloat( PARAM_ACCEPTANCE_RATE, 0.0 ) * 100.0f;
    
    // 2017-01-15 can't allow this to happen that we don't have the percent.
 
    
    
    
    // TODO: make this obsolete, convert to actual file to take this data from
    // auto tmp_no_init_freqs_as_parameters = Parameters::getInt( PARAM_ALLELE_NO_INIT_FREQS, 0 );
    //auto tmp_no_init_freqs_as_parameters = zparams.GetInt( PARAM_ALLELE_NO_INIT_FREQS, 1 );
    //_no_init_freqs_as_parameters = ( tmp_no_init_freqs_as_parameters == 1 );
}
*/

bool CMulator::InitInitialAlleleFreqs( ZParams zparams )
{
    FLOAT_TYPE sanity_allele_init_freq_sum = 0.0;
    
    //FLOAT_TYPE highest_freq = -1.0;
    
    
    for (auto current_allele_num=0; current_allele_num<_num_alleles; ++current_allele_num) {
        
        
        std::string current_allele_freq_str = fits_constants::PARAM_ALLELE_INIT_FREQ + std::to_string(current_allele_num);
        
        //std::cout << "parameter " << current_allele_freq_str << std::endl;
        if ( !zparams.IsParameterFound(current_allele_freq_str) ) {
            //std::cout << "not found " << current_allele_freq_str << std::endl;
            return false;
        }

        
        _allele_init_freqs[current_allele_num] = zparams.GetDouble(current_allele_freq_str);
        
        if ( _allele_init_freqs[current_allele_num] < 0.0f || _allele_init_freqs[current_allele_num] > 1.0f ) {
            std::string tmp_str = "Error: invalid frequency value (" + std::to_string(_allele_init_freqs[current_allele_num]) + ") for allele " + std::to_string(current_allele_num);
            throw tmp_str;
        }
        
        sanity_allele_init_freq_sum += _allele_init_freqs[current_allele_num];
        
        //_sim_data[0][current_allele_num] = _allele_init_freqs[current_allele_num];
        _all_simulated_data(0, current_allele_num) = _allele_init_freqs[current_allele_num];
        _observed_simulated_data(0, current_allele_num) = _allele_init_freqs[current_allele_num];
        
    } // end fitness and init freq
    
    if ( sanity_allele_init_freq_sum != 1.0 ) {
        std::string tmp_str = "Error: allele initial frequencies do not sum up to 1 (" + std::to_string(sanity_allele_init_freq_sum) + ")";
        throw tmp_str;
    }
    return true;
    /*
    if ( std::fabs( 1.0 - sanity_allele_init_freq_sum ) > _epsilon_float_compare ) {
        std::cerr << _name_of_run << ": Allele initial frequencies don't sum up to 1.0: " << sanity_allele_init_freq_sum << std::endl;
        throw( _name_of_run + ": Allele initial frequencies don't sum up to 1.0: " + std::to_string(sanity_allele_init_freq_sum) );
    }
     */
}


bool CMulator::InitFitnessValues( ZParams zparams )
{
    
    for (auto current_allele_num=0; current_allele_num<_num_alleles; current_allele_num++) {
        
        std::string current_allele_fitness_str = fits_constants::PARAM_ALLELE_FITNESS + std::to_string(current_allele_num);
        //std::string current_allele_min_fitness = PARAM_ALLELE_MIN_FITNESS + std::to_string(current_allele_num);
        //std::string current_allele_max_fitness = PARAM_ALLELE_MAX_FITNESS + std::to_string(current_allele_num);
        
        if ( !zparams.IsParameterFound( current_allele_fitness_str ) ) {
            return false;
        }
        
        try {
            _allele_fitness[current_allele_num] = zparams.GetFloat(current_allele_fitness_str);
        }
        catch (...) {
            return false;
        }
        
        
        if ( _allele_fitness[current_allele_num] < 0 ) {
            //std::cerr << "Missing fitness value for allele " << current_allele_num << std::endl;
            std::string tmp_str = "Error: Negative value " + std::to_string(_allele_fitness[current_allele_num]) + " set for allele " + std::to_string(current_allele_num);
            throw tmp_str;
        }

    }
            
    return true;
}


bool CMulator::InitFitnessRange( ZParams zparams )
{
    for (auto current_allele_num=0; current_allele_num<_num_alleles; current_allele_num++) {
        
        std::string current_allele_min_fitness = fits_constants::PARAM_ALLELE_MIN_FITNESS + std::to_string(current_allele_num);
        std::string current_allele_max_fitness = fits_constants::PARAM_ALLELE_MAX_FITNESS + std::to_string(current_allele_num);
        
        // if missing, use default. if wrongly specify throw an error
        if ( !zparams.IsParameterFound( current_allele_min_fitness ) ) {
            // _allele_min_fitness[current_allele_num] = fits_constants::ALLELE_FITNESS_DEFAULT_MIN;
            return false;
        }
        else {
            _allele_min_fitness[current_allele_num] = zparams.GetDouble( current_allele_min_fitness );
        }
        
        if ( !zparams.IsParameterFound( current_allele_max_fitness ) ) {
            // _allele_max_fitness[current_allele_num] = fits_constants::ALLELE_FITNESS_DEFAULT_MAX;
            return false;
        }
        else {
            _allele_max_fitness[current_allele_num] = zparams.GetDouble( current_allele_max_fitness );
        }
    
        
        // now verify the values are legit
        if ( _allele_max_fitness[current_allele_num] < _allele_min_fitness[current_allele_num] ) {
            std::string tmp_str = "Error: maximum fitness for allele " + std::to_string(current_allele_num) +
            " is smaller than its minimum (" + std::to_string( _allele_max_fitness[current_allele_num] ) +
            "<" + std::to_string( _allele_min_fitness[current_allele_num] ) + ")";
            throw tmp_str;
        }
        
        if ( _allele_min_fitness[current_allele_num] < 0 || _allele_max_fitness[current_allele_num] < 0 ) {
            std::string tmp_str = "Error: range for allele " + std::to_string(current_allele_num) +
            " contains negative values";
            throw tmp_str;
        }
    }
    
    return true;
}


bool CMulator::InitPopulationSize( ZParams zparams )
{
    if (!zparams.IsParameterFound( fits_constants::PARAM_POPULATION_SIZE ) ) {
        _N = -1;
        _N0 = _N;
        return false;
    }
    
    _N = zparams.GetInt( fits_constants::PARAM_POPULATION_SIZE );
    _N0 = _N;
    
    return true;
}


bool CMulator::InitPopulationSizeRange( ZParams zparams )
{
    try {
        _Nlog_min = zparams.GetDouble( fits_constants::PARAM_MIN_LOG_POPSIZE );
        _Nlog_max = zparams.GetDouble( fits_constants::PARAM_MAX_LOG_POPSIZE );
    }
    catch (...) {
        return false;
    }
    
    return true;
}




bool CMulator::InitSampleSize( ZParams zparams )
{
    try {
        _sample_size = zparams.GetInt( fits_constants::PARAM_SAMPLE_SIZE );
    }
    catch (...) {
        return false;
    }
    
    if ( _sample_size < 0 ) {
        std::string tmp_str = "Error: sample size must be positive (" + std::to_string(_sample_size) + ")";
        throw tmp_str;
    }
    
    return true;
}


bool CMulator::InitRepeats( ZParams zparams )
{
    try {
        _repeats = zparams.GetInt(fits_constants::PARAM_SIM_REPEATS);
    }
    catch(...) {
        _repeats = fits_constants::PARAM_SIM_REPEATS_DEFAULT;
    }
    
    if ( _repeats <= 0 ) {
        std::string tmp_str = "Error: Number of repoeats must be positive (" + std::to_string(_sample_size) + ")";
        throw tmp_str;
    }
    
    return true;
}


void CMulator::InitMemberVariables( ZParams zparams )
{
    
    // probably the most critical parameter to load - sized of vectors and matrices depend on it
    try {
        _num_alleles = zparams.GetInt(fits_constants::PARAM_NUM_ALLELES);
        _available_num_alleles = true;
    }
    catch (...) {
        // we have to have the number of alleles for initialization of
        // some data objects
        std::string tmp_str = "Error: number of alleles must be specified in the parameters file.";
        throw tmp_str;
    }
    
    
    try {
        _num_generations = zparams.GetInt(fits_constants::PARAM_NUM_GENERATIONS);
        _available_num_generations = true;
    }
    catch (...) {
        _available_num_generations = false;
    }
    
    
    try {
            _generation_shift = zparams.GetInt(fits_constants::PARAM_GENERATION_SHIFT);
    }
    catch (...) {
        // TODO: say that we don't ahave shift.. not that it matters
    }
    
    // the following require number of alleles and generations
    // so once they are defined we can (and must) initialize them
    try {
        _mutation_rates_matrix.resize( _num_alleles, _num_alleles );
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
    }
    catch (...) {
        // the most probably reason for failure to allocate space is _num_alleles < 0,
        // so adding this information should help fix the problem
        std::string tmp_str = "Error while setting size of data vectors/matrices. Number of alleles = " + std::to_string(_num_alleles);
        throw tmp_str;
    }
    
    
    // using default value so always true
    InitRepeats(zparams);
    
    
    // if a certain parameter is missing, its range (used for inference) must be available
    _available_mutrate = InitMutationRates(zparams);
    
    _available_fitness = InitFitnessValues(zparams);
    
    _available_popsize = InitPopulationSize(zparams);
    
    
    _available_fitness_range = InitFitnessRange(zparams);
    _available_mutrate_range = InitMutatioRateRange(zparams);
    _available_popsize_range = InitPopulationSizeRange(zparams);
    _available_init_freqs = InitInitialAlleleFreqs(zparams);
    
    _use_observed_data = InitSampleSize(zparams);
    

    
    _available_bottleneck = false;
    _bottleneck_size = zparams.GetInt( fits_constants::PARAM_BOTTLENECK_SIZE, fits_constants::PARAM_BOTTLENECK_NOT_DEFINED );
    _bottleneck_interval = zparams.GetInt( fits_constants::PARAM_BOTTLENECK_INTERVAL, fits_constants::PARAM_BOTTLENECK_NOT_DEFINED );
    if ( _bottleneck_size > 0 && _bottleneck_interval < 0) {
        std::string tmp_str = "Error: if a bottleneck is defined, its interval must also be defined.";
        throw tmp_str;
    }
    if ( _bottleneck_size > 0 && _bottleneck_interval > 0) {
        _available_bottleneck = true;
    }
    
    
    _logistic_growth = ( zparams.GetInt( fits_constants::PARAM_LOGISTIC_GROWTH, fits_constants::PARAM_DEFAULT_VAL_INT ) > 0 );
    if ( _logistic_growth ) {
        try {
            _logistic_growth_K = zparams.GetInt( fits_constants::PARAM_LOGISTIC_GROWTH_K );
            _logistic_growth_r = zparams.GetDouble( fits_constants::PARAM_LOGISTIC_GROWTH_r );
            _logistic_growth_t = 0;
        }
        catch (...) {
            std::string tmp_str = "Error: logistic growth is set, but not all parameters (K, r) are found.";
            throw tmp_str;
        }
    }
    

    
    // Not Mandatory, General Parameters
    // I fear that it would cause lots of errors, we'll see
    // auto tmp_default_epsilon = 2.0f * std::numeric_limits<float>::epsilon();
    //_epsilon_float_compare = zparams.GetFloat( PARAM_EPSILON_FLOAT_COMPARE, tmp_default_epsilon );
    // _epsilon_float_compare = zparams.GetDouble( fits_constants::PARAM_EPSILON_FLOAT_COMPARE, tmp_default_epsilon );
    
    
    // seed already initialized, replace only if manual seed is given
    auto original_seed = _time_for_seeding;
    _time_for_seeding = zparams.GetUnsignedInt(fits_constants::PARAM_MANUAL_SEED, original_seed);
    
    
    // I want to make sure these are initialized
    // only possible after basic data was read
    _wt_allele_index = -1;
    
    
    
    
    try {
        InitInitialAlleleFreqs(zparams);
        _available_init_freqs = true;
    }
    catch (...) {
        _available_init_freqs = false;
    }
    
    
    // Internal State
    _current_generation = 0;
    
    _initialized_with_parameters = true;
}

