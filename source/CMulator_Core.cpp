//
//  CMulator_Core.cpp
//  fits
//
//  Created by Tal Zinger on 15/11/2016.
//  Copyright © 2016 Stern Lab. All rights reserved.
//

#include "CMulator.h"


void CMulator::PerformChecksBeforeEvolution()
{
    if (!_initialized_with_parameters) {
        std::cerr << " EvolveToGeneration: object not initialized with parameters." << std::endl;
        throw " EvolveToGeneration: object not initialized with parameters.";
    }
    
    if ( _allele_init_freqs.empty() ) {
        std::cerr << " EvolveToGeneration: initial frequencies not provided." << std::endl;
        throw " EvolveToGeneration: initial frequencies not provided.";
    }
    
    /*
    if ( _sim_data.empty() ) {
        std::cerr <<	" EvolveToGeneration: sim_data uninitialized." << std::endl;
        throw " EvolveToGeneration: sim_data uninitialized.";
    }
    */
    
    if ( _all_simulated_data.size1() < 1 || _all_simulated_data.size2() < 1 ) {
        std::cerr <<	" EvolveToGeneration: sim_data uninitialized." << std::endl;
        throw " EvolveToGeneration: _all_simulated_data uninitialized.";
    }
    
    if ( _allele_fitness.empty() ) {
        if ( _allele_max_fitness.empty() || _allele_min_fitness.empty() ) {
            std::cerr <<	" EvolveToGeneration: fitness data empty but so are min or max." << std::endl;
            throw " EvolveToGeneration: fitness data empty but so are min and max.";
        }
    }
}

std::vector<FLOAT_TYPE> CMulator::Freqs2Binomial2Freqs( const std::vector<FLOAT_TYPE> &freqs_vec, int popsize )
{
    if ( freqs_vec.size() != _num_alleles ) {
        std::cerr << "Freqs2Binomial2Freqs: mismatch in size, freqs_vec is " <<
        freqs_vec.size() << " while _num_alleles is " << _num_alleles << std::endl;
        
        throw "Freqs2Binomial2Freqs: size mismatch between vector of values to number of alleles.";
    }
    
    auto currentSum = 0.0f;
    
    std::vector<FLOAT_TYPE> newvals(freqs_vec);
    
    // perform binomial sampling, keep in mind actually selected individuals may be !=N
    for (auto allele_i = 0; allele_i<_num_alleles; ++allele_i) {
        
        auto currentP = freqs_vec[allele_i];
        
        boost::random::binomial_distribution<int> local_bin_distrib( popsize, currentP );
        
        newvals[allele_i] = static_cast<FLOAT_TYPE>(local_bin_distrib(_boost_gen));
        
        currentSum += newvals[allele_i];
    }
    
    // Normalize frequencies
    for ( auto allele_i = 0; allele_i<_num_alleles; ++allele_i ) {
        newvals[allele_i] = newvals[allele_i] / currentSum;
    }

    return newvals;
}


int CMulator::EvolveToGeneration( int target_generation )
{
    
    // this distribution is unique - it takes a closed range [inclusive]
    //if ( _debug_mode ) {
     //   std::cout << "Initializing local int distribution for alt N to max generation " << target_generation-1 << std::endl;
    //}
    
    boost::random::uniform_int_distribution<int> local_int_distrib(1, target_generation-1);
    
    if ( _alt_N > 0 && _alt_generation < 0) {
        _alt_generation = local_int_distrib(_boost_gen);
    }
    
    PerformChecksBeforeEvolution();
    
    // just to inform the calling function what was the initial generation
    int init_generation = _current_generation;
    
    //std::cout << " current generation " << _current_generation << std::endl;
    //std::vector<FLOAT_TYPE> newvals(_num_alleles);
    //FLOAT_TYPE currentSum = 0.0;
    
    // should we revert to original population size
    bool alt_popsize_flag = false;
    
    //if ( _debug_mode ) {
    //    std::cout << "Starting evolution for loop" << std::endl;
    //}
    
    
    // copy current frequencies to a convinient vector
    //std::vector<FLOAT_TYPE> previous_generation_freqs(_num_alleles);
    //std::vector<FLOAT_TYPE> alt_previous_generation_freqs(_num_alleles);
    
    // generation 0 is given, evolve until desired generation
    for (_current_generation=1; _current_generation<=target_generation; _current_generation++) {
        
        boost::numeric::ublas::matrix_row<MATRIX_TYPE> prev_gen_row( _all_simulated_data, _current_generation-1 );
        /*
        if ( _debug_mode ) {
            std::cout << "previous generation: ";
            for ( auto val : prev_gen_row ) std::cout << val << "\t";
            std::cout << std::endl;
        }
        */
        
        if ( _alt_generation == _current_generation ) {
            _old_N = GetPopulationSize();
            SetPopulationSize(_alt_N);
            alt_popsize_flag = true;
            /*
            if (_debug_mode) {
                std::cout << "Generation "
                << _current_generation
                << ": changed to alt popsize "
                << _N << " instead of "
                << _old_N << std::endl;
            }
             */
        }
        
        if ( alt_popsize_flag ) {
            SetPopulationSize(_old_N);
            alt_popsize_flag = false;
            
            /*
            if (_debug_mode) {
                std::cout << "Generation "
                << _current_generation
                << ": reverted popsize to "
                << _N << std::endl;
            }
             */
        }
        
        
        // normalize fitness
        
        //if ( _debug_mode ) {
        //    std::cout << "Normalizing fitness." << std::endl;
        //}
        
        std::vector<FLOAT_TYPE> allele_fitness_bar(_num_alleles, 0.0f);
        std::vector<FLOAT_TYPE> allele_fitness_adjusted(_num_alleles, 0.0f);
        
        std::transform( _allele_fitness.cbegin(), _allele_fitness.cend(),
                       prev_gen_row.cbegin(),
                       allele_fitness_bar.begin(),
                       []( FLOAT_TYPE fitness, FLOAT_TYPE freq ) {
                           return fitness*freq;
                       } );
        
        auto wbar = std::accumulate( allele_fitness_bar.cbegin(),
                                     allele_fitness_bar.cend(), 0.0f );
        
        
        // normalize fitness
        std::transform( _allele_fitness.cbegin(), _allele_fitness.cend(),
                        allele_fitness_adjusted.begin(),
                        [wbar]( FLOAT_TYPE f){ return f / wbar; } );
        
        
        /*
        if ( _debug_mode ) {
            std::cout << "fitness values: ";
            for ( auto val : _allele_fitness ) std::cout << val << "\t";
            std::cout << std::endl;
            
            std::cout << "fitness values bar: ";
            for ( auto val : allele_fitness_bar ) std::cout << val << "\t";
            std::cout << std::endl;
            
            std::cout << "adjusted fitnesss: ";
            for ( auto val : allele_fitness_adjusted ) std::cout << val << "\t";
            std::cout << std::endl;
            std::cout << "wbar = " << wbar << std::endl;
            std::cout << "Multiplications" << std::endl;
        }
         */
        // TODO: replace these nested for loops with matrix multiplication
        // vector of last generation X mutation_rate_matrix -> tmp1_vector
        // tmp1_vector X fitness_vector -> new_generation
        std::vector<FLOAT_TYPE> vec_after_deterministic(_num_alleles, 0.0f);
        for ( auto currenly_evolving_allele=0; currenly_evolving_allele<_num_alleles; currenly_evolving_allele++ ) {
            
            for ( auto currently_other_allele=0; currently_other_allele<_num_alleles; currently_other_allele++ ) {
                
                auto using_matrix_sim_data =
                _mutation_rates_matrix(currently_other_allele,currenly_evolving_allele) *
                prev_gen_row[currently_other_allele] *
                allele_fitness_adjusted[currently_other_allele];
                
                vec_after_deterministic[currenly_evolving_allele] += using_matrix_sim_data;
            }
        }
        
        //if ( _debug_mode ) {
        //    std::cout << "Normalize probailities" << std::endl;
        //}

        // normalize probabilities
        FLOAT_TYPE sumall = std::accumulate( vec_after_deterministic.begin(),
                                             vec_after_deterministic.end(), 0.0f );
        
        
        std::transform( vec_after_deterministic.begin(),
                        vec_after_deterministic.end(),
                        vec_after_deterministic.begin(),
                        [sumall]( FLOAT_TYPE f ) { return f/sumall; });
        
        // logistic growth
        if (_logistic_growth) {
            _N = static_cast<int>( _logistic_growth_K *(_N0/(_N0+(exp(-_logistic_growth_r*static_cast<FLOAT_TYPE>(_logistic_growth_t))*(_logistic_growth_K-_N0)))) );
            _logistic_growth_t++;
            
        }
        
        //if ( _debug_mode ) {
        //    std::cout << "Stochastic step - Freqs2Binomial2Freqs; population size = " << _N << std::endl;
        //}
        
        auto vec_after_stochastic = Freqs2Binomial2Freqs(vec_after_deterministic, _N);
        
        
        // observed data is before bottleneck
        if ( _sample_size>0 ) {
            
            //if ( _debug_mode ) {
            //    std::cout << "Observation - Freqs2Binomial2Freqs; population size = " << _sample_size << std::endl;
            //}
            
            auto vec_observed = Freqs2Binomial2Freqs( vec_after_deterministic, _sample_size );
            
            for (auto allele_i = 0; allele_i<_num_alleles; ++allele_i) {
                _observed_simulated_data( _current_generation, allele_i ) = vec_observed[allele_i];
            }
        }
        
        
        // apply bottleneck if relevant
        if ( (_bottleneck_interval) > 0 &&
            (_current_generation % _bottleneck_interval == 0) ) {
            
            _logistic_growth_t = 0;
            _N0 = _bottleneck_size;	// this should be like that as this is used when logistic growth is in use.
            
            // like stochastic step but with other population size
            // brackets remains from multithreaded period where mutex was used via lock guard (has to be scoped)
        
            //if ( _debug_mode ) {
            //    std::cout << "Bottleneck step - Freqs2Binomial2Freqs; population size = " << _bottleneck_size << std::endl;
            //}
            
            vec_after_stochastic = Freqs2Binomial2Freqs(vec_after_deterministic, _bottleneck_size);
            
        } // bottleneck
        
        boost::numeric::ublas::matrix_row<MATRIX_TYPE> current_gen_row( _all_simulated_data, _current_generation );
        std::copy( vec_after_stochastic.cbegin(), vec_after_stochastic.cend(), current_gen_row.begin() );
        
        
        /*
        for (auto allele_i = 0; allele_i<_num_alleles; ++allele_i) {
            _all_simulated_data( _current_generation, allele_i ) = vec_after_stochastic[allele_i];
        }
        */
        
        
    } // evolution for loop
    
    return init_generation;
}


int CMulator::EvolveAllGenerations()
{
    if (!_initialized_with_parameters) {
        throw " EvolveAllGenerations: object not initialized with parameters.";
    }
    
    return EvolveToGeneration(_num_generations);
}

