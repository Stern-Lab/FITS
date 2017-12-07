//
//  PriorSampler.hpp
//  fits
//
//  Created by Tal Zinger on 14/11/2016.
//  Copyright © 2016 Stern Lab. All rights reserved.
//

#ifndef PriorSampler_hpp
#define PriorSampler_hpp

#include <iostream>
#include <vector>
#include <chrono>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/lognormal_distribution.hpp>

#include <limits>

#include "CMulator.h"
#include "ZParams.h"
#include "fits_constants.h"

enum PriorDistributionType {
    UNIFORM,
    LOG_UNIFORM,
    FITNESS_COMPOSITE,
    SMOOTHED_COMPOSITE,
    LETHAL_UNIFORM,
    INTERVAL_UNIFORM,
    CUSTOM_FILE,
    UNDEFINED
};

template<class C>
class PriorSampler {
public:
    
    int GetRndNumber() {
     
        return _rnd_gen();
    }
    
    void SetManualCategoryProportions(ZParams zparams)
    {
        bool at_least_one_found = false;
        
        auto manual_lth = zparams.GetFloat( fits_constants::PARAM_PRIOR_FRACTION_LTH, fits_constants::PARAM_PRIOR_FRACTION_DEFAULT );
        auto manual_del = zparams.GetFloat( fits_constants::PARAM_PRIOR_FRACTION_DEL, fits_constants::PARAM_PRIOR_FRACTION_DEFAULT );
        auto manual_neu = zparams.GetFloat( fits_constants::PARAM_PRIOR_FRACTION_NEU, fits_constants::PARAM_PRIOR_FRACTION_DEFAULT );
        auto manual_adv = zparams.GetFloat( fits_constants::PARAM_PRIOR_FRACTION_ADV, fits_constants::PARAM_PRIOR_FRACTION_DEFAULT );
        
        at_least_one_found = at_least_one_found || (manual_lth >= 0.0);
        at_least_one_found = at_least_one_found || (manual_del >= 0.0);
        at_least_one_found = at_least_one_found || (manual_neu >= 0.0);
        at_least_one_found = at_least_one_found || (manual_adv >= 0.0);
        
        if (at_least_one_found) {
            auto all_found = (manual_del >= 0.0) && (manual_neu >= 0.0) && (manual_adv >= 0.0);
            
            if ( !all_found) {
                std::cerr << "Prior fraction LTH : " << ( manual_lth >= 0 ? "FOUND" : "MISSING" ) << std::endl;
                std::cerr << "Prior fraction DEL : " << ( manual_del >= 0 ? "FOUND" : "MISSING" ) << std::endl;
                std::cerr << "Prior fraction NEU : " << ( manual_neu >= 0 ? "FOUND" : "MISSING" ) << std::endl;
                std::cerr << "Prior fraction ADV : " << ( manual_adv >= 0 ? "FOUND" : "MISSING" ) << std::endl;
                
                throw "Coomposite prior: At least one but not all manual fractions are given.";
            }
            
            auto tmp_sum = manual_lth + manual_del + manual_neu + manual_adv;
            
            if (tmp_sum - 1.0f > std::numeric_limits<float>::epsilon() ) {
                std::cerr << "Prior fractions LTH/DEL/NEU/ADV do not sum up to 1.0 (" << tmp_sum << ")" << std::endl;
                throw "Prior fractions LTH/DEL/NEU/ADV do not sum up to 1.0";
            }
            
            _fitness_lth_prob = manual_lth;
            _fitness_del_prob = manual_del;
            _fitness_neu_prob = manual_neu;
            _fitness_adv_prob = manual_adv;
        }
    }
                                 
private:
    std::vector<C> _min_vector;
    std::vector<C> _max_vector;
    PriorDistributionType _distrib_type;
    
    // this is meant to be used with a posterior distribution
    // loaded as a new prior for this iteration
    MATRIX_TYPE _custom_distribution;
    
    
    // once these were generated for each sample. moved it here to make things faster
    //boost::numeric::ublas::matrix_column<MATRIX_TYPE> _fitness_composite_prob_col;
    //boost::numeric::ublas::matrix_column<MATRIX_TYPE> _fitness_composite_val_col;
    //boost::random::discrete_distribution<int> _fitness_composite_dist_disc;
    
    boost::mt19937 _rnd_gen;
    unsigned int _rnd_seed;
    
    MATRIX_TYPE _distrib_matrix;
    
    float _fitness_lth_prob;
    float _fitness_del_prob;
    float _fitness_neu_prob;
    float _fitness_adv_prob;
    
    const float FITNESS_NEU_VAL = 1.0;
    
    void InitializeDistribMatrix() {
        _distrib_matrix(0,0) = 0.1f;     _distrib_matrix(0,1) = 0.00f;
        _distrib_matrix(1,0) = 0.005f;   _distrib_matrix(1,1) = 0.05f;
        _distrib_matrix(2,0) = 0.005f;   _distrib_matrix(2,1) = 0.10f;
        _distrib_matrix(3,0) = 0.01f;    _distrib_matrix(3,1) = 0.15f;
        _distrib_matrix(4,0) = 0.01f;    _distrib_matrix(4,1) = 0.20f;
        _distrib_matrix(5,0) = 0.015f;   _distrib_matrix(5,1) = 0.25f;
        _distrib_matrix(6,0) = 0.015f;   _distrib_matrix(6,1) = 0.30f;
        _distrib_matrix(7,0) = 0.02f;    _distrib_matrix(7,1) = 0.35f;
        _distrib_matrix(8,0) = 0.02f;    _distrib_matrix(8,1) = 0.40f;
        _distrib_matrix(9,0) = 0.025f;   _distrib_matrix(9,1) = 0.45f;
        _distrib_matrix(10,0) = 0.025f;  _distrib_matrix(10,1) = 0.50f;
        _distrib_matrix(11,0) = 0.03f;   _distrib_matrix(11,1) = 0.55f;
        _distrib_matrix(12,0) = 0.03f;   _distrib_matrix(12,1) = 0.60f;
        _distrib_matrix(13,0) = 0.035f;  _distrib_matrix(13,1) = 0.65f;
        _distrib_matrix(14,0) = 0.04f;   _distrib_matrix(14,1) = 0.70f;
        _distrib_matrix(15,0) = 0.045f;  _distrib_matrix(15,1) = 0.75f;
        _distrib_matrix(16,0) = 0.05f;   _distrib_matrix(16,1) = 0.80f;
        _distrib_matrix(17,0) = 0.065f;  _distrib_matrix(17,1) = 0.85f;
        _distrib_matrix(18,0) = 0.07f;   _distrib_matrix(18,1) = 0.90f;
        _distrib_matrix(19,0) = 0.085f;  _distrib_matrix(19,1) = 0.95f;
        
        _distrib_matrix(20,0) = 0.1;     _distrib_matrix(20,1) = 1.00f;
        _distrib_matrix(21,0) = 0.05;    _distrib_matrix(21,1) = 1.05f;
        _distrib_matrix(22,0) = 0.03;    _distrib_matrix(22,1) = 1.10f;
        _distrib_matrix(23,0) = 0.02;    _distrib_matrix(23,1) = 1.15f;
        _distrib_matrix(24,0) = 0.02;    _distrib_matrix(24,1) = 1.20f;
        _distrib_matrix(25,0) = 0.01;    _distrib_matrix(25,1) = 1.25f;
        _distrib_matrix(26,0) = 0.01;    _distrib_matrix(26,1) = 1.30f;
        _distrib_matrix(27,0) = 0.01;    _distrib_matrix(27,1) = 1.35f;
        _distrib_matrix(28,0) = 0.01;    _distrib_matrix(28,1) = 1.40f;
        _distrib_matrix(29,0) = 0.01;    _distrib_matrix(29,1) = 1.45f;
        _distrib_matrix(30,0) = 0.005;   _distrib_matrix(30,1) = 1.50f;
        _distrib_matrix(31,0) = 0.005;   _distrib_matrix(31,1) = 1.55f;
        _distrib_matrix(32,0) = 0.005;   _distrib_matrix(32,1) = 1.60f;
        _distrib_matrix(33,0) = 0.005;   _distrib_matrix(33,1) = 1.65f;
        _distrib_matrix(34,0) = 0.005;   _distrib_matrix(34,1) = 1.70f;
        _distrib_matrix(35,0) = 0.001;   _distrib_matrix(35,1) = 1.75f;
        _distrib_matrix(36,0) = 0.001;   _distrib_matrix(36,1) = 1.80f;
        _distrib_matrix(37,0) = 0.001;   _distrib_matrix(37,1) = 1.85f;
        _distrib_matrix(38,0) = 0.001;   _distrib_matrix(38,1) = 1.90f;
        _distrib_matrix(39,0) = 0.001;   _distrib_matrix(39,1) = 1.95f;
    }
    
public:
    PriorSampler() :
    _min_vector(0),
    _max_vector(0),
    //_rnd_gen(0),
    _distrib_matrix(40,2),
    _fitness_lth_prob(fits_constants::FITNESS_LTH_PROB),
    _fitness_del_prob(fits_constants::FITNESS_DEL_PROB),
    _fitness_neu_prob(fits_constants::FITNESS_NEU_PROB),
    _fitness_adv_prob(fits_constants::FITNESS_ADV_PROB),
    _distrib_type(PriorDistributionType::UNIFORM)
    {
        _rnd_seed = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
        _rnd_gen.seed(_rnd_seed);
        InitializeDistribMatrix();
    }
    
    PriorSampler( std::vector<C> min_values, std::vector<C> max_values, PriorDistributionType distrib_type ) :
    //_min_vector(min_values),
    //_max_vector(max_values),
    _distrib_type(distrib_type),
    //_rnd_gen(0),
    _distrib_matrix(40,2),
    _fitness_lth_prob(fits_constants::FITNESS_LTH_PROB),
    _fitness_del_prob(fits_constants::FITNESS_DEL_PROB),
    _fitness_neu_prob(fits_constants::FITNESS_NEU_PROB),
    _fitness_adv_prob(fits_constants::FITNESS_ADV_PROB)
    {
        _min_vector = min_values;
        _max_vector = max_values;
        
        if ( _min_vector.size() == 0 || _max_vector.size() == 0 ) {
            std::cerr << "No minimum or no maximum values given for prior sampling" << std::endl;
            throw "No minimum or no maximum values given for prior sampling";
        }
        
        //std::cout << "prior init...";
        _rnd_seed = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
        _rnd_gen.seed(_rnd_seed);
        InitializeDistribMatrix();
        
        //_fitness_composite_dist_disc( _fitness_composite_prob_col );
        //std::cout << "prior init done" << std::endl;
    }
    
    PriorSampler( boost::numeric::ublas::matrix<C> min_matrix,
                  boost::numeric::ublas::matrix<C> max_matrix,
                  PriorDistributionType distrib_type ) :
    _min_vector(0),
    _max_vector(0),
    _distrib_type(distrib_type),
    //_rnd_gen(0),
    _fitness_lth_prob(fits_constants::FITNESS_LTH_PROB),
    _fitness_del_prob(fits_constants::FITNESS_DEL_PROB),
    _fitness_neu_prob(fits_constants::FITNESS_NEU_PROB),
    _fitness_adv_prob(fits_constants::FITNESS_ADV_PROB)
    {
        //std::cout << "prior init123...";
        
        if ( _min_vector.size() != _max_vector.size() ) {
            std::cerr << "Prior sampler: size of matrices do not match: min (" <<
            min_matrix.size1() << "," << min_matrix.size2() << ") vs max (" <<
            max_matrix.size1() << "," << max_matrix.size2() << ")" << std::endl;
            throw "Prior sampler: size of matrices do not match.";
        }
        
        for ( auto row=0; row<min_matrix.size1(); ++row ) {
            for ( auto col=0; col<min_matrix.size2(); ++col ) {
                _min_vector.push_back( min_matrix(row,col) );
                _max_vector.push_back( max_matrix(row,col) );
            }
        }
    }
    
    PriorSampler( PriorSampler<C> &original ) :
    _min_vector(original._min_vector),
    _max_vector(original._max_vector),
    _rnd_gen(0),
    _distrib_type(PriorDistributionType::UNIFORM),
    _fitness_lth_prob(original._fitness_lth_prob),
    _fitness_del_prob(original._fitness_del_prob),
    _fitness_neu_prob(original._fitness_neu_prob),
    _fitness_adv_prob(original._fitness_adv_prob)
    {
        _rnd_seed = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
        _rnd_gen.seed(_rnd_seed);
    }
    
    // TODO: write special method for matrix/vector/single sampling
    std::vector<std::vector<C>> SamplePrior( std::size_t num_samples )
    {
        std::vector<std::vector<C>> return_vector;
        return_vector.reserve(num_samples);
        
        for ( auto current_sample=0; current_sample<num_samples; ++current_sample ) {
            
            std::vector<C> tmp_vector;
            tmp_vector.reserve(_min_vector.size());
            
            for ( auto current_item=0; current_item<_min_vector.size(); ++current_item ) {
                tmp_vector.push_back( SampleFromDistribution(_min_vector[current_item], _max_vector[current_item]) );
            }
            
            if ( tmp_vector.size() == 0 ) {
                std::cerr << "sample number " << current_sample << " vector is empty" << std::endl;
            }
            //return_vector.push_back(std::move(tmp_vector));
            return_vector.push_back(tmp_vector);
        }
        
        return return_vector;
    }
    
private:
    C SampleFromDistribution( C min, C max )
    {
        C tmp_val = 0.0f;
        
        if ( min == max ) {
            return min;
        }
        
        if ( min > max ) {
            std::cerr << "Sampling from prior, but " << min << " is larger than " << max << std::endl;
            throw "min>max";
        }
        switch (_distrib_type) {
                
            case PriorDistributionType::UNIFORM: {
                boost::random::uniform_real_distribution<FLOAT_TYPE> _uniform_dist(min, max);
                tmp_val = _uniform_dist(_rnd_gen);
                
                break;
            }
                
            case PriorDistributionType::INTERVAL_UNIFORM: {
                FLOAT_TYPE fitness_interval = 0.001f;
                auto diff = max - min;
                auto num_intervals = static_cast<int>(std::round(diff / fitness_interval));
                boost::random::uniform_int_distribution<> dist(0, num_intervals);
                
                tmp_val = min + (dist(_rnd_gen)-2) * fitness_interval;
                if ( tmp_val < min ) tmp_val = min;
                break;
            }
            
            case PriorDistributionType::LOG_UNIFORM: {
                boost::random::uniform_01<> dist;
                auto power = min + (max-min) * dist(_rnd_gen);
                tmp_val = std::pow(10.0, power);
                break;
            }
                
            case PriorDistributionType::LETHAL_UNIFORM: {
                FLOAT_TYPE my_lethal_prob = 0.1;
                boost::random::discrete_distribution<int> dist_disc( { my_lethal_prob, 1.0f-my_lethal_prob } );
                boost::random::uniform_real_distribution<FLOAT_TYPE> dist_uniform(min, max);
                
                auto chosen_category = dist_disc(_rnd_gen);
                
                if ( chosen_category == 0 ) {
                    //std::cout << "lethal chosen" << std::endl;
                    tmp_val = 0.0f;
                }
                else {
                    //std::cout << "non-lethal chosen" << std::endl;
                    tmp_val = dist_uniform(_rnd_gen);
                }
                break;
            }
                
            case PriorDistributionType::SMOOTHED_COMPOSITE:
            case PriorDistributionType::FITNESS_COMPOSITE: {
                
                // Moving these to be member variables causes some weird call for boost's uniform int distribution
                // with inexsiting matrix - does not compile and hard to find the cause here since boost's headers
                // are blamed by the compiler.
                boost::numeric::ublas::matrix_column<MATRIX_TYPE> fitness_composite_prob_col(_distrib_matrix, 0);
                boost::numeric::ublas::matrix_column<MATRIX_TYPE> fitness_composite_val_col(_distrib_matrix, 1);
                
                
                boost::random::discrete_distribution<int> fitness_composite_dist_disc( fitness_composite_prob_col );
                
                
                do {
                    auto tmp_idx = fitness_composite_dist_disc(_rnd_gen);
                    
                    if ( _distrib_type == SMOOTHED_COMPOSITE ) {
                        
                        // the min is the last item - nothing to smooth
                        if ( tmp_idx == fitness_composite_val_col.size() - 1 ) {
                            tmp_val = fitness_composite_val_col(tmp_idx);
                        }
                        else {
                            auto tmp_min = fitness_composite_val_col(tmp_idx);
                            auto tmp_max = fitness_composite_val_col(tmp_idx+1);
                            
                            if ( tmp_min == tmp_max ) {
                                tmp_val = tmp_min;
                            }
                            else {
                                //std::cout << "uniform between: " << tmp_min << " " << tmp_max << std::endl;
                                boost::random::uniform_real_distribution<FLOAT_TYPE> smoothed_distrib( tmp_min, tmp_max );
                                tmp_val = smoothed_distrib(_rnd_gen);
                            }
                        }
                    }
                    else {
                        // no smoothing
                        tmp_val = fitness_composite_val_col(tmp_idx);
                    }
                    
                } while (tmp_val < min || tmp_val > max );
                
                break;
            }
            default:
                std::cerr << "Undefined prior distribution." << std::endl;
                throw "Undefined prior distribution.";
        }
        
        return tmp_val;
    }
    
};

#endif /* PriorSampler_hpp */
