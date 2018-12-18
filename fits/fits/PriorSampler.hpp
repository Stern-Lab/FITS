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

#ifndef PriorSampler_hpp
#define PriorSampler_hpp

#include <iostream>
#include <vector>
#include <chrono>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/random/mersenne_twister.hpp>
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
    FITNESS_COMPOSITE,
    SMOOTHED_COMPOSITE,
    FITNESS_LOGNORMAL,
    CUSTOM_FILE,
    UNDEFINED
};

template<class C>
class PriorSampler {
public:
    int GetRndNumber() {
        return _rnd_gen();
    }
    
private:
    std::vector<C> _min_vector;
    std::vector<C> _max_vector;
    PriorDistributionType _distrib_type;
    
    FLOAT_TYPE _log_normal_sigma, _log_normal_mu, _log_normal_lethal;
    
    // this is meant to be used with a posterior distribution
    // loaded as a new prior for this iteration
    MATRIX_TYPE _custom_distribution;
    
    boost::mt19937 _rnd_gen;
    unsigned int _rnd_seed;
    
    MATRIX_TYPE _distrib_matrix;
    
    const FLOAT_TYPE FITNESS_NEU_VAL = 1.0;
    
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
    
    void LoadDistribMatrix( std::string filename )
    {
        throw "PriorSampler LoadDistribMatrix unimplemented.";
    }
    
public:
    PriorSampler() :
    _min_vector(0),
    _max_vector(0),
    _distrib_matrix(40,2),
    _distrib_type(PriorDistributionType::UNIFORM),
    _log_normal_mu(-0.248f), _log_normal_sigma(0.149f), _log_normal_lethal(0.045f)
    {
        _rnd_seed = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
        _rnd_gen.seed(_rnd_seed);
        InitializeDistribMatrix();
    }
    
    PriorSampler( std::string filename ) :
    _min_vector(0),
    _max_vector(0),
    _distrib_matrix(1,1),
    _distrib_type(PriorDistributionType::CUSTOM_FILE),
    _log_normal_mu(0), _log_normal_sigma(0), _log_normal_lethal(0)
    {
        _rnd_seed = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
        _rnd_gen.seed(_rnd_seed);
        LoadDistribMatrix();
    }
    
    PriorSampler( std::vector<C> min_values, std::vector<C> max_values, PriorDistributionType distrib_type ) :
    _log_normal_mu(-0.248f), _log_normal_sigma(0.149f), _log_normal_lethal(0.045f),
    _distrib_type(distrib_type),
    _distrib_matrix(40,2)
    {
        _min_vector = min_values;
        _max_vector = max_values;
        
        if ( _min_vector.size() == 0 || _max_vector.size() == 0 ) {
            std::cerr << "No minimum or no maximum values given for prior sampling" << std::endl;
            throw "No minimum or no maximum values given for prior sampling";
        }
        
        _rnd_seed = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
        _rnd_gen.seed(_rnd_seed);
        InitializeDistribMatrix();
    }
    
    PriorSampler( boost::numeric::ublas::matrix<C> min_matrix,
                  boost::numeric::ublas::matrix<C> max_matrix,
                  PriorDistributionType distrib_type ) :
    _log_normal_mu(-0.248f), _log_normal_sigma(0.149f), _log_normal_lethal(0.045f),
    _min_vector(0),
    _max_vector(0),
    _distrib_type(distrib_type)
    {
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
    _log_normal_mu(original._log_normal_mu), _log_normal_sigma(original._log_normal_sigma), _log_normal_lethal(original._log_normal_lethal),
    _min_vector(original._min_vector),
    _max_vector(original._max_vector),
    _rnd_gen(0),
    _distrib_type(PriorDistributionType::UNIFORM)
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
                
            case PriorDistributionType::FITNESS_LOGNORMAL: {
                
                // values and algorithm inspired by Bons et al. doi: 10.1093/ve/vey029
                //FLOAT_TYPE sigma =
                //FLOAT_TYPE mu = -0.248f;
                //FLOAT_TYPE lethal_fraction = 0.045f;
                
                // not using the uniform01 for consistency with other distributions, e.g. () oprator
                boost::random::uniform_real_distribution<FLOAT_TYPE> lethal_distrib(0.0f, 1.0f);
                
                auto tmp_lethal_prob = lethal_distrib(_rnd_gen);
                if ( tmp_lethal_prob <= _log_normal_lethal ) {
                    tmp_val = 0.0f;
                }
                else {
                    boost::random::lognormal_distribution<FLOAT_TYPE> lognorm_distrib( _log_normal_mu, _log_normal_sigma );
                    do {
                        tmp_val = lognorm_distrib(_rnd_gen);
                    } while ( tmp_val < min || tmp_val > max );
                }
                break;
            }
            case PriorDistributionType::UNIFORM: {
                boost::random::uniform_real_distribution<FLOAT_TYPE> _uniform_dist(min, max);
                tmp_val = _uniform_dist(_rnd_gen);
                
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
