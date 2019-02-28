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

#include "ResultsStats.hpp"

void ResultsStats::CalculateStatsMutation()
{
    CMulator local_sim(_zparams);
    
    _num_alleles = _result_vector[0].fitness_values.size();
    
    _num_results = _result_vector.size();
    
    boost::accumulators::accumulator_set<
    FLOAT_TYPE,
    boost::accumulators::stats<
    boost::accumulators::tag::variance,
    boost::accumulators::tag::mean,
    boost::accumulators::tag::median,
    boost::accumulators::tag::min,
    boost::accumulators::tag::max> > acc_distance;
    
    boost::numeric::ublas::matrix< boost::accumulators::accumulator_set<
    FLOAT_TYPE,
    boost::accumulators::stats<
    boost::accumulators::tag::variance,
    boost::accumulators::tag::mean,
    boost::accumulators::tag::median,
    boost::accumulators::tag::min,
    boost::accumulators::tag::max> > > acc_matrix(_num_alleles,_num_alleles);
    
    boost::numeric::ublas::matrix< std::vector<FLOAT_TYPE> > raw_rates_matrix(_num_alleles,_num_alleles);
    
    prior_matrix.resize(_num_alleles,_num_alleles);
    
    min_mutation_rates.resize(_num_alleles, _num_alleles);
    max_mutation_rates.resize(_num_alleles, _num_alleles);
    mean_mutation_rates.resize(_num_alleles, _num_alleles);
    median_mutation_rates.resize(_num_alleles, _num_alleles);
    // normalized_median_mutation_rates.resize(_num_alleles, _num_alleles);
    levenes_pval_matrix.resize(_num_alleles,_num_alleles);
    mad_mutation_rates.resize(_num_alleles, _num_alleles);
    
    for ( auto sim_result : _result_vector ) {
        
        acc_distance(sim_result.distance_from_actual);
        
        for ( auto row=0; row<_num_alleles; ++row ) {
            
            for ( auto col=0; col<_num_alleles; ++col ) {
                
                acc_matrix(row,col)( sim_result.mutation_rates(row,col) );
                raw_rates_matrix(row,col).push_back(sim_result.mutation_rates(row,col));
            }
        }
    }
    
    for ( auto row=0; row<_num_alleles; ++row ) {
        
        for ( auto col=0; col<_num_alleles; ++col ) {
            
            min_mutation_rates(row,col) = boost::accumulators::min( acc_matrix(row,col) );
            max_mutation_rates(row,col) = boost::accumulators::max( acc_matrix(row,col) );
            
            mean_mutation_rates(row,col) = boost::accumulators::mean( acc_matrix(row,col) );
            
            
            // median_mutation_rates(row,col) = GetMedian(median_matrix(row,col));
            median_mutation_rates(row,col) = boost::accumulators::median( acc_matrix(row,col) );
        }
    }
    
    
    for ( auto row=0; row<_num_alleles; ++row ) {
        
        for ( auto col=0; col<_num_alleles; ++col ) {
            
            boost::accumulators::accumulator_set<
            FLOAT_TYPE,
            boost::accumulators::stats<
            boost::accumulators::tag::median > > acc_distance_for_mad;
            
            for ( auto current_rate : raw_rates_matrix(row, col) ) {
                
                acc_distance_for_mad( std::fabs( current_rate - median_mutation_rates(row,col) ) );
            }
            
            mad_mutation_rates(row, col) = boost::accumulators::median(acc_distance_for_mad);
        }
    }
    
    _distance_min = boost::accumulators::min(acc_distance);
    _distance_max = boost::accumulators::max(acc_distance);
    _distance_mean = boost::accumulators::mean(acc_distance);
    _distance_sd = std::sqrt(boost::accumulators::variance(acc_distance));
    
    /*
    for ( auto row=0; row<_num_alleles; ++row ) {
        
        FLOAT_TYPE tmp_sum = 0.0f;
        
        for ( auto col=0; col<_num_alleles; ++col ) {
            
            tmp_sum += median_mutation_rates(row,col);
        }
        
        for ( auto col=0; col<_num_alleles; ++col ) {
            
            normalized_median_mutation_rates(row,col) = median_mutation_rates(row,col) / tmp_sum;
        }
    }
    */
    
    // levene's pval
    auto min_prior_mutation_rates = local_sim.GetMinMutationRateMatrix();
    auto max_prior_mutation_rates = local_sim.GetMaxMutationRateMatrix();
    
    PriorSampler<int> sampler( min_prior_mutation_rates,
                              max_prior_mutation_rates,
                              PriorDistributionType::UNIFORM);
    //std::cout << "min rates: " << min_mutation_rates << std::endl;
    //std::cout << "max rates: " << max_mutation_rates << std::endl;
    // auto mutrate_vector_list = sampler.SamplePrior( _num_results );
    
    if ( _zparams.GetInt( "Debug", 0 ) > 0 ) {
        std::cout << "building prior" << std::endl;
    }
    
    
    // for (auto current_mutrate_vector : mutrate_vector_list) {
    for ( auto current_mutrate_vector : _prior_distrib ) {
        for ( auto i=0; i<current_mutrate_vector.size(); ++i ) {
            
            auto row = i / _num_alleles;
            auto col = i % _num_alleles;
            
            prior_matrix(row,col).push_back( std::pow( 10, current_mutrate_vector[i] ) );
            
            if ( _zparams.GetInt( "Debug", 0 ) > 0 ) {
                std::cout << std::pow( 10, current_mutrate_vector[i] ) << fits_constants::FILE_FIELD_DELIMITER;
            }
            
        }
    }
    
    for ( auto row=0; row<_num_alleles; ++row ) {
        
        for ( auto col=0; col<_num_alleles; ++col ) {
            
            if ( row==col) {
                levenes_pval_matrix(row,col) = -1.0f;
            }
            else {
                levenes_pval_matrix(row,col) = LevenesTest2( raw_rates_matrix(row,col),
                                                            prior_matrix(row,col) );
            }
            
        }
    }
    //std::cout << "mutrate matrix:"<< std::endl;
    //std::cout << levenes_pval_matrix;
}


std::string ResultsStats::GetSummaryMutRate( bool table_only )
{
    std::stringstream ss;
    
    if ( table_only ) {
        
        // header
        ss << "from" << fits_constants::FILE_FIELD_DELIMITER;
        ss << "to" << fits_constants::FILE_FIELD_DELIMITER;
        ss << "median" << fits_constants::FILE_FIELD_DELIMITER;
        ss << "MAD" << fits_constants::FILE_FIELD_DELIMITER;
        ss << "min" << fits_constants::FILE_FIELD_DELIMITER;
        ss << "max" << fits_constants::FILE_FIELD_DELIMITER;
        ss << "pval" << std::endl;
        
        // fill data
        for ( auto row=0; row<_num_alleles; ++row ) {
            
            for (auto col=0; col<_num_alleles; ++col ) {
                
                if ( row == col ) {
                    continue;
                }
                
                auto tmp_median = median_mutation_rates(row,col);
                auto tmp_mad = mad_mutation_rates(row,col);
                auto tmp_min = min_mutation_rates(row,col);
                auto tmp_max = max_mutation_rates(row,col);
                auto tmp_pval = levenes_pval_matrix(row,col);
                
                ss << row << fits_constants::FILE_FIELD_DELIMITER;
                ss << col << fits_constants::FILE_FIELD_DELIMITER;
                
                // median
                if ( levenes_pval_matrix(row,col) < fits_constants::LEVENES_SIGNIFICANCE ) {
                    ss << boost::format("%-.2e") % tmp_median << fits_constants::FILE_FIELD_DELIMITER;
                }
                else {
                    ss << boost::format("*%-.2e") % tmp_median << fits_constants::FILE_FIELD_DELIMITER;
                }
                
                // the rest
                ss << boost::format("%-.2e") % tmp_mad << fits_constants::FILE_FIELD_DELIMITER;
                ss << boost::format("%-.2e") % tmp_min << fits_constants::FILE_FIELD_DELIMITER;
                ss << boost::format("%-.2e") % tmp_max << fits_constants::FILE_FIELD_DELIMITER;
                ss << boost::format("%-.2e") % tmp_pval;
                
                if ( col < _num_alleles-1 ) {
                    ss << fits_constants::FILE_FIELD_DELIMITER;
                }
                
                if ( row < _num_alleles-1 || col < _num_alleles-1 ) {
                    ss << std::endl;
                }
            }
        } // row
        
        return ss.str();
    } // table only
    
    
    ss << "Mutation Rate Report" << std::endl;
    ss << "===============" << std::endl;
    ss << GetSummaryHeader();
    
    
    
    //auto current_time_raw = std::chrono::system_clock::now();
    //auto current_time = std::chrono::system_clock::to_time_t(current_time_raw);
    //auto current_time_final = *std::localtime(&current_time);
    //ss << std::put_time(&current_time_final, "%F %T") << std::endl;
    
    //std::cout << "Simulation results used for calculations: " << _num_results << std::endl;
    
    //ss << "=====================" << std::endl;
    
    //if ( _rejection_threshold > 0.0f ) {
    //   ss << "Rejection threshold set to " << boost::format("%-10d") % _rejection_threshold << std::endl;
    //}
    
    
    
    ss << "Population size (N) is " << _zparams.GetInt(fits_constants::PARAM_POPULATION_SIZE, -1) << std::endl;
    
    auto tmp_scaling_str = _zparams.GetString( fits_constants::PARAM_SCALING,
                                              fits_constants::PARAM_SCALING_DEFAULT );
    
    if ( tmp_scaling_str.compare(fits_constants::PARAM_SCALING_OFF) == 0 ) {
        // ss << "Data has not been scaled" << std::endl;
    }
    else {
        ss << "Data was scaled using " << tmp_scaling_str << std::endl;
    }
    
    
    if ( _single_mutrate_used ) {
        ss << "Used a single mutation rate." << std::endl;
    }
    else {
        ss << "Used individual mutation rates." << std::endl;
    }
    
    ss << "Distance metric: " << _distance_metric << std::endl;
    
    if ( _zparams.GetInt(fits_constants::PARAM_SAMPLE_SIZE, 0) > 0 ) {
        ss << " (sampled " << _zparams.GetInt(fits_constants::PARAM_SAMPLE_SIZE, 0) << ")" << std::endl;
    }
    
    if ( _single_mutrate_inferred ) {
        ss << "Inferred a single mutation rate." << std::endl;
    }
    
    ss << "====================" << std::endl;
    
    
    
    // header
    ss << boost::format("%-6s") % "from";
    ss << boost::format("%-6s") % "to";
    ss << boost::format("%-12s") % "median";
    ss << boost::format("%-12s") % "MAD";
    ss << boost::format("%-12s") % "min";
    ss << boost::format("%-12s") % "max";
    ss << boost::format("%-12s") % "pval" << std::endl;
    
    // fill data
    for ( auto row=0; row<_num_alleles; ++row ) {
        
        for (auto col=0; col<_num_alleles; ++col ) {
            
            if ( row == col ) {
                continue;
            }
            
            auto tmp_median = median_mutation_rates(row,col);
            auto tmp_mad = mad_mutation_rates(row,col);
            auto tmp_min = min_mutation_rates(row,col);
            auto tmp_max = max_mutation_rates(row,col);
            auto tmp_pval = levenes_pval_matrix(row,col);
            
            
            ss << boost::format("%-6s") % row;
            ss << boost::format("%-6s") % col;
            
            // median
            if ( levenes_pval_matrix(row,col) < fits_constants::LEVENES_SIGNIFICANCE ) {
                ss << boost::format("%-12.2e") % tmp_median;
            }
            else {
                ss << boost::format("*%-11.2e") % tmp_median;
            }
            
            // the rest
            ss << boost::format("%-12.2e") % tmp_mad;
            ss << boost::format("%-12.2e") % tmp_min;
            ss << boost::format("%-12.2e") % tmp_max;
            ss << boost::format("%-12.2e") % tmp_pval;
            ss << std::endl;
        }
    } // row
    
    
    if ( _zparams.GetInt( fits_constants::PARAM_DUMP_PARAMETERS, 0) > 1 ) {
        ss << _zparams.GetAllParameters() << std::endl;
    }
    
    
    return ss.str();
    
}
