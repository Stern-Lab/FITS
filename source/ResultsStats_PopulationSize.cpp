//
//  ResultsStats_PopulationSize.cpp
//  fits
//
//  Created by Tal Zinger on 16/11/2016.
//  Copyright © 2016 Stern Lab. All rights reserved.
//

#include "ResultsStats.hpp"


void ResultsStats::CalculateStatsPopulationSize(const std::vector<SimulationResult>& result_vector)
{
    _num_results = result_vector.size();
    
    boost::accumulators::accumulator_set<
    FLOAT_TYPE,
    boost::accumulators::stats<
    boost::accumulators::tag::variance,
    boost::accumulators::tag::mean,
    boost::accumulators::tag::min,
    boost::accumulators::tag::max> > acc_distance;
    
    boost::accumulators::accumulator_set<
    FLOAT_TYPE,
    boost::accumulators::stats<
    boost::accumulators::tag::variance,
    boost::accumulators::tag::mean,
    boost::accumulators::tag::min,
    boost::accumulators::tag::max> > acc_population_size;
    

    std::vector<int> popsize_storage;
    for ( auto sim_result : result_vector ) {
        acc_distance(sim_result.distance_from_actual);
        
        acc_population_size(sim_result.N);
        
        //std::cout << "N=" << sim_result.N << std::endl;
        popsize_storage.push_back(sim_result.N);
    }
    
    _pop_min = boost::accumulators::min(acc_population_size);
    _pop_max = boost::accumulators::max(acc_population_size);
    _pop_mean = boost::accumulators::mean(acc_population_size);
    _pop_sd = std::sqrt(boost::accumulators::variance(acc_population_size));
    
    _pop_median = GetMedian(popsize_storage);
    
    _distance_min = boost::accumulators::min(acc_distance);
    _distance_max = boost::accumulators::max(acc_distance);
    _distance_mean = boost::accumulators::mean(acc_distance);
    _distance_sd = std::sqrt(boost::accumulators::variance(acc_distance));
    
    // levene's test for population size
    levenes_pval.resize(1, -1.0f);
    PriorDistributionType prior_type = PriorDistributionType::UNIFORM;
    
    std::vector<int> minN {_zparams.GetInt(fits_constants::PARAM_MIN_LOG_POPSIZE)};
    std::vector<int> maxN {_zparams.GetInt(fits_constants::PARAM_MAX_LOG_POPSIZE)};
    
    PriorSampler<int> sampler( minN, maxN, PriorDistributionType::UNIFORM );
    
    auto popsize_vector_list = sampler.SamplePrior(_num_results);
    std::vector<int> prior_vec_int;
    for ( auto current_vec : popsize_vector_list ) {
        prior_vec_int.push_back( std::pow( 10, current_vec[0] ) );
    }
    
    //auto prior_vec_int = sampler.SamplePrior(_num_results)[1];
    
    /*
    std::cout << "prior size=" << prior_vec_int.size() << std::endl;
    for ( auto current_n : prior_vec_int ) {
        std::cout << current_n << "\t";
    }
    std::cout << std::endl;
    
    std::cout << "posteiror size=" << popsize_storage.size() << std::endl;
    for ( auto current_n : popsize_storage ) {
        std::cout << current_n << "\t";
    }
    std::cout << std::endl;
    */
    
    std::vector<float> posterior_vec;
    std::vector<float> prior_vec_float; // this provides conversion
    
    prior_vec_float.resize( prior_vec_int.size() );
    posterior_vec.resize( prior_vec_int.size() );
    
    std::copy( prior_vec_int.cbegin(), prior_vec_int.cend(), prior_vec_float.begin() );
    std::copy( popsize_storage.cbegin(), popsize_storage.cend(), posterior_vec.begin() );
    
    /*
    std::transform( prior_vec_int.cbegin(), prior_vec_int.cend(),
                   prior_vec_float.begin(),
                   []( int val ) {
                       return static_cast<float>(val);
                   } );
     */
    /*
    std::transform( popsize_storage.cbegin(), popsize_storage.cend(),
                   posterior_vec.begin(),
                   []( int val ) {
                       return static_cast<float>(val);
                   } );
     */
    
    /*
    std::cout << "prior2 size=" << prior_vec_float.size() << std::endl;
    for ( auto current_n : prior_vec_float ) {
        std::cout << current_n << "\t";
    }
    std::cout << std::endl;
    
    std::cout << "posteiror2 size=" << posterior_vec.size() << std::endl;
    for ( auto current_n : posterior_vec ) {
        std::cout << current_n << "\t";
    }
    std::cout << std::endl;
     */
    
    levenes_pval[0] = LevenesTest2(posterior_vec, prior_vec_float );
    
    //std::cout << "levenes p = " << levenes_pval[0] <<  std::endl;
}


std::string ResultsStats::GetSummaryPopSize()
{
    std::stringstream ss;
    
    ss.imbue(std::locale(fits_constants::used_locale));
    
    
    ss << "Population Size Report" << std::endl;
    ss << GetSummaryHeader();
    
    //ss << "FITS v"<< fits_constants::current_version_str << std::endl;
    
    //auto current_time_raw = std::chrono::system_clock::now();
    //auto current_time = std::chrono::system_clock::to_time_t(current_time_raw);
    //auto current_time_final = *std::localtime(&current_time);
    //ss << std::put_time(&current_time_final, "%F %T") << std::endl;
    
    //std::cout << "Simulation results used for calculations: " << _num_results << std::endl;
    
    //ss << "=======================" << std::endl;
    
    //if ( _rejection_threshold > 0.0f ) {
      //  ss << "Rejection threshold set to " << boost::format("%-10d") % _rejection_threshold << std::endl;
    //}
    
    
    if ( _single_mutrate_used ) {
        ss << "Used a single mutation rate." << std::endl;
    }
    
    ss << "====================" << std::endl;
    
    ss << boost::format("%-12s") % "median";
    ss << boost::format("%-12s") % "mean";
    ss << boost::format("%-12s") % "low";
    ss << boost::format("%-12s") % "high";
    ss << boost::format("%-10s") % "minldist";
    ss << boost::format("%-10s") % "maxldist";
    ss << boost::format("%-10s") % "Levene's p";
    ss << std::endl;
    
    ss << boost::format("%-12.2e") % _pop_median;
    ss << boost::format("%-12.2e") % _pop_mean;
    ss << boost::format("%-12.2e") % _pop_min;
    ss << boost::format("%-12.2e") % _pop_max;
    ss << boost::format("%-10.3d") % _distance_min;
    ss << boost::format("%-10.3d") % _distance_max;
    ss << boost::format("%-10.3d") % levenes_pval[0];
    
    ss << std::endl;
    
    if ( _zparams.GetInt( fits_constants::PARAM_DUMP_PARAMETERS, 0) > 1 ) {
        ss << _zparams.GetAllParameters() << std::endl;
    }
    
    return ss.str();
}
