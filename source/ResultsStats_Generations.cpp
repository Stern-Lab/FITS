//
//  ResultsStats_Generations.cpp
//  fits
//
//  Created by Tal Zinger on 2-6-2017.
//  Copyright © 2016 Stern Lab. All rights reserved.
//

#include "ResultsStats.hpp"


void ResultsStats::CalculateStatsGenerations(const std::vector<SimulationResult>& result_vector)
{
    std::cout << "CalculateStatsGenerations starting..." << std::endl;
    
    _num_results = result_vector.size();
    _num_timepoints = result_vector[0].actual_generations.size();
    
    //std::cout << "checkpoint 1" << std::endl;
    
    boost::accumulators::accumulator_set<
    FLOAT_TYPE,
    boost::accumulators::stats<
    // boost::accumulators::tag::median,
    boost::accumulators::tag::variance,
    boost::accumulators::tag::mean,
    boost::accumulators::tag::min,
    boost::accumulators::tag::max> > acc_distance;
    
    boost::accumulators::accumulator_set<
    int,
    boost::accumulators::stats<
    // boost::accumulators::tag::median,
    boost::accumulators::tag::variance,
    boost::accumulators::tag::mean,
    boost::accumulators::tag::min,
    boost::accumulators::tag::max> > acc_generation_interval;
    
    std::vector< boost::accumulators::accumulator_set<
    int,
    boost::accumulators::stats<
    boost::accumulators::tag::variance,
    boost::accumulators::tag::mean,
    boost::accumulators::tag::min,
    boost::accumulators::tag::max> > > acc_vec_generations(_num_timepoints);

    _num_generations = result_vector[0].actual_generations.size();
    _first_generation = result_vector[0].generation_shift;
    
    
    std::vector<int> generation_interval_vec;
    for ( auto sim_result : result_vector ) {
        acc_distance(sim_result.distance_from_actual);
        acc_generation_interval(sim_result.generation_interval);
        
        //std::cout << "interval " << sim_result.generation_interval << std::endl;
        generation_interval_vec.push_back(sim_result.generation_interval);
        //for ( auto time_point=0; time_point<_num_timepoints; ++time_point ) {
        //    acc_vec_generations[time_point](sim_result.actual_generations[time_point]);
        //}
    }
    
    
    //_inferred_generations.resize( acc_vec_generations.size());
    //for ( auto i=0; i<acc_vec_generations.size(); ++i ) {
        // MEDIAN
        //_inferred_generations[i] = boost::accumulators::median(acc_vec_generations[i]);
    //}
    _gen_interval_min = boost::accumulators::min(acc_generation_interval);
    _gen_interval_max = boost::accumulators::max(acc_generation_interval);
    _gen_interval_mean = boost::accumulators::mean(acc_generation_interval);
    _gen_interval_sd = std::sqrt(boost::accumulators::variance(acc_generation_interval));
    // MEDIAN
    _gen_interval_median = GetMedian(generation_interval_vec);
    
    _distance_min = boost::accumulators::min(acc_distance);
    _distance_max = boost::accumulators::max(acc_distance);
    _distance_mean = boost::accumulators::mean(acc_distance);
    _distance_sd = std::sqrt(boost::accumulators::variance(acc_distance));
    
    
    std::cout << "CalculateStatsGenerations ended." << std::endl;
}


std::string ResultsStats::GetSummaryGenerations()
{
    std::cout << "GetSummaryGenerations starting..." << std::endl;
    std::stringstream ss;
    
    ss.imbue(std::locale(fits_constants::used_locale));
    
    ss << "Generation Interval Report" << std::endl;
    
    ss << GetSummaryHeader();
    //ss << "FITS v"<< fits_constants::current_version_str << std::endl;
    
    //auto current_time_raw = std::chrono::system_clock::now();
    //auto current_time = std::chrono::system_clock::to_time_t(current_time_raw);
    //auto current_time_final = *std::localtime(&current_time);
    //ss << std::put_time(&current_time_final, "%F %T") << std::endl;
    
    //std::cout << "Simulation results used for calculations: " << _num_results << std::endl;
    
    //ss << "=======================" << std::endl;
    
    //if ( _rejection_threshold > 0.0f ) {
    //    ss << "Rejection threshold set to " << boost::format("%-10d") % _rejection_threshold << std::endl;
   // }
    
    if ( _single_mutrate_used ) {
        ss << "Used a single mutation rate." << std::endl;
    }
    
    
    ss << boost::format("%-12s") % "median";
    ss << boost::format("%-12s") % "mean";
    ss << boost::format("%-12s") % "low";
    ss << boost::format("%-12s") % "high";
    ss << boost::format("%-10s") % "minldist";
    ss << boost::format("%-10s") % "maxldist";
    ss << std::endl;
    
    ss << boost::format("%-12.2e") % _gen_interval_median;
    ss << boost::format("%-12.2e") % _gen_interval_mean;
    ss << boost::format("%-12.2e") % _gen_interval_min;
    ss << boost::format("%-12.2e") % _gen_interval_max;
    ss << boost::format("%-10.2e") % _distance_min;
    ss << boost::format("%-10.2e") % _distance_max;
    ss << std::endl;
    
    //ss << "Inferred generations: ";
    //ss << _inferred_generations[0];
    //for ( auto idx=1; idx<_num_timepoints; ++idx ) {
    //    ss << "," << _inferred_generations[idx];
   // }
    //ss << std::endl;
    
    std::cout << "GetSummaryGenerations ended." << std::endl;
    
    return ss.str();
}
