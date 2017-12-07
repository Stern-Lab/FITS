//
//  ReportStats.cpp
//  fits
//
//  Created by Tal Zinger on 16/11/2016.
//  Copyright © 2016 Stern Lab. All rights reserved.
//

#include "ResultsStats.hpp"

ResultsStats::ResultsStats(ZParams zparams)
:
_running_time_sec(0),
_num_results(0),
_num_alleles(0),
_num_generations(0),
allele_mean_fitness(_num_alleles),
allele_sd_fitness(_num_alleles),
allele_median_fitness(_num_alleles),
allele_min_fitness(_num_alleles, -1.0),
allele_max_fitness(_num_alleles, -1.0),
allele_min_95percentile_fitness(_num_alleles, -1.0),
allele_max_95percentile_fitness(_num_alleles, -1.0),
allele_pval(_num_alleles, 0),
lethal_counter(_num_alleles, 0),
deleterious_counter(_num_alleles, 0),
neutral_counter(_num_alleles, 0),
advantageous_counter(_num_alleles, 0),
lethal_percent(_num_alleles, -1),
deleterious_percent(_num_alleles, -1),
neutral_percent(_num_alleles, -1),
advantageous_percent(_num_alleles, -1),
allele_category(_num_alleles, AlleleCategory::Undefined),
_rejection_threshold(0.0f),
_results_count(0),
_allele_Nu(_num_alleles, 0.0f),
_single_mutrate_inferred(false),
_single_mutrate_used(false),
_first_generation(0),
_zparams(zparams),
levenes_pval(_num_alleles, -1.0f)
{
    _prior_distrib = "";
}

void ResultsStats::SetRejectionThreshold(FLOAT_TYPE new_val)
{
    _rejection_threshold = new_val;
}

FLOAT_TYPE ResultsStats::GetRejectionThreshold()
{
    return _rejection_threshold;
}

void ResultsStats::SetResultsCount( std::size_t count )
{
    _results_count = count;
}

std::string ResultsStats::GetPrintCommonHeaderStr()
{
    std::string tmp_str = "";
    
    return tmp_str;
}

void ResultsStats::WritePriorDistribToFile( const std::vector<std::vector<FLOAT_TYPE>>& prior_distrib, std::string filename )
{
    std::ofstream outfile(filename, std::ofstream::out | std::ofstream::trunc);
    
    if (!outfile.is_open()) {
        std::cerr << "unable to open file for writing: " << filename << std::endl;
        throw "unable to open file for writing: " + filename;
    }
    

    // 2016-09-04 adding sanity check to report problem
    if ( prior_distrib.empty() ) {
        std::cerr << "Error in writing prior: result_vector is empty. Filename: " << filename << std::endl;
        std::string my_err = "";
        my_err = "Error in writing prior: result_vector is empty. Filename: ";
        my_err += filename;
        throw my_err.c_str();
    }
 
    
    for ( auto current_vec : prior_distrib ) {
        
        for ( auto current_val : current_vec ) {
            
            outfile << current_val << "\t";
        }
        outfile << std::endl;
    }
    
    outfile.close();
}

void ResultsStats::WritePriorDistribToFile( const std::vector<std::vector<int>>& prior_distrib, std::string filename )
{
    std::ofstream outfile(filename, std::ofstream::out | std::ofstream::trunc);
    
    if (!outfile.is_open()) {
        std::cerr << "unable to open file for writing: " << filename << std::endl;
        throw "unable to open file for writing: " + filename;
    }
    
    
    // 2016-09-04 adding sanity check to report problem
    if ( prior_distrib.empty() ) {
        std::cerr << "Error in writing prior: result_vector is empty. Filename: " << filename << std::endl;
        std::string my_err = "";
        my_err = "Error in writing prior: result_vector is empty. Filename: ";
        my_err += filename;
        throw my_err.c_str();
    }
    
    
    for ( auto current_vec : prior_distrib ) {
        
        for ( auto current_val : current_vec ) {
            
            outfile << current_val << "\t";
        }
        outfile << std::endl;
    }
    
    outfile.close();
}

// Leven'es test for two groups
// Used to test whether the prior is different than the posterior
// Code is based on the Wikipedia article:
// https://en.wikipedia.org/wiki/Levene's_test
//
// returns the pval of the comparison
float ResultsStats::LevenesTest2( std::vector<float> group1, std::vector<float> group2 )
{
    //using namespace boost::accumulators;
    
    float Ytilde1 = GetMedian(group1);
    float Ytilde2 = GetMedian(group2);
    
    std::vector<float> Z1j(group1.size());
    std::vector<float> Z2j(group2.size());
    
    // std::cout << "median group1=" << Ytilde1 << " group2=" << Ytilde2 << std::endl;
    
    // group 1
    // std::cout << "group 1: " << std::endl;
    for ( auto j=0; j<group1.size(); ++j ) {
        Z1j[j] = std::fabs(group1[j] - Ytilde1);
        // std::cout << "\t|" << group1[j] << " - " << Ytilde1 << "| = " << Z1j[j] << std::endl;
    }
    //std::cout << std::endl;
    
    // group 2
    // std::cout << "group 2: " << std::endl;
    for ( auto j=0; j<group2.size(); ++j ) {
        Z2j[j] = std::fabs(group2[j] - Ytilde2);
        // std::cout << "\t|" << group2[j] << " - " << Ytilde2 << "| = " << Z2j[j] << std::endl;
    }
    // std::cout << std::endl;
    
    boost::accumulators::accumulator_set<float,
    boost::accumulators::stats<
    boost::accumulators::tag::mean>> accZ0, accZ1, accZ2;
    
    // std::cout << "z1j: ";
    for ( auto val : Z1j ) {
        accZ0(val);
        accZ1(val);
        // std::cout << val << " ";
    }
    // std::cout << std::endl;
    
    // std::cout << "z2j: ";
    for ( auto val : Z2j ) {
        accZ0(val);
        accZ2(val);
        // std::cout << val << " ";
    }
    // std::cout << std::endl;
    
    // this should always be mean, not median
    auto Z0 = boost::accumulators::mean(accZ0);
    auto Z1 = boost::accumulators::mean(accZ1);
    auto Z2 = boost::accumulators::mean(accZ2);
    
    auto N1 = static_cast<float>( group1.size() );
    auto N2 = static_cast<float>( group2.size() );
    auto N = N1 + N2;
    auto k = 2.0f; // prior & posterior
    
    //auto alpha = 0.05;
    auto alpha = 0.01f;
    
    //auto numerator = static_cast<float>( (N-k)*( N1*(Z1-Z0)*(Z1-Z0) + N2*(Z2-Z0)*(Z2-Z0) ) );
    auto numerator = static_cast<float>( (N-k)*( N1*std::pow(Z1-Z0, 2) + N2*std::pow(Z2-Z0, 2) ) );
    
    // sum for each group
    float sum_diff_group1 = 0.0f;
    for (auto z1j_val : Z1j) {
        sum_diff_group1 += std::pow( (z1j_val-Z1), 2.0f );
    }
    
    float sum_diff_group2 = 0.0f;
    for (auto z2j_val : Z2j) {
        sum_diff_group2 += std::pow( (z2j_val-Z2), 2.0f );
    }
    
    // final result
    auto denominator = static_cast<float>(k-1.0f)*( sum_diff_group1 + sum_diff_group2 );
    
    auto W = numerator / denominator;
    
    // quantile of the F distribution is F( alpha=0.05, k-1 degrees, N-k degrees )
    auto numerator_df = k-1;
    auto denominator_df = N-k;
    
    auto f_dist = boost::math::fisher_f( numerator_df, denominator_df );
    
    // std::cout << "W=" << W << std::endl;
    
    auto pval = boost::math::cdf( boost::math::complement(f_dist, W ));
    
    return pval;
}


std::string ResultsStats::GetSummaryHeader()
{
     std::stringstream ss;
    
    ss << "FITS v"<< fits_constants::current_version_str << std::endl;
    
    auto current_time_raw = std::chrono::system_clock::now();
    auto current_time = std::chrono::system_clock::to_time_t(current_time_raw);
    auto current_time_final = *std::localtime(&current_time);
    ss << std::put_time(&current_time_final, "%F %T") << std::endl;
    
    ss << "Simulation results used for calculations: " << _num_results << std::endl;
    
    if ( _running_time_sec > 0 ) {
        ss << "Running time " << _running_time_sec << " seconds" << std::endl;
    }
    
    ss << "=====================" << std::endl;
    
    if ( _rejection_threshold > 0.0f ) {
        ss << "Rejection threshold set to " << boost::format("%-10d") % _rejection_threshold << std::endl;
    }
    
    return ss.str();
}

FLOAT_TYPE ResultsStats::GetMedian( std::vector<FLOAT_TYPE> vec )
{
    auto median_idx = std::floor(vec.size() / 2);
    
    std::nth_element( vec.begin(), vec.begin()+median_idx, vec.end() );
    
    if ( vec.size() % 2 == 0 ) {
        return ( vec[median_idx] + vec[median_idx-1] ) / 2;
    }
    
    return vec[median_idx];
}

int ResultsStats::GetMedian( std::vector<int> vec )
{
    auto median_idx = std::floor(vec.size() / 2);
    
    
    //std::cout << " median idx=" << median_idx << std::endl;
    
    std::nth_element( vec.begin(), vec.begin()+median_idx, vec.end() );
    
    if ( vec.size() % 2 == 0 ) {
        //std::cout << " median value=" << ( vec[median_idx] + vec[median_idx-1] ) / 2 << std::endl;
        return ( vec[median_idx] + vec[median_idx-1] ) / 2;
    }
    
    //std::cout << " median value=" << vec[median_idx] << std::endl;
    
    return vec[median_idx];
}
