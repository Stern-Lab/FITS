//
//  ResultsStats.hpp
//  fits
//
//  Created by Tal Zinger on 16/11/2016.
//  Copyright © 2016 Stern Lab. All rights reserved.
//

#ifndef ReportStats_hpp
#define ReportStats_hpp

#include <vector>
#include <string>
#include <sstream>
#include <locale>
#include <cmath>

// for time reporting
#include <chrono>
#include <iomanip>

// for Levene's test
#include <boost/math/distributions.hpp>

#include <boost/math/distributions/fisher_f.hpp>

#include <boost/accumulators/numeric/functional.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/framework/features.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/format.hpp>
#include <boost/locale.hpp>

#include "PriorSampler.hpp"
#include "CMulator.h"
#include "SimulationResult.hpp"
#include "ZParams.h"
#include "fits_constants.h"

enum AlleleCategory {
    // cannot be assigned with fitness value
    WT,
    Undefined,
    
    // pass relaxed threshold
    Possible_lethal,
    Possible_deleterious,
    Possible_neutral,
    Possible_advantageous,
    
    // pass strict threshold
    Lethal,
    Deleterious,
    Neutral,
    Adventageous
};


class ResultsStats {
    
    const FLOAT_TYPE FITNESS_LETHAL_MIN = 0.0;
    const FLOAT_TYPE FITNESS_LETHAL_MAX = 0.5;
    const FLOAT_TYPE FITNESS_DELETERIOUS_MIN = 0.5;
    const FLOAT_TYPE FITNESS_DELETERIOUS_MAX = 0.95;
    const FLOAT_TYPE FITNESS_NEUTRAL_MIN = 0.95;
    const FLOAT_TYPE FITNESS_NEUTRAL_MAX = 1.05;
    const FLOAT_TYPE FITNESS_ADVANTAGEOUS_MIN = 1.05;
    const FLOAT_TYPE FITNESS_ADVANTAGEOUS_MAX = 2.0;
    const FLOAT_TYPE FITNESS_NEUTRAL = 1.0f;
    
    const int CATEGORY_INCLUSION_RELAXED_THRESHOLD = 50;
    const int CATEGORY_INCLUSION_STRICT_THRESHOLD = 95;

    std::size_t _running_time_sec;
    
    std::string GetSummaryHeader();
    
public:
    ResultsStats(ZParams zparams);
    
    FLOAT_TYPE GetMedian( std::vector<FLOAT_TYPE> vec );
    int GetMedian( std::vector<int> vec );
    
    void SetRunningTimeSec(std::size_t duration) { _running_time_sec = duration; }
    std::size_t GetRunningTimeSec() { return _running_time_sec; }

    void SetResultsCount( std::size_t count );
    
    std::string _prior_distrib;
    void SetPriorDistrib( std::string prior ) { _prior_distrib=prior; }
    
    bool _single_mutrate_inferred;
    void SetSingleMutrateInferred( bool is_single_inferred ) { _single_mutrate_inferred=is_single_inferred; }
    
    bool _single_mutrate_used;
    void SetSingleMutrateUsed(  bool is_single_used) { _single_mutrate_used=is_single_used; }
    
    std::size_t _results_count;
    
    // These will calculate and store in internal variables
    // will later be private

    void CalculateStatsFitness(const std::vector<SimulationResult>& result_vector);
    void CalculateStatsPopulationSize(const std::vector<SimulationResult>& result_vector);
    void CalculateStatsMutation(const std::vector<SimulationResult>& result_vector);
    void CalculateStatsGenerations(const std::vector<SimulationResult>& result_vector);
     
    std::string GetSummaryFitness();
    std::string GetSummaryPopSize();
    std::string GetSummaryMutRate();
    std::string GetSummaryGenerations();
    
    void SetRejectionThreshold(FLOAT_TYPE new_val);
    FLOAT_TYPE GetRejectionThreshold();
    
    // Accumulators
    double _pop_min;
    double _pop_max;
    double _pop_mean;
    double _pop_sd;
    double _pop_median;
    

    float _levenes_w;
    
    FLOAT_TYPE _distance_min;
    FLOAT_TYPE _distance_max;
    FLOAT_TYPE _distance_mean;
    FLOAT_TYPE _distance_sd;
    
    int _gen_interval_min;
    int _gen_interval_max;
    int _gen_interval_mean;
    int _gen_interval_sd;
    int _gen_interval_median;
    int _gen_interval_mad;
    std::size_t _first_generation;
    std::size_t _num_generations;
    std::vector<int> _inferred_generations;
    int _num_timepoints;
    
    std::string AlleleCategory2String( AlleleCategory category );
    
    /* general data */
    std::size_t _num_results;
    std::size_t _num_alleles;
    
    /* Different strategies for ABC rejection */
    FLOAT_TYPE _rejection_threshold;
    

    /* fitness inference */
    std::vector<FLOAT_TYPE> allele_mean_fitness;
    std::vector<FLOAT_TYPE> allele_sd_fitness;
    std::vector<FLOAT_TYPE> allele_median_fitness;
    
    std::vector<FLOAT_TYPE> allele_min_fitness;
    std::vector<FLOAT_TYPE> allele_max_fitness;
    
    std::vector<FLOAT_TYPE> allele_min_95percentile_fitness;
    std::vector<FLOAT_TYPE> allele_max_95percentile_fitness;
    
    // measure of disperssion
    std::vector<double> allele_MAD;
    
    
    // P(w<1|data)
    std::vector<FLOAT_TYPE> allele_pval;
    
    std::vector<unsigned int> lethal_counter;
    std::vector<unsigned int> deleterious_counter;
    std::vector<unsigned int> neutral_counter;
    std::vector<unsigned int> advantageous_counter;
    
    std::vector<unsigned int> lethal_percent;
    std::vector<unsigned int> deleterious_percent;
    std::vector<unsigned int> neutral_percent;
    std::vector<unsigned int> advantageous_percent;
    
    std::vector<float> levenes_pval;
    boost::numeric::ublas::matrix<FLOAT_TYPE> levenes_pval_matrix;
    
    boost::numeric::ublas::matrix< std::vector<FLOAT_TYPE> > prior_matrix;
    
    std::vector<AlleleCategory> allele_category;
    
    /* mutation rate inference */
    boost::numeric::ublas::matrix<FLOAT_TYPE> min_mutation_rates;
    boost::numeric::ublas::matrix<FLOAT_TYPE> max_mutation_rates;
    boost::numeric::ublas::matrix<FLOAT_TYPE> mean_mutation_rates;
    boost::numeric::ublas::matrix<FLOAT_TYPE> median_mutation_rates;
    boost::numeric::ublas::matrix<FLOAT_TYPE> normalized_median_mutation_rates;
    
    /* population size inference */
    int max_population_size;
    int min_population_size;
    int mean_population_size;
    int median_population_size;
    int normalized_median_population_size;
    
    void WriteStringToFile( std::string filename, std::string str );
    
    std::string GetMutrateDistrib(const std::vector<SimulationResult>& result_vector);
    std::string GetFitnessDistrib(const std::vector<SimulationResult>& result_vector);
    std::string GetPopsizeDistrib(const std::vector<SimulationResult>& result_vector);
    
    void WriteFitnessDistribToFile(const std::vector<SimulationResult>& result_vector, std::string filename);
    void WriteMutRateDistribToFile(const std::vector<SimulationResult>& result_vector, std::string filename);
    void WritePopSizeDistribToFile(const std::vector<SimulationResult>& result_vector, std::string filename);
    void WriteGenerationsDistribToFile(const std::vector<SimulationResult>& result_vector, std::string filename);
    
    void WritePriorDistribToFile( const std::vector<std::vector<FLOAT_TYPE>>& prior_distrib, std::string filename );
    
    void WritePriorDistribToFile( const std::vector<std::vector<int>>& prior_distrib, std::string filename );
    
    std::vector<FLOAT_TYPE> GetMinFitnessVec() { return allele_min_fitness; }
    std::vector<FLOAT_TYPE> GetMaxFitnessVec() { return allele_max_fitness; }
    
    float LevenesTest2( std::vector<float> group1, std::vector<float> group2 );
    
private: // here for init-order - require the public vars to be initialize
    
    
    
    std::string GetPrintCommonHeaderStr();
    std::vector<FLOAT_TYPE> _allele_Nu; // for allele i u=f(wt->i)
    ZParams _zparams;
};


#endif /* ReportStats_hpp */
