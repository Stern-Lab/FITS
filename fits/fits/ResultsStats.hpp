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
#include <boost/accumulators/statistics/median.hpp>
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
    
    std::string _distance_metric;
    
    std::string GetSummaryHeader();
    
    bool _is_multi_position;
    
    PRIOR_DISTRIB _prior_distrib;
    PriorDistributionType _prior_type;
    std::vector<SimulationResult> _result_vector;
    
public:
    ResultsStats( ZParams zparams, PriorDistributionType prior_type, const PRIOR_DISTRIB &prior_distrib, const std::vector<SimulationResult>& result_vector );
    
    FLOAT_TYPE GetMedian( std::vector<FLOAT_TYPE> vec );
    int GetMedian( std::vector<int> vec );
    
    void SetRunningTimeSec(std::size_t duration) { _running_time_sec = duration; }
    std::size_t GetRunningTimeSec() { return _running_time_sec; }

    void SetResultsCount( std::size_t count );
    
    void SetMultiPosition( bool is_multi_position ) { _is_multi_position = is_multi_position; }
    bool GetMultiPosition() { return _is_multi_position; }
    
    //std::string _prior_type_str;
    
    void SetPriorType( PriorDistributionType prior_type ) { _prior_type = prior_type; }
    void SetPriorDistrib( PRIOR_DISTRIB prior_distrib ) { _prior_distrib = prior_distrib; }
    
    bool _single_mutrate_inferred;
    void SetSingleMutrateInferred( bool is_single_inferred ) { _single_mutrate_inferred=is_single_inferred; }
    
    bool _single_mutrate_used;
    void SetSingleMutrateUsed(  bool is_single_used) { _single_mutrate_used=is_single_used; }
    
    std::size_t _results_count;
    
    // These will calculate and store in internal variables
    // will later be private

    /*
    void CalculateStatsFitness(const std::vector<SimulationResult>& result_vector);
    void CalculateStatsPopulationSize(const std::vector<SimulationResult>& result_vector);
    void CalculateStatsMutation(const std::vector<SimulationResult>& result_vector);
    void CalculateStatsGenerations(const std::vector<SimulationResult>& result_vector);
     */
    
    void CalculateStatsFitness();
    void CalculateStatsPopulationSize();
    void CalculateStatsMutation();
     
    std::string GetSummaryFitness( bool table_only = false );
    std::string GetSummaryPopSize( bool table_only = false );
    std::string GetSummaryMutRate( bool table_only = false );
    
    void SetRejectionThreshold(FLOAT_TYPE new_val);
    FLOAT_TYPE GetRejectionThreshold();
    
    // Accumulators
    double _pop_min;
    double _pop_max;
    double _pop_mean;
    double _pop_sd;
    double _pop_median;
    
    FLOAT_TYPE _levenes_w;
    
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
    std::vector<unsigned int> unassigned_counter;
    
    std::vector<double> lethal_percent;
    std::vector<double> deleterious_percent;
    std::vector<double> neutral_percent;
    std::vector<double> advantageous_percent;
    std::vector<double> unassigned_percent;
    
    std::vector<FLOAT_TYPE> levenes_pval;
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
    
    // std::string GetMutrateDistrib(const std::vector<SimulationResult>& result_vector);
    // std::string GetFitnessDistrib(const std::vector<SimulationResult>& result_vector);
    // std::string GetPopsizeDistrib(const std::vector<SimulationResult>& result_vector);
    
    // void WriteFitnessDistribToFile(const std::vector<SimulationResult>& result_vector, std::string filename);
    // void WriteMutRateDistribToFile(const std::vector<SimulationResult>& result_vector, std::string filename);
    // void WritePopSizeDistribToFile(const std::vector<SimulationResult>& result_vector, std::string filename);
    
    void WritePosterior( bool is_multi_position, FactorToInfer factor, const std::vector<SimulationResult>& accepted_results_vec, const std::vector<SimulationResult>& all_results_vec, std::string filename );
    
    void WritePriorDistribToFile( FactorToInfer factor_to_infer, const PRIOR_DISTRIB& prior_distrib, std::string filename );
    
    // void WritePriorDistribToFile( const std::vector<std::vector<int>>& prior_distrib, std::string filename );
    
    std::vector<FLOAT_TYPE> GetMinFitnessVec() { return allele_min_fitness; }
    std::vector<FLOAT_TYPE> GetMaxFitnessVec() { return allele_max_fitness; }
    
    FLOAT_TYPE LevenesTest2( std::vector<FLOAT_TYPE> group1, std::vector<FLOAT_TYPE> group2 );
    
    
    
private: // here for init-order - require the public vars to be initialize
    
    std::vector<FLOAT_TYPE> DownsampleVector( const std::vector<FLOAT_TYPE> &vec, std::size_t sample_size );
    std::string GetPrintCommonHeaderStr();
    std::vector<FLOAT_TYPE> _allele_Nu; // for allele i u=f(wt->i)
    ZParams _zparams;
};


#endif /* ReportStats_hpp */
