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

/* Many Constants used throughout FITS */

#ifndef fits_constants_h
#define fits_constants_h

#include <string>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/odeint/util/ublas_wrapper.hpp>

// 20170611 because some diatcnes turn 0 trying to go to double
//#define FLOAT_TYPE float
#define FLOAT_TYPE double
#define MATRIX_TYPE boost::numeric::ublas::matrix<FLOAT_TYPE>
#define INT_MATRIX boost::numeric::ublas::matrix<int>

namespace fits_constants {
    
    
    /* GENERAL */
    const std::string current_version_str = "1.1";
    const FLOAT_TYPE LEVENES_SIGNIFICANCE = 0.05;
    
    /* COMMAND LINE ARGUMANTS */
    const std::string ARG_INFER_FITNESS = "-fitness";
    const std::string ARG_INFER_MUTATION = "-mutation";
    const std::string ARG_INFER_POPSIZE = "-popsize";
    
    
    /* PRIOR DISTRIBUTION
       define fitness category index, to be used with the discrete distribution */
    const int FITNESS_LTH_IDX = 0;
    const int FITNESS_DEL_IDX = 1;
    const int FITNESSS_NEU_IDX = 2;
    const int FITNESS_ADV_IDX = 3;
    

    /* PARAMETERS */
    
    /* General */
    const std::string PARAM_INFERENCE_ARG = "infer"; // will be used instead of command line args
    const std::string PARAM_INFERENCE_ARG_FITNESS = "fitness";
    const std::string PARAM_INFERENCE_ARG_MUTATION = "mutation";
    const std::string PARAM_INFERENCE_ARG_POPSIZE = "popsize";
    const std::string PARAM_INFERENCE_ARG_EMPTY = "";
    
    const std::string PARAM_DUMP_PARAMETERS = "_dump_all_parameteres";
    
    /* Population Size */
    const std::string PARAM_POPULATION_SIZE = "_N";
    const std::string PARAM_SAMPLE_SIZE = "_sample_size";
    const int PARAM_SAMPLE_SIZE_DEFAULT = 0;
    
    const std::string PARAM_MIN_LOG_POPSIZE = "_Nlog_min";
    const std::string PARAM_MAX_LOG_POPSIZE = "_Nlog_max";
    
    const std::string PARAM_ALT_POPULATION_SIZE = "_alt_N";
    const std::string PARAM_ALT_GENERATION = "_alt_generation";
    
    const std::string PARAM_BOTTLENECK_SIZE = "_bottleneck_size";
    const std::string PARAM_BOTTLENECK_INTERVAL = "_bottleneck_interval";
    
    const std::string PARAM_REJECTION_THRESHOLD = "_rejection_threshold";
    
    
    /* Number of generations */
    const std::string PARAM_NUM_GENERATIONS = "_num_generations";
    
    const std::string PARAM_MIN_FIXED_INT_GENERATION = "_min_fixed_generation_interval";
    const std::string PARAM_MAX_FIXED_INT_GENERATION = "_max_fixed_generation_interval";
    const std::string PARAM_MIN_INT_GENERATION = "_min_generation_interval";
    const std::string PARAM_MAX_INT_GENERATION = "_max_generation_interval";
    
    const std::string PARAM_NUM_ALLELES = "_num_alleles";
    
    const std::string PARAM_ALLELE_INIT_FREQ = "_allele_init_freq";
    //const std::string PARAM_ALLELE_NO_INIT_FREQS = "_no_init_freqs_as_parameters"; // assume it will be provided later
    
    
    /* Mutation Rates */
    const std::string PARAM_MUTATION_RATE = "_mutation_rate";
    const std::string PARAM_MIN_LOG_MUTATION_RATE = "_min_log_mutation_rate";
    const std::string PARAM_MAX_LOG_MUTATION_RATE = "_max_log_mutation_rate";
    
    const std::string PARAM_SINGLE_MUTATION_RATE = "_single_mutation_rate";
    const std::string PARAM_MIN_LOG_SINGLE_MUTATION_RATE = "_min_log_single_mutation_rate";
    const std::string PARAM_MAX_LOG_SINGLE_MUTATION_RATE = "_max_log_single_mutation_rate";
    
    const std::string PARAM_ALLELE_FITNESS = "_allele_fitness";
    const std::string PARAM_ALLELE_MAX_FITNESS = "_allele_max_fitness";
    const std::string PARAM_ALLELE_MIN_FITNESS = "_allele_min_fitness";
    
    const std::string PARAM_SIM_ID = "_name_of_run";
    const std::string PARAM_SIM_REPEATS = "_num_repeats";
    
    const std::string PARAM_MANUAL_SEED = "_manual_seed";
    
    const std::string PARAM_LOGISTIC_GROWTH = "_logistic_growth";
    const std::string PARAM_LOGISTIC_GROWTH_K = "_logistic_growth_K";
    const std::string PARAM_LOGISTIC_GROWTH_r = "_logistic_growth_r";
    const FLOAT_TYPE ALLELE_FITNESS_DEFAULT_MIN = 0.0;
    const FLOAT_TYPE ALLELE_FITNESS_DEFAULT_MAX = 2.0;
    
    const std::string PARAM_EPSILON_FLOAT_COMPARE = "_epsilon_float_compare";
    
    const FLOAT_TYPE PARAM_EPSILON_SIM_DEFAULT = 0.001f; // default range within actual data to accept simulations
    
    const int PARAM_DEFAULT_VAL_INT = -1;
    const FLOAT_TYPE PARAM_DEFAULT_VAL_FLOAT = -1.0f;
    const FLOAT_TYPE PARAM_DEFAULT_NEGATIVE_FLOAT = -1.0f;
    const std::string PARAM_DEFAULT_VAL_STRING = "NA";
    
    const std::string PARAM_ACCEPTANCE_RATE = "_acceptance_rate";
    const FLOAT_TYPE ACCEPTANCE_RATE_DEFAULT = 0.01f;
    
    const std::string PARAM_ACCEPTANCE_LIMIT = "_acceptance_limit";
    const std::size_t ACCEPTANCE_LIMIT_DEFAULT = 0;
    
    //const std::string PARAM_SKIP_STOCHASTIC_STEP = "_skip_stochastic_step";
    //const int PARAM_SKIP_STOCHASTIC_STEP_DEFAULT = 0;
    
    // prior distributions
    const std::string PARAM_PRIOR_DISTRIB = "_prior_distribution";
    const std::string PARAM_PRIOR_DISTRIB_UNIFORM = "uniform";
    const std::string PARAM_PRIOR_DISTRIB_LOGNORMAL = "log_normal";
    const std::string PARAM_PRIOR_DISTRIB_COMPOSITE = "fitness_composite";
    const std::string PARAM_PRIOR_DISTRIB_SMOOTHED_COMPOSITE = "smoothed_composite";
    const std::string PARAM_PRIOR_DISTRIB_DEFAULT = "uniform";
    
    // composite prior for fitness inference
    // if one is manually chosen - all rest must be also
    const std::string PARAM_PRIOR_FRACTION_LTH = "_prior_fraction_lth";
    const std::string PARAM_PRIOR_FRACTION_DEL = "_prior_fraction_del";
    const std::string PARAM_PRIOR_FRACTION_NEU = "_prior_fraction_neu";
    const std::string PARAM_PRIOR_FRACTION_ADV = "_prior_fraction_adv";
    const FLOAT_TYPE PARAM_PRIOR_FRACTION_DEFAULT = -1.0f;
    
    const std::string PARAM_DISTANCE = "distance_metric";
    const std::string PARAM_DISTANCE_L1 = "L1";
    const std::string PARAM_DISTANCE_L2 = "L2";
    
    const std::string PARAM_SCALING = "scaling";
    const std::string PARAM_SCALING_MAD = "mad";
    const std::string PARAM_SCALING_SD = "sd";
    const std::string PARAM_SCALING_DEFAULT = "off";
    const std::string PARAM_SCALING_OFF = "off";
    
    const std::string PARAM_COVERAGE_SWITCH = "_coverage";
    const int PARAM_DEFAULT_COVERAGE_SWITCH = 0;
    /* END Constants */

}
#endif /* fits_constants_h */
