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
#define PRIOR_DISTRIB_VECTOR std::vector< std::vector<FLOAT_TYPE> >
#define PRIOR_DISTRIB_MATRIX std::vector< std::vector<FLOAT_TYPE> >

enum FactorToInfer {
    Fitness,
    MutationRate,
    PopulationSize
};


namespace fits_constants {
    
    /* GENERAL */
    const std::string current_version_str = "1.3.2";
    const FLOAT_TYPE LEVENES_SIGNIFICANCE = 0.05f;
    const std::string FILE_FIELD_DELIMITER = "\t";
    
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
    const std::string PARAM_DUMP_PARAMETERS = "dump_all_parameteres";
    
    /* Population Size */
    const std::string PARAM_POPULATION_SIZE = "N";
    const std::string PARAM_SAMPLE_SIZE = "sample_size";
    const int PARAM_SAMPLE_SIZE_DEFAULT = 0;
    
    const std::string PARAM_MIN_LOG_POPSIZE = "Nlog_min";
    const std::string PARAM_MAX_LOG_POPSIZE = "Nlog_max";
    
    //const std::string PARAM_ALT_POPULATION_SIZE = "alt_N";
    //const std::string PARAM_ALT_GENERATION = "alt_generation";
    
    const std::string PARAM_BOTTLENECK_SIZE = "bottleneck_size";
    const std::string PARAM_BOTTLENECK_INTERVAL = "bottleneck_interval"; // if size is defined, interval must be too
    const int PARAM_BOTTLENECK_NOT_DEFINED = -1;
    
    const std::string PARAM_REJECTION_THRESHOLD = "rejection_threshold";
    
    const std::string PARAM_WT_ALLELE = "wt_allele";
    
    
    // const std::string PARAM_IMEMDIATE_REJECTION = "immediate_rejection";
    
    const std::string PARAM_GENERATION_SHIFT = "generation_shift";
    
    /* Number of generations */
    const std::string PARAM_NUM_GENERATIONS = "num_generations";
    
    
    const std::string PARAM_NUM_ALLELES = "num_alleles";
    const int PARAM_NUM_ALLELES_DEFAULT = -1;
    
    const std::string PARAM_ALLELE_INIT_FREQ = "init_freq_allele";
    
    
    /* Mutation Rates */
    const std::string PARAM_MUTATION_RATE = "mutation_rate";
    const std::string PARAM_MIN_LOG_MUTATION_RATE = "min_log_mutation_rate";
    const std::string PARAM_MAX_LOG_MUTATION_RATE = "max_log_mutation_rate";
    
    const std::string PARAM_SINGLE_MUTATION_RATE = "single_mutation_rate";
    const std::string PARAM_MIN_LOG_SINGLE_MUTATION_RATE = "min_log_single_mutation_rate";
    const std::string PARAM_MAX_LOG_SINGLE_MUTATION_RATE = "max_log_single_mutation_rate";
    
    const std::string PARAM_ALLELE_FITNESS = "fitness_allele";
    const std::string PARAM_ALLELE_MAX_FITNESS = "max_fitness_allele";
    const std::string PARAM_ALLELE_MIN_FITNESS = "min_fitness_allele";
    
    //const std::string PARAM_SIM_ID = "name_of_run";
    const std::string PARAM_SIM_REPEATS = "num_samples_from_prior";
    const std::size_t PARAM_SIM_REPEATS_DEFAULT = 100000;
    
    const std::string PARAM_MANUAL_SEED = "manual_seed";
    
    const std::string PARAM_LOGISTIC_GROWTH = "logistic_growth";
    const std::string PARAM_LOGISTIC_GROWTH_K = "logistic_growth_K";
    const std::string PARAM_LOGISTIC_GROWTH_r = "logistic_growth_r";
    const FLOAT_TYPE ALLELE_FITNESS_DEFAULT_MIN = 0.0f;
    const FLOAT_TYPE ALLELE_FITNESS_DEFAULT_MAX = 2.0f;
    const FLOAT_TYPE ALLELE_FITNESS_WT = 1.0f;
    
    // const std::string PARAM_EPSILON_FLOAT_COMPARE = "epsilon_float_compare";
    
    const FLOAT_TYPE PARAM_EPSILON_SIM_DEFAULT = 0.001f; // default range within actual data to accept simulations
    
    const FLOAT_TYPE EPSILON_FLOAT_TOLERANCE = 0.0001;
    
    const int PARAM_DEFAULT_VAL_INT = -1;
    const FLOAT_TYPE PARAM_DEFAULT_VAL_FLOAT = -1.0f;
    const FLOAT_TYPE PARAM_DEFAULT_NEGATIVE_FLOAT = -1.0f;
    const std::string PARAM_DEFAULT_VAL_STRING = "NA";
    
    const std::string PARAM_ACCEPTANCE_RATE = "acceptance_rate";
    const FLOAT_TYPE ACCEPTANCE_RATE_DEFAULT = 0.01f;
    
    const std::string PARAM_ACCEPTANCE_LIMIT = "acceptance_limit";
    const std::size_t ACCEPTANCE_LIMIT_DEFAULT = 0;
    
    // const std::string PARAM_SKIP_STOCHASTIC_STEP = "skip_stochastic_step";
    // const int PARAM_SKIP_STOCHASTIC_STEP_DEFAULT = 0;
    
    // prior distributions
    const std::string PARAM_PRIOR_FILENAME = "prior_file";
    
    const std::string PARAM_FITNESS_PRIOR_DISTRIB = "fitness_prior";
    
    const std::string PARAM_PRIOR_DISTRIB_UNIFORM = "uniform";
    const std::string PARAM_PRIOR_DISTRIB_LOGNORMAL = "log_normal";
    const std::string PARAM_PRIOR_DISTRIB_LOGNORMAL_SIGMA = "log_normal_sigma";
    const FLOAT_TYPE PARAM_PRIOR_DISTRIB_LOGNORMAL_SIGMA_DEFAULT = 0.149f;
    const std::string PARAM_PRIOR_DISTRIB_LOGNORMAL_MU = "log_normal_mu";
    const FLOAT_TYPE PARAM_PRIOR_DISTRIB_LOGNORMAL_MU_DEFAULT = -0.248f;
    const std::string PARAM_PRIOR_DISTRIB_LOGNORMAL_LETHAL = "log_normal_lethal";
    const FLOAT_TYPE PARAM_PRIOR_DISTRIB_LOGNORMAL_LETHAL_DEFAULT = 0.045f;
    
    const std::string PARAM_PRIOR_DISTRIB_COMPOSITE = "fitness_composite";
    const std::string PARAM_PRIOR_DISTRIB_SMOOTHED_COMPOSITE = "smoothed_composite";
    const std::string PARAM_PRIOR_DISTRIB_FITNESS_DEFAULT = "smoothed_composite";
    
    // composite prior for fitness inference
    // if one is manually chosen - all rest must be also
    // const std::string PARAM_PRIOR_FRACTION_LTH = "prior_fraction_lth";
    // const std::string PARAM_PRIOR_FRACTION_DEL = "prior_fraction_del";
    // const std::string PARAM_PRIOR_FRACTION_NEU = "prior_fraction_neu";
    // const std::string PARAM_PRIOR_FRACTION_ADV = "prior_fraction_adv";
    // const FLOAT_TYPE PARAM_PRIOR_FRACTION_DEFAULT = -1.0f;
    
    const std::string PARAM_DISTANCE = "distance_metric";
    const std::string PARAM_DISTANCE_L1 = "L1";
    const std::string PARAM_DISTANCE_L2 = "L2";
    
    const std::string PARAM_SCALING = "scaling";
    const std::string PARAM_SCALING_MAD = "mad";
    const std::string PARAM_SCALING_SD = "sd";
    const std::string PARAM_SCALING_DEFAULT = "off";
    const std::string PARAM_SCALING_OFF = "off";
    
    const std::string PARAM_VERBOSE_SWITCH = "verbose";
    const int PARAM_VERBOSE_SWITCH_ON = 1;
    const int PARAM_VERBOSE_SWITCH_OFF = 0;
    
    // const std::string PARAM_COVERAGE_SWITCH = "coverage";
    //const int PARAM_DEFAULT_COVERAGE_SWITCH = 0;
    /* END Constants */

}
#endif /* fits_constants_h */
