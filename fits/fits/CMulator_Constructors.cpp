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

/****************/
/* Constructors */
/****************/

#include "CMulator.h"
using namespace fits_constants;

// Default ctor - marked as uninitialized
CMulator::CMulator() :
_all_simulated_data(),
_observed_simulated_data(),
_initialized_with_parameters(false),
_debug_mode(false),
_acceptance_rate(ACCEPTANCE_RATE_DEFAULT),
_use_observed_data(false),
_use_single_mutation_rate(false)
{
    _time_for_seeding = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    _boost_gen.seed(_time_for_seeding);
}


// Constructor using parameters object
CMulator::CMulator( const ZParams zparams ) :
_initialized_with_parameters(false),
_all_simulated_data(),
_observed_simulated_data(),
_debug_mode(false),
_use_observed_data(false),
_use_single_mutation_rate(false)
{
    InitMemberVariables(zparams);
    
    _time_for_seeding = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    _boost_gen.seed(_time_for_seeding);
}


// Copy Constructor
CMulator::CMulator( const CMulator &original ) :
_debug_mode(original._debug_mode),
_N(original._N),
_N0(original._N0),
_num_generations(original._num_generations),
_num_alleles(original._num_alleles),
_bottleneck_interval(original._bottleneck_interval),
_bottleneck_size(original._bottleneck_size),
_generation_shift(original._generation_shift),
_logistic_growth(original._logistic_growth),
_logistic_growth_K(original._logistic_growth_K),
_logistic_growth_t(original._logistic_growth_t),
_logistic_growth_r(original._logistic_growth_r),
_repeats(original._repeats),
_allele_fitness(original._allele_fitness),
_allele_min_fitness(original._allele_min_fitness),
_wt_allele_index(original._wt_allele_index),
_current_generation(original._current_generation),
_allele_max_fitness(original._allele_max_fitness),
_initialized_with_parameters(original._initialized_with_parameters),
//_epsilon_float_compare(original._epsilon_float_compare),
_allele_init_freqs(original._allele_init_freqs),
_acceptance_rate(original._acceptance_rate),
_use_observed_data(original._use_observed_data),
_use_single_mutation_rate(original._use_single_mutation_rate),
_all_simulated_data(original._all_simulated_data),
_observed_simulated_data(original._observed_simulated_data),
_mutation_rates_matrix(original._mutation_rates_matrix)
//_name_of_run(original._name_of_run)
{
    // do not copy seed so random numbers won't repeat
    _time_for_seeding = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    _boost_gen.seed(_time_for_seeding);
}

