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


#include "CMulator.h"
using namespace fits_constants;


std::vector<FLOAT_TYPE> CMulator::GetAlleleFitnessValues() const
{
    return _allele_fitness;
}


int CMulator::GetNumOfGenerations() const
{
    return _num_generations;
}


void CMulator::SetNumOfGeneration(int generations)
{
    if ( generations == _num_generations ) {
        return;
    }
    
    if ( !IsValid_Generation(generations) ) {
        std::string tmp_str = "SetNumOfGeneration: invalid value (" + std::to_string(generations) + ").";
        throw tmp_str;
    }
    
    _num_generations = generations;
    
    Reset_Soft();
}

/*
std::string CMulator::GetSimUID() const
{
    return _name_of_run;
}
*/

/*
void CMulator::SetSimUID( std::string new_sim_uid )
{
    _name_of_run = new_sim_uid;
}
*/

/*
void CMulator::SetBottleneckSize( int size )
{
    if ( !IsValid_BottleneckSize(size) ) {
        throw " SetBottleneckSize: invalid value (" + std::to_string(size) + ").";
    }
    
    _bottleneck_size = size;
}
*/

int CMulator::GetBottleneckSize() const
{
    return _bottleneck_size;
}

/*
void CMulator::SetAlleleFitnessValue( int allele, FLOAT_TYPE fitness )
{
    if ( !IsValid_Fitness(fitness) ) {
        throw " SetAlleleFitnessValue: invalid fitness value (" + std::to_string(fitness) + ").";
    }
    
    if ( !IsValid_Allele(allele) ) {
        throw " SetAlleleFitnessValue: invalid allele (" + std::to_string(allele) + ").";
    }
    
    _allele_fitness[allele] = fitness;
}
*/

void CMulator::SetAlleleInitFreqs( std::vector<FLOAT_TYPE> freqs )
{
    float current_freq_sum = 0;
    
    for ( auto tmp_idx=0; tmp_idx<freqs.size(); ++tmp_idx ) {
        if ( freqs[tmp_idx] < 0 ) {
            
            std::string tmp_str = "Negative frequency value (" + std::to_string( freqs[tmp_idx] ) +
            " for allele " + std::to_string(tmp_idx);
            
            throw tmp_str;
        }
        
        current_freq_sum += freqs[tmp_idx];
    }
    
    if ( current_freq_sum != 1.0f ) {
        std::string tmp_str = "Setting initial frequencies for alleles, but they don't sum up to 1 (" +
        std::to_string(current_freq_sum) + ")";
        
        throw tmp_str;
    }
    
    _allele_init_freqs.resize( freqs.size() );
    _allele_init_freqs = freqs;
    
    // for consistency reasons, we must reset the simulated data
    // this will also update the init freqs
    Reset_Soft();
}

std::vector<FLOAT_TYPE> CMulator::GetAlleleInitFreqs() const
{
    return _allele_init_freqs;
}


std::vector<FLOAT_TYPE> CMulator::GetAlleleMaxFitnessValues() const
{
    return _allele_max_fitness;
}


std::vector<FLOAT_TYPE> CMulator::GetAlleleMinFitnessValues() const
{
    return _allele_min_fitness;
}


FLOAT_TYPE CMulator::GetAlleleFreq( int generation, int allele ) const
{
    if ( !IsValid_Generation(generation) ) {
        std::string tmp_str = "GetAlleleFreq: illegal generation " + std::to_string(generation);
        throw tmp_str;
    }
    
    if ( !IsValid_Allele(allele) ) {
        std::string tmp_str = "GetAlleleFreq: illegal allele " + std::to_string(allele);
        throw tmp_str;
    }
    
    if (_use_observed_data) {
        return _observed_simulated_data(generation, allele);
    }
    
    return _all_simulated_data(generation, allele);
}


std::vector<FLOAT_TYPE> CMulator::GetRawFrequencyData()
{
    std::vector<FLOAT_TYPE> tmp_freqs;
    
    for ( auto cur_generation=0; cur_generation<_num_generations; ++cur_generation ) {
        
        for (int cur_allele=0; cur_allele<_num_alleles; ++cur_allele) {
            
            if (_use_observed_data) {
                tmp_freqs.push_back( _observed_simulated_data(cur_generation, cur_allele) );
            }
            else {
                tmp_freqs.push_back( _all_simulated_data(cur_generation, cur_allele) );
            }
        }
    }
    
    return tmp_freqs;
}


std::vector<FLOAT_TYPE> CMulator::GetRawFrequencyData( std::vector<int> selected_generations )
{
    std::vector<FLOAT_TYPE> tmp_freqs;
    
    // sort so the results would be given in a predictable order
    std::sort( selected_generations.begin(), selected_generations.end() );
    
    for ( auto cur_generation : selected_generations ) {
        
        auto tmp_cur_generation = cur_generation - _generation_shift;
        
        for ( auto cur_allele=0; cur_allele<_num_alleles; cur_allele++ ) {
            
            if ( !IsValid_Generation(cur_generation) ) {
                std::string tmp_str = "GetRawFrequencyData: generation out of range (" +
                std::to_string(tmp_cur_generation) + "). number of generations is " + std::to_string(_num_generations);
                throw tmp_str;
            }
            
            if (_use_observed_data) {
                tmp_freqs.push_back( _observed_simulated_data(cur_generation, cur_allele) );
            }
            else {
                tmp_freqs.push_back( _all_simulated_data(cur_generation, cur_allele) );
            }
        }
    }
    
    return tmp_freqs;
}


FLOAT_TYPE CMulator::GetMutationRate( int from, int to ) const
{
    if ( from < 0 || from >= _num_alleles ) {
        std::string tmp_str = " GetMutationRate: from allele out of range (" + std::to_string(from) + ").";
        throw tmp_str;
    }
    
    if ( to < 0 || to >= _num_alleles ) {
        std::string tmp_str = " GetMutationRate: from allele out of range (" + std::to_string(from) + ").";
        throw tmp_str;
    }
    
    return _mutation_rates_matrix(from,to);
}

/*
void CMulator::SetMutationRate( int from, int to, FLOAT_TYPE rate )
{
    if ( from < 0 || from >= _num_alleles ) {
        throw " SetMutationRate: from allele out of range (" + std::to_string(from) + ").";
    }
    
    if ( to < 0 || to >= _num_alleles ) {
        throw " SetMutationRate: from allele out of range (" + std::to_string(from) + ").";
    }
    
    if ( rate < 0 || rate > 1.0f ) {
        throw " SetMutationRate: rate of mutation not valid (" + std::to_string(rate) + ").";
    }
    
    _mutation_rates_matrix(from,to) = rate;
}
*/

void CMulator::SetAlleleInitFreq( int allele, FLOAT_TYPE freq )
{
    if ( !IsValid_Allele(allele) ) {
        std::string tmp_str = " SetAlleleInitFreq: allele out of range (" + std::to_string(allele) + ").";
        throw tmp_str;
    }
    
    if ( !IsValid_Frequency(freq) && freq != PARAM_DEFAULT_VAL_FLOAT ) {
        std::string tmp_str = " SetAlleleInitFreq: frequency out of range (" + std::to_string(freq) + ").";
        throw tmp_str;
    }
    
    _allele_init_freqs[allele] = freq;
    
    Reset_Soft();
}


void CMulator::SetAlleleMinFitness( int allele, FLOAT_TYPE fitness )
{
    if ( !IsValid_Allele(allele) ) {
        std::string tmp_str = " SetAlleleMinFitness: allele out of range (" + std::to_string(allele) + ").";
        throw tmp_str;
    }
    
    if ( !IsValid_Fitness(fitness  && fitness != PARAM_DEFAULT_VAL_FLOAT ) ) {
        std::string tmp_str = " SetAlleleMinFitness: fitness invalid (" + std::to_string(fitness) + ").";
        throw tmp_str;
    }
    
    _allele_min_fitness[allele] = fitness;
}


void CMulator::SetAlleleMaxFitness( int allele, FLOAT_TYPE fitness )
{
    if ( !IsValid_Allele(allele) ) {
        std::string tmp_str = " SetAlleleMaxFitness: allele out of range (" + std::to_string(allele) + ").";
        throw tmp_str;
    }
    
    if ( !IsValid_Fitness(fitness)  && fitness != PARAM_DEFAULT_VAL_FLOAT ) {
        std::string tmp_str = " SetAlleleMaxFitness: fitness invalid (" + std::to_string(fitness) + ").";
        throw tmp_str;
    }
    
    _allele_max_fitness[allele] = fitness;
}

/*
void CMulator::SetBottleneckInterval( int interval )
{
    if ( !IsValid_BottleneckInterval(interval) ) {
        throw " SetBottleneckInterval: invalid value (" + std::to_string(interval) + ").";
    }
    
    _bottleneck_interval = interval;
}
*/

int CMulator::GetBottleneckInterval() const
{
    return _bottleneck_interval;
}


int CMulator::GetPopulationSize() const
{
    return _N;
}


int CMulator::GetGenerationShift() const
{
    if (!_initialized_with_parameters) {
        std::string tmp_str = "GetGenerationShift: object not initialized with parameters.";
        throw tmp_str;
    }
    
    return _generation_shift;
}

/*
void CMulator::SetGenerationShift(int shift)
{
    if ( shift < 0 ) {
        throw " SetGenerationShift: invalid value (" + std::to_string(shift) + ").";
    }
    
    _generation_shift = shift;
}
*/

int CMulator::GetAlleleNumber() const
{
    if (!_initialized_with_parameters) {
        std::string tmp_str = " GetAlleleNumber: object not initialized with parameters.";
        throw tmp_str;
    }
    
    return _num_alleles;
}


FLOAT_TYPE CMulator::GetInitAlleleFreq(int allele) const
{
    if (!_initialized_with_parameters) {
        std::string tmp_str = "GetInitAlleleFreq: object not initialized with parameters.";
        throw tmp_str;
    }
    
    return _allele_init_freqs[allele];
}


int CMulator::GetWTAllele() const
{
    return _wt_allele_index;
}


void CMulator::SetWTAllele(int wt_index, FLOAT_TYPE min_fitness_nonwt, FLOAT_TYPE max_fitness_nonwt)
{
    for (int i = 0; i < _num_alleles; ++i) {
        SetAlleleMaxFitness(i, max_fitness_nonwt);
        SetAlleleMinFitness(i, min_fitness_nonwt);
    }
    
    SetAlleleMaxFitness(wt_index, 1.0f);
    SetAlleleMinFitness(wt_index, 1.0f);
    
    _wt_allele_index = wt_index;
}

void CMulator::SetWTAllele(int wt_index)
{
    SetAlleleMaxFitness(wt_index, 1.0f);
    SetAlleleMinFitness(wt_index, 1.0f);
    
    _wt_allele_index = wt_index;
}


void CMulator::SetFitnessValues( const std::vector<FLOAT_TYPE> &_given_allele_fitness )
{
    if ( _given_allele_fitness.size() != _num_alleles ) {
        std::string tmp_str = " SetFitnessValues: num of alleles is " +
        std::to_string(_num_alleles) +
        " but trying to set " +
        std::to_string(_given_allele_fitness.size()) + " values ";
        
        throw tmp_str;
    }
    
    // test that the fitness values are valid
    for ( auto tmp_idx=0; tmp_idx<_given_allele_fitness.size(); ++tmp_idx ) {
        
        auto sampled_fitness = _given_allele_fitness[tmp_idx];
        
        if ( !IsValid_Fitness(sampled_fitness) ) {
            std::string tmp_str = "Invalid fitness value sampled from prior for allele" +
            std::to_string(tmp_idx) + ": " + std::to_string(sampled_fitness);
            
            throw tmp_str;
        }
    }
    
    _allele_fitness = _given_allele_fitness;
    _available_fitness = true;
}


unsigned int CMulator::GetRandomSeed() const
{
    return _time_for_seeding;
}


void CMulator::SetRandomSeed(unsigned int new_seed)
{
    _time_for_seeding = new_seed;
}


void CMulator::SetMutationRateMatrix( const MATRIX_TYPE& new_mutation_matrix )
{
    
    // not using matrix row because it cannot accept const matrix
    for ( auto current_row=0; current_row<new_mutation_matrix.size1(); ++current_row ) {
        
        float sum_mutation_rates = 0;
        
        for ( auto current_col=0; current_col<new_mutation_matrix.size2(); ++current_col ) {
            
            auto tmp_mutrate = new_mutation_matrix(current_row, current_col);
            
            if ( tmp_mutrate < 0 ) {
                std::string tmp_str = "Negative mutation rate from allele " +
                std::to_string(current_row) + " to allele " + std::to_string(current_col);
                
                throw tmp_str;
            }
            
            sum_mutation_rates += tmp_mutrate;
        }
        
        if ( std::fabs(1.0 - sum_mutation_rates) > 2.0f * std::numeric_limits<float>::epsilon() ) {
            std::string tmp_str = "Sum of mutation rates from allele " + std::to_string(current_row) + " is not 1 (" + std::to_string(sum_mutation_rates) + ")";
            throw tmp_str;
        }
    }
    
    _mutation_rates_matrix = new_mutation_matrix;
}
    

MATRIX_TYPE CMulator::GetMutationRateMatrix() const
{
    return _mutation_rates_matrix;
}

INT_MATRIX CMulator::GetMinMutationRateMatrix() const
{
    return _min_log_mutation_rate_matrix;
}

INT_MATRIX CMulator::GetMaxMutationRateMatrix() const
{
    return _max_log_mutation_rate_matrix;
}


std::size_t CMulator::GetRepeats() const
{
    return _repeats;
}

void CMulator::SetPopulationSize( int N )
{
    if ( !IsValid_PopulationSize(N) ) {
        std::string tmp_str = "Invalid population size: " + std::to_string(N);
        throw tmp_str;
    }
    
    _N = N;
    _available_popsize = true;
}
