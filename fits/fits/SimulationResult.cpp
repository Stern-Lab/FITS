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

#include "SimulationResult.hpp"

SimulationResult::SimulationResult()
: pos(-1),
N(0),
wt_index(-1),
sim_id(""),
distance_metric(""),
distance_from_actual(-1.0),
actual_generations(),
generation_shift(0),
prior_sample_index(0),
fitness_values(),
sim_data_matrix(),
mutation_rates(),
generation_interval(0),
num_generations(0)
{}


SimulationResult::SimulationResult(const SimulationResult &original)
: pos(original.pos),
N(original.N),
wt_index(original.wt_index),
sim_id(original.sim_id),
distance_metric(original.distance_metric),
distance_from_actual(original.distance_from_actual),
actual_generations(original.actual_generations),
generation_shift(original.generation_shift),
prior_sample_index(original.prior_sample_index),
fitness_values(original.fitness_values),
sim_data_matrix(original.sim_data_matrix),
mutation_rates(original.mutation_rates),
generation_interval(original.generation_interval),
num_generations(original.num_generations)
{}


SimulationResult::SimulationResult(const CMulator& sim_object)
: pos(-1),
sim_id(sim_object.GetSimUID()),
distance_from_actual(-1.0f),
distance_metric(""),
fitness_values(sim_object.GetAlleleFitnessValues()),
wt_index(sim_object.GetWTAllele()),
N(sim_object.GetPopulationSize()),
mutation_rates(sim_object.GetMutationRateMatrix()),
prior_sample_index(0),
sim_data_matrix(sim_object.GetAllOutputAsMatrix()),
generation_shift(sim_object.GetGenerationShift()),
actual_generations(),
num_generations(sim_object.GetNumOfGenerations())
{}

SimulationResult::SimulationResult(const CMulator& sim_object, std::vector<int> actual_gens)
: pos(-1),
sim_id(sim_object.GetSimUID()),
distance_from_actual(-1.0f),
distance_metric(""),
fitness_values(sim_object.GetAlleleFitnessValues()),
wt_index(sim_object.GetWTAllele()),
N(sim_object.GetPopulationSize()),
mutation_rates(sim_object.GetMutationRateMatrix()),
prior_sample_index(0),
generation_shift(sim_object.GetGenerationShift()),
actual_generations(actual_gens),
num_generations(sim_object.GetNumOfGenerations())
{
    if ( actual_generations.empty() ) {
        std::cerr << "Actual Generation list is empty" << std::endl;
        throw "Actual Generation list is empty";
    }
    
    sim_data_matrix = sim_object.GetAllOutputAsMatrix( actual_generations );
}


// move constructor
// creates a minimal object that will be destructed with minimal effort
SimulationResult::SimulationResult( SimulationResult&& other ) noexcept
{
    swap(other);
}


// Operators
// the less-than is the only one used for comparisons
bool SimulationResult::operator<(const SimulationResult& result) const
{
    return distance_from_actual < result.distance_from_actual;
}


// Swap functions
void SimulationResult::swap( SimulationResult& other )
{
    std::swap(pos, other.pos);
    std::swap(N, other.N);
    std::swap(wt_index, other.wt_index);
    sim_id.swap(other.sim_id);
    
    std::swap(distance_metric,other.distance_metric),
    std::swap(distance_from_actual, other.distance_from_actual);
    actual_generations.swap(other.actual_generations);
    std::swap(generation_shift, other.generation_shift);
    
    std::swap(prior_sample_index, other.prior_sample_index);

    fitness_values.swap(other.fitness_values);
    sim_data_matrix.swap(other.sim_data_matrix);
    mutation_rates.swap(other.mutation_rates);
    
    std::swap(generation_interval, other.generation_interval);
    std::swap(num_generations, other.num_generations);
}


void SimulationResult::swap(SimulationResult& res1, SimulationResult& res2)
{
    std::cerr << std::endl << " Swap with references is used but not implemented!" << std::endl;
    throw " Swap with references is used but not implemented!";
}


void SimulationResult::swap(SimulationResult *res1, SimulationResult *res2)
{
    std::cerr << std::endl << " Swap *res is used but not implemented!" << std::endl;
    throw " Swap *res is used but not implemented!";
}


SimulationResult& SimulationResult::operator=(SimulationResult other)
{
    swap(other);
    return *this;
}


std::vector<FLOAT_TYPE> SimulationResult::GetSDForEachAllele()
{
    
    std::vector<FLOAT_TYPE> allele_sd_vec(sim_data_matrix.size2(), 0.0f);
    
    for ( auto current_allele=0; current_allele<sim_data_matrix.size2(); ++current_allele ) {
        
        boost::numeric::ublas::matrix_column<MATRIX_TYPE> current_col( sim_data_matrix, current_allele );
        
        boost::accumulators::accumulator_set<
        FLOAT_TYPE,
        boost::accumulators::stats<
        //boost::accumulators::tag::median,
        boost::accumulators::tag::variance,
        boost::accumulators::tag::mean,
        boost::accumulators::tag::min,
        boost::accumulators::tag::max> > allele_accumulator;
        
        //std::cout << "Accumulating allele " << current_allele << std::endl;
        for ( auto current_freq : current_col ) {
            allele_accumulator( current_freq );
            std::cout << current_freq << std::endl;
        }
        
        //current_col = current_col * 2.0f;
        
        allele_sd_vec[current_allele] = std::sqrt( boost::accumulators::variance(allele_accumulator) );
        
        //std::cout << " sd is " << allele_sd_vec[current_allele] << std::endl;
    }
    
    
    return allele_sd_vec;
}
