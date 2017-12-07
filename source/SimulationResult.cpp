//
//  SimulationResult.cpp
//  fits
//
//  Created by Tal Zinger on 15/11/2016.
//  Copyright © 2016 Stern Lab. All rights reserved.
//

#include "SimulationResult.hpp"


/*
 std::string sim_id;
 FLOAT_TYPE distance_from_actual;
 
 std::vector<FLOAT_TYPE> fitness_values;
 
 MATRIX_TYPE mutation_rates;
 int N;
 
 int wt_index;
 
 // need to save memory, reaching over 13GB 2016-04-20
 // reduced number of simulation so checking this again
 std::string raw_data;
 */


SimulationResult::SimulationResult()
: N(0),
wt_index(-1),
sim_id(""),
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

/*
SimulationResult::SimulationResult(std::string id, FLOAT_TYPE distance)
: sim_id(id),
distance_from_actual(distance),
fitness_values(),
wt_index(-1),
N(0),
mutation_rates(),
prior_sample_index(0),
sim_data_matrix(),
generation_shift(0),
actual_generations()
{}
*/

SimulationResult::SimulationResult(const SimulationResult &original)
: N(original.N),
wt_index(original.wt_index),
sim_id(original.sim_id),
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
: sim_id(sim_object.GetSimUID()),
distance_from_actual(-1.0f),
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
: sim_id(sim_object.GetSimUID()),
distance_from_actual(-1.0f),
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
    std::swap(N, other.N);
    std::swap(wt_index, other.wt_index);
    sim_id.swap(other.sim_id);
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
    
    /*
    std::swap(res1.N, res2.N);
    std::swap(res1.wt_index, res2.wt_index);
    std::swap(res1.distance_from_actual, res2.distance_from_actual);
    std::swap(res1.prior_sample_index, res2.prior_sample_index );
    std::swap(res1.generation_shift, res2.generation_shift );
    
    
    res1.sim_id.swap(res2.sim_id);
    res1.fitness_values.swap(res2.fitness_values);
    res1.sim_data_matrix.swap(res2.sim_data_matrix);
    res1.actual_generations.swap(res2.actual_generations);
     */
}


void SimulationResult::swap(SimulationResult *res1, SimulationResult *res2)
{
    std::cerr << std::endl << " Swap *res is used but not implemented!" << std::endl;
    throw " Swap *res is used but not implemented!";
    
    /*
    std::swap(res1->N, res2->N);
    std::swap(res1->wt_index, res2->wt_index);
    std::swap(res1->distance_from_actual, res2->distance_from_actual);
    std::swap(res1->prior_sample_index, res2->prior_sample_index );
    std::swap(res1->generation_shift, res2->generation_shift );
    
    
    res1->sim_id.swap(res2->sim_id);
    res1->fitness_values.swap(res2->fitness_values);
    res1->sim_data_matrix.swap(res2->sim_data_matrix);
    res1->actual_generations.swap(res2->actual_generations);
     */
}


// hack I learned that unifies assignment for copy and move
// other is passed by value, so it's a copy
// swapping will essentialy move the data, leaving other with garbage
// either way, "other" will be destructed, leaving us with the data
// and the original rvalue is untouched
SimulationResult& SimulationResult::operator=(SimulationResult other)
{
    swap(other);
    return *this;
}

// TODO: move this to result stats - it needs to go through all 100 best results or something
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


