//
//  SimulationResult.hpp
//  fits
//
//  Created by Tal Zinger on 15/11/2016.
//  Copyright © 2016 Stern Lab. All rights reserved.
//

#ifndef SimulationResult_hpp
#define SimulationResult_hpp

#include <vector>
#include <string>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <boost/accumulators/numeric/functional.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/framework/features.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/moment.hpp>
//#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/format.hpp>

#include "CMulator.h"

struct SimulationResult {
    
    // general
    int N;
    int wt_index;
    std::string sim_id;
    FLOAT_TYPE distance_from_actual;
    std::vector<int> actual_generations;

    int generation_shift;
    
    // in order to uniquely identify the sample from the prior used to simulate this result
    std::size_t prior_sample_index;
    
    // fitness
    std::vector<FLOAT_TYPE> fitness_values;
    
    // raw data
    MATRIX_TYPE sim_data_matrix;
    
    // mutation rate
    MATRIX_TYPE mutation_rates;
    
    // generations
    int generation_interval;
    int num_generations;
    
    
    // constructors
    SimulationResult();
    //SimulationResult(std::string id, FLOAT_TYPE distance);
    SimulationResult(const SimulationResult &original); // copy constructor
    SimulationResult(SimulationResult&& other) noexcept; // move constructor
    SimulationResult(const CMulator& sim_object);
    
    SimulationResult(const CMulator& sim_object, std::vector<int> actual_gens);
    
    // swap for sorting
    void swap(SimulationResult& other);
    
    // not sure these swaps are needed
    void swap(SimulationResult& res1, SimulationResult& res2);
    void swap(SimulationResult *res1, SimulationResult *res2);
    
    // for sorting so the rest are redundant
    bool operator<(const SimulationResult& result) const;
    
    SimulationResult& operator=(SimulationResult other);
    
    // this is used when the simulation result is taken as pseudo-data
    std::vector<FLOAT_TYPE> GetSDForEachAllele();
};

#endif /* SimulationResult_hpp */
