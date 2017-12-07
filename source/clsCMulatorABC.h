//
//  clsCMulatorABC.hpp
//
//  Created by Tal Zinger
//  Copyright (c) 2017 Tal Zinger. All rights reserved.
//

#ifndef clsCMulatorABC_hpp
#define clsCMulatorABC_hpp

#include <future>
#include <vector>
#include <chrono>
#include <algorithm>

// performance
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/stats.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/odeint/util/ublas_wrapper.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>

#include <cmath> //for pow
#include <iomanip> // for put_time
#include <string>
#include <locale>

#include "ZParams.h"
#include "CMulator.h"
#include "fits_constants.h"

#include "SimulationResult.hpp"
#include "PriorSampler.hpp"
#include "ActualDataFile.hpp"

enum FactorToInfer {
    Fitness, MutationRate, PopulationSize, Generations,
    None // just to catch extreme cases
};


class clsCMulatorABC 
{
private:
    
    double _total_running_time;
    
    boost::mt19937 _boost_gen;
    
    const FLOAT_TYPE THRESHOLD_RESET_VALUE = -1.0;
    
    std::size_t _repeats;
    // std::size_t _percent_to_keep;
    std::size_t _sims_to_keep;
    
    FactorToInfer _factor_to_infer;
    
    unsigned int _num_alleles;
    
	FLOAT_TYPE _rejection_threshold;
    bool _use_rejection_threshold; 
    
    ZParams _zparams;
    
    PriorDistributionType _prior_type;

	std::vector<SimulationResult> _simulation_result_vector;
    
    ActualDataFile _actual_data_file;
    
    std::vector<ActualDataEntry> _actual_data_vector;
    std::vector<FLOAT_TYPE> _actual_data_raw_freqs;
    std::vector<int> _selected_actual_generations;
    

    
    // stores samples from the prior
    std::vector< std::vector<FLOAT_TYPE> > _float_prior_archive;
    std::vector< std::vector<int> > _int_prior_archive;
    
public:
    
    clsCMulatorABC( ZParams sim_params, ActualDataFile actual_data_file );

    FLOAT_TYPE GetMedian( std::vector<FLOAT_TYPE> vec );
    
    std::vector<int> GetUniqueIndexSet( int num_items );

    double GetTotalRunningTime() { return _total_running_time; }
    
    FLOAT_TYPE ResetRejectionThreshold();
	FLOAT_TYPE SetRejectionThreshold(FLOAT_TYPE new_threshold);
	FLOAT_TYPE GetRejectionThreshold();
    
    std::size_t GetNumberOfKeptResults();
    
    void SetImmediateRejection(bool new_val);

    std::size_t GetRepeats();
    
    std::vector<SimulationResult> GetResultsVector(bool only_accepted_results=false);
    
    void RunABCInference( FactorToInfer factor, std::size_t number_of_batches );
    std::vector<SimulationResult> RunFitnessInferenceBatch( std::size_t num_simulations );
    std::vector<SimulationResult> RunPopulationSizeInferenceBatch( std::size_t num_simulations );
    std::vector<SimulationResult> RunMutationInferenceBatch( std::size_t num_simulations );
    std::vector<SimulationResult> RunGenerationInferenceBatch( std::size_t num_simulations );

    std::vector<FLOAT_TYPE> GetSDPerAllele( std::size_t start_idx, std::size_t end_idx );
    std::vector<FLOAT_TYPE> GetMADPerAllele( std::size_t start_idx, std::size_t end_idx );
    void DivideEachAllele( std::size_t start_idx, std::size_t end_idx, std::vector<FLOAT_TYPE> value_vector );
    
    
    //void ScaleSimulationResults();
    void CalculateResultsDistances( FLOAT_TYPE scaling_factor = 1.0f );
    
    std::vector< std::vector<FLOAT_TYPE>> GetPriorFloat();
    std::vector< std::vector<int>> GetPriorInt();
    
    std::string GetPriorFloatAsString();
    std::string GetPriorIntAsString();
    
    void DoCoverageCook();
    
    std::pair<bool,FLOAT_TYPE> DoLevenesTest();
    void DoCoverageTest();
    
    MATRIX_TYPE GetAlleleCoveragePvals() const;
    // using unsigend int and not size_t to allow easier cast to float
    std::vector<unsigned int> CoverageSingleDatasetFitness( std::size_t dataset_idx, std::size_t start_idx, std::size_t end_idx ) const;
    
    // std::vector<SimulationResult> RunGenerationInferenceBatch( std::size_t num_simulations );
    
    // TODO: make default with vector of 1.0 for scaling
    FLOAT_TYPE GetDistanceSimActual( const MATRIX_TYPE &actual_data, const MATRIX_TYPE &sim_data, const std::vector<FLOAT_TYPE> &scaling_vector );
    
    void WriteStringToFile( std::string filename, std::string str );
    
    void WriteSimDataToFile( std::string filename, SimulationResult &res );
};

#endif
