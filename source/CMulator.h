//
//  CMulator.h
//  simFourAllelesV1
//
//  Created by Tal Zinger on 6/9/15.
//  Copyright (c) 2015 Tal Zinger. All rights reserved.
//


//#define _SCL_SECURE_NO_WARNINGS

//#ifndef __simFourAllelesV1__CMulator__
//#define __simFourAllelesV1__CMulator__
#ifndef __CMulator__
#define __CMulator__

#include <cmath>
#include <chrono>
#include <random>
#include <iostream>
#include <string>
#include <numeric>

#include <boost/numeric/ublas/io.hpp>
#include <boost/random/binomial_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>


#include "ZParams.h"
#include "fits_constants.h"






class CMulator {
public:
    
    
private:
    
    bool _debug_mode;
    
    // No need to create a member of this object.
    // create this object only when initializing and then use only member variables, not parameters
    // ZParams _zparams;
	
	/* Basic Parameters For Simulation */
	// TODO: turn all int to unsigend int
	int _N;		// current population size
	int _N0; 	// founding population size
    
    int _alt_N; // alternative population size to introduce noise to simulation
    int _old_N; // the original value, to revert from altenative value
    int _alt_generation; // if given, in this generation the popsize will alter; otherwise it would be random.
    
	int _num_generations;
	int _num_alleles;
	int _bottleneck_interval;
	int _bottleneck_size;
    
    int _sample_size;
	int _generation_shift;
    
    // 2017-02-26 changed _uid to _name_of_run
	std::string _name_of_run;
    
	//bool _no_init_freqs_as_parameters; // init freqs will be taken later
	bool _logistic_growth;
	int _logistic_growth_K;		// population capacity
	int _logistic_growth_t;		// time; couples to number of generations but may be reset following bottleneck
	FLOAT_TYPE _logistic_growth_r;	// growth rate
	int _repeats;

	int _wt_allele_index;
	//bool _infer_wt_automatically;		// don't expect min and max values for each alleles. infer the wt, assign 1 to min and max, and assign arbitrary min and max for rest of alleles

    //bool _infer_mutation_rate;
    
    FLOAT_TYPE _acceptance_rate;
    
	
    /* Fitness */
	//FLOAT_TYPE _fitness_increment;
	std::vector<FLOAT_TYPE> _allele_fitness;
	std::vector<FLOAT_TYPE> _allele_init_freqs;
	std::vector<FLOAT_TYPE> _allele_max_fitness;
	std::vector<FLOAT_TYPE> _allele_min_fitness;
	
    // 2016-11-15
    // What data is available to the simulator - the reset would have to be supplied before actually running simulation
    bool _available_mutrate;
    bool _available_mutrate_range;
    bool _available_fitness;
    bool _available_fitness_range;
    bool _available_popsize;
    bool _available_popsize_range;
    bool _available_essential;
    bool _available_inference;
    bool _available_init_freqs;
    bool _available_generation_range;
	
	// store a single simulation - column for allele, row for generation
	// TODO: maybe switch to matrix data type to enable matrix multiplications?
	//typedef boost::multi_array<FLOAT_TYPE,2> allele_dataArray;
	//typedef allele_dataArray::index _allele_index;
	//allele_dataArray _sim_data;
    
    /* flags */
    bool _use_observed_data; // flag - when observed is used instead of all for user
    bool _use_single_mutation_rate;
    
    //boost::numeric::ublas::matrix<FLOAT_TYPE> _all_simulated_data;
    //boost::numeric::ublas::matrix<FLOAT_TYPE> _observed_simulated_data;
	
    MATRIX_TYPE _all_simulated_data;
    MATRIX_TYPE _observed_simulated_data;
    
	/* Mutation rates */
	//typedef boost::multi_array<FLOAT_TYPE,2> allele_mutationArray;
	//typedef allele_mutationArray::index _mutation_rate_index;
	//allele_mutationArray _mutation_rates;
    
    MATRIX_TYPE _mutation_rates_matrix;
    INT_MATRIX _min_log_mutation_rate_matrix;
    INT_MATRIX _max_log_mutation_rate_matrix;
	/* END Basic Parameters For Simulation */
	
    unsigned int _Nlog_min;
    unsigned int _Nlog_max;
	
	/* Runtime Parameters */
	int _current_generation;
	
	/* END Runtime Parameters */
	
	
	
	/* Assisting Functions */
	
    void WarnAgainstDeprecatedParameters(ZParams zparams);
    
    void PerformChecksBeforeEvolution();
    void InitBasicVariables( ZParams zparams );
    void InitInitialAlleleFreqs( ZParams zparams );
    
    void InitGeneralInferenceVariables( ZParams zparams );
    
    void InitPopulationSizeInference( ZParams zparams );
    
    // boolean flag to know whether we should be expecting the mutation rate at all (or the range was given)
    void InitMutationRates( ZParams zparams, bool available_mutrate_range );
    
    void InitMutationInferenceVariables( ZParams zparams );
    
    void InitFitnessValues( ZParams zparams );
    void InitFitnessInference( ZParams zparams );
    void InitGenerationInference( ZParams zparams );
    
    std::vector<FLOAT_TYPE> Freqs2Binomial2Freqs( const std::vector<FLOAT_TYPE> &freqs_vec, int popsize );
    
	/* END Assisting Functions */
	
	
	/* Technical Parameters */
	// for stochastic step
	unsigned long long _my_random_seed;
    boost::mt19937 _boost_gen;

	unsigned int _time_for_seeding;
	
	//int _first_sim_in_batch;
	
	// if the object was created with no parameter file, make sure it is later initialized with parameters
	bool _initialized_with_parameters;
	
	// to mitigate floating point rounding errors
	FLOAT_TYPE _epsilon_float_compare;
	/* END Technical Parameters */
	
	
	
	/* Internal Sanity Checks */
	bool IsValid_Allele( int allele ) const { return ( allele >= 0 && allele < _num_alleles ); }
    bool IsValid_NumAlleles( int num_alleles ) const { return num_alleles >= 2; }
	
    bool IsValid_Generation( int generation ) const { return ( generation >=0  ); }
    
	bool IsValid_Frequency( FLOAT_TYPE freq ) const { return freq >= 0 && freq <= 1.0; }
	bool IsValid_Fitness( FLOAT_TYPE fitness ) const { return fitness >= 0; }
	bool IsValid_PopulationSize( int N ) const { return N>0; }
	//bool IsValid_SamplingSize( int sample_size ) const { return sample_size > 0 && sample_size <=_N; }
	bool IsValid_BottleneckInterval( int interval ) const { return interval < _num_generations; } // non-positive treated as no bottleneck
	bool IsValid_BottleneckSize( int size ) const { return size > 0 && size <= _N; }
	/* End Internal Sanity Checks */
	
// set functions
// used to be public
    // but many of them require reset to make the simulated data in sync with the parameters
    // I don't like this and don't like the static param, library
    
public:
    bool IsObservedDataUsed() { return _use_observed_data; }
    bool IsSingleMutrateUsed() { return _use_single_mutation_rate; }
    
    // using the available data, is the simulator able to do that
    bool IsAbleToInferFitness() const {
        if (!_available_mutrate) std::cerr << "mutrate not available" << std::endl;
        if (!_available_popsize) std::cerr << "popsize not available" << std::endl;
        if (!_available_essential) std::cerr << "essential not available" << std::endl;
        if (!_available_fitness_range) std::cerr << "fitness range not available" << std::endl;
        return _available_mutrate && _available_popsize && _available_essential && _available_fitness_range;
    }
    bool IsAbleToInferGeneration() const {
        if (!_available_mutrate) std::cerr << "mutrate not available" << std::endl;
        if (!_available_popsize) std::cerr << "popsize not available" << std::endl;
        if (!_available_fitness) std::cerr << "fitness not available" << std::endl;
        if (!_available_essential) std::cerr << "essential not available" << std::endl;
        if (!_available_generation_range) std::cerr << "generation range not available" << std::endl;
        return _available_mutrate && _available_fitness && _available_popsize && _available_essential && _available_generation_range;
    }
    bool IsAbleToInferPopulationSize() const {
        if (!_available_mutrate) std::cerr << "mutrate not available" << std::endl;
        if (!_available_fitness) std::cerr << "fitness not available" << std::endl;
        if (!_available_essential) std::cerr << "essential not available" << std::endl;
        if (!_available_popsize_range) std::cerr << "popsize range not available" << std::endl;
        return _available_mutrate && _available_fitness && _available_essential && _available_popsize_range;
    }
    bool IsAbleToInferMutationRate() const {
        if (!_available_mutrate_range) std::cerr << "mutrate range not available" << std::endl;
        if (!_available_popsize) std::cerr << "popsize not available" << std::endl;
        if (!_available_essential) std::cerr << "essential not available" << std::endl;
        if (!_available_fitness) std::cerr << "fitness not available" << std::endl;
        return _available_fitness && _available_popsize && _available_essential && _available_mutrate_range;
    }
    bool IsAbleToSimulate() const {
        if (!_available_mutrate) std::cerr << "mutrate rates not available" << std::endl;
        if (!_available_popsize) std::cerr << "popsize not available" << std::endl;
        if (!_available_essential) std::cerr << "essential not available" << std::endl;
        if (!_available_fitness) std::cerr << "fitness not available" << std::endl;
        if (!_available_init_freqs) std::cerr << "initial frequencies not available" << std::endl;
        if ( _num_generations < 1 ) std::cerr << "zero generations to simulate" << std::endl;
        
        return _available_mutrate && _available_popsize &&_available_essential && _available_fitness && _available_init_freqs && ( _num_generations >= 1 );
    }
private:
    // Past Set functions
    
    void SetAlleleNumber( int alleles );
    void SetBottleneckInterval( int interval );
    void SetBottleneckSize( int size );
    
    void SetSamplingSize( int sampling_size );
    void SetSkipStochasticStep(bool skip);
    
    
    void SetRepeats( int repeats );
    void SetTopPercentToKeep(unsigned int new_percent);
    void SetTopSimsToKeep(unsigned int new_sims);
    
    void EnableLogisticGrowth();
    void EnableLogisticGrowth( int k, FLOAT_TYPE r );
    void DisableLogisticGrowth();
    bool GetLogisticGrowth( int &k, FLOAT_TYPE &r ) const;
    
    void SetAlleleMinFitness( int allele, FLOAT_TYPE fitness );
    void SetAlleleMaxFitness( int allele, FLOAT_TYPE fitness );
    void SetAlleleFitnessValue( int allele, FLOAT_TYPE fitness );
    
    
    void SetMutationRate( int from, int to, FLOAT_TYPE rate );
    
public:
	/* Constructors */
	CMulator();
	//CMulator( std::string param_filename );
    
	CMulator( const CMulator &original );
    CMulator( const ZParams &zparams );
    
    //TODO: move ctor
	/* END Constructors */
	
    
	/* Get/Set Functions */
    
    // functions relying on actual data
    void SetGenerationShift( int shift );
    void SetAlleleInitFreq( int allele, FLOAT_TYPE freq );
    void SetNumOfGeneration( int generations );
    void SetWTAllele( int wt_index, FLOAT_TYPE min_fitness_nonwt, FLOAT_TYPE max_fitness_nonwt );
    void SetWTAllele( int wt_index);
    void SetFitnessValues( const std::vector<FLOAT_TYPE> &_allele_fitness );
    
    void InitMemberVariables( ZParams zparams );
    
	// to be used only when the object is still uninitialized (default ctor used)
    
	void SetSimUID( std::string new_name_of_run );
	int GetPopulationSize() const;
	void SetPopulationSize( int N );
    
	int GetNumOfGenerations() const;
	
	int GetAlleleNumber() const;
	
	int GetBottleneckInterval() const;
		
	int GetBottleneckSize() const;
	
	int GetSamplingSize() const;
	
	int GetGenerationShift() const;
	
	std::string GetSimUID() const;

	bool GetSkipStochasticStep() const;
	
    // to avoid direct access to array - allow to raise exceptions
    //FLOAT_TYPE GetSimulatedFrequency( int generation, int allele ) const;
    //void SetSimulatedFrequency(  int generation, int allele, FLOAT_TYPE frequency );

	unsigned int GetTopPercentToKeep() const;
	
	unsigned int GetTopSimsToKeep() const;
	
	// this shouldn't be here.. but rather in a wrapping class
	// but I can't afford any other thread calling parameters class
	int GetRepeats() const;
	
	FLOAT_TYPE GetInitAlleleFreq(int allele) const;
	
	std::vector<FLOAT_TYPE> GetAlleleInitFreqs() const;
	void SetAlleleInitFreqs( std::vector<FLOAT_TYPE> freqs );

	std::string GetRFileName() const;
	void SetRFileName( const std::string& filename );
	
    std::vector<FLOAT_TYPE> GetAlleleFitnessValues() const;
	std::vector<FLOAT_TYPE> GetAlleleMaxFitnessValues() const;
	std::vector<FLOAT_TYPE> GetAlleleMinFitnessValues() const;

	int GetWTAllele() const;
	
	FLOAT_TYPE GetMutationRate( int from, int to ) const;
    
    void SetMutationRateMatrix( MATRIX_TYPE new_mutation_matrix );
    MATRIX_TYPE GetMutationRateMatrix() const;
    INT_MATRIX GetMinMutationRateMatrix() const;
    INT_MATRIX GetMaxMutationRateMatrix() const;
    
	FLOAT_TYPE GetAlleleFreq( int generation, int allele ) const;
	/* END Get/Set Functions */
	
	
	/* Other Information Retreival */
	unsigned int GetRandomSeed() const;
	void SetRandomSeed(unsigned int new_seed);
	/* END Other Information  Retrieval */
	
	
	/* Operational Functions */
	int EvolveToGeneration( int target_generation );
	int EvolveAllGenerations();
	
    
	/* Output Functions */
    std::string GetAllOutputAsTextForR( bool header = true ) const;
    std::string GetAllOutputAsText( bool header = true, std::string delimiter = "\t" ) const;
	//std::string GetAllOutputAsCSV() const;
    
    MATRIX_TYPE GetAllOutputAsMatrix() const;
    MATRIX_TYPE GetAllOutputAsMatrix( std::vector<int> actual_generations ) const;
    

    std::vector<FLOAT_TYPE> GetRawFrequencyData();
	
	// export only subset of generation, if actual data is not continuous (i.e missing generations)
    std::vector<FLOAT_TYPE> GetRawFrequencyData( std::vector<int> selected_generations );
	/* END Output Functions */
	
	/* Reset */
	// delete simulated data, set current generation to 0.
	void Reset_Soft();
	void Reset_Soft( std::string new_name_of_run );
	/* END Reset */	
	
};

#endif /* defined(__CMulator__) */
