#include <iostream>
#include <chrono>
#include <fstream>
#include <iomanip>

#include <regex> // for printing only executable name

#include <boost/format.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random/bernoulli_distribution.hpp>

#include "CMulator.h"
#include "clsCMulatorABC.h"

#include "PriorSampler.hpp"
#include "SimulationResult.hpp"
#include "ResultsStats.hpp"
#include "ActualDataFile.hpp"

#include "fits_constants.h"


void TestMutationRates()
{
    boost::numeric::ublas::matrix<FLOAT_TYPE> min(2, 2);
    boost::numeric::ublas::matrix<FLOAT_TYPE> max(2, 2);
    
    for (auto i = 0; i < min.size1(); ++i) {
        for (auto j = 0; j < min.size2(); ++j) {
            
            min(i, j) = 0;
            max(i, j) = 0.01;
        }
    }
}
void test_parameters( std::string filename )
{
    ZParams my_zparams(filename, true);
    
    CMulator sim(my_zparams);
    
    // todo: use flag to dump the object as it is created, i don't retain the params along with member variables
    //std::cout << sim.GetAllZParams() << std::endl;
}


void test_range()
{
    //PriorSampler<float>::PriorDistributionType dist_type = PriorSampler<float>::PriorDistributionType::UNIFORM;
    //template class PriorSampler<float>;
    std::vector<FLOAT_TYPE> min {0.0};
    std::vector<FLOAT_TYPE> max {2.0};
    
    std::vector<unsigned int> min_int { 40, 40, 40 };
    std::vector<unsigned int> max_int { 100, 100, 100 };
    
    //PriorSampler<float> sampler( min, max, PriorDistributionType::SMOOTHED_COMPOSITE );
    PriorSampler<FLOAT_TYPE> sampler( min, max, PriorDistributionType::FITNESS_COMPOSITE );
    
    auto res_vec = sampler.SamplePrior(100000);
    
    for ( auto vec : res_vec ) {
        for ( auto val : vec ) {
            
            std::cout << val << fits_constants::FILE_FIELD_DELIMITER;
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    
}

void test_actualdata( std::string filename )
{
    ActualDataFile datafile;
    datafile.LoadActualData(filename);
    
    auto vec = datafile.GetActualFrequencies();
    
    std::cout << std::endl;
    for ( auto val : vec ) {
        std::cout << val << fits_constants::FILE_FIELD_DELIMITER;
    }
    std::cout << std::endl;
    
}

