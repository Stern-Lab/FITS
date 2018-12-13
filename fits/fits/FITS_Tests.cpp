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

void TestActualData()
{
    std::string single_pos_filename = "/Users/talzinger/Nobackup/Test/Test20181203/actual_data_test1.txt";
    std::string multi_pos_filename = "/Users/talzinger/Nobackup/Test/Test20181203/actual_data_test2.txt";
    
    ActualDataFile data_file;
    
    
    try {
        //data_file.LoadActualData(single_pos_filename);
        data_file.LoadActualData(multi_pos_filename);
        
        //auto data_matrix = data_file.GetActualFreqsAsMatrix(6);
        //std::cout << data_matrix << std::endl;
        
        auto vec = data_file.GetPositionNumbers();
        for ( auto val : vec ) {
            
            std::cout << "position " << val << std::endl;
            
            //auto vec2 = data_file.GetActualGenerations(true, val);
            auto val2 = data_file.GetNumberOfAlleles(val);
            std::cout << val2 << " ";
            //for (auto val2 : vec2) {
            //    std::cout << val2 << " ";
            //}
            
        }
        std::cout << std::endl;
    }
    catch (const char* str) {
        std::cerr << "Error const char*: " << str << std::endl;
    }
    catch (const std::string str) {
        std::cerr << "Error str: " << str << std::endl;
    }
    
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

