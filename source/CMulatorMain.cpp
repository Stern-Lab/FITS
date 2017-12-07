// CMulatorMain.cpp : Defines the entry point for the console application.
//

#define DEBUG_VERBOSE

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

int InferABC( FactorToInfer factor,
                std::string param_filename, std::string actual_data_filename,
                std::string posterior_output_filename, std::string summary_output_filename,
                std::string prior_output_filename )
{
    std::cout << "Parameter file: " << param_filename << std::endl;
    std::cout << "Actual data file: " << actual_data_filename << std::endl;
    std::cout << "Posterior distribution file: " << posterior_output_filename << std::endl;
    std::cout << "Summary file: " << summary_output_filename << std::endl;
    
    if ( prior_output_filename.compare("") != 0 ) {
        std::cout << "Prior distribution file: " << prior_output_filename << std::endl;
    }
 
    
    ActualDataFile actual_data_file;
    try {
        std::cout << "Reading actual data... ";
        actual_data_file.LoadActualData(actual_data_filename);
        std::cout << "Done" << std::endl;
        
        //std::cout << "num of alleles: " << actual_data_file.GetNumberOfAlleles() << std::endl;
        //std::cout << "sd for actual data: " << std::endl;;
        //auto resvec = actual_data_file.GetSDPerAllele();
        
        //for ( auto val : resvec ) {
        //    std::cout << "\t" << val;
       // }
        //std::cout << std::endl;
        //return 0;
    }
    catch (std::exception& e) {
        std::cerr << "Exception while loading actual data: " << e.what() << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "Unknown exception while loading actual data." << std::endl;
        return 1;
    }
    
    
    ZParams my_zparams;
    try {
        std::cout << "Reading parameters... ";
        my_zparams.ReadParameters(param_filename, true);
        std::cout << "Done" << std::endl;
    }
    catch (std::exception& e) {
        std::cerr << "Exception while loading parameters: " << e.what() << std::endl;
        return 1;
    }
    catch (std::string str_exp) {
        std::cerr << "Exception while loading parameters: " << str_exp << std::endl;
        return 1;
    }
    catch (const char* str_exp) {
        std::cerr << "Exception while loading parameters: " << str_exp << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "Unknown exception while loading parameters." << std::endl;
        return 1;
    }
    
    
    clsCMulatorABC abc_object_sim( my_zparams, actual_data_file );
    try {
        std::cout << "Starting ABC." << std::endl;
        
        //abc_object_sim.GetUniqueIndexSet(10);
        abc_object_sim.SetImmediateRejection(false);
        //std::cout << "checkpoint alpha" << std::endl;
        auto num_batches = my_zparams.GetInt(fits_constants::PARAM_SIM_REPEATS);
        if ( num_batches > 10 ) num_batches = 10;
        abc_object_sim.RunABCInference(factor, num_batches);
        
        
        if (my_zparams.GetInt(fits_constants::PARAM_COVERAGE_SWITCH, fits_constants::PARAM_DEFAULT_COVERAGE_SWITCH) > 0) {
            std::cout << "Starting coverage." << std::endl;
            //abc_object_sim.DoCoverageTest();
            
        }
        
        //auto start_coverage_time = std::chrono::high_resolution_clock::now();
        //DoCoverageTest();
        //auto end_coverage_time = std::chrono::high_resolution_clock::now();
        //auto coverage_elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_coverage_time - start_coverage_time);
        //total_running_time = static_cast<double>(coverage_elapsed_ms.count()) / 1000.0;
        //std::cout << "Coverage running time: " << total_running_time << " seconds" << std::endl;
        
        std::cout << "ABC Finished." << std::endl;
    }
    catch (std::exception& e) {
        std::cerr << "exception in runnig abc: " << e.what() << std::endl;
        return 1;
    }
    catch (const char *str) {
        std::cerr << "exception in runnig abc: " << str << std::endl;
        return 1;
    }
    catch (std::string str) {
        std::cerr << "exception in runnig abc: " << str << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "unkown exception while runnig abc." << std::endl;
        return 1;
    }
    
    
    
    std::cout << std::endl << "Writing posterior distribution and summary files... " << std::endl << std::endl;
    try {
        ResultsStats result_stats(my_zparams);
        
        // std::cout << " getting results" << std::endl;
        auto accepted_results_vector = abc_object_sim.GetResultsVector(true);
        
        result_stats.SetRejectionThreshold( abc_object_sim.GetRejectionThreshold() );
        
        auto single_mutation_rate = my_zparams.GetFloat( fits_constants::PARAM_SINGLE_MUTATION_RATE, 0.0f );
        result_stats.SetSingleMutrateUsed( single_mutation_rate != 0.0f );
        
        /*
        for ( auto sim_counter=0; sim_counter<accepted_results_vector.size(); ++sim_counter ) {
            std::string tmp_filename;
            tmp_filename += "sim_out_gens" +  std::to_string(accepted_results_vector[sim_counter].num_generations) + "_" + std::to_string(sim_counter) + ".txt";
            abc_object_sim.WriteSimDataToFile(tmp_filename, accepted_results_vector[sim_counter]);
        }
        */
        
        std::cout << "Summary:" << std::endl;
        
        switch (factor) {
            case Fitness: {
                
                result_stats.SetPriorDistrib( my_zparams.GetString(fits_constants::PARAM_PRIOR_DISTRIB, fits_constants::PARAM_PRIOR_DISTRIB_DEFAULT ) );
                
                result_stats.CalculateStatsFitness(accepted_results_vector);
                
                std::string tmp_summary_str = result_stats.GetSummaryFitness();
                
                
                std::cout << tmp_summary_str << std::endl;
                
                
                result_stats.WriteFitnessDistribToFile(accepted_results_vector, posterior_output_filename);
                
                
                result_stats.WriteStringToFile(summary_output_filename, tmp_summary_str);
                
                if ( prior_output_filename.compare("") != 0 ) {
                    std::cout << "Prior... " << std::endl;
                    
                    auto tmp_prior = abc_object_sim.GetPriorFloat();
                    
                    result_stats.WritePriorDistribToFile(tmp_prior, prior_output_filename);
                }
                
                break;
            }
                
            case PopulationSize: {
                result_stats.CalculateStatsPopulationSize(accepted_results_vector);
                std::string tmp_summary_str = result_stats.GetSummaryPopSize();
                
                std::cout << tmp_summary_str << std::endl;
                
                auto tmp_prior = abc_object_sim.GetPriorInt();
                result_stats.WritePriorDistribToFile(tmp_prior, prior_output_filename);
                
                result_stats.WritePopSizeDistribToFile(accepted_results_vector, posterior_output_filename);
                result_stats.WriteStringToFile(summary_output_filename, tmp_summary_str);
                
                if ( prior_output_filename.compare("") != 0 ) {
                    std::cout << "Prior... " << std::endl;
                    
                    auto tmp_prior = abc_object_sim.GetPriorInt();
                    
                    result_stats.WritePriorDistribToFile(tmp_prior, prior_output_filename);
                }
                break;
            }
                
            case MutationRate: {
                auto infer_single_mutrate = my_zparams.GetInt(fits_constants::PARAM_MIN_LOG_SINGLE_MUTATION_RATE, 0);
                
                result_stats.SetSingleMutrateInferred( infer_single_mutrate != 0 );
                result_stats.CalculateStatsMutation(accepted_results_vector);
                std::string tmp_summary_str = result_stats.GetSummaryMutRate();
                
                std::cout << tmp_summary_str << std::endl;
                
                result_stats.WriteMutRateDistribToFile(accepted_results_vector, posterior_output_filename);
                result_stats.WriteStringToFile(summary_output_filename, tmp_summary_str);
                
                if ( prior_output_filename.compare("") != 0 ) {
                    std::cout << "Prior... " << std::endl;
                    
                    auto tmp_prior = abc_object_sim.GetPriorFloat();
                    
                    result_stats.WritePriorDistribToFile(tmp_prior, prior_output_filename);
                }
                break;
            }
                
            case Generations: {
                result_stats.CalculateStatsGenerations(accepted_results_vector);
                std::string tmp_summary_str = result_stats.GetSummaryGenerations();
                
                std::cout << tmp_summary_str << std::endl;
                
                result_stats.WriteGenerationsDistribToFile(accepted_results_vector, posterior_output_filename);
                result_stats.WriteStringToFile(summary_output_filename, tmp_summary_str);
                
                auto tmp_prior = abc_object_sim.GetPriorInt();
                result_stats.WritePriorDistribToFile(tmp_prior, prior_output_filename);
                break;
            }
        }
    }
    catch (std::exception& e) {
        std::cerr << "Exception caught while attempting to report stats: " << e.what() << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "Unknown exception while attempting to report stats." << std::endl;
        return 1;
    }
    
    
    return 0;
}


int RunSingleSimulation(std::string param_filename, std::string output_filename)
{
    // simulator object
    CMulator sim_object;
    ZParams my_zparams;
    
    try {
        my_zparams.ReadParameters(param_filename, true);
        
        sim_object.InitMemberVariables(my_zparams);
        
        if ( !sim_object.IsAbleToSimulate() ) {
            throw "Missing parameters for simulation.";
        }
        
        //sim_object.ReadParametersFromFile(param_filename);
        // sim_object.InitializeFromParamFile(param_filename);
    }
    catch (const char* txt) {
        std::cerr << "Exception caught while initializing simulator object:" << std::endl;
        std::cerr << txt << std::endl;
        return 1;
    }
    catch (std::string txt) {
        std::cerr << "Exception caught while initializing simulator object:" << std::endl;
        std::cerr << txt << std::endl;
        return 1;
    }
    catch (std::exception& e) {
        std::cerr << "Exception caught while initializing simulator object:" << std::endl;
        std::cerr << e.what() << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "Unknown exception while initializing simulator object." << std::endl;
        throw;
        return 1;
    }
    
    std::cout << "Simulator object initialized." << std::endl;
    std::cout << "Random seed set to " << sim_object.GetRandomSeed() << std::endl;
    std::cout << "Running simulation..." << std::endl;
    try {
        sim_object.EvolveAllGenerations();
    }
    catch (const char* txt) {
        std::cerr << "Exception: " << txt << std::endl;
    }
    std::cout << "Done." << std::endl;
    
    //std::string sim_output = sim_object.GetAllOutputAsText();
    std::string sim_output = sim_object.GetAllOutputAsTextForR(true);
    
    // also output to screen
    // TODO: make this via parameter
    std::cout << sim_output << std::endl;
    
    std::cout << std::endl;
    
    //std::cout << sim_object.GetAllOutputAsMatrix() << std::endl;
    
    std::ofstream outfile(output_filename, std::ofstream::out | std::ofstream::trunc);
    
    if (!outfile.is_open()) {
        std::cerr << "unable to open file for writing: " << output_filename << std::endl;
        return 1;
    }
    std::cout << "Writing to file...";
    outfile << sim_output;
    outfile.close();
    std::cout << " Done." << std::endl;
    
    return 0;
}


void TestMutationRates()
{
    boost::numeric::ublas::matrix<float> min(2, 2);
    boost::numeric::ublas::matrix<float> max(2, 2);
    
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
    std::vector<float> min {0.0};
    std::vector<float> max {2.0};
    
    std::vector<unsigned int> min_int { 40, 40, 40 };
    std::vector<unsigned int> max_int { 100, 100, 100 };
    
    //PriorSampler<float> sampler( min, max, PriorDistributionType::SMOOTHED_COMPOSITE );
    PriorSampler<float> sampler( min, max, PriorDistributionType::FITNESS_COMPOSITE );
    
    auto res_vec = sampler.SamplePrior(100000);
    
    for ( auto vec : res_vec ) {
        for ( auto val : vec ) {
            
            std::cout << val << "\t";
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
        std::cout << val << "\t";
    }
    std::cout << std::endl;
    
}

void print_syntaxes(std::string exec_name)
{
    /*
    std::regex rx{"\\S*[\\/\\\\](fits\\S*)"}; // (name) (value) pairs, ignore comments with #
    
    // match regex
    std::smatch matches; // matched strings go here
    
    if (std::regex_search(exec_name , matches, rx, std::regex_constants::match_any )) {
        // 0 is the whole line
        
        for ( auto a : matches ) {
            std::cout << a << std::endl;
        }
        //exec_name = matches[ ];
    }
*/
    // tired of regex and want this to look pretty
    exec_name = "fits ";
    std::cout << "\t" << exec_name << fits_constants::ARG_INFER_FITNESS
    << " <param_file> <actual_data_file> <posterior_file> <summary_file> (optional: <prior_file>)" << std::endl << std::endl;
    
    std::cout << "\t" << exec_name << fits_constants::ARG_INFER_MUTATION
    << " <param_file> <actual_data_file> <posterior_file> <summary_file> (optional: <prior_file>)" << std::endl << std::endl;
    
    std::cout << "\t" << exec_name << fits_constants::ARG_INFER_POPSIZE
    << " <param_file> <actual_data_file> <posterior_file> <summary_file> (optional: <prior_file>)" << std::endl << std::endl;
    
    /*
     std::cout << "\t" << exec_name << fits_constants::ARG_INFER_GENERATION
    << " <param_file> <actual_data_file> <posterior_file> <summary_file> (optional: <prior_file>)" << std::endl << std::endl;
     */
    
    std::cout << "\t" << exec_name << "-simulate <param_file> <output_file>" << std::endl << std::endl;
}

void print_welcome()
{
    std::cout << std::endl << "================================================";
    std::cout << std::endl << "                  FITS v" << fits_constants::current_version_str;
    std::cout << std::endl << "    Flexible Inference from Time-Series data    ";
    std::cout << std::endl << "         (c) Tal Zinger, Stern Lab, TAU         ";
    std::cout << std::endl << "================================================";
    std::cout << std::endl;
    
    
    
}

bool IsInferenceRun( std::string first_argument )
{
    return ( first_argument.compare( fits_constants::ARG_INFER_FITNESS ) == 0
            || first_argument.compare( fits_constants::ARG_INFER_MUTATION ) == 0
            || first_argument.compare( fits_constants::ARG_INFER_POPSIZE ) == 0
            || first_argument.compare( fits_constants::ARG_INFER_GENERATION ) == 0 );
}

int main(int argc, char* argv[])
{
    
    if (argc <= 1) {
        print_welcome();
        print_syntaxes(argv[0]);
        return 1;
    }
    
    std::string tmp_first_param = argv[1];
    
    // This would make numbers nicer
    std::cout.imbue(std::locale(fits_constants::used_locale));
    
    std::string param_filename = "";
    std::string actual_data_filename = "";
    std::string posterior_output_filename = "";
    std::string summary_output_filename = "";
    std::string prior_output_filename = "";
    
    if ( IsInferenceRun(tmp_first_param) ) {
        
        // prior is optional
        if (argc != 6 && argc != 7) {
            std::cout << "illegal number of arguments (" << argc-1 << "). Syntax is:" << std::endl;
            print_syntaxes(argv[0]);
            return 1;
        }
        
        param_filename = argv[2];
        actual_data_filename = argv[3];
        posterior_output_filename = argv[4];
        summary_output_filename = argv[5];
        
        if ( argc == 7 ) {
            prior_output_filename = argv[6];
        }
    }
    
    if ( tmp_first_param.compare( fits_constants::ARG_INFER_FITNESS ) == 0 ) {
        
        print_welcome();
        
        std::cout << "Inferring fitness" << std::endl;
        
        return InferABC( FactorToInfer::Fitness,
                         param_filename, actual_data_filename,
                         posterior_output_filename, summary_output_filename,
                         prior_output_filename);
    }
    
    if ( tmp_first_param.compare( fits_constants::ARG_INFER_MUTATION ) == 0 ) {
        
        print_welcome();
        
        std::cout << "Inferring mutation rate" << std::endl;
    
        return InferABC( FactorToInfer::MutationRate,
                        param_filename, actual_data_filename,
                        posterior_output_filename, summary_output_filename,
                        prior_output_filename);
    }
    
    if ( tmp_first_param.compare( fits_constants::ARG_INFER_POPSIZE ) == 0 ) {
        
        print_welcome();
        
        std::cout << "Inferring population size" << std::endl;
        
        return InferABC( FactorToInfer::PopulationSize,
                        param_filename, actual_data_filename,
                        posterior_output_filename, summary_output_filename,
                        prior_output_filename );
    }
    
    if ( tmp_first_param.compare( fits_constants::ARG_INFER_GENERATION ) == 0 ) {
        
        print_welcome();
        
        std::cout << "Inferring generation interval" << std::endl;
        
        return InferABC( FactorToInfer::Generations,
                        param_filename, actual_data_filename,
                        posterior_output_filename, summary_output_filename,
                        prior_output_filename );
    }
    
    if (tmp_first_param == "-simulate") {
        
        print_welcome();
        
        if (argc != 4) {
            std::cout << "illegal number of arguments (" << argc-1 << "). Syntax is:" << std::endl;
            print_syntaxes(argv[0]);
            return 1;
        }
        
        std::cout << "Running a single simulation" << std::endl;
        
        param_filename = argv[2];
        std::string output_filename = argv[3];
        
        return RunSingleSimulation(param_filename, output_filename);
    }
    
    if (tmp_first_param == "-test_range") {
        try {
            test_range();
        }
        catch( const char* txt ) {
            std::cerr << "exception in test_range: " << txt << std::endl;
        }
        
        return 1;
    }
 
    if (tmp_first_param == "-test_actual") {
        
        test_actualdata( argv[2] );
        return 1;
    }
    
    if (tmp_first_param == "-levene") {
        
        ZParams myparams;
        std::string strparams = argv[2];
        myparams.ReadParameters(strparams, false);
        ResultsStats result_stats(myparams);
        
        //std::vector<float> vec1 {1.0f, 1.4f, 1.6f, 4.6f, 9.5f, 3.6f};
        //std::vector<float> vec2 {1.0f, 1.4f, 1.6f, 1.6f, 1.5f, 3.6f};
        
        std::vector<float> vec1 {1,3,5,6,8};
        std::vector<float> vec2 {1,3,6,7,8};
        
        auto res = result_stats.LevenesTest2(vec1, vec2);
        std::cout << "levenes output==" << res << std::endl;
        return 1;
    }
    
    
        
    // we should never reach here...
    print_welcome();
    std::cout << "Invalid suntax (option chosen was " << tmp_first_param << "). Use only the following:" << std::endl;
    print_syntaxes(argv[0]);
    return 1;
}

