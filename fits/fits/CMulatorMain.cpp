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
    
    
    ZParams my_zparams;
    std::cout << "Reading parameters... ";
    try {
        my_zparams.ReadParameters(param_filename, true);
        std::cout << "Done." << std::endl;
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
    
    
    ActualDataFile actual_data_file;
    std::cout << "Reading actual data... ";
    try {
        actual_data_file.LoadActualData(actual_data_filename);
        
        auto positions_detected = actual_data_file.GetNumberOfPositions();
        
        if (positions_detected>1) {
            std::cout << "Done - Multiple positions detected (" << positions_detected << ")." << std::endl;
        }
        else {
            std::cout << "Done." << std::endl << std::endl;
        }
        
    }
    catch (std::exception& e) {
        std::cerr << "Exception while loading actual data: " << e.what() << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "Unknown exception while loading actual data." << std::endl;
        return 1;
    }
    
    
    std::cout << "Running simulations:" << std::endl;
    FLOAT_TYPE tmp_rejection_threshold = -1.0f;
    std::size_t running_time_sec = 0;
    std::vector<SimulationResult> accepted_results_vector(0);
    std::vector< std::vector<FLOAT_TYPE>> used_prior_distrib;
    auto results_to_accept = -1;
    
    try {
        
        auto positions_detected = actual_data_file.GetNumberOfPositions();
        
        if (positions_detected>1) {
            
            auto actual_positions_vec = actual_data_file.GetPositionNumbers();
            
            for ( auto current_position_num : actual_positions_vec ) {
                
                auto current_position_data = actual_data_file.GetPosition(current_position_num);
                std::cout << "-- Position " << current_position_num << " --" << std::endl;
                
                /* ----------------- */
                /*  Run simulations  */
                /* ----------------- */
                clsCMulatorABC abc_object_sim( my_zparams, current_position_data );
                
                results_to_accept = abc_object_sim.GetNumberOfKeptResults() * positions_detected;
                
                abc_object_sim.SetImmediateRejection(false);
                
                auto num_batches = my_zparams.GetInt(fits_constants::PARAM_SIM_REPEATS);
                
                if ( num_batches > 10 ) {
                    num_batches = 10;
                }
                
                abc_object_sim.RunABCInference(factor, num_batches);
                running_time_sec += abc_object_sim.GetSunningTimeSec();
                tmp_rejection_threshold = abc_object_sim.GetRejectionThreshold();
                
                /* --------------- */
                /*  Store results  */
                /* --------------- */
                std::cout << "Gathering results... ";
                auto tmp_accepted_results_vec = abc_object_sim.GetResultsVector(false);
                for ( auto tmp_result : tmp_accepted_results_vec ) {
                    accepted_results_vector.push_back(tmp_result);
                }
                std::cout << "Done." << std::endl;
                
                if ( prior_output_filename.compare("") != 0 ) {
                    auto tmp_prior = abc_object_sim.GetPriorFloat();
                    
                    for ( auto current_vec : tmp_prior ) {
                        used_prior_distrib.push_back(current_vec);
                    }
                }
            }
            
            
            std::cout << "Finished running simulations." << std::endl;
            
            std::cout << "Sorting data from multiple positions... ";
            std::nth_element(accepted_results_vector.begin(),
                             accepted_results_vector.begin() + results_to_accept,
                             accepted_results_vector.end());
            
            accepted_results_vector.erase( accepted_results_vector.begin() + results_to_accept, accepted_results_vector.end() );
            std::cout << "Done." << std::endl;
        }
        else {
            
            // single position
            
            /* ----------------- */
            /*  Run simulations  */
            /* ----------------- */
            clsCMulatorABC abc_object_sim( my_zparams, actual_data_file.GetFirstPosition() );
            
            abc_object_sim.SetImmediateRejection(false);
            
            auto num_batches = my_zparams.GetInt(fits_constants::PARAM_SIM_REPEATS);
            
            if ( num_batches > 10 ) {
                num_batches = 10;
            }
            
            abc_object_sim.RunABCInference(factor, num_batches);
            running_time_sec += abc_object_sim.GetSunningTimeSec();
            tmp_rejection_threshold = abc_object_sim.GetRejectionThreshold();
            
            /* --------------- */
            /*  Store results  */
            /* --------------- */
            std::cout << "Gathering results... ";
            auto tmp_accepted_results_vec = abc_object_sim.GetResultsVector(true);
            for ( auto tmp_result : tmp_accepted_results_vec ) {
                accepted_results_vector.push_back(tmp_result);
            }
            std::cout << "Done." << std::endl;
            
            if ( prior_output_filename.compare("") != 0 ) {
                auto tmp_prior = abc_object_sim.GetPriorFloat();
                
                for ( auto current_vec : tmp_prior ) {
                    used_prior_distrib.push_back(current_vec);
                }
            }
            std::cout << "Finished running simulations." << std::endl << std::endl;
        }
        
        
        /* --------------------- */
        /*  Process simulations  */
        /* --------------------- */
        std::cout << "Processing results (" <<  accepted_results_vector.size() << " simulations):" << std::endl;
        
        ResultsStats result_stats(my_zparams);
        
        // result_stats.SetRejectionThreshold( abc_object_sim.GetRejectionThreshold() );
        result_stats.SetRejectionThreshold( tmp_rejection_threshold );
        
        auto single_mutation_rate = my_zparams.GetFloat( fits_constants::PARAM_SINGLE_MUTATION_RATE, 0.0f );
        result_stats.SetSingleMutrateUsed( single_mutation_rate != 0.0f );
        
        // result_stats.SetRunningTimeSec( abc_object_sim.GetTotalRunningTimeSec() );
        result_stats.SetRunningTimeSec( running_time_sec );
        
        
        
        switch (factor) {
            case Fitness: {
                
                try {
                    result_stats.SetPriorDistrib( my_zparams.GetString(fits_constants::PARAM_PRIOR_DISTRIB, fits_constants::PARAM_PRIOR_DISTRIB_DEFAULT ) );
                    result_stats.CalculateStatsFitness(accepted_results_vector);
                }
                catch (std::string str) {
                    std::cerr << "Exception caught: " << str << std::endl;
                    return 1;
                }
                catch (const char* str) {
                    std::cerr << "Exception caught: " << str << std::endl;
                    return 1;
                }
                catch (std::exception& e) {
                    std::cerr << "Exception caught: " << e.what() << std::endl;
                    return 1;
                }
                catch (...) {
                    throw "Error while calculating stats.";
                }
                
                
                try {
                    std::string tmp_summary_str = result_stats.GetSummaryFitness();
                    
                    std::cout << "Writing posterior... ";
                    result_stats.WriteFitnessDistribToFile(accepted_results_vector, posterior_output_filename);
                    std::cout << "Done." << std::endl;
                
                    std::cout << "Writing summary... ";
                    result_stats.WriteStringToFile(summary_output_filename, tmp_summary_str);
                    std::cout << "Done." << std::endl;
                    
                    if ( prior_output_filename.compare("") != 0 ) {
                        std::cout << "Writing prior... ";
                        
                        // auto tmp_prior = abc_object_sim.GetPriorFloat();
                        // result_stats.WritePriorDistribToFile(tmp_prior, prior_output_filename);
                        
                        result_stats.WritePriorDistribToFile(used_prior_distrib, prior_output_filename);
                        std::cout << "Done." << std::endl;
                    }
                    
                    std::cout << std::endl << "Summary:" << std::endl;
                    std::cout << tmp_summary_str << std::endl;
                }
                catch (std::string str) {
                    std::cerr << "Exception caught: " << str << std::endl;
                    return 1;
                }
                catch (const char* str) {
                    std::cerr << "Exception caught: " << str << std::endl;
                    return 1;
                }
                catch (std::exception& e) {
                    std::cerr << "Exception caught: " << e.what() << std::endl;
                    return 1;
                }
                catch (...) {
                    throw "Error writing results files.";
                }
                
                break;
            }
                
                
            case PopulationSize: {
                
                try {
                    result_stats.SetPriorDistrib( my_zparams.GetString(fits_constants::PARAM_PRIOR_DISTRIB, fits_constants::PARAM_PRIOR_DISTRIB_UNIFORM ) );
                    //result_stats.CalculateStatsFitness(accepted_results_vector);
                    result_stats.CalculateStatsPopulationSize(accepted_results_vector);
                }
                catch (std::string str) {
                    std::cerr << "Exception caught: " << str << std::endl;
                    return 1;
                }
                catch (const char* str) {
                    std::cerr << "Exception caught: " << str << std::endl;
                    return 1;
                }
                catch (std::exception& e) {
                    std::cerr << "Exception caught: " << e.what() << std::endl;
                    return 1;
                }
                catch (...) {
                    throw "Error while calculating stats.";
                }
                
                
                
                try {
                    std::string tmp_summary_str = result_stats.GetSummaryPopSize();
                    
                    
                    std::cout << "Writing posterior... ";
                    result_stats.WritePopSizeDistribToFile(accepted_results_vector, posterior_output_filename);
                    std::cout << "Done." << std::endl;
                
                    std::cout << "Writing summary... ";
                    result_stats.WriteStringToFile(summary_output_filename, tmp_summary_str);
                    std::cout << "Done." << std::endl;
                    
                    if ( prior_output_filename.compare("") != 0 ) {
                        std::cout << "Writing prior... ";
                        
                        // auto tmp_prior = abc_object_sim.GetPriorFloat();
                        // result_stats.WritePriorDistribToFile(tmp_prior, prior_output_filename);
                        
                        result_stats.WritePriorDistribToFile(used_prior_distrib, prior_output_filename);
                        std::cout << "Done." << std::endl;
                    }
                    
                    std::cout << std::endl << "Summary:" << std::endl;
                    std::cout << tmp_summary_str << std::endl;
                }
                catch (const char* str) {
                    std::cerr << "Exception caught: " << str << std::endl;
                    return 1;
                }
                catch (std::string str) {
                    std::cerr << "Exception caught: " << str << std::endl;
                    return 1;
                }
                catch (std::exception& e) {
                    std::cerr << "Exception caught: " << e.what() << std::endl;
                    return 1;
                }
                catch (...) {
                    throw "Error writing results files.";
                }
                
                break;
            }
                
            case MutationRate: {
                
                auto infer_single_mutrate = my_zparams.GetInt(fits_constants::PARAM_MIN_LOG_SINGLE_MUTATION_RATE, 0);
                
                try {
                    result_stats.SetSingleMutrateInferred( infer_single_mutrate != 0 );
                    result_stats.SetPriorDistrib( my_zparams.GetString(fits_constants::PARAM_PRIOR_DISTRIB, fits_constants::PARAM_PRIOR_DISTRIB_UNIFORM ) );
                    result_stats.CalculateStatsMutation(accepted_results_vector);
                }
                catch (std::string str) {
                    std::cerr << "Exception caught: " << str << std::endl;
                    return 1;
                }
                catch (const char* str) {
                    std::cerr << "Exception caught: " << str << std::endl;
                    return 1;
                }
                catch (std::exception& e) {
                    std::cerr << "Exception caught: " << e.what() << std::endl;
                    return 1;
                }
                catch (...) {
                    throw "Error while calculating stats.";
                }
                
                
                try {
                    std::string tmp_summary_str = result_stats.GetSummaryMutRate();
                    
                    std::cout << "Writing posterior... ";
                    result_stats.WriteMutRateDistribToFile(accepted_results_vector, posterior_output_filename);
                    std::cout << "Done." << std::endl;
                    
                    std::cout << "Writing summary... ";
                    result_stats.WriteStringToFile(summary_output_filename, tmp_summary_str);
                    std::cout << "Done." << std::endl;
                    
                    if ( prior_output_filename.compare("") != 0 ) {
                        std::cout << "Writing prior... ";
                        
                        // auto tmp_prior = abc_object_sim.GetPriorFloat();
                        // result_stats.WritePriorDistribToFile(tmp_prior, prior_output_filename);
                        
                        result_stats.WritePriorDistribToFile(used_prior_distrib, prior_output_filename);
                        std::cout << "Done." << std::endl;
                    }
                    
                    std::cout << std::endl << "Summary:" << std::endl;
                    std::cout << tmp_summary_str << std::endl;
                }
                catch (const char* str) {
                    std::cerr << "Exception caught: " << str << std::endl;
                    return 1;
                }
                catch (std::string str) {
                    std::cerr << "Exception caught: " << str << std::endl;
                    return 1;
                }
                catch (std::exception& e) {
                    std::cerr << "Exception caught: " << e.what() << std::endl;
                    return 1;
                }
                catch (...) {
                    throw "Error writing results files.";
                }
                
                break;
            }
        }
    }
    catch (const char* str) {
        std::cerr << "Exception caught while attempting to report stats: " << str << std::endl;
        return 1;
    }
    catch (std::string str) {
        std::cerr << "Exception caught while attempting to report stats: " << str << std::endl;
        return 1;
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




void print_syntaxes(std::string exec_name)
{
    exec_name = "fits ";
    std::cout << "\t" << exec_name << fits_constants::ARG_INFER_FITNESS
    << " <param_file> <actual_data_file> <posterior_file> <summary_file> (optional: <prior_file>)" << std::endl << std::endl;
    
    std::cout << "\t" << exec_name << fits_constants::ARG_INFER_MUTATION
    << " <param_file> <actual_data_file> <posterior_file> <summary_file> (optional: <prior_file>)" << std::endl << std::endl;
    
    std::cout << "\t" << exec_name << fits_constants::ARG_INFER_POPSIZE
    << " <param_file> <actual_data_file> <posterior_file> <summary_file> (optional: <prior_file>)" << std::endl << std::endl;
    
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
            || first_argument.compare( fits_constants::ARG_INFER_POPSIZE ) == 0 );
}

int main(int argc, char* argv[])
{
    if (argc <= 1) {
        print_welcome();
        print_syntaxes(argv[0]);
        return 1;
    }
    
    std::string tmp_first_param = argv[1];
    
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
    
    
    // we should never reach here...
    print_welcome();
    std::cout << "Invalid suntax (option chosen was " << tmp_first_param << "). Use only the following:" << std::endl;
    print_syntaxes(argv[0]);
    return 1;
}

