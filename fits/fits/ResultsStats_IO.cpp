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

#include "ResultsStats.hpp"

/*
std::string ResultsStats::GetFitnessDistrib(const std::vector<SimulationResult>& result_vector)
{
    std::stringstream ss;
    
    // header
    //ss << "sim_id" << "\t" << "distance";
    ss << "distance";
    for (auto i = 0; i < result_vector[0].fitness_values.size(); i++) {
        ss << fits_constants::FILE_FIELD_DELIMITER << "allele" << i;
    }
    ss << std::endl;
    
    
    // 2016-09-04 adding sanity check to report problem
    if ( result_vector.empty() ) {
        throw "Error in GetFitnessDistrib: result_vector is empty.";
    }
    
    for ( auto tmp_entry : result_vector ) {
        
        // ss << tmp_entry.sim_id << "\t" << boost::format("%-10.3d") % tmp_entry.distance_from_actual;
        //ss << boost::format("%-10.3d") % tmp_entry.distance_from_actual;
        ss << tmp_entry.distance_from_actual;
        
        for ( auto tmpval : tmp_entry.fitness_values ) {
            ss << fits_constants::FILE_FIELD_DELIMITER << tmpval;
        }
        
        ss << std::endl;
    }
    
    return ss.str();
}
*/

/*
void ResultsStats::WriteFitnessDistribToFile(const std::vector<SimulationResult>& result_vector, std::string filename)
{
    std::ofstream outfile(filename, std::ofstream::out | std::ofstream::trunc);
    
    if (!outfile.is_open()) {
        std::cerr << "unable to open file for writing: " << filename << std::endl;
        throw "unable to open file for writing: " + filename;
    }
    
    // header
    outfile << fits_constants::FILE_FIELD_DELIMITER << "distance";
    for (auto i = 0; i < result_vector[0].fitness_values.size(); i++) {
        outfile << fits_constants::FILE_FIELD_DELIMITER << "allele" << i;
    }
    outfile << std::endl;
    
    
    // 2016-09-04 adding sanity check to report problem
    if ( result_vector.empty() ) {
        std::cerr << "Error in WriteToFile: result_vector is empty. Filename: " << filename << std::endl;
        std::string my_err = "";
        my_err = "Error in WriteToFile: result_vector is empty. Filename: ";
        my_err += filename;
        throw my_err.c_str();
    }
    
    for ( auto tmp_entry : result_vector ) {
        
        
        //outfile << tmp_entry.sim_id << "\t" << tmp_entry.distance_from_actual;
        //outfile << fits_constants::FILE_FIELD_DELIMITER << boost::format("%-10.3d") % tmp_entry.distance_from_actual;
        outfile << fits_constants::FILE_FIELD_DELIMITER << tmp_entry.distance_from_actual;
        
        for ( auto tmpval : tmp_entry.fitness_values ) {
            outfile << fits_constants::FILE_FIELD_DELIMITER << tmpval;
        }
        
        outfile << std::endl;
    }
    
    outfile.close();
}
*/

/*
std::string ResultsStats::GetMutrateDistrib(const std::vector<SimulationResult>& result_vector)
{
    std::stringstream ss;
    
    // header
    //ss << "sim_id" << fits_constants::FILE_FIELD_DELIMITER << "distance";
    ss << fits_constants::FILE_FIELD_DELIMITER << "distance";
    
    for (auto i = 0; i < result_vector[0].fitness_values.size(); ++i) {
        for (auto j = 0; j < result_vector[0].fitness_values.size(); ++j) {
            ss << fits_constants::FILE_FIELD_DELIMITER << "allele" << i << "_" << j;
        }
    }
    ss << std::endl;
    
    
    // 2016-09-04 adding sanity check to report problem
    if ( result_vector.empty() ) {
        std::cerr << "Error in GetMutrateDistribc: result_vector is empty." << std::endl;
        std::string my_err = "";
        my_err = "Error in GetMutrateDistribc: result_vector is empty. ";
        throw my_err.c_str();
    }
    
    for ( auto tmp_entry : result_vector ) {
        
        // ss << tmp_entry.sim_id << "\t" << tmp_entry.distance_from_actual;
        ss << fits_constants::FILE_FIELD_DELIMITER << tmp_entry.distance_from_actual;
        
        for (auto row = 0; row < result_vector[0].mutation_rates.size1(); ++row) {
            for (auto col = 0; col < result_vector[0].mutation_rates.size2(); ++col) {
                
                
                //ss << fits_constants::FILE_FIELD_DELIMITER << boost::format("%.3d") % tmp_entry.mutation_rates(row,col);
                //ss << fits_constants::FILE_FIELD_DELIMITER << boost::format("*%.2e") % tmp_entry.mutation_rates(row,col);
                ss << fits_constants::FILE_FIELD_DELIMITER << tmp_entry.mutation_rates(row,col);
            }
        }
        ss << std::endl;
    }
    
    return ss.str();
}
*/

/*
void ResultsStats::WriteMutRateDistribToFile(const std::vector<SimulationResult>& result_vector, std::string filename)
{
    std::ofstream outfile(filename, std::ofstream::out | std::ofstream::trunc);
    
    if (!outfile.is_open()) {
        std::cerr << "unable to open file for writing: " << filename << std::endl;
        throw "unable to open file for writing: " + filename;
    }
    
    // header
    outfile << "distance";
    
    for (auto i = 0; i < result_vector[0].fitness_values.size(); ++i) {
        for (auto j = 0; j < result_vector[0].fitness_values.size(); ++j) {
            outfile << fits_constants::FILE_FIELD_DELIMITER << "allele" << i << "_" << j;
        }
    }
    outfile << std::endl;
    
    
    // 2016-09-04 adding sanity check to report problem
    if ( result_vector.empty() ) {
        std::cerr << "Error in WriteToFile: result_vector is empty. Filename: " << filename << std::endl;
        std::string my_err = "";
        my_err = "Error in WriteToFile: result_vector is empty. Filename: ";
        my_err += filename;
        throw my_err.c_str();
    }
    
    for ( auto tmp_entry : result_vector ) {
        
        outfile << tmp_entry.distance_from_actual;
        
        for (auto row = 0; row < result_vector[0].mutation_rates.size1(); ++row) {
            for (auto col = 0; col < result_vector[0].mutation_rates.size2(); ++col) {
                
                //outfile << fits_constants::FILE_FIELD_DELIMITER << boost::format("%.2e") % tmp_entry.mutation_rates(row,col);
                outfile << fits_constants::FILE_FIELD_DELIMITER << tmp_entry.mutation_rates(row,col);
            }
        }
        outfile << std::endl;
    }
    
    outfile.close();
}
*/

/*
std::string ResultsStats::GetPopsizeDistrib(const std::vector<SimulationResult>& result_vector)
{
    std::stringstream ss;
    
    
    // header
    ss << "distance" << fits_constants::FILE_FIELD_DELIMITER << "N" << std::endl;
    
    
    // 2016-09-04 adding sanity check to report problem
    if ( result_vector.empty() ) {
        std::cerr << "Error in GetPopsizeDistrib: result_vector is empty." << std::endl;
        std::string my_err = "";
        my_err = "Error in GetPopsizeDistrib: result_vector is empty. Filename: ";
        throw my_err.c_str();
    }
    
    for ( auto tmp_entry : result_vector ) {
        ss << tmp_entry.distance_from_actual << fits_constants::FILE_FIELD_DELIMITER
        //<< boost::format("%e") % tmp_entry.N << std::endl;
        << tmp_entry.N << std::endl;
    }
    
    return ss.str();
}
*/

/*
void ResultsStats::WritePopSizeDistribToFile(const std::vector<SimulationResult>& result_vector, std::string filename)
{
    std::ofstream outfile(filename, std::ofstream::out | std::ofstream::trunc);
    
    if (!outfile.is_open()) {
        std::cerr << "unable to open file for writing: " << filename << std::endl;
        throw "unable to open file for writing: " + filename;
    }
    
    // header
    outfile << "distance" << fits_constants::FILE_FIELD_DELIMITER << "N" << std::endl;
    
    
    // 2016-09-04 adding sanity check to report problem
    if ( result_vector.empty() ) {
        std::cerr << "Error in WriteToFile: result_vector is empty. Filename: " << filename << std::endl;
        std::string my_err = "";
        my_err = "Error in WriteToFile: result_vector is empty. Filename: ";
        my_err += filename;
        throw my_err.c_str();
    }
    
    for ( auto tmp_entry : result_vector ) {
        outfile << tmp_entry.distance_from_actual << fits_constants::FILE_FIELD_DELIMITER
                // << boost::format("%e") % tmp_entry.N << std::endl;
        << tmp_entry.N << std::endl;
    }
    
    outfile.close();
}
*/


void ResultsStats::WriteStringToFile( std::string filename, std::string str )
{
    std::ofstream outfile(filename, std::ofstream::out | std::ofstream::trunc);
    
    if (!outfile.is_open()) {
        std::cerr << "unable to open file for writing: " << filename << std::endl;
        throw "unable to open file for writing: " + filename;
    }
    
    outfile << str;
    
    outfile.close();
}

std::string ResultsStats::GetPosterior( bool is_multi_position, FactorToInfer factor, const std::vector<SimulationResult>& accepted_results_vec, const std::vector<SimulationResult>& all_results_vec )
{
    if ( accepted_results_vec.empty() ) {
        std::string tmp_str = "Error in WriteMultiPositionPosterior: accepted results vector is empty.";
        throw tmp_str;
    }
    
    /* 2019-01-03 No longer keeping per-position posterior to save memory (8>GB for 10 positions)
    if ( all_results_vec.empty() && is_multi_position ) {
        std::string tmp_str = "Error in WriteMultiPositionPosterior: all-results vector is empty.";
        throw tmp_str;
    }
    */
    
    std::stringstream ss;
    
    // header
    ss << "position" << fits_constants::FILE_FIELD_DELIMITER << "distance";
    switch (factor) {
        case Fitness: {
            for (auto i = 0; i < accepted_results_vec[0].fitness_values.size(); ++i) {
                ss << fits_constants::FILE_FIELD_DELIMITER << "allele" << i;
            }
            break;
        }
            
        case MutationRate: {
            for (auto i = 0; i < accepted_results_vec[0].fitness_values.size(); ++i) {
                for (auto j = 0; j < accepted_results_vec[0].fitness_values.size(); ++j) {
                    ss << fits_constants::FILE_FIELD_DELIMITER << "allele" << i << "_" << j;
                }
            }
            break;
        }
            
        case PopulationSize: {
            ss << fits_constants::FILE_FIELD_DELIMITER << "N";
            break;
        }
    }
    ss << std::endl;
    
    
    // in any case accepted results will be the accepted results..
    for ( auto tmp_entry : accepted_results_vec ) {
        
        // position of -1 to mark the "total" or "global" position
        ss << "N/A" << fits_constants::FILE_FIELD_DELIMITER;
        
        if ( is_multi_position ) {
            // note this is the sum of distances
            ss << tmp_entry.sum_distance;
        }
        else {
            ss << tmp_entry.distance_from_actual;
        }
        
        
        switch (factor) {
            case Fitness: {
                for ( auto tmpval : tmp_entry.fitness_values ) {
                    ss << fits_constants::FILE_FIELD_DELIMITER << tmpval;
                }
                break;
            }
                
            case MutationRate: {
                for (auto row = 0; row < tmp_entry.mutation_rates.size1(); ++row) {
                    for (auto col = 0; col < tmp_entry.mutation_rates.size2(); ++col) {
                        ss << fits_constants::FILE_FIELD_DELIMITER << tmp_entry.mutation_rates(row,col);
                    }
                }
                break;
            }
                
            case PopulationSize: {
                ss << fits_constants::FILE_FIELD_DELIMITER << tmp_entry.N << std::endl;
                break;
            }
        }
        ss << std::endl;
    } // end of accepted results
    
    
    return ss.str();
    
    /* 2019-01-03 No longer keeping per-position posterior to save memory (8>GB for 10 positions)
    
    // only single position, we're done writing
    if ( !is_multi_position ) {
        return ss.str();
    }
    
    
    // now write per-position data
    for ( auto tmp_entry : all_results_vec ) {
        
        ss << tmp_entry.pos << fits_constants::FILE_FIELD_DELIMITER;
        
        ss << tmp_entry.distance_from_actual;
        
        switch (factor) {
            case Fitness: {
                for ( auto tmpval : tmp_entry.fitness_values ) {
                    ss << fits_constants::FILE_FIELD_DELIMITER << tmpval;
                }
                break;
            }
                
            case MutationRate: {
                for (auto row = 0; row < tmp_entry.mutation_rates.size1(); ++row) {
                    for (auto col = 0; col < tmp_entry.mutation_rates.size2(); ++col) {
                        ss << fits_constants::FILE_FIELD_DELIMITER << tmp_entry.mutation_rates(row,col);
                    }
                }
                break;
            }
                
            case PopulationSize: {
                ss << fits_constants::FILE_FIELD_DELIMITER << tmp_entry.N << std::endl;
                break;
            }
        }
        ss << std::endl;
    } // end of all results
    */
    // return ss.str();
}


void ResultsStats::WritePosterior( bool is_multi_position, FactorToInfer factor, const std::vector<SimulationResult>& accepted_results_vec, const std::vector<SimulationResult>& all_results_vec, std::string filename )
{
    std::ofstream outfile(filename, std::ofstream::out | std::ofstream::trunc);
    
    if (!outfile.is_open()) {
        std::string tmp_str = "Error in WriteMultiPositionPosterior: unable to open file " + filename;
        throw tmp_str;
    }
    
    try {
        auto str_posterior = GetPosterior( is_multi_position, factor, accepted_results_vec, all_results_vec );
        outfile << str_posterior;
    }
    catch ( std::string str ) {
        outfile.close();
        throw str;
    }
    catch ( ... ) {
        outfile.close();
        std::string tmp_str = "Unknown error while getting posterior";
        throw tmp_str;
    }
    
    outfile.close();
}

void ResultsStats::WritePriorDistribToFile( FactorToInfer factor_to_infer, const PRIOR_DISTRIB_VECTOR& prior_distrib, std::string filename )
{
    std::ofstream outfile(filename, std::ofstream::out | std::ofstream::trunc);
    
    if (!outfile.is_open()) {
        std::string tmp_str = "Error in WritePriorDistribToFile: unable to open file " + filename;
        throw tmp_str;
    }
    
    try {
        auto str_prior = GetPrior( factor_to_infer, prior_distrib );
        outfile << str_prior;
    }
    catch ( std::string str ) {
        outfile.close();
        throw str;
    }
    catch ( ... ) {
        outfile.close();
        std::string tmp_str = "Unknown error while getting prior";
        throw tmp_str;
    }
    
    outfile.close();
}

std::string ResultsStats::GetPrior( FactorToInfer factor_to_infer, const PRIOR_DISTRIB_VECTOR& prior_distrib )
// void ResultsStats::WritePriorDistribToFile( FactorToInfer factor_to_infer, const PRIOR_DISTRIB_VECTOR& prior_distrib, std::string filename )
{
    std::stringstream ss;
    
    if ( prior_distrib.empty() ) {
        std::string tmp_str = "Error in getting prior: result_vector is empty.";
        throw tmp_str;
    }
    
    // print title according to the prior type
    switch (factor_to_infer) {
        case Fitness: {
            
            for (auto i = 0; i < prior_distrib[0].size(); ++i) {
                ss << "allele" << i;
                
                if ( i< prior_distrib[0].size()-1 ) {
                    ss << fits_constants::FILE_FIELD_DELIMITER;
                }
            }
            
            
            ss << std::endl;
            break;
        }
            
        case MutationRate: {
            
            auto num_alleles = prior_distrib[0].size() / 2;
            
            for (auto i = 0; i < prior_distrib[0].size(); ++i) {
                auto from_allele = i / num_alleles;
                auto to_allele = i % num_alleles;
                
                ss << "allele" << from_allele << "_" << to_allele;
                
                if ( i< prior_distrib[0].size()-1 ) {
                    ss << fits_constants::FILE_FIELD_DELIMITER;
                }
            }
            ss << std::endl;
            break;
        }
            
        case PopulationSize: {
            
            ss << "N" << std::endl;
            break;
        }
    }
    
    for ( auto current_vec : prior_distrib ) {
        
        for ( auto i=0; i<current_vec.size(); ++i ) {
            
            auto current_val = current_vec[i];
            ss << current_val;
            
            if ( i< current_vec.size()-1 ) {
                ss << fits_constants::FILE_FIELD_DELIMITER;
            }
        }
        ss << std::endl;
    }
    
    return ss.str();
}
