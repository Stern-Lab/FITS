#include "clsCMulatorABC.h"

void clsCMulatorABC::WriteStringToFile( std::string filename, std::string str )
{
    std::ofstream outfile(filename, std::ofstream::out | std::ofstream::trunc);
    
    if (!outfile.is_open()) {
        std::cerr << "unable to open file for writing: " << filename << std::endl;
        throw "unable to open file for writing: " + filename;
    }
    
    outfile << str;
    
    outfile.close();
}

void clsCMulatorABC::WriteSimDataToFile( std::string filename, SimulationResult &res )
{
    std::ofstream outfile(filename, std::ofstream::out | std::ofstream::trunc);
    
    auto tmp_matrix = res.sim_data_matrix;
    
    for ( auto gen=0; gen<tmp_matrix.size1(); ++gen ) {
        for ( auto allele=0; allele<tmp_matrix.size2(); ++allele ) {
            
            outfile << std::to_string( tmp_matrix(gen,allele) ) << "\t";
        }
        outfile << std::endl;
    }
    
    outfile.close();
}

FLOAT_TYPE clsCMulatorABC::GetMedian( std::vector<FLOAT_TYPE> vec )
{
    auto median_idx = std::floor(vec.size() / 2);
    
    std::nth_element( vec.begin(), vec.begin()+median_idx, vec.end() );
    
    if ( vec.size() % 2 == 0 ) {
        return ( vec[median_idx] + vec[median_idx-1] ) / 2;
    }
    
    return vec[median_idx];
}
