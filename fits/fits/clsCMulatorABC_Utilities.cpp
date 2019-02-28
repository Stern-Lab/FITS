/*
    FITS - Flexible Inference from Time-Series data
    (c) 2016-2019 by Tal Zinger
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

#include "clsCMulatorABC.h"

void clsCMulatorABC::WriteStringToFile( std::string filename, std::string str )
{
    std::ofstream outfile(filename, std::ofstream::out | std::ofstream::trunc);
    
    if (!outfile.is_open()) {
        std::string tmp_str = "unable to open file for writing: " + filename;
        throw tmp_str;
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
            
            outfile << std::to_string( tmp_matrix(gen,allele) ) << fits_constants::FILE_FIELD_DELIMITER;
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
