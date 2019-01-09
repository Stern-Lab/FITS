//
//  PriorFile.hpp
//  fits
//
//  Created by Tal Zinger on 22/12/2018.
//  Copyright © 2018 Stern Lab. All rights reserved.
//

#ifndef PriorFile_hpp
#define PriorFile_hpp

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/algorithm/string.hpp> // for splitting
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>


#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>

#include "fits_constants.h"



class PriorFile {

    boost::mt19937 _rnd_gen;
    unsigned int _rnd_seed;
    
    PRIOR_DISTRIB_VECTOR _prior_distribution;
    std::vector<std::string> _header_titles;
    
    bool _is_initialized;
    
public:
    PriorFile();
    
    void ReadPriorFromFile( std::string filename );
    
    PRIOR_DISTRIB_VECTOR GetPriorAsVector();
    
    std::size_t GetPriorSize() { return _prior_distribution.size(); }
    
    PRIOR_DISTRIB_VECTOR ResamplePriorAsVector( std::size_t num_samples, bool overwrite_internal );
    
    // MATRIX_TYPE GetPriorMatrix( std::size_t rows, std::size_t cols );
    
    const int PRIOR_FILE_COLUMN_PROBABILITY = 0;
    const int PRIOR_FILE_COLUMN_VALUE = 1;
    
    bool IsInitialized() { return _is_initialized; }
};

#endif /* PriorFile_hpp */
