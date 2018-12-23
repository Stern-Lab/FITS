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

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "fits_constants.h"


struct PriorFileEntry {
    double probability;
    double value;
    
    bool operator<( const PriorFileEntry& other ) const
    {
        return value < other.value;
    }
};


class PriorFile {
  
    std::vector<PriorFileEntry> _prior_entry_vec;
    
public:
    PriorFile( std::string filename );
    
    std::vector<PriorFileEntry> ReadPriorFromFile( std::string filename );
    
    MATRIX_TYPE GetPriorMatrix();
    
    const int PRIOR_FILE_COLUMN_PROBABILITY = 0;
    const int PRIOR_FILE_COLUMN_VALUE = 1;
};

#endif /* PriorFile_hpp */
