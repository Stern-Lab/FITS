//
//  PriorFile.cpp
//  fits
//
//  Created by Tal Zinger on 22/12/2018.
//  Copyright © 2018 Stern Lab. All rights reserved.
//

#include "PriorFile.hpp"

/*
MATRIX_TYPE GetPriorMatrix( std::size_t rows, std::size_t cols )
{
    auto expected_size = rows * cols;
    
    throw "unimplemented";
}
*/


void PriorFile::ReadPriorFromFile( std::string filename )
{
    std::ifstream infile(filename);
    
    if (!infile.is_open()) {
        std::string tmp_str = "error while opening prior file: " + filename;
        throw tmp_str;
    }
    
    std::string tmp_line;
    bool is_first_line = true;
    std::size_t expected_num_columns = 0;
    std::size_t current_line_num = 0;
    
    _prior_distribution.clear();
    
    while ( std::getline(infile, tmp_line) ) {
        
        ++current_line_num;
        
        if (is_first_line) {
            
            std::vector<std::string> line_fields;
            try {
                boost::split(line_fields, tmp_line, boost::is_any_of( fits_constants::FILE_FIELD_DELIMITER ));
                expected_num_columns = line_fields.size();
                
                for ( auto& tmp_field : line_fields ) {
                    boost::trim( tmp_field );
                    _header_titles.push_back( tmp_field );
                }
            }
            catch (...) {
                std::string tmp_str = "Error while parsing prior file to columns. Header line:\n" + tmp_line;
                
                std::cerr << tmp_str << std::endl;
                
                throw tmp_str;
            }
            
            is_first_line = false;
            continue;
        }
        
        
        std::vector<FLOAT_TYPE> tmp_prior_sample;
        std::vector<std::string> line_fields;
        try {
            boost::split(line_fields, tmp_line, boost::is_any_of( fits_constants::FILE_FIELD_DELIMITER ));
            
            if ( line_fields.size() != expected_num_columns ) {
                std::string tmp_str = "line " + std::to_string(current_line_num) +
                " contains " + std::to_string(line_fields.size()) +
                " columns, expected " + std::to_string(expected_num_columns) + ".";
                
                std::cerr << tmp_str << std::endl;
                
                throw tmp_str;
            }
            
            for ( auto tmp_field : line_fields ) {
                boost::trim( tmp_field );
                
                if ( tmp_field.empty() ) {
                    continue;
                }
                
                auto tmp_val = boost::lexical_cast<double>( tmp_field );
                tmp_prior_sample.push_back( tmp_val );
            }
        }
        catch (...) {
            std::string tmp_str = "Error while parsing prior file to columns. Line from file:\n" + tmp_line;
            std::cerr << tmp_str << std::endl;
            throw tmp_str;
        }
        
        _prior_distribution.push_back( tmp_prior_sample );
    }
    
    
    if ( _prior_distribution.empty() ) {
        std::string tmp_str = "No entries found in prior file.";
        throw tmp_str;
    }
    
    _is_initialized = true;
}


PriorFile::PriorFile()
{
    _is_initialized = false;
    _header_titles.clear();
    _prior_distribution.clear();
    
    _rnd_seed = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    _rnd_gen.seed(_rnd_seed);
}


PRIOR_DISTRIB_VECTOR PriorFile::GetPriorAsVector()
{
    if (!_is_initialized) {
        std::string tmp_str = "Prior uninitialized - file not loaded";
        throw tmp_str;
    }
    
    return _prior_distribution;
}


PRIOR_DISTRIB_VECTOR PriorFile::ResamplePriorAsVector( std::size_t num_samples, bool overwrite_internal )
{
    PRIOR_DISTRIB_VECTOR resampled_prior;
    
    if (!_is_initialized) {
        std::string tmp_str = "Prior uninitialized - file not loaded";
        throw tmp_str;
    }
    
    if ( num_samples < 1 ) {
        std::string tmp_str = "Prior resampled with invalid number of times: " + std::to_string(num_samples);
        throw tmp_str;
    }
    
    boost::random::uniform_int_distribution<std::size_t> sample_idx_distribution( 0, _prior_distribution.size() );
    for ( std::size_t sample_counter=0; sample_counter<num_samples; ++sample_counter ) {
        
        auto rnd_idx = sample_idx_distribution(_rnd_gen);
        resampled_prior.push_back( _prior_distribution[rnd_idx] );
    }
    
    if ( overwrite_internal ) {
        _prior_distribution = resampled_prior;
    }
    
    return resampled_prior;
}
