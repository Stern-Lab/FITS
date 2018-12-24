//
//  PriorFile.cpp
//  fits
//
//  Created by Tal Zinger on 22/12/2018.
//  Copyright © 2018 Stern Lab. All rights reserved.
//

#include "PriorFile.hpp"


PriorFile::PriorFile( std::string filename ) :
_prior_entry_vec(0)
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
    
    std::vector<PriorFileEntry> prior_entry_vec;
    
    while ( std::getline(infile, tmp_line) ) {
        
        ++current_line_num;
        
        if (is_first_line) {
            
            std::vector<std::string> line_fields;
            try {
                boost::split(line_fields, tmp_line, boost::is_any_of( fits_constants::FILE_FIELD_DELIMITER ));
                expected_num_columns = line_fields.size();
            }
            catch (...) {
                std::string tmp_str = "Error while parsing prior file to columns. Header line:\n" + tmp_line;
                throw tmp_str;
            }
            
            if ( line_fields.size() != 2 ) {
                std::string tmp_str = "Error while reading data file (line " + std::to_string(current_line_num) +
                "): line contains " + std::to_string(line_fields.size())  + " columns, expected " + std::to_string(expected_num_columns) + ".\n";
                
                throw tmp_str;
            }
            
            is_first_line = false;
            continue;
        }
        
        
        std::vector<std::string> line_fields;
        try {
            boost::split(line_fields, tmp_line, boost::is_any_of( fits_constants::FILE_FIELD_DELIMITER ));
        }
        catch (...) {
            std::string tmp_str = "Error while parsing prior file to columns. Line from file:\n" + tmp_line;
            throw tmp_str;
        }
        
        if ( line_fields.size() != 2 ) {
            std::string tmp_str = "Error while reading data file (line " + std::to_string(current_line_num) +
            "): line contains " + std::to_string(line_fields.size())  + " columns, expected " + std::to_string(expected_num_columns) + ".\n";
            
            throw tmp_str;
        }
        
        
        PriorFileEntry tmp_prior_entry;
        try {
            boost::trim( line_fields[PRIOR_FILE_COLUMN_VALUE] );
            boost::trim( line_fields[PRIOR_FILE_COLUMN_PROBABILITY] );
        }
        catch (...) {
            std::string tmp_str = "Error while parsing prior file while trimming (line " + std::to_string(current_line_num) + ")";
            throw tmp_str;
        }
        
        
        try {
            tmp_prior_entry.value = boost::lexical_cast<double>( line_fields[PRIOR_FILE_COLUMN_VALUE] );
        }
        catch (...) {
            std::string tmp_str = "Error while loading prior: line " + std::to_string(current_line_num)+ " does not contain a valid value.";
            throw tmp_str;
        }
        
        try {
            tmp_prior_entry.probability = boost::lexical_cast<double>( line_fields[PRIOR_FILE_COLUMN_PROBABILITY] );
            
            if ( tmp_prior_entry.probability < 0 || tmp_prior_entry.probability > 1 ) {
                std::string tmp_str = "Probability must be between 0 and 1";
                throw tmp_str;
            }
        }
        catch ( std::string str ) {
            std::string tmp_str = "Error while loading prior: line " + std::to_string(current_line_num)+ " does not contain a valid probability: " + str;
            throw tmp_str;
        }
        catch (...) {
            std::string tmp_str = "Error while loading prior: line " + std::to_string(current_line_num)+ " does not contain a valid probability.";
            throw tmp_str;
        }
        
        prior_entry_vec.push_back( tmp_prior_entry );
    }
    
    std::sort( prior_entry_vec.begin(), prior_entry_vec.end() );
    
    if ( prior_entry_vec.empty() ) {
        std::string tmp_str = "No entries found in prior file.";
        throw tmp_str;
    }
    
    for ( auto idx=1; idx<prior_entry_vec.size(); ++idx ) {
        if ( prior_entry_vec[idx].value == prior_entry_vec[idx-1].value ) {
            std::string tmp_str = "Duplicate probabilities in prior file for value " +
            std::to_string(prior_entry_vec[idx].value) +
            "(" + std::to_string(prior_entry_vec[idx-1].probability) +
            " and " + std::to_string(prior_entry_vec[idx-1].probability) + ")";
            
            throw tmp_str;
        }
    }
    
    
    
}
