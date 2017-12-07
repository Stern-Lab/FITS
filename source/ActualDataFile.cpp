//
//  ActualDataFile.cpp
//  fits
//
//  Created by Tal Zinger on 16/11/2016.
//  Copyright © 2016 Stern Lab. All rights reserved.
//

#include "ActualDataFile.hpp"


ActualDataFile::ActualDataFile() :
_actual_data(0),
_is_initialized(false),
_actual_generations(0),
_actual_frequencies(0),
_init_frequencies(0),
_wt_index(-1),
_num_alleles(-1)
{}


ActualDataFile::ActualDataFile( const ActualDataFile& other ) :
_actual_data(other._actual_data),
_is_initialized(other._is_initialized),
_actual_generations(other._actual_generations),
_actual_frequencies(other._actual_frequencies),
_init_frequencies(other._init_frequencies),
_wt_index(other._wt_index),
_num_alleles(other._num_alleles)
{}


// based on code taken from the GUI
// 2016-04-12
void ActualDataFile::LoadActualData( std::string filename )
{
    _actual_data.clear();
    
    _actual_generations.clear();
    _actual_frequencies.clear();
    _init_frequencies.clear();
    _wt_index = -1;
    _num_alleles = -1;
    
    std::ifstream infile(filename);
    
    if (!infile.is_open()) {
        std::cerr << "error while opening actual data file: " << filename << std::endl;
        throw "error while opening actual data file";
    }
    
    std::string tmp_line;
    bool is_first_line = true;
    
    
    // todo: check if newline needs to be normalized
    while (std::getline(infile, tmp_line)) {
        
        if (is_first_line) {
            is_first_line = false;
            continue;
        }
        
        
        std::vector<std::string> line_fields;
        try {
            boost::split(line_fields, tmp_line, boost::is_any_of("\t"));
            
        }
        catch (std::exception& e) {
            std::cerr << "Error while parsing actual data file to columns. Line: " << std::endl << tmp_line << std::endl;
            throw e;
        }
        catch (...) {
            std::cerr << "Unknown error while parsing actual data file to columns. Line: " << std::endl << tmp_line << std::endl;
            std::cerr << "NOT TERMINATING" << std::endl;
        }
        
        
        if (line_fields.size() <=1) {
            std::cout << "skipping line" << std::endl;
            continue;
        }
        
        // 2016-11-3 had enough of useleff information being read so not forcing all columns to be present
        // except generation, allele, frequency.
        if (line_fields.size() < ACTUAL_DATA_COLUMNS) {
            std::cerr << "incompatible actual data file with " << line_fields.size() << " columns, expected " << ACTUAL_DATA_COLUMNS << std::endl;
            throw "incompatible actual data file.";
        }
        
        
        ActualDataEntry tmp_data_entry;
        
        try {
            boost::trim(line_fields[ACTUAL_DATA_COLUMN_GENERATION]);
            boost::trim(line_fields[ACTUAL_DATA_COLUMN_ALLELE]);
            boost::trim(line_fields[ACTUAL_DATA_COLUMN_FREQ]);
        }
        catch (...) {
            std::cerr << "Error while parsing actual data file (trim)." << std::endl;
            throw "Error while pasring actual data file (trim).";
        }
        
        try {
            tmp_data_entry.gen = boost::lexical_cast<int>(line_fields[ACTUAL_DATA_COLUMN_GENERATION]);
            tmp_data_entry.base = boost::lexical_cast<int>(line_fields[ACTUAL_DATA_COLUMN_ALLELE]);
            tmp_data_entry.freq = boost::lexical_cast<FLOAT_TYPE>(line_fields[ACTUAL_DATA_COLUMN_FREQ]);
            tmp_data_entry.pos = -1;
            tmp_data_entry.ref = -1;
            tmp_data_entry.read_count = -1;
        }
        catch (...) {
            std::cerr << "Error while parsing actual data file (cast)." << std::endl;
            throw "Error while pasring actual data file (cast).";
        }
        
        
        try {
            _actual_data.push_back(tmp_data_entry);
        }
        catch (...) {
            std::cerr << "Error while adding actual data entry" << std::endl;
            throw "Error while adding actual data entry";
        }
    }
    
    try {
        std::sort(_actual_data.begin(), _actual_data.end());
    }
    catch (...) {
        std::cerr << "Error while sorting actual data." << std::endl;
        throw "Error while sorting actual data.";
    }
    
    infile.close();
    
    std::cout << "Done." << std::endl;
    
    _is_initialized = true;
}


std::vector<int> ActualDataFile::GetActualGenerations(bool only_unique)
{
    if (!_is_initialized) {
        std::cerr << "GetActualGenerations - object uninitialized" << std::endl;
    }
    
    if ( _actual_generations.size() > 0 ) {
        return _actual_generations;
    }
    
    
    std::vector<int> tmp_gen_vec;
    
    for (auto current_entry : _actual_data) {
        tmp_gen_vec.push_back(current_entry.gen);
    }
    
    // keep only unique values
    // sort because unique works on censecutive repeats
    // erase because no deletions occur but only running over existing data
    // no need to sort - done when loading
    //std::sort(tmp_gen_vec.begin(), tmp_gen_vec.end());
    
    if (only_unique) {
        auto last = std::unique(tmp_gen_vec.begin(), tmp_gen_vec.end());
        tmp_gen_vec.erase(last, tmp_gen_vec.end());
    }
    
    _actual_generations = tmp_gen_vec;
    return tmp_gen_vec;
}


// get vector of the frequencies of the alleles in the first generation
// assumes this is sorted
std::vector<FLOAT_TYPE> ActualDataFile::GetInitFreqs()
{
    if (!_is_initialized) {
        std::cerr << "GetInitFreqs - object uninitialized" << std::endl;
    }
    
    if ( _init_frequencies.size() > 0 ) {
        return _init_frequencies;
    }
    
    
    auto first_generation = _actual_data[0].gen;
    std::vector<FLOAT_TYPE> tmp_freqs(0);
    
    for (auto current_entry : _actual_data) {
        if (current_entry.gen == first_generation) {
            tmp_freqs.push_back(current_entry.freq);
        }
    }
    
    _init_frequencies = tmp_freqs;
    return tmp_freqs;
}

std::vector<FLOAT_TYPE> ActualDataFile::GetActualFrequencies()
{
    if (!_is_initialized) {
        std::cerr << "GetActualFrequencies - object uninitialized" << std::endl;
    }
    
    if ( _actual_frequencies.size() > 0 ) {
        return _actual_frequencies;
    }
    std::vector<FLOAT_TYPE> tmp_freq_vec;
    
    
    for (auto current_entry : _actual_data) {
        tmp_freq_vec.push_back(current_entry.freq);
    }
    
    _actual_frequencies = tmp_freq_vec;
    return tmp_freq_vec;
}


int ActualDataFile::GetFirstGeneration()
{
    if (!_is_initialized) {
        std::cerr << "GetFirstGeneration - object uninitialized" << std::endl;
    }
    
    // not need to sort, it's done when loaded
    //std::sort(_actual_data.begin(), _actual_data.end());
    
    auto ret_gen = _actual_data[0].gen;
    return ret_gen;
}

int ActualDataFile::GetLastGeneration()
{
    if (!_is_initialized) {
        std::cerr << "GetLastGeneration - object uninitialized" << std::endl;
    }
    
    // not need to sort, it's done when loaded
    //std::sort(_actual_data.begin(), _actual_data.end());
    
    auto ret_gen = _actual_data[ _actual_data.size() - 1 ].gen;
    return ret_gen;
}


int ActualDataFile::GetNumberOfAlleles()
{
    if ( _actual_data.empty() ) {
        throw "Get number of alleles - actual data vector is empty.";
    }
    
    if ( _actual_data.size() <= 1 ) {
        throw "Get number of alleles - actual data vector contains only 1 entry.";
    }
    
    if ( _num_alleles > -1 ) {
        return _num_alleles;
    }
    
    // not need to sort, it's done when loaded
    // std::sort(_actual_data.begin(), _actual_data.end());
    
    auto max_known_allele = 0;
    auto allele_ascending = true;
    
    while ( max_known_allele < _actual_data.size() && allele_ascending ) {
        ++max_known_allele;
        allele_ascending =
        _actual_data[max_known_allele].base > _actual_data[max_known_allele-1].base;
    
    }
    
    _num_alleles = max_known_allele;
    return max_known_allele;
}


MATRIX_TYPE ActualDataFile::GetActualFreqsAsMatrix()
{
    auto num_alleles = GetNumberOfAlleles();
    auto actual_generations = GetActualGenerations();
    auto num_generations = actual_generations.size();
    
    MATRIX_TYPE freq_matrix( num_generations, num_alleles );
    
    auto actual_freqs = GetActualFrequencies();
    
    for ( auto i=0; i<actual_freqs.size(); ++i ) {
        
        auto current_row = i / num_alleles;
        auto current_col = i %  num_alleles;
        
        freq_matrix(current_row, current_col) = actual_freqs[i];
    }
    
    return freq_matrix;
}


int ActualDataFile::GetWTIndex()
{
    if ( _wt_index > -1 ) {
        return _wt_index;
    }
    
    auto first_generation = _actual_data[0].gen;
    
    auto current_wt_idx = -1;
    FLOAT_TYPE current_wt_freq = -1.0f;
    
    for (auto current_entry : _actual_data) {
        if ( current_entry.gen == first_generation ) {
            if (current_entry.freq > current_wt_freq) {
                current_wt_idx = current_entry.base;
                current_wt_freq = current_entry.freq;
            }
        }
    }
    
    _wt_index = current_wt_idx;
    return current_wt_idx;
}
