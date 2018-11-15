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


/* This class deals with actual/observed frequency data */

#ifndef ActualDataFile_hpp
#define ActualDataFile_hpp


#include <string>
#include <vector>

#include <boost/algorithm/string.hpp> // for splitting

#include <boost/numeric/ublas/matrix.hpp>

#include <boost/accumulators/numeric/functional.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/framework/features.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/format.hpp>

#include "CMulator.h"

struct ActualDataEntry {
    int pos;
    int gen;
    int base;
    FLOAT_TYPE freq;
    int ref;
    int read_count;
    
    static constexpr int EMPTY_INT = -1;
    static constexpr FLOAT_TYPE EMPTY_FLOAT = -1.0;
    
    ActualDataEntry()
    : pos(EMPTY_INT), gen(EMPTY_INT), base(EMPTY_INT), freq(EMPTY_FLOAT), ref(EMPTY_INT), read_count(EMPTY_INT) {}
    
    ActualDataEntry( int position, int generation, int base_num, FLOAT_TYPE frequency, int reference=-1, int reads=1000 )
    : pos(position), gen(generation), base(base_num), freq(frequency), ref(reference), read_count(reads) {}
    
    ActualDataEntry( const ActualDataEntry& original )
    : pos(original.pos), gen(original.gen), base(original.base), freq(original.freq), ref(original.ref), read_count(original.read_count) {}
    
    const bool SameBase( const ActualDataEntry& other ) { return base==other.base; }
    const bool SameGeneration( const ActualDataEntry &other ) { return gen==other.gen; }
    const bool SameRef( const ActualDataEntry& other ) { return ref==other.ref; }
    const bool SameReadCount( const ActualDataEntry &other ) { return read_count==other.read_count; }
    
    // expected usage - for sorting in ascending generations, sor comparison with simulated data
    bool operator<( const ActualDataEntry& other ) const;
    
    ActualDataEntry& operator=(ActualDataEntry other);
    void swap(ActualDataEntry& other);
    
    
    bool AllDataFilled() { return (pos > EMPTY_INT && gen > EMPTY_INT && base > EMPTY_INT && freq > EMPTY_INT && ref > EMPTY_FLOAT && read_count > EMPTY_INT); }
    
};


class ActualDataFile {
    
    bool _is_initialized;
    
public:
    
    ActualDataFile();
    ActualDataFile( const ActualDataFile& other );
    
    void LoadActualData( std::string filename );
    
    int GetNumberOfAlleles();
    int GetWTIndex();
    
    std::vector<int> GetActualGenerations(bool only_unique = true);
    std::vector<FLOAT_TYPE> GetActualFrequencies();
    MATRIX_TYPE GetActualFreqsAsMatrix();
    
    int GetFirstGeneration();
    int GetLastGeneration();
    
    std::vector<FLOAT_TYPE> GetInitFreqs();
    
    //std::vector<FLOAT_TYPE> GetSDPerAllele();

    // fields for actual data
    const int ACTUAL_DATA_EMPTY_CELL = -1;
    const int ACTUAL_DATA_COLUMN_GENERATION = 0;
    const int ACTUAL_DATA_COLUMN_ALLELE = 1;
    const int ACTUAL_DATA_COLUMN_FREQ = 2;
    const int ACTUAL_DATA_COLUMNS = 3;

private:
    std::vector<ActualDataEntry> _actual_data;
    
    // to avoid searching the vector each and every
    std::vector<int> _actual_generations;
    std::vector<FLOAT_TYPE> _actual_frequencies;
    std::vector<FLOAT_TYPE> _init_frequencies;
    int _wt_index;
    int _num_alleles;
};

#endif /* ActualDataFile_hpp */
