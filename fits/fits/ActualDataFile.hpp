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
    //int base;
    int allele;
    FLOAT_TYPE freq;
    int ref;
    int read_count;
    
    static constexpr int EMPTY_INT = -1;
    static constexpr FLOAT_TYPE EMPTY_FLOAT = -1.0f;
    
    ActualDataEntry()
    : pos(EMPTY_INT), gen(EMPTY_INT), allele(EMPTY_INT), freq(EMPTY_FLOAT), ref(EMPTY_INT), read_count(EMPTY_INT) {}
    
    ActualDataEntry( int position, int generation, int allele_num, FLOAT_TYPE frequency, int reference=-1, int reads=1000 )
    : pos(position), gen(generation), allele(allele_num), freq(frequency), ref(reference), read_count(reads) {}
    
    ActualDataEntry( const ActualDataEntry& original )
    : pos(original.pos), gen(original.gen), allele(original.allele), freq(original.freq), ref(original.ref), read_count(original.read_count) {}
    
    const bool SameAllele( const ActualDataEntry& other ) { return allele==other.allele; }
    const bool SameGeneration( const ActualDataEntry &other ) { return gen==other.gen; }
    const bool SameRef( const ActualDataEntry& other ) { return ref==other.ref; }
    const bool SameReadCount( const ActualDataEntry &other ) { return read_count==other.read_count; }
    
    // expected usage - for sorting in ascending generations, sor comparison with simulated data
    bool operator<( const ActualDataEntry& other ) const;
    
    ActualDataEntry& operator=(ActualDataEntry other);
    void swap(ActualDataEntry& other);
    
    
    bool AllDataFilled() { return ( pos > EMPTY_INT && gen > EMPTY_INT && allele > EMPTY_INT && freq > EMPTY_INT && ref > EMPTY_FLOAT && read_count > EMPTY_INT); }
    
};


struct ActualDataPositionData {
    
    static const int NO_POSITION_SPECIFIED = -1;
    
    ActualDataPositionData();
    ActualDataPositionData( const ActualDataPositionData& other );
    
    std::vector<ActualDataEntry> _actual_data;
    
    int _position;
    int _wt_index;
    int _num_alleles;
    
    // to avoid searching the vector each and every time
    std::vector<int> _actual_generations;
    std::vector<FLOAT_TYPE> _actual_frequencies;
    std::vector<FLOAT_TYPE> _init_frequencies;
    
    void Clear();
    
    MATRIX_TYPE GetActualFreqsAsMatrix();
    
    int GetFirstGeneration();
    int GetLastGeneration();
    
    int GetNumberOfAlleles();
    int GetWTIndex();
    
    std::vector<int> GetActualGenerations( bool only_unique = true);
    std::vector<FLOAT_TYPE> GetActualFrequencies();
    std::vector<FLOAT_TYPE> GetInitFreqs();
    
    bool operator<( const ActualDataPositionData& other ) const;
};

class ActualDataFile {
    
    bool _is_initialized;
    bool _multi_positions;
    
    std::vector<ActualDataPositionData> _position_data;
    
    static const int NO_POSITION_SPECIFIED = -1;
    
    // fields for actual data
    const int ACTUAL_DATA_EMPTY_CELL = -1;
    const int ACTUAL_DATA_COLUMN_GENERATION = 0;
    const int ACTUAL_DATA_COLUMN_ALLELE = 1;
    const int ACTUAL_DATA_COLUMN_FREQ = 2;
    const int ACTUAL_DATA_COLUMN_POSITION = 3;
    const int ACTUAL_DATA_MINIMAL_COLS = 3;
    
    void ValidateMultiPosition(int position);
    
    std::vector<ActualDataEntry> DataFileToEntries( std::string filename );
    
public:
    
    ActualDataFile();
    ActualDataFile( const ActualDataFile& other );
    
    void LoadActualData( std::string filename );
    
    //Cannot be different between positions
    int GetNumberOfAlleles();
    //int GetNumberOfAlleles( int position = NO_POSITION_SPECIFIED );
    int GetWTIndex( int position = NO_POSITION_SPECIFIED );
    
    MATRIX_TYPE GetActualFreqsAsMatrix( int position = NO_POSITION_SPECIFIED );
    
    int GetFirstGeneration( int position = NO_POSITION_SPECIFIED );
    int GetLastGeneration( int position = NO_POSITION_SPECIFIED );
    
    ActualDataPositionData GetFirstPosition();
    ActualDataPositionData GetPosition( int position );
    int GetNumberOfPositions();
    std::vector<int> GetPositionNumbers();
    
    
    //std::vector<FLOAT_TYPE> GetSDPerAllele();
    std::vector<int> GetActualGenerations( bool only_unique = true, int position = NO_POSITION_SPECIFIED );
    std::vector<FLOAT_TYPE> GetActualFrequencies( int position = NO_POSITION_SPECIFIED );
    std::vector<FLOAT_TYPE> GetInitFreqs( int position = NO_POSITION_SPECIFIED );
};

#endif /* ActualDataFile_hpp */
