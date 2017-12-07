//
//  ActualDataFile_Entry.cpp
//  fits
//
//  Created by Tal Zinger on 16/11/2016.
//  Copyright © 2016 Stern Lab. All rights reserved.
//

#include "ActualDataFile.hpp"


void ActualDataEntry::swap(ActualDataEntry& other)
{
    std::swap(pos, other.pos);
    std::swap(gen, other.gen);
    std::swap(base, other.base);
    std::swap(freq, other.freq);
    std::swap(ref, other.ref);
    std::swap(read_count, other.read_count);
}


// used for sorting only
bool ActualDataEntry::operator<( const ActualDataEntry& other ) const
{
    if (gen < other.gen) {
        return true;
    }
    
    if (gen > other.gen) {
        return false;
    }
    
    // we reached here, so gen == other.gen
    return base < other.base;
}


ActualDataEntry& ActualDataEntry::operator=(ActualDataEntry other)
{
    swap(other);
    return *this;
}


