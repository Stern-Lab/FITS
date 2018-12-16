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

#include "ActualDataFile.hpp"


void ActualDataEntry::swap(ActualDataEntry& other)
{
    std::swap(gen, other.gen);
    std::swap(allele, other.allele);
    std::swap(freq, other.freq);
    std::swap(ref, other.ref);
    std::swap(read_count, other.read_count);
}


// used for sorting only
bool ActualDataEntry::operator<( const ActualDataEntry& other ) const
{
    // first, sort by position
    //if (pos < other.pos) {
     //   return true;
   // }
    //if (pos > other.pos) {
     //   return false;
    //}
    
    if (gen < other.gen) {
        return true;
    }
    if (gen > other.gen) {
        return false;
    }
    
    // same position, same generation - sort by allele
    return allele < other.allele;
}


ActualDataEntry& ActualDataEntry::operator=(ActualDataEntry other)
{
    swap(other);
    return *this;
}


