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

