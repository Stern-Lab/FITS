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

/* This is a simple class for input and output of parameters */

#ifndef ZParams_h
#define ZParams_h

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <regex>
#include <boost/lexical_cast.hpp>

/*
 I assume that reading the parameters will occur only once, so no optimization is
 necessary. Also, I assume this will handle a few dozens parameters at most,
 so no optimization is necessary. Only readability and usability and rubustness.
 */

/* TODO: maybe convert from string to numeric through other than atoi and co. */

class ZParams {
    
public:
    ZParams();
    ZParams( std::string filename, bool ReadOnly );
    ZParams( const ZParams& other );
    
    
    /* Utility */
    void Clear();
    bool IsEmpty() const;
    std::string GetAllParameters() const;
    
    
    /* Data Input */
    // TODO: make all string const by ref?
    void ReadParameters( const std::string filename, bool ReadOnly );
    
    void AddParameter( const std::string paramName, const int value );
    void AddParameter( const std::string paramName, const double value );
    void AddParameter( const std::string paramName, const float value );
    void AddParameter( const std::string paramName, const std::string value );
    void AddParameter( const std::string paramName, const bool value );
    
    void UpdateParameter( const std::string paramName, const int value );
    void UpdateParameter( const std::string paramName, const double value );
    void UpdateParameter( const std::string paramName, const float value );
    void UpdateParameter( const std::string paramName, const std::string value );
    void UpdateParameter( const std::string paramName, const bool value );
    
    bool IsParameterFound( const std::string paramName );
    
    /* Data Output */
    int GetInt(const std::string paramName, const int defaultValue) const;
    float GetFloat(const std::string paramName, const float defaultValue) const;
    double GetDouble(const std::string paramName, const double defaultValue) const;
    std::string GetString(const std::string paramName, const std::string defaultValue) const;
    
    // TODO: unsigned for int, others
    // Note that conversion from singned to unsigned is well defined, so this is actually unnecessary.
    unsigned int GetUnsignedInt(const std::string paramName, const unsigned int defaultValue) const;
    unsigned long GetUnsignedLong(const std::string paramName, const unsigned long defaultValue) const;
    
    
    /* These equivalent functions will throw exceptions instead of giving default values */
    int GetInt(const std::string paramName) const;
    unsigned int GetUnsignedInt(const std::string paramName) const;
    unsigned long GetUnsignedLong(const std::string paramName) const;
    float GetFloat(const std::string paramName) const;
    double GetDouble(const std::string paramName) const;
    std::string GetString(const std::string paramName) const;
    
private:
    
    std::map<std::string, std::string> _param_map;
    bool _is_initialized;
    bool _read_only; // won't accept updates
    
    void WarnIfDecimal( std::string param_name, std::string param_val ) const;
    void WarnIfInteger( std::string param_name, std::string param_val ) const;
    
};
#endif /* ZParams_h */
