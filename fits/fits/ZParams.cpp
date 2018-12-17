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

#include "ZParams.h"


// ctor
ZParams::ZParams() :
_param_map(),
_is_initialized(false),
_read_only(false)
{}


// copy ctor
ZParams::ZParams( const ZParams &other ) :
_param_map(other._param_map),
_is_initialized(other._is_initialized),
_read_only(other._read_only)
{}

void ZParams::Clear()
{
    if ( _read_only ) {
        std::string err_msg = "ZParams: Read only. Cannot clear.";
        throw err_msg.c_str();
    }
    
    _param_map.clear();
}

void ZParams::ReadParameters( const std::string filename, bool ReadOnly = false )
{
    _read_only = ReadOnly;
    
    std::ifstream f ( filename );
    
    if ( !f.is_open() ) {
        std::string err_msg = "ZParams: error while opening file " + filename + "\n";
        throw err_msg.c_str();
    }
    
    
    std::regex parameter_regex{"(^\\S+)\\s(\\S+)"}; // (name) (value) pairs, ignore comments with #
    std::regex comment_regex{"^[#]"}; // (name) (value) pairs, ignore comments with #
    std::string line;
    
    std::size_t line_count = 0;
    std::size_t parameter_count = 0;
    
    while(getline(f, line)) {
        
        ++line_count;
        
        // match regex
        std::smatch matches; // matched strings go here
        
        // empty line
        if ( line.empty() ) {
            continue;
        }
        
        // a comment is found
        if ( std::regex_search(line , matches, comment_regex) ) {
            continue;
        }
        
        // parameter line
        if ( std::regex_search(line , matches, parameter_regex) ) {
            // 0 is the whole line
            std::string param_name = matches[1];
            std::string param_val_str = matches[2];
            // std::cout << "name: " << param_name << " val: " << param_val_str << std::endl;
            AddParameter(param_name, param_val_str);
            ++parameter_count;
            continue;
        }
        
        std::string err_msg = "ZParams: error while reading file " + filename + ": the following line (" +  std::to_string(line_count) + ") is not a valid parameter or comment:\n" + line + "\n";
        throw err_msg.c_str();
    }
    
    if (f.bad()) {
        std::string err_msg = "ZParams: error while reading file " + filename + "\n";
        throw err_msg.c_str();
    }
    
    if ( parameter_count == 0 ) {
        std::string err_msg = "ZParams: error while reading file " + filename + " - No valid parameters found.\n";
        throw err_msg.c_str();
    }
    
    _is_initialized = true;
}


ZParams::ZParams( const std::string filename, bool ReadOnly ) :
_param_map(),
_is_initialized(false),
_read_only(ReadOnly)
{
    ReadParameters( filename, ReadOnly );
}


void ZParams::AddParameter( const std::string paramName, const std::string value )
{
    if ( _read_only && _is_initialized ) {
        std::string err_msg = "ZParams: Read only. Cannot add parameters.";
        throw err_msg.c_str();
    }
    
    std::pair<std::string, std::string> tmp_pair(paramName, value);
    
    try {
        _param_map.insert( tmp_pair );
    }
    catch (...) {
        std::string err_msg = "ZParams: Error while adding parameter ";
        err_msg += paramName;
        throw err_msg.c_str();
    }
}


void ZParams::AddParameter( const std::string paramName, const int value )
{
    std::string tmp_str = std::to_string(value);
    
    AddParameter(paramName, tmp_str);
}


void ZParams::AddParameter( const std::string paramName, const double value )
{
    std::string tmp_str = std::to_string(value);
    
    AddParameter(paramName, tmp_str);
}


void ZParams::AddParameter( const std::string paramName, const float value )
{
    std::string tmp_str = std::to_string(value);
    
    AddParameter(paramName, tmp_str);
}


void ZParams::AddParameter( const std::string paramName, const bool value )
{
    std::string tmp_str = std::to_string(value);
    
    AddParameter(paramName, tmp_str);
}


unsigned int ZParams::GetUnsignedInt(const std::string paramName) const
{
    std::string retreived_str = GetString(paramName);
    
    auto tmp_value = boost::lexical_cast<unsigned int>(retreived_str);
    
    return tmp_value;
}


unsigned int ZParams::GetUnsignedInt(const std::string paramName, const unsigned int defaultValue) const
{
    unsigned int retreived_value = 0;
    
    try {
        retreived_value = GetUnsignedInt(paramName);
    }
    catch (...) {
        return defaultValue;
    }
    
    return retreived_value;
}


int ZParams::GetInt(const std::string paramName) const
{
    std::string retreived_str = GetString(paramName);
    
    std::size_t decimal_pos = retreived_str.find(".");
    if ( decimal_pos != std::string::npos ) {
        std::cerr << "Warning: parameter is treated as integer, but found to be float ("
        << paramName << "=" << retreived_str << std::endl;
    }
    
    auto tmp_value = boost::lexical_cast<int>(retreived_str);
    
    return tmp_value;
}


int ZParams::GetInt(const std::string paramName, const int defaultValue) const
{
    int retreived_value = 0;
    
    try {
        retreived_value = GetInt(paramName);
    }
    catch (...) {
        return defaultValue;
    }
    
    return retreived_value;
}


unsigned long ZParams::GetUnsignedLong(const std::string paramName) const
{
    std::string retreived_str = GetString(paramName);
    
    std::size_t decimal_pos = retreived_str.find(".");
    if ( decimal_pos != std::string::npos ) {
        std::cerr << "Warning: parameter is treated as integer, but found to be float ("
        << paramName << "=" << retreived_str << std::endl;
    }
    
    auto tmp_value = boost::lexical_cast<unsigned long>(retreived_str);;
    
    return tmp_value;
}


unsigned long ZParams::GetUnsignedLong(const std::string paramName, const unsigned long defaultValue) const
{
    unsigned long retreived_value = 0;
    
    try {
        retreived_value = GetUnsignedLong(paramName);
    }
    catch (...) {
        return defaultValue;
    }
    
    return retreived_value;
}


void ZParams::UpdateParameter( const std::string paramName, const std::string value )
{
    if ( _read_only && _is_initialized ) {
        std::string err_msg = "ZParams: Read only. Cannot update parameters.";
        err_msg += paramName;
        throw err_msg.c_str();
    }
    
    try {
        _param_map.at( paramName ) = value;
    }
    catch (...) {
        AddParameter( paramName,  value );
    }
}


std::string ZParams::GetString(const std::string paramName) const
{
    
    if (!_is_initialized) {
        std::string err_msg = "ZParams uninitialized. Error while attempting to get ";
        err_msg += paramName;
        throw err_msg.c_str();
    }
    
    std::string retreived_value;
    
    // the original exception doesn't say what parameter is missing
    try {
        retreived_value = _param_map.at(paramName);
    }
    catch (...) {
        std::string err_msg = "ZParams: Error while attempting to get ";
        err_msg += paramName;
        throw err_msg.c_str();
    }
    
    return retreived_value;
}

std::string ZParams::GetString(const std::string paramName, const std::string defaultValue) const
{
    std::string retreived_value = "";
    
    try {
        retreived_value = GetString(paramName);
    }
    catch (...) {
        return defaultValue;
    }
    
    return retreived_value;
}


float ZParams::GetFloat(const std::string paramName) const
{
    std::string retreived_str = GetString(paramName);
    
    auto tmp_value = boost::lexical_cast<float>(retreived_str);
    
    return tmp_value;
}


float ZParams::GetFloat(const std::string paramName, const float defaultValue) const
{
    float retreived_value = 0.0;
    
    try {
        retreived_value = GetFloat(paramName);
    }
    catch (...) {
        return defaultValue;
    }
    
    return retreived_value;
}


double ZParams::GetDouble(const std::string paramName) const
{
    std::string retreived_str = GetString(paramName);
    
    auto tmp_value = boost::lexical_cast<double>(retreived_str);
    
    return tmp_value;
}


double ZParams::GetDouble(const std::string paramName, const double defaultValue) const
{
    double retreived_value = 0.0f;
    
    try {
        retreived_value = GetDouble(paramName);
    }
    catch (...) {
        return defaultValue;
    }
    
    return retreived_value;
}


std::string ZParams::GetAllParameters() const
{
    std::string tmp_str = "name\t\tvalue\n";
    
    for ( auto current_pair : _param_map ) {
        tmp_str += current_pair.first + "\t\t" + current_pair.second + "\n";
    }
    
    return tmp_str;
}


bool ZParams::IsEmpty() const
{
    return _param_map.empty();
}
