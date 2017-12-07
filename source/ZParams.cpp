//
//  ZParams.cpp
//  ZParams
//
//  Created by Tal Zinger on 07/09/2016.
//  Copyright © 2016 Stern Lab. All rights reserved.
//

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
    
    
    std::regex rx{"(^[^#]\\S+)\\s(\\S+)"}; // (name) (value) pairs, ignore comments with #
    std::string line;
    
    while(getline(f, line)) {
        // match regex
        std::smatch matches; // matched strings go here
        if (std::regex_search(line , matches, rx)) {
            // 0 is the whole line
            std::string param_name = matches[1];
            std::string param_val_str = matches[2];
            // std::cout << "name: " << param_name << " val: " << param_val_str << std::endl;
            AddParameter(param_name, param_val_str);
        }
    }
    
    if (f.bad()) {
        std::string err_msg = "ZParams: error while reading file " + filename + "\n";
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
    double retreived_value = 0.0;
    
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
