#ifndef __NO_DEFAULTS_H
#define __NO_DEFAULTS_H

/// \file nodefaults.h
/// \brief Implements NO_DEFAULTS


/// Disables default copy constructor and assignment operators

/// From http://home.pb.net/~tglenn/CPPNOTES.html. Inherit from this
/// class in order to inhibit the automatic generation of the default
/// copy constructor and the default assignment operator of the
/// derived class.
class NO_DEFAULTS {
public:
    NO_DEFAULTS() {};
private:
    // hide these - DO NOT IMPLEMENT!
    NO_DEFAULTS(const NO_DEFAULTS&);
    NO_DEFAULTS& operator=(const NO_DEFAULTS&);
};

#endif

