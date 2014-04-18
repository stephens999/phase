/*
 * Test for error conditioning in program
 * 
 * $Id: errcheck.hpp,v 1.1 2001/04/07 07:18:10 nali Exp $
 */

#ifndef ERRCHECK_H
#define ERRCHECK_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
using namespace::std;

// Open file to read
inline void assure ( std::ifstream & istr,
		     const std::string & filename = "") 
{
    using namespace std;
    if ( ! istr ) {
	cerr << "Could not open file " << filename 
	     << " to read." << endl;
	exit(1);
    }
}

// Open file to write
inline void assure ( std::ofstream & ostr,
		     const std::string & filename = "") 
{
    using namespace std;
    if ( ! ostr ) {
	cerr << "Could not open file " << filename
	     << " to write." << endl;
	exit(1);
    }
}

#endif  // ERRCHECK_H

// {{{ Log
// 
// $Log: errcheck.hpp,v $
// Revision 1.1  2001/04/07 07:18:10  nali
// Revised to concentrate on implementing SD method.
// Discarded support for EM algorithm, etc..
//
// }}}
