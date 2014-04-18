/*
 * $Id: utility.hpp,v 1.7 2003/06/14 00:24:06 stephens Exp $
 */
#ifndef UTILITY_CC_H
#define UTILITY_CC_H

#include <string>
#include <map>
#include <vector>
//#include "crand/ranlib.h"
using namespace::std;

int proc_args ( int, char **, 
		std::map<string, string> &, 
		std::map<string, int> &, 
		std::map<string, double> &,
	        int &, int &, int &);

extern "C" void init_genrand(unsigned long s);
extern "C" double genrand_real2(void);

double ranf();
// generate a random integer according to a user-defined density
int rint2 ( const std::vector<double> & , double psum = -1.0 );
double rgamma(double n,double lambda);
void rdirichlet(const double * a, const int k, double * b);
double rnorm(double,double);
double logdnorm(double,double,double);
void rperm(vector<int> &, int);

#endif

// {{{ Log
//
// $Log: utility.hpp,v $
// Revision 1.7  2003/06/14 00:24:06  stephens
// Adding files, and committing a lot of changes
// that have resulted in version 2.0 of PHASE
//
// Revision 1.6  2002/02/27 18:56:45  stephens
// Commiting the current source, which is essentially that released as
// PHASE 1.0. The number of quadrature points is reduced to 2. Some code has
// been added to cope with recombination, but this is not included in the release
// and is still under development.
//
// Revision 1.5  2001/10/21 18:33:30  stephens
// fixed bug in GibbsUpdate
//
// Revision 1.3  2001/10/12 23:49:27  stephens
// Various updates, particularly cosmetic, but some bug fixes.
// this version tested on MS and SNP data without missing alleles
// gives very similar answers to the original phase.
//
// Revision 1.2  2001/04/24 19:36:37  nali
// rint2 is enhanced so that it is not necessary to provide the sum of the
// (not normalized) probabilities.
//
// Revision 1.1  2001/04/17 22:08:44  nali
//
// Major revisement in overall structure. Created new class ClassPop and
// almost all global functions now became member functions of ClassPop, most
// of them private.
//
// "mult" removed in update_phase_NR. No other changes in terms of algorithm.
// Haven't check the results yet.
//
// proc_args() is moved to utility.cpp, which also defines a couple of other
// global functions.
//
// }}}
