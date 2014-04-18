/*
 * Define a class for multidimensional array Q
 *
 * $Id: arrayQ.hpp,v 1.9 2003/06/14 00:24:04 stephens Exp $
 */

#ifndef ARRAY_Q_H
#define ARRAY_Q_H

#include "constants.hpp"
#include <string>
#include <vector>
#include <cmath>
using namespace::std;

class ArrayQ {
    double **** array;
    int nchr;             // Number of chromosomes, range of 2nd index
    int nalleles; // Number of alleles
                  // range of the 4th and 5th indices
    // Priave copy constructor, this disable pass by value.
    // We only need one Q anyway.

    // Boundary checking when debugging
    // Return 0 when no error
    int check_bound (int, int, int, int) const;
    // Calculate Q
    void CalcQ ( double theta, double delta );    
public:
  ArrayQ ();  
  ArrayQ (const ArrayQ &); // copy constructor
  const ArrayQ & operator= ( const ArrayQ & );

  ArrayQ ( char LociType,          // Loci Type
             int NInd,                 // Number of Individuals
             double theta ,
	     double delta);
    
  ~ArrayQ ();
  // Return const reference so that changine the elements this way
  // is not allowed.
  const double & operator () (int j, int k, int l, int m)
        const {
#ifdef CHECK_BOUNDARY
      assert ( check_bound (j, k, l, m) == 0 );
#endif
      return array[j][k][l][m];
    }
};

// prob of moving from copying nfrom,cfrom to nt,cto
// when there are nc chromosomes
inline double TransitionProb ( int nc,
                               int nfrom,
                               int cfrom,
                               int nto,
                               int cto,
                               double rho )
{    
    return (1.0 / nc ) * ( 1 - exp( -rho/nc ) )+
        ( nfrom == nto ) * ( cfrom == cto ) * exp( -rho/nc );
}

inline double TransitionProb ( int nc,
                               int nfrom,
                               int cfrom,
			       int tfrom,
                               int nto,
                               int cto,
			       int tto,
                               double rho )
{    
    return (WEIGHTS[tto] / (SS*nc) ) * ( 1 - exp( -rho/nc ) )+
        ( nfrom == nto ) * ( cfrom == cto ) * ( tfrom == tto ) * exp( -rho/nc );
}


inline double TransitionProb ( int nc,
			       int nfrom,
                               int cfrom,
			       int tfrom,
                               int nto,
                               int cto,
			       int tto,
                               double rho,
			       double expsave )
{    
    return (WEIGHTS[tto] / (SS*nc) ) * ( 1 - expsave )+
        ( nfrom == nto ) * ( cfrom == cto ) * ( tfrom == tto ) * expsave;
}




/** Probability of observing an allelic type given it is copied from 
    another allele allele.
*/
inline double PrHitTarg ( int r, int n, int t, int from, int to, 
                          const vector<ArrayQ *> & Q ) 
{
  //cout << "r=" << r << ", n=" << n << ", t=" << t << ", from=" << from << ", to=" << to << ", Q= " << (*Q[r])(n, t, from, to);
    return (*Q[r])(n, t, from, to);

}

inline double PrHitTarg ( int n, int from, int to, double theta ) 
{
  if(from==to)
    return (theta/(n+theta)*0.5+(n/(n+theta)));
  else
    return theta/(n+theta)*0.5;

}
inline double FuzzyPrHitTarg ( int r, int n, int t, float from, int to, 
                          const vector<ArrayQ *> & Q ) 
{
  //cout << "r=" << r << ", n=" << n << ", t=" << t << ", from=" << from << ", to=" << to << ", Q= " << (*Q[r])(n, t, from, to);
    return from * PrHitTarg(r,n,t,1,to,Q) + (1-from) * PrHitTarg(r,n,t,0,to,Q);

}

inline double FuzzyPrHitTarg ( int n, float from, int to, double theta ) 
{
  return  from * PrHitTarg(n,1,to,theta) + (1-from) * PrHitTarg(n,0,to,theta);
}

inline double PrHitTargNoQuad ( int n, int from, int to, double theta ) 
{
  if(from==to)
    return (theta/(n+theta)*0.5+(n/(n+theta)));
  else
    return theta/(n+theta)*0.5;

}




#endif // ARRAY_Q_H

// {{{ Log
// 
// $Log: arrayQ.hpp,v $
// Revision 1.9  2003/06/14 00:24:04  stephens
// Adding files, and committing a lot of changes
// that have resulted in version 2.0 of PHASE
//
// Revision 1.8  2002/02/27 18:56:44  stephens
// Commiting the current source, which is essentially that released as
// PHASE 1.0. The number of quadrature points is reduced to 2. Some code has
// been added to cope with recombination, but this is not included in the release
// and is still under development.
//
// Revision 1.7  2001/10/12 23:49:26  stephens
// Various updates, particularly cosmetic, but some bug fixes.
// this version tested on MS and SNP data without missing alleles
// gives very similar answers to the original phase.
//
// Revision 1.6  2001/05/21 20:31:09  stephens
// added conditions #ifdef CHECK_BOUNDARY
// to avoid always checking bounds of ArrayQ.
//
// Revision 1.5  2001/05/21 20:17:15  nali
// No real changes
//
// Revision 1.4  2001/04/20 00:32:45  nali
// Put reading loci types into the contructor of ClassPop
//
// Revision 1.3  2001/04/19 19:50:08  nali
// Minor changes to the other files.
//
// Revision 1.2  2001/04/17 22:08:44  nali
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
// Revision 1.1  2001/04/09 16:25:50  nali
// A class for multidimensional array Q.
//
// }}}
