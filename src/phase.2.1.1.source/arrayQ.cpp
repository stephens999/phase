
/* 
 * Implementation of ArrayQ
 *
 * $Id: arrayQ.cpp,v 1.10 2003/06/14 00:24:04 stephens Exp $
 */

#include "arrayQ.hpp"
#include "tnt/tnt.h"
#include "tnt/cmat.h"
#include <iostream>
#include <cassert>
#include <cmath>
#include <iomanip>

using namespace std;


// default constructor
ArrayQ::ArrayQ ():
  nchr(0),
  nalleles(),
  array(NULL)
{  
}

// copy constructor
ArrayQ::ArrayQ (const ArrayQ & q2) :
  nchr(q2.nchr),
  nalleles(q2.nalleles)
{
// Allocate memory for array.

  array = new double *** [ nchr ];
  for (int n = 0; n < nchr; ++n) {
    //for (int n = nchr-4; n < nchr; ++n) {
    array[n] = new double ** [ SS ];
    for (int s = 0; s < SS; ++s) {
      array[n][s] = new double * [ nalleles ];
      for (int k = 0; k < nalleles; ++k) {
	array[n][s][k] = new double [ nalleles ];
	for(int j=0; j< nalleles; j++)
	  array[n][s][k][j] = q2.array[n][s][k][j];
      }
    }
  }
}

ArrayQ::ArrayQ ( const char LociType,          // Loci Type
                 int NInd,                         // Number of Individuals
                 const double theta,
		 const double delta) : nchr(2*NInd)
{
    nalleles = 2;
    if ( LociType == 'M' ) nalleles = KMAX;
    // Allocate memory for array.
    //cout << "Allocating Memory for Q" << endl;
    //cout << "nchr = " << nchr << endl;

    array = new double *** [ nchr ];
    // For microsatelites, the number of alleles is
    // set to be the maximum number allowed.
    for (int n = 0; n < nchr; ++n) {
      //for (int n = nchr-4; n < nchr; ++n) {
      array[n] = new double ** [ SS ];
      for (int s = 0; s < SS; ++s) {
	array[n][s] = new double * [ nalleles ];
	for (int k = 0; k < nalleles; ++k) {
	  array[n][s][k] = new double [ nalleles ];
	}
      }
    }

    //cout << "Memory allocation for Q done" << endl;
    // Calculate Q
    CalcQ ( theta, delta );
}

const ArrayQ & ArrayQ::operator=(const ArrayQ & rhs)
{
  
  if(array){
    //cout << "Deleting Old Q" << endl;
    //cout << "NLoci : " << nloci << endl;    
    //for (int n = nchr-4; n < nchr; ++n) {
    for (int n = 0; n < nchr; ++n) {
      for (int s = 0; s < SS; ++s) {
	for (int k = 0; k < nalleles; ++k) {
	  delete [] array[n][s][k];
	}
	delete [] array[n][s];
      }
      delete [] array[n];
    }
    delete [] array;
    //cout << "Done Deleting old Q" << endl;
  }

  if(this != &rhs){    
    nchr = rhs.nchr;
    nalleles = rhs.nalleles;

    //cout << "Now reserving memory for new Q" << endl;    
    //cout << "nchr = " << nchr << endl;
    
    array = new double *** [ nchr ];
    for (int n = 0; n < nchr; ++n) {
      //for (int n = nchr-4; n < nchr; ++n) {
      array[n] = new double ** [ SS ];
      for (int s = 0; s < SS; ++s) {
	array[n][s] = new double * [ nalleles ];
	for (int k = 0; k < nalleles; ++k) {
	  array[n][s][k] = new double [ nalleles ];
	  for(int j=0; j< nalleles; j++)
	    array[n][s][k][j] = rhs.array[n][s][k][j];
	}
      }
    }
  }  
  return *this;
}

ArrayQ::~ArrayQ ( )
{
  //cout << "DELETING Q" << endl;
  //cout << "nloci = " << nloci << endl;
  //cout << "nchr = " << nchr << endl;
  //for (int n = nchr-4; n < nchr; ++n) {
  for (int n = 0; n < nchr; ++n) {
    for (int s = 0; s < SS; ++s) {
      for (int k = 0; k < nalleles; ++k) {
	delete [] array[n][s][k];
      }
      delete [] array[n][s];
    }
    delete [] array[n];
  }
  delete [] array;
}

// Private functions
int ArrayQ::check_bound (int j, int k, int l, int m) const
{
  if ( j >= nchr || j < 1 ) {
    cerr << "Second index of Q out of range! " << j 
	 << "\nwhere as max index = " << nchr - 1 << endl;
    return 2;
  } 
  if ( k >= SS || k < 0 ) {
    cerr << "Third index of Q out of range! " << k 
	 << "\nwhere as max index = " << SS - 1 << endl;
    return 3;
  }
  if ( l >= nalleles || l < 0 ) {
    cerr << "Fourth index of Q out of range! " << l 
	 << "\nwhere as max index = " << nalleles - 1 << endl;
    return 4;
  }
  if ( m >= nalleles || m < 0 ) {
    cerr << "Fifth index of Q out of range! " << m 
	 << "\nwhere as max index = " << nalleles - 1 << endl;
    return 5;
  }
  return 0;
}

void ArrayQ::CalcQ (double theta, double delta)
{
  //cout << "Computing the matrix Q; please wait" << endl;
  //cout << "theta = " << theta << endl;
  //cout << "delta = " << delta << endl;
 
  static bool firsttime = true;
  static TNT::Matrix<double> PSNP (2, 2, "0.0 1.0 1.0 0.0"); 
  static vector <TNT::Matrix<double> > PSNpower ( NPOWER, PSNP );

  // Set up mutation transition matrix for SNPs

  if(firsttime){  
    // Set up powers of PSNP
    PSNpower[0] = TNT::Matrix<double> (2, 2, "1.0 0.0 0.0 1.0");
    for (int m = 1; m < NPOWER; ++m) {
      PSNpower[m] = matmult ( PSNpower[m-1], PSNP );
    }
    firsttime = false;
  }
  
  // }}}

  // Temporary variables use in calculation
  TNT::Matrix<double> Qtemp;
  double lambda = 0.0;
  double ptemp  = 0.0;
  // Q[n] corresponds to calculating Pr(n+1th chr | first n)
  // i.e. Pr(chr n | chr 0,...,n-1)
  //for (int n = nchr-4; n < nchr; ++n) {
  for (int n = 0; n < nchr; ++n) {
    // cerr << "Chromosome " <<  setw(6) << n << " done\033[A" << endl; 
    for (int t = 0; t < SS; ++t) {
      if(n)
	lambda = theta * TIMES[t] / n;
      else
	lambda = 1e100;
      if ( nalleles == 2 ) { // SNPs
	Qtemp = TNT::Matrix<double> (2, 2, 0.0);
	for (int m = 0; m < NPOWER; ++m) {
	  ptemp = pow(lambda, m) * exp(-lambda) / 
	    exp(lgamma(m+1));
	  Qtemp += PSNpower[m] * ptemp;
	}
      } else { // MSs
	
	//set up mutation matrix
	TNT::Matrix<double> PMS (KMAX, KMAX, delta/KMAX);
	for (int k = 1; k < KMAX - 1; ++k) {
	  PMS[k][k-1] += (1-delta)*0.5;
	  PMS[k][k+1] += (1-delta)*0.5;
	}
	PMS[0][1] +=1.0*(1-delta);
	PMS[KMAX-1][KMAX-2] += 1.0*(1-delta);
	
	// Set up powers of PMS
	vector <TNT::Matrix<double> > PMSpower ( NPOWER );
	PMSpower[0] = TNT::Matrix<double> (KMAX, KMAX, 0.0);
	for (int i = 0; i < KMAX; ++i) PMSpower[0][i][i] = 1.0;
	for (int m = 1; m < NPOWER; ++m) {
	  PMSpower[m] = matmult ( PMSpower[m-1], PMS );
	}
		  		  
	Qtemp = TNT::Matrix<double> (KMAX, KMAX, 0.0);
	for (int m = 0; m < NPOWER; ++m) {
	  ptemp = pow(lambda, m) * exp(-lambda) / 
	    exp(lgamma(m+1));
	  Qtemp += PMSpower[m] * ptemp;
	}
      }                
      for (int i = 0; i < nalleles; ++i) {
	for (int j = 0; j < nalleles; ++j) {
	  array[n][t][i][j] = Qtemp[i][j];
#ifdef DEBUG
	  cout << n << "," << t << "," << i << "," << j << ":" << Qtemp[i][j] << endl;
#endif
	}
      }
    }
  }
  //cout << "Done Computing Q" << endl;
}
                        
// {{{ Log
//
// $Log: arrayQ.cpp,v $
// Revision 1.10  2003/06/14 00:24:04  stephens
// Adding files, and committing a lot of changes
// that have resulted in version 2.0 of PHASE
//
// Revision 1.9  2001/10/12 23:49:26  stephens
// Various updates, particularly cosmetic, but some bug fixes.
// this version tested on MS and SNP data without missing alleles
// gives very similar answers to the original phase.
//
// Revision 1.8  2001/10/09 23:03:24  stephens
// Modified to deal with very large files - put in checks for underflow
// and overflow, and changed computation of Q, DiffProb etc to only compute those elements that are necessary.
//
// Revision 1.7  2001/06/19 06:00:21  nali
// No change, testing cvs
//
// Revision 1.6  2001/05/23 02:16:55  stephens
// Corrected some bugs in computation of CC. Reduced number of
// loci we update by one when the total number of unknowns is <=5.
// Confirmed program was giving vaguely sensible resutls, and is
// competitive on the speed test example (36secs vs 31 for old version)
//
// Revision 1.5  2001/05/21 20:17:15  nali
// No real changes
//
// Revision 1.4  2001/04/19 19:50:08  nali
// Minor changes to the other files.
//
// Revision 1.3  2001/04/17 22:08:43  nali
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
// Revision 1.2  2001/04/11 20:29:30  stephens
// Fixed bugs in calculation of Q:
// - replaced mult_element with matmult for finding powers of P
// - replaced definition of P[0] with "1.0 0.0 0.0 1.0"; previously
// this was given in 2 strings, which apparently caused an error.
// - revised iteration so that calculation of powers of PSN
//  started with 1, and not 0
//
// Revision 1.1  2001/04/09 16:25:50  nali
// A class for multidimensional array Q.
//
// }}}
