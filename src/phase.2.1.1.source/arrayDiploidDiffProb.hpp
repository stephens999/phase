/*
 * Define a class for DiploidDiffProb, which holds the probability of
 * seeing i SNP differences at r loci in times t0,t1
 *
 */

#ifndef ARRAY_DIPLOIDDIFFPROB_H
#define ARRAY_DIPLOIDDIFFPROB_H

#include "constants.hpp"
#include "arrayQ.hpp"
#include <string>
#include <vector>
#include <cmath>
using namespace::std;

class ArrayDiploidDiffProb {

    double ***** array; // [n][t0][t1][r][i] : prob of i differences, to the power r 
     
    int Nind;             // Number of individuals, range of 1st index
    int nloci;            // Number of loci, range of 3rd index

    // Privave copy constructor, this disable pass by value.
    // We only need one DiploidDiffProb anyway.
    ArrayDiploidDiffProb ( const ArrayDiploidDiffProb & ) {}
 
    // Boundary checking when debugging
    // Return 0 when no error
    int check_bound (int, int, int, int, int) const;

    // Calculate DiploidDiffProb
    void CalcDiploidDiffProb (const vector<ArrayQ *> & Q, const string & LociType);    

public:

  // constructor and destructor
    ArrayDiploidDiffProb ( const string & LociType,          // Loci Type
             int NInd,                 // Number of Individuals
             vector<ArrayQ *> & Q );
    ~ArrayDiploidDiffProb ();
    // Return const reference so that changine the elements this way
    // is not allowed.
    const double & operator () (int i, int j, int k, int l, int m)
        const {
#ifdef CHECK_BOUNDARY
      assert ( check_bound (i, j, k, l, m) == 0 );
#endif
     return array[i][j][k][l][m];
    }
};


#endif // ARRAY_DIPLOIDDIFFPROB_H

