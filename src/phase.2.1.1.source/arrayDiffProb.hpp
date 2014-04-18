/*
 * Define a class for DiffProb, which holds the probability of
 * seeing n SNP differences (or similarities) in time t
 * (=Q^n)
 *
 */

#ifndef ARRAY_DIFFPROB_H
#define ARRAY_DIFFPROB_H

#include "constants.hpp"
#include "arrayQ.hpp"
#include <string>
#include <vector>
#include <cmath>
using namespace::std;

class ArrayDiffProb {
    double **** array; // [nchr][t][ndiff][0 or 1] (0 or 1 indicates difference or sames)
     
    int nchr;             // Number of chromosomes, range of 1st index
    int nloci;            // Number of loci, range of 3rd index

  // Privave copy constructor, this disable pass by value.
    // We only need one DiffProb anyway.
    ArrayDiffProb ( const ArrayDiffProb & ) {}
    // Boundary checking when debugging
    // Return 0 when no error
    int check_bound (int, int, int, int) const;
    // Calculate DiffProb
    void CalcDiffProb (const vector<ArrayQ *> & Q, const string & LociType);    
public:

  // constructor and destructor
    ArrayDiffProb ( const string & LociType,          // Loci Type
             int NInd,                 // Number of Individuals
             vector<ArrayQ *> & Q );
    ~ArrayDiffProb ();
    // Return const reference so that changine the elements this way
    // is not allowed.
    const double & operator () (int i, int j, int k, int l)
        const {
#ifdef CHECK_BOUNDARY
      assert ( check_bound (i, j, k, l) == 0 );
#endif
     return array[i][j][k][l];
    }
};


#endif // ARRAY_DIFFPROB_H

