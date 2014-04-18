#ifndef ARRAY_DIFFCOUNT_H
#define ARRAY_DIFFCOUNT_H

#include "constants.hpp"
#include "arrayCC.hpp"
#include "arrayDiffProb.hpp"
#include "indnode.hpp"
//#include "crand/ranlib.h"
#include "utility.hpp"

#include <cmath>
#include <vector>
using namespace::std;

class ArrayDiffCount {
    int N;                        // total number of individuals
    int **** ndiff;
 
    // copy constructor disabled.
    // Can only declare a vector of pointers to DiffCount
  //ArrayDiffCount ( const ArrayDiffCount & ) {
  //  cout << "DiffCount Copy constructor called" << endl;
  //   exit(1);
  //  }
public:
    ArrayDiffCount ( int = 0 );
    ArrayDiffCount ( const ArrayDiffCount & );
    ~ArrayDiffCount ();
  
  const ArrayDiffCount & operator=(const ArrayDiffCount &);
    void resize (int);

    /** Direct access function to the elements of DiffCount 
     */
    int operator() ( int, int, int, int ) const ;
     
    // Compute DiffCount for id at list of given positions
    void compute   ( const std::vector<CIndividual> &,
		     const std::vector<int> &);
  
  // update DiffCount after changing alleles at locus in individual id
    void Update (int id, 
		 const std::vector<CIndividual> &,
		 int locus, int, int);
  // update DiffCount after flipping phase at locus in individual id
  void Update (int id, 
		 const std::vector<CIndividual> &,
		 int locus);

    double calc_prob(int, int,
		     const ArrayDiffProb &);

    double CombineProb ( int, int, const ArrayDiffProb &, const ArrayCC &);

    void set_element(int,int,int,int,int);
    int get_size() const;
};


inline void ArrayDiffCount::set_element(int n0, int c0, int n1, int c1, int value)
{
  ndiff[n0][c0][n1][c1]= value;
}

inline int ArrayDiffCount::get_size() const {
    return N;
}

inline int ArrayDiffCount::operator() (int n0, int c0, int n1, int c1) const
{
#ifdef DEBUG
    assert ( n0 >=0 && n0 < N);
    assert ( c0 >= 0 && c0 < 2 );
    assert ( n1 >= 0 && n1 < N );
    assert ( c1  >= 0 && c1 < 2);
#endif 
   return ndiff [n0][c0][n1][c1];
}


//
// resymettrize the diff count data on all individuals, 
// after updating the phase of individual id
//
void Resymmetrize(ArrayDiffCount & DiffCount, int id);


#endif // ARRAY_DIFFCOUNT_H

