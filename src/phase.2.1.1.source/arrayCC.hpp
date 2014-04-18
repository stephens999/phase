#ifndef ARRAY_CC_H
#define ARRAY_CC_H

#include "constants.hpp"
#include "arrayQ.hpp"
#include "arrayDiffProb.hpp"
#include "indnode.hpp"
//#include "crand/ranlib.h"
#include "utility.hpp"

#include <cmath>
#include <vector>
using namespace::std;

class ArrayCC {
    int N;                        // total number of individuals
    double ***** veccc;
    
    // copy constructor disabled.
    // Can only declare a vector of pointers to CC
    // ArrayCC ( const ArrayCC & ) {
    //   cout << "CC Copy constructor called" << endl;
    // }
public:
    ArrayCC ( int = 0 );
    ~ArrayCC ();
    ArrayCC ( const ArrayCC & );
  
  const ArrayCC & operator= (const ArrayCC &);

    void resize ( int );

    /** Direct access function to the elements of CC 
     */
    double operator() ( int, int, int, int, int ) const ;
  
    // Compute CC for id at known phase positions
    // void compute   ( int id, 
    //                  const ArrayQ &, 
    //                  const std::vector<CIndividual> &);

    // Compute CC for id at list of given phase positions
    void compute   ( int id, 
                     const vector<ArrayQ *> &, 
                     const std::vector<CIndividual> &,
		     const std::vector<int> &, 
		     const std::string &,
		     const ArrayDiffProb &);
  
    void Update (int id, 
		 const vector<ArrayQ *> &, 
		 const std::vector<CIndividual> &,
		 int locus, int, int);

    void Update (int id, 
	       const vector<ArrayQ *> &, 
	       const std::vector<CIndividual> &,
		 int locus);

    double calc_prob(int id);

    void set_element (int, int, int, int, int, double);

    int get_size() const;
};

inline int ArrayCC::get_size() const {
    return N;
}

inline double ArrayCC::operator() (int n0, int c0, int n1, int c1, int t) const
{
#ifdef DEBUG
    assert ( n0 >= 0 && n0 < N );
    assert ( c0 >= 0 && c0 < 2 );
    assert (n1 >=0 && n1 < N);
    assert (c1 >= 0 && c1 < 2 );
    assert ( t >= 0 && t < SS );
#endif
    return veccc [n0][c0][n1][c1][t];
}

inline void ArrayCC::set_element(int n0, int c0, int n1, int c1, int t, double value)
{
  veccc[n0][c0][n1][c1][t]=value;
}

double CombineProb ( int id, const ArrayCC &, const ArrayCC &);
void Resymmetrize ( ArrayCC & CC, int id);

#endif // ARRAY_CC_H

