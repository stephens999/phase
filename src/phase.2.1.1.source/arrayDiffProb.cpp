
/* 
 * Implementation of ArrayDiffProb
 *
 */

#include "arrayDiffProb.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace std;

ArrayDiffProb::ArrayDiffProb ( const string & LociType,          // Loci Type
                 int NInd,                         // Number of Individuals
                 vector<ArrayQ *> & Q ) : nchr(2*NInd)
{
    nloci = LociType.size();
  
    // Allocate memory for array.
    cout << "Allocating memory for DiffProb" << endl;

    array = new double *** [ nchr ];

    // currently not storing DiffProb for n< nchr-2, to save space
    for (int n = nchr-2; n < nchr; ++n) {
        array[n] = new double ** [ SS ];        
        for (int t = 0; t < SS; ++t) {
            array[n][t] = new double * [ nloci+1 ];
	    for (int r = 0; r < (nloci+1); ++r) {
	      array[n][t][r] = new double [ 2 ];
	    }
	}
    }
    
    // Calculate DiffProb
    CalcDiffProb ( Q, LociType );

}

ArrayDiffProb::~ArrayDiffProb ( )
{  
  for (int n = nchr-2; n < nchr; ++n) {
    for (int s = 0; s < SS; ++s) {
      for (int r = 0; r < (nloci+1); ++r) {
	delete [] array[n][s][r];
      }
      delete [] array[n][s];
    }
    delete [] array[n];
  }
  delete [] array;
}

// Private functions
int ArrayDiffProb::check_bound (int i, int j, int k, int l) const
{
    if ( i >= nchr || i < 1 ) {
        cerr << "First index of DiffProb out of range!" << i 
             << "\nwhere as max index = " << nchr - 1 << endl;
        return 1;
    }
    if ( j >= SS || j < 0 ) {
        cerr << "Second index of DiffProb out of range! " << j 
             << "\nwhere as max index = " << SS - 1 << endl;
        return 2;
    } 
    if ( k >= (nloci+1) || k < 0 ) {
        cerr << "Third index of DiffProb out of range! " << k 
             << "\nwhere as max index = " << nloci - 1 << endl;
        return 3;
    }
    if ( l >= 2 || l < 0 ) {
        cerr << "Fourth index of DiffProb out of range! " << l 
             << "\nwhere as max index = " << 1 << endl;
        return 4;
    }
    //cerr << i << " " << j << " " << k << " " << l << endl;
    return 0;
}

void ArrayDiffProb::CalcDiffProb ( const vector<ArrayQ *> & Q, const string & LociType )
{
  
  int firstSNP= 0; //first SNP locus;
  while((LociType[firstSNP]!='S') && (firstSNP<(LociType.size()-1)))
    firstSNP++;
 
  cout << "computing DiffProb; please wait" << endl;
  for(int n=nchr-2; n < nchr; n++)
    for(int t=0; t<SS ; t++)
      for(int r=0; r< (nloci+1); r++)
	{
	  array[n][t][r][0]=exp(r*log((*Q[firstSNP])(n,t,0,0))-10*log((*Q[firstSNP])(n,0,0,0))); // the minus bit is just a constant there
	  array[n][t][r][1]=exp(r*log((*Q[firstSNP])(n,t,0,1))-10*log((*Q[firstSNP])(n,0,0,0))); // to reduce problems of underflow
	}


}
                        
