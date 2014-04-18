
/* 
 * Implementation of ArrayDiploidDiffProb
 * $Id: arrayDiploidDiffProb.cpp,v 1.5 2003/06/14 00:24:04 stephens Exp $
 */

#include "arrayDiploidDiffProb.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace std;

// array[n][t0][t1][r][i] holds the probability of i SNP differences
// at r loci (actually, to the power r), if the (n+1)st individual copied
// at times t0 and t1. It is essentially a look-up table
// used to avoid a lot of unnecessary multiplications in computation of arrayFF.
// i=0:  00 copies 00
// i=1: 00 copies 01
// i=2: 00 copies 11
// i=3: 01 copies 01

ArrayDiploidDiffProb::ArrayDiploidDiffProb ( const string & LociType,          // Loci Type
					     int N,                         // Number of Individuals
					     vector<ArrayQ *> & Q ) : Nind(N)
{
    nloci = LociType.size();
  
    // Allocate memory for array.
    cout << "Allocating memory for ArrayDiploidDiffProb" << endl;

    array = new double **** [ Nind ];
    int n=Nind-1;
    //    for (int n = 0; n < Nind; ++n) {
        array[n] = new double *** [ SS ];        
        for (int t0 = 0; t0 < SS; ++t0) {
            array[n][t0] = new double ** [ SS ];
	    for(int t1 =0; t1 < SS; ++t1) {
	      array[n][t0][t1] = new double * [nloci+1];
	      for (int r = 0; r < (nloci+1); ++r) {
		array[n][t0][t1][r] = new double [ 4 ];
		
	      }
	    }
	}
	//}
    // Calculate DiffProb
    CalcDiploidDiffProb ( Q, LociType );

}

ArrayDiploidDiffProb::~ArrayDiploidDiffProb ( )
{  
  for (int n = Nind-1; n < Nind; ++n) {
    for (int t0 = 0; t0 < SS; ++t0) {
       for (int t1 = 0; t1 < SS; ++t1) {    
	 for (int r = 0; r < (nloci+1); ++r) {
	   delete [] array[n][t0][t1][r];
	 }
	 delete [] array[n][t0][t1];
       }
       delete [] array[n][t0];
    }
    delete [] array[n];
  }
  delete [] array;
}

// Private functions
int ArrayDiploidDiffProb::check_bound (int i, int j, int k, int l, int m) const
{
    if ( i >= Nind || i < 1 ) {
        cerr << "First index of DiploidDiffProb out of range!" << i 
             << "\nwhere as max index = " << Nind - 1 << endl;
        return 1;
    } 
    if ( j >= SS || j < 0 ) {
        cerr << "Second index of DiploidDiffProb out of range! " << j 
             << "\nwhere as max index = " << SS - 1 << endl;
        return 2;
    }
    if ( k >= SS || k < 0 ) {
        cerr << "Third index of DiploidDiffProb out of range! " << k 
             << "\nwhere as max index = " << SS - 1 << endl;
        return 3;
    } 
    if ( l >= (nloci+1) || l < 0 ) {
        cerr << "Fourth index of DiploidDiffProb out of range! " << l 
             << "\nwhere as max index = " << nloci << endl;
        return 4;
    }
    if ( m >= 4 || m < 0 ) {
      cerr << "Fifth index of DiploidDiffProb out of range! " << m 
	   << "\nwhere as max index = " << 4 << endl;
        return 5;
    }
    return 0;
}

void ArrayDiploidDiffProb::CalcDiploidDiffProb ( const vector<ArrayQ *> & Q, const string & LociType )
{
  
  cout << "Computing DiploidDiffProb" << endl;

  int firstSNP= 0; //first SNP locus;
  while((LociType[firstSNP]!='S') && (firstSNP<(LociType.size()-1)))
    firstSNP++;
 
  //  for(int n=1; n < Nind; n++)
  int n= Nind-1;
    for(int t0=0; t0<SS ; t0++)
      for(int t1=0; t1<SS ; t1++)
	for(int r=0; r< (nloci+1); r++)
	  {
	    int nchr=n+n;
	    array[n][t0][t1][r][0]=
	      exp(r*log((*Q[firstSNP])(nchr,t0,0,0)*(*Q[firstSNP])(nchr+1,t1,0,0) + 
			(*Q[firstSNP])(nchr,t1,0,0)*(*Q[firstSNP])(nchr+1,t0,0,0) ));

	    array[n][t0][t1][r][1]=
	      exp(r*log( (*Q[firstSNP])(nchr,t0,0,1)*(*Q[firstSNP])(nchr+1,t1,0,0) + 
			 (*Q[firstSNP])(nchr,t1,0,1)*(*Q[firstSNP])(nchr+1,t0,0,0)));

	    array[n][t0][t1][r][2]=
	      exp(r*log( (*Q[firstSNP])(nchr,t0,0,1)*(*Q[firstSNP])(nchr+1,t1,0,1) + 
			 (*Q[firstSNP])(nchr,t1,0,1)*(*Q[firstSNP])(nchr+1,t0,0,1))); 

	    array[n][t0][t1][r][3]=
	      exp(r*log( (*Q[firstSNP])(nchr,t0,0,0)*(*Q[firstSNP])(nchr+1,t1,0,0) + 
			 (*Q[firstSNP])(nchr,t1,0,1)*(*Q[firstSNP])(nchr+1,t0,0,1)));

	  }


}
                        
