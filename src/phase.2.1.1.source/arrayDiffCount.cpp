#include "arrayDiffCount.hpp"
#include "utility.hpp"

using namespace std;

// Constructor
ArrayDiffCount::ArrayDiffCount ( int size ) : N ( size )
{
    ndiff = new int *** [ N ];

    for( int n0 =0; n0 < N; n0++ ){
      ndiff[n0] = new int ** [ 2 ];
      for(int c0 = 0; c0 < 2; c0++ ){
	ndiff[n0][c0] = new int * [ N ];
	for(int n1 = 0; n1 < N; n1++) {
            ndiff[n0][c0][n1] = new int [2];
            ndiff[n0][c0][n1][0] = ndiff[n0][c0][n1][1] = 0;
        }
      }
    }
}
       
// Copy Constructor
ArrayDiffCount::ArrayDiffCount ( const ArrayDiffCount & dc2 ) :
    N ( dc2.get_size() )
{ 
    ndiff = new int *** [ N ];

    for( int n0 =0; n0 < N; n0++ ){
      ndiff[n0] = new int ** [ 2 ];
      for(int c0 = 0; c0 < 2; c0++ ){
	ndiff[n0][c0] = new int * [ N ];
	for(int n1 = 0; n1 < N; n1++) {
            ndiff[n0][c0][n1] = new int [2];
            ndiff[n0][c0][n1][0] = dc2(n0,c0,n1,0);
	    ndiff[n0][c0][n1][1] = dc2(n0,c0,n1,1);
        }
      }
    }
}
       
// Destructor
ArrayDiffCount::~ArrayDiffCount ()
{
  //  cout << "DiffCount destructor called" << endl;
  for( int n0 = 0; n0 < N; n0++){
    for( int c0 = 0; c0 < 2; c0++){
      for( int n1 = 0; n1 < N; n1++){
	delete [] ndiff[n0][c0][n1];
      }
      delete [] ndiff[n0][c0];
    }
    delete [] ndiff [n0];
  }
  delete [] ndiff;

}

// Assignment
const ArrayDiffCount & ArrayDiffCount::operator= ( const ArrayDiffCount & rhs )
{ 

  if(ndiff){
    for( int n0 = 0; n0 < N; n0++){
      for( int c0 = 0; c0 < 2; c0++){
	for( int n1 = 0; n1 < N; n1++){
	  delete [] ndiff[n0][c0][n1];
	}
	delete [] ndiff[n0][c0];
      }
      delete [] ndiff [n0];
    }
    delete [] ndiff;
  }
  
  if(this != &rhs){
    N = rhs.get_size();
    
    this->ndiff = new int *** [ N ];
    
    for( int n0 =0; n0 < N; n0++ ){
      this->ndiff[n0] = new int ** [ 2 ];
      for(int c0 = 0; c0 < 2; c0++ ){
	this->ndiff[n0][c0] = new int * [ N ];
	for(int n1 = 0; n1 < N; n1++) {
	  this->ndiff[n0][c0][n1] = new int [2];
	  this->ndiff[n0][c0][n1][0] = rhs(n0,c0,n1,0);
	  this->ndiff[n0][c0][n1][1] = rhs(n0,c0,n1,1);
	}
      }
    }
  }
  return *this;
}
       


void ArrayDiffCount::resize ( int size )
{
  //  cout << "DiffCount resizer called with size " << size << endl;
  for( int n0 = 0; n0 < N; n0++){
    for( int c0 = 0; c0 < 2; c0++){
      for( int n1 = 0; n1 < N; n1++){
	delete [] ndiff[n0][c0][n1];
      }
      delete [] ndiff[n0][c0];
    }
    delete [] ndiff [n0];
  }
  delete [] ndiff;
  
  N=size;
  ndiff = new int *** [ N ];

  for( int n0 =0; n0 < N; n0++ ){
    ndiff[n0] = new int ** [ 2 ];
    for(int c0 = 0; c0 < 2; c0++ ){
      ndiff[n0][c0] = new int * [ N ];
      for(int n1 = 0; n1 < N; n1++) {
	ndiff[n0][c0][n1] = new int [2];
	ndiff[n0][c0][n1][0] = ndiff[n0][c0][n1][1] =0;
      }
    }
  }
}



/** Compute DiffCount
    DiffCount[n0][c0][n1][c1] is the number of differences at positions in uselist
    (usually all SNP positions) between the individual n0's chromosome c0, and ind n1, chr c1.
*/
void ArrayDiffCount::compute ( const vector<CIndividual> & pop, 
			       const vector<int> & uselist)
{
  for( int n0 = 0; n0 < N; n0++){
    for (int c0 = 0; c0 < 2;  ++c0) {
      for (int n1 = 0; n1 < N;  ++n1) {
	for (int c1 = 0; c1 < 2;  ++c1) {
	  ndiff[n0][c0][n1][c1] = NDiff(pop,n0,c0,n1,c1,uselist);
	}
      }
    }
  }
}

//
// update the matrix DiffCount, 
// for the new alleles now at locus
// in individual id
//
void ArrayDiffCount::Update( int id, 
			     const vector<CIndividual> & pop,
			     int locus, int oldtarg0, int oldtarg1)
{ 
#ifdef DEBUG
  cout << "Calling Update with id=" << id << " and locus = " << locus << endl;
#endif
  
  int targ0 = pop[id].get_haplotype (0, locus);
  int targ1 = pop[id].get_haplotype (1, locus);
  int from;
  int adjust;

  //if targ0==targ1, no update needed; otherwise...
  if((targ0 != oldtarg0) || (targ1 != oldtarg1))
    {
      for(int n1=0; n1 < N ; n1++){
	if(n1!=id){
	  for(int c1=0; c1<2 ; c1++){
	    from=pop[n1].get_haplotype(c1,locus);
	    //cout << "Ndiff=" << ndiff[id][0][n1][c1] << "," << ndiff[id][1][n1][c1] << endl;
	    //cout << from << " " << targ0 << oldtarg0 << targ1 << oldtarg1 << endl;
	    //adjust = (from!=targ0) - (from!=targ1) ;
	    ndiff[id][0][n1][c1] += ((from!=targ0) - (from!=oldtarg0));
	    //cout << ndiff[id][0][n1][c1] << endl;
	    ndiff[id][1][n1][c1] += ((from!=targ1) - (from!=oldtarg1));
	    
	    //cout << ndiff[id][1][n1][c1] << endl;
	  }
	}
      }
    }

}
//
// update after flipping of phase at locus in id
//
void ArrayDiffCount::Update( int id, 
			     const vector<CIndividual> & pop,
			     int locus)
{ 
  Update(id, pop, locus, pop[id].get_haplotype(1,locus), pop[id].get_haplotype(0,locus));
}


//
// compute the prob of two chromosomes in individual id based on just DiffCount
// differences (NSNP is total number of snps, so NSNP-DiffCount
// gives number of positions that are the same)
//
double ArrayDiffCount::calc_prob(int id, int NSNP, const ArrayDiffProb & DiffProb)
{
#ifdef DEBUG
  cout << "ArrayDiffCount::calc_prob called with id =" << id << " and NSNP =" << NSNP << endl;
#endif
  double firstprob = 0.0; // probability of first chromosome (c=0)
  double secondprob = 0.0; // probability of second chromosome (c=1)
  int nchr= N + N -2;
  int nchrplus1 = nchr + 1;

  for(int t0=0;t0<SS;++t0){
    double firstterm=0;
    double secondterm=0;
    for (int n0 = 0; n0 < N; ++n0) {
      for (int c0 = 0; c0 < 2; ++c0) {
	firstterm += DiffProb(nchr, t0, ndiff[id][0][n0][c0], 1) * 
	  DiffProb(nchr, t0, NSNP - ndiff[id][0][n0][c0], 0);

	secondterm += DiffProb(nchrplus1, t0, ndiff[id][1][n0][c0], 1) * 
	  DiffProb(nchrplus1, t0, NSNP - ndiff[id][1][n0][c0], 0);
      }
    }

    // correct for fact that neither can copy itself, and first cannot copy second
    firstterm -= DiffProb(nchr, t0, ndiff[id][0][id][0], 1) * 
      DiffProb(nchr, t0, NSNP - ndiff[id][0][id][0], 0);
    firstterm -= DiffProb(nchr, t0, ndiff[id][1][id][0], 1) * 
      DiffProb(nchr, t0, NSNP - ndiff[id][1][id][0], 0 );
    secondterm -= DiffProb(nchrplus1, t0, ndiff[id][1][id][1], 1) * 
      DiffProb(nchrplus1, t0, NSNP - ndiff[id][1][id][1], 0 );

    firstprob+= WEIGHTS[t0]*firstterm;
    secondprob+= WEIGHTS[t0]*secondterm;
  }
  double result=firstprob*secondprob;
 
  return result;

}

//
// compute probability based on both DiffCount and CC
// for combining SNPs and micsats
// 
double ArrayDiffCount::CombineProb ( int id, int NSNP, const ArrayDiffProb & DiffProb, const ArrayCC & CC )
{
  double firstprob = 0.0; // probability of first chromosome (c=0)
  double secondprob = 0.0; // probability of second chromosome (c=1)
  int nchr= N + N -2;
  int nchrplus1 = nchr + 1;

  for(int t0=0;t0<SS;++t0){
    double firstterm=0;
    double secondterm=0;
    for (int n0 = 0; n0 < N; ++n0) {
      for (int c0 = 0; c0 < 2; ++c0) {
	if(n0 !=id){
	  firstterm += DiffProb(nchr, t0, ndiff[id][0][n0][c0], 1) * 
	    DiffProb(nchr, t0, NSNP - ndiff[id][0][n0][c0], 0) * CC(id,0,n0,c0,t0);
	}
	if(n0 !=id || c0 !=1){
	  secondterm += DiffProb(nchrplus1, t0, ndiff[id][1][n0][c0], 1) * 
	    DiffProb(nchrplus1, t0, NSNP - ndiff[id][1][n0][c0], 0) * CC(id,1,n0,c0,t0);
	}
      }
    }

    // correct for fact that neither can copy itself, and first cannot copy second
    //  firstterm -= DiffProb(nchr, t0, ndiff[0][id][0], 1) * 
    //  DiffProb(nchr, t0, NSNP - ndiff[0][id][0], 0) * CC(id,0,t0,0);
    //firstterm -= DiffProb(nchr, t0, ndiff[1][id][0], 1) * 
    //  DiffProb(nchr, t0, NSNP - ndiff[1][id][0], 0 ) * CC(id,1,t0,0);
    //secondterm -= DiffProb(nchrplus1, t0, ndiff[1][id][1], 1) * 
    //  DiffProb(nchrplus1, t0, NSNP - ndiff[1][id][1], 0 ) * CC(id,1,t0,1);

    firstprob+= WEIGHTS[t0]*firstterm;
    secondprob+= WEIGHTS[t0]*secondterm;
  }
  double result=firstprob*secondprob;
  
  if(log(result)<-1000){
    cout << "Warning: potential underflow problem" << result << endl;
  }
  if(log(result)>1000){
    cout << "Warning: potential overflow problem" << result << endl;
  }
  return result; 

}

//
// resymettrize the diff count data on all individuals, 
// after updating the phase at id
//
void Resymmetrize(ArrayDiffCount & DiffCount, int id)
{
   for(int n0 = 0; n0 < DiffCount.get_size(); n0++){
      if(n0!=id){
	DiffCount.set_element(n0,0,id,0,DiffCount(id,0,n0,0));
	DiffCount.set_element(n0,0,id,1,DiffCount(id,1,n0,0));
	DiffCount.set_element(n0,1,id,0,DiffCount(id,0,n0,1));
	DiffCount.set_element(n0,1,id,1,DiffCount(id,1,n0,1));
      }
    }
}
