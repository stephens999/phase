#include "arrayCC.hpp"
#include "utility.hpp"

using namespace std;

// Constructor
ArrayCC::ArrayCC ( int size ) : N ( size )
{
  veccc = new double **** [N];
  for (int n0 = 0; n0 < N; ++n0) {
    veccc[n0] = new double *** [2];
    for (int c0 = 0; c0 < 2; ++c0) {
      veccc[n0][c0] = new double ** [N];
      for (int n1 = 0; n1 < N; ++n1) {
	veccc[n0][c0][n1] = new double * [2];
	for (int c1 = 0; c1 < 2; ++c1) {
	  veccc[n0][c0][n1][c1] = new double [SS];
	  for (int t = 0; t < SS; ++t) {
	    // Initialize veccc[i][j][k][l]
	    veccc[n0][c0][n1][c1][t] = 1.0;
	  }
        }
      }
    }
  }
}

// Copy Constructor
ArrayCC::ArrayCC ( const  ArrayCC & CC2 ) : N ( CC2.get_size() )
{
  veccc = new double **** [N];
  for (int n0 = 0; n0 < N; ++n0) {
    veccc[n0] = new double *** [2];
    for (int c0 = 0; c0 < 2; ++c0) {
      veccc[n0][c0] = new double ** [N];
      for (int n1 = 0; n1 < N; ++n1) {
	veccc[n0][c0][n1] = new double * [2];
	for (int c1 = 0; c1 < 2; ++c1) {
	  veccc[n0][c0][n1][c1] = new double [SS];
	  for (int t = 0; t < SS; ++t) {
	    // Initialize veccc[i][j][k][l]
	    veccc[n0][c0][n1][c1][t] = CC2(n0,c0,n1,c1,t);
	  }
        }
      }
    }
  }
}

const ArrayCC & ArrayCC::operator=(const ArrayCC & rhs)
{
  
  if(veccc){
    for (int n0 = 0; n0 < N; ++n0) {
      for (int c0 = 0; c0 < 2; ++c0) {
	for (int n1 = 0; n1 < N; ++n1) {
	  for (int c1 = 0; c1 < 2; ++c1) {
	    delete [] veccc[n0][c0][n1][c1];
	  }
	  delete [] veccc[n0][c0][n1];
	}
	delete [] veccc[n0][c0];
      }
      delete [] veccc[n0];
    }
    delete [] veccc;
  }

 if(this != &rhs){
   N = rhs.get_size();
    this->veccc = new double **** [N];
    for (int n0 = 0; n0 < N; ++n0) {
      this->veccc[n0] = new double *** [2];
      for (int c0 = 0; c0 < 2; ++c0) {
	this->veccc[n0][c0] = new double ** [N];
	for (int n1 = 0; n1 < N; ++n1) {
	  this->veccc[n0][c0][n1] = new double * [2];
	  for (int c1 = 0; c1 < 2; ++c1) {
	    this->veccc[n0][c0][n1][c1] = new double [SS];
	    for (int t = 0; t < SS; ++t) {
	      // Initialize veccc[i][j][k][l]
	      this->veccc[n0][c0][n1][c1][t] = rhs(n0,c0,n1,c1,t);
	    }
	  }
	}
      }
    }
  }
  return *this;
}

// Destructor
ArrayCC::~ArrayCC ()
{
  for (int n0 = 0; n0 < N; ++n0) {
    for (int c0 = 0; c0 < 2; ++c0) {
      for (int n1 = 0; n1 < N; ++n1) {
	for (int c1 = 0; c1 < 2; ++c1) {
	  delete [] veccc[n0][c0][n1][c1];
	}
	delete [] veccc[n0][c0][n1];
      }
      delete [] veccc[n0][c0];
    }
    delete [] veccc[n0];
  }
  delete [] veccc;
}

void ArrayCC::resize ( int size )
{
  for (int n0 = 0; n0 < N; ++n0) {
    for (int c0 = 0; c0 < 2; ++c0) {
      for (int n1 = 0; n1 < N; ++n1) {
	for (int c1 = 0; c1 < 2; ++c1) {
	  delete [] veccc[n0][c0][n1][c1];
	}
	delete [] veccc[n0][c0][n1];
      }
      delete [] veccc[n0][c0];
    }
    delete [] veccc[n0];
  }
  delete [] veccc;

  N = size;

  veccc = new double **** [N];
  for (int n0 = 0; n0 < N; ++n0) {
    veccc[n0] = new double *** [2];
    for (int c0 = 0; c0 < 2; ++c0) {
      veccc[n0][c0] = new double ** [N];
      for (int n1 = 0; n1 < N; ++n1) {
	veccc[n0][c0][n1] = new double * [2];
	for (int c1 = 0; c1 < 2; ++c1) {
	  veccc[n0][c0][n1][c1] = new double [SS];
	  for (int t = 0; t < SS; ++t) {
	    // Initialize veccc[i][j][k][l]
	    veccc[n0][c0][n1][c1][t] = 1.0;
	  }
        }
      }
    }
  }  

}


/** Compute CC
    CC[n0][c0][t0][c] is the probability of the known parts of
    chromosome c of individual id, given it was copied from chromosome c0 of
    individual n0, with time elapse t0.
*/

//  void ArrayCC::compute ( int id,
//                          const ArrayQ & Q,
//                          const vector<CIndividual> & pop)
//  {
//      int nloci = pop[id].get_nloci();
//      // create array containing loci with known phase
//      vector <int> knownvec;
//      for ( int locus = 0; locus < nloci; ++locus ) {
//          if ( !pop[id].is_unknown(locus) )
//              knownvec.push_back(locus);
//      }
//      compute ( id, Q, pop, knownvec );
//  }

/** Compute CC

    CC[n0][c0][n1][c1][t] is the probability of
    chromosome c0 of individual n0, given it was copied from chromosome c1 of
    individual n1, with time elapse t, at positions in uselist.
*/

void ArrayCC::compute ( int id,
                        const vector<ArrayQ *> & Q,
                        const vector<CIndividual> & pop,
                        const vector<int> & uselist, 
			const std::string & loci_type,
			const ArrayDiffProb & DiffProb)
{
    int from = 0;
    int targ0 = 0;
    int targ1= 0;
    int nchr = N + N - 2;
    int nchrplus1 = nchr + 1;


#ifdef TREAT_SNPS_AS_MS

    for ( vector<int>::const_iterator u = uselist.begin();
	  u != uselist.end(); ++u ) {
      // compute CC based only on loci in uselist	         
      targ0 = pop[id].get_haplotype (0, *u);
      targ1 = pop[id].get_haplotype (1, *u);
      for (int n1 = 0; n1 < N;  ++n1) {
	for (int c1 = 0; c1 < 2;  ++c1) {
	  from = pop[n1].get_haplotype (c1, *u);
	  for (int t = 0; t < SS; ++t) {
	    if(u==uselist.begin()){
	      if (n1 == id ) {
	    /* The first chromosome (c=0) is not allowed to copy 
	       from itself (n1==id)*/
		veccc[id][0][n1][c1][t] = 0.0;
	    
	    /* The second chromosome (c==1) of id is allowed
	       to copy from the first one (c1==0), but not 
	       from n1 = id and c1 = 1 */
		if(c1==1){
		  veccc[id][1][n1][c1][t]=0.0;
		}
		else
		  veccc[id][1][n1][c1][t] =1;
	      }
	      else{
		veccc[id][0][n1][c1][t] = 1;
		veccc[id][1][n1][c1][t] = 1;
	      }
	    }    	   
   	    veccc[id][0][n1][c1][t] *= 
	      PrHitTarg(*u, nchr , t, from, targ0, Q);
	    veccc[id][1][n1][c1][t] *= 
	      PrHitTarg(*u, nchrplus1 , t, from, targ1, Q);	
	  }  
	}
      }
    }

#endif //TREAT_SNPS_AS_MS

#ifndef TREAT_SNPS_AS_MS
   
    for (int n1 = 0; n1 < N;  ++n1) {
      for (int c1 = 0; c1 < 2;  ++c1) {
	// initialise veccc
	if (n1 == id ) {
	    /* The first chromosome (c=0) is not allowed to copy 
	       from itself (n1==id)*/
	  for(int t=0; t<SS; t++){
	    veccc[id][0][n1][c1][t] = 0.0;
	    
	    /* The second chromosome (c==1) of id is allowed
	       to copy from the first one (c1==0), but not 
	       from n1 = id and c1 = 1 */
	    if(c1==1)
	      veccc[id][1][n1][c1][t]=0.0;
	    else
	      veccc[id][1][n1][c1][t]=1;
	  }
	}
	else{
	  for(int t=0; t<SS; t++){
	    veccc[id][0][n1][c1][t]=1;
	    veccc[id][1][n1][c1][t]=1;
	  }
	}

	int diff0=0; // count numbers of different SNPS and same SNPS
	int diff1=0; // for chromosomes 0 and 1 of id. 
	int same0=0;
	int same1=0;
	for( vector<int>::const_iterator u = uselist.begin(); u != uselist.end(); ++u ) {
	  from = pop[n1].get_haplotype (c1, *u);
	  targ0 = pop[id].get_haplotype (0, *u);
	  targ1 = pop[id].get_haplotype (1, *u);
	  if(loci_type[*u]!='S')
	    {
	      for(int t=0; t<SS; t++){
		veccc[id][0][n1][c1][t] *= 
		  PrHitTarg(*u, nchr , t, from, targ0, Q);
		veccc[id][1][n1][c1][t] *= 
		  PrHitTarg(*u, nchrplus1 , t, from, targ1, Q);	
	      }
	    }
	  else
	    {
	      diff0 += (from!=targ0);
	      same0 += (from==targ0);
	      diff1 += (from!=targ1);
	      same1 += (from==targ1);
	    }
	  //cout << nchr << " " << diff0 << " " << same0 << " " << diff1 << " " << same1 << endl;
	}
	for(int t=0; t<SS; t++){	  
	  veccc[id][0][n1][c1][t] *= DiffProb(nchr,t,same0,0) * DiffProb(nchr,t,diff0,1);
	  veccc[id][1][n1][c1][t] *= DiffProb(nchrplus1,t,same1,0) * DiffProb(nchrplus1,t,diff1,1);
	}
      }
    }   
#endif // TREAT_SNPS_AS_MS


}


// update the matrix CC for the new phase now at locus
// (involves dividing by the prob of the old type, and
// multiplying by the prob of the new type)
//
void ArrayCC::Update( int id,
                        const vector<ArrayQ *> & Q,
                        const vector<CIndividual> & pop,
                        int locus, int oldtarg0, int oldtarg1)
{
  int from = 0;
  int nchr = N + N - 2;
  int nchrplus1 = nchr + 1;
  
  int targ0 = pop[id].get_haplotype (0, locus);
  int targ1 = pop[id].get_haplotype (1, locus);

  //if targ0=oldtarg0 and targ1=oldtarg1, no update needed; otherwise...
  if((targ0 != oldtarg0) || (targ1 != oldtarg1))
    {
      for(int n1=0;n1<N;n1++){
	for(int c1=0; c1<2 ; c1++){
	  from = pop[n1].get_haplotype (c1, locus);
	  for(int t=0; t<SS; t++){
	    veccc[id][0][n1][c1][t] /= PrHitTarg(locus, nchr , t, from, oldtarg0, Q);
	    veccc[id][0][n1][c1][t] *= PrHitTarg(locus, nchr , t, from, targ0, Q);
	    veccc[id][1][n1][c1][t] /= PrHitTarg(locus, nchrplus1 , t, from, oldtarg1, Q);
	    veccc[id][1][n1][c1][t] *= PrHitTarg(locus, nchrplus1 , t, from, targ1, Q);
	  }
	}
      }
    }
}

//
// Update after a phase flip
//
void ArrayCC::Update( int id,
                        const vector<ArrayQ *> & Q,
                        const vector<CIndividual> & pop,
                        int locus)
{
  Update(id,Q,pop,locus,pop[id].get_haplotype(1,locus),pop[id].get_haplotype(0,locus));
}

//
// compute the prob of the two chromosomes based on CC
//
double ArrayCC::calc_prob(int id)
{
  double firstprob = 0.0; // probability of first chromosome (c=0)
  double secondprob = 0.0; // probability of second chromosome (c=1)
  for(int t=0;t<SS;++t){
    double firstterm=0;
    double secondterm=0;
    for (int n1 = 0; n1 < N; ++n1) {
      for (int c1 = 0; c1 < 2; ++c1) {
	firstterm += veccc[id][0][n1][c1][t];
	secondterm +=veccc[id][1][n1][c1][t];
      }
    }
    firstprob+= WEIGHTS[t]*firstterm;
    secondprob+= WEIGHTS[t]*secondterm;
  }
  return firstprob*secondprob;
}


double CombineProb ( int id, const ArrayCC & CC, const ArrayCC & BB )
{
    int N = CC.get_size();
    double firstprob = 0.0; // probability of first chromosome (c=0)
    double secondprob = 0.0; // probability of second chromosome (c=1)
    for(int t=0;t<SS;++t){
      double firstterm=0;
      double secondterm=0;
      for (int n1 = 0; n1 < N; ++n1) {
        for (int c1 = 0; c1 < 2; ++c1) {	  
	  firstterm += CC(id,0,n1,c1,t) * BB(id,0,n1,c1,t);
	  secondterm += CC(id,1,n1,c1,t) * BB(id,1,n1,c1,t);
	}
      }    
      firstprob+= WEIGHTS[t]*firstterm;
      secondprob+= WEIGHTS[t]*secondterm;
    }
    return firstprob*secondprob;
}

//
// make sure prob of copying n1 from n2 = prob copy n2 from n1.
// called after updating CC for id.
//
void Resymmetrize( ArrayCC & CC, int id)
{
  for(int n1 = 0; n1 < CC.get_size(); n1++){
      if(n1!=id){
	for(int t=0;t<SS;t++)
	  {
	    CC.set_element(n1,0,id,0,t,CC(id,0,n1,0,t));
	    CC.set_element(n1,0,id,1,t,CC(id,1,n1,0,t));
	    CC.set_element(n1,1,id,0,t,CC(id,0,n1,1,t));
	    CC.set_element(n1,1,id,1,t,CC(id,1,n1,1,t));
	  }
      }
  }
}
