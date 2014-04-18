/*
 * $Id: arrayFF.cpp,v 1.12 2003/06/14 00:24:04 stephens Exp $
 */

#include "arrayFF.hpp"
#include "utility.hpp"

using namespace std;

// Constructor
ArrayFF::ArrayFF ( int size ) :
    N ( size ),
    D6 ( SS ),
    D5 ( SS * D6 ),
    D4 ( 2 * D5 ),
    D3 ( N * D4 ),
    D2 ( 2 * D3 ),
    D1 ( N * D2 ),
    vecff ( D1, 0.0 ),
    total_sum ( 0.0 )
{     
}

// Destructor
ArrayFF::~ArrayFF () 
{
}

// This bit is wrong: based on old order of FF, no longer correct
//double & ArrayFF::operator() (int n0, int c0, int t0, 
//                              int n1, int c1, int t1)
//{
//#ifdef CHECK_BOUNDARY
//    assert ( n0 >= 0 && n0 < N );
//    assert ( c0 >= 0 && c0 < 2 );
//    assert ( t0 >= 0 && t0 < SS );
//    assert ( n1 >= 0 && n1 < N );
//    assert ( c1 >= 0 && c1 < 2 );
//    assert ( t1 >= 0 && t1 < SS );
//#endif // CHECK_BOUNDARY    
//    return vecff [ n0 * D2 + c0 * D3 + t0 * D4 
//                 + n1 * D5 + c1 * D6 + t1     ];
//}

// Interface, sample n0-t1
void ArrayFF::sample_froms ( int id,
                             int * pn0, 
                             int * pc0, 
                             int * pt0, 
                             int * pn1, 
                             int * pc1,
                             int * pt1,
                             const vector<ArrayQ *> & Q, const ArrayDiffProb & DiffProb, 
                             const ArrayDiploidDiffProb & DiploidDiffProb,
			     const vector<CIndividual> & pop,
			     const string & loci_type)
{
    // Compute FF first
    compute ( id, Q, DiffProb, DiploidDiffProb, pop, loci_type );

    // Now choose n0 -- t1 at random, according to FF
    int temp = rint2 ( vecff, total_sum );
   
    //cout << "Sum:" << total_sum << endl;
    //for(int i=0; i< vecff.size(); i++)
    //  cout << vecff[i] << " ";
    //cout << endl;

    int ptr = 0;
    for (int n0 = 0; n0 < N; ++n0) {
      for (int c0 = 0; c0 < 2; ++c0) {
	for (int n1 = 0; n1 < N; ++n1) {
	  for (int c1 = 0; c1 < 2; ++c1) {
	    for (int t0 = 0; t0 < SS; ++t0) {
	      for (int t1 = 0; t1 < SS; ++t1) {
		if ( ptr == temp ) {
		  *pn0 = n0;
		  *pn1 = n1;
		  *pc0 = c0;
		  *pc1 = c1;
		  *pt0 = t0;
		  *pt1 = t1;
		  goto ENDLOOP;
		} else {
		  ++ptr;
		}
	      }
	    }
	  }
	}
      }
    }
 ENDLOOP:
    return;
}

// Private functions

/** Compute FF 

    FF[n0][c0][t0][n1][c1][t1] is the probability that the two
    chromosomes of individual id are copied from chromosome c0 of
    individual n0, with time elapse t0, and from chromosome c1 of n1,
    with time t1, given the observed genotypes of id.
*/

void ArrayFF::compute ( int id,
                        const vector<ArrayQ *> & Q, const ArrayDiffProb & DiffProb, 
			const ArrayDiploidDiffProb & DiploidDiffProb,
                        const vector<CIndividual> & pop, const std::string & loci_type )
{
    int ptr = 0;
    int from0 = 0;
    int from1 = 0;
    int targ0 = 0;
    int targ1 = 0;
    ArrayCC CC(N);
    int nchr=N+N-2;
    int nchrplus1 = N+N-1;
    int Nindminus1 = N-1;
	
    // static int ctype[2][2][2][2]={0,1,1,2,1,3,3,1,1,3,3,1,2,1,1,0};    
    // look-up table for what kind of diploid copy is being made
    // (for SNPs only)
    //
    // from0 from1 targ0 targ1 ctype
    // 0     0     0     0     0
    // 0     0     0     1     1
    // 0     0     1     0     1
    // 0     0     1     1     2
    // 0     1     0     0     1
    // 0     1     0     1     3
    // 0     1     1     0     3
    // 0     1     1     1     1
    // 1     0     0     0     1
    // 1     0     0     1     3
    // 1     0     1     0     3
    // 1     0     1     1     1
    // 1     1     0     0     2
    // 1     1     0     1     1
    // 1     1     1     0     1
    // 1     1     1     1     0


    total_sum=0;
     
    CC.compute(id,Q,pop,pop[id].get_known_pos(),loci_type,DiffProb);
    int num_unknown= pop[id].numunknown();

    vector<int> diff(4);
    
    int r,ctype;

    for (int n0 = 0; n0 < N;  ++n0) {        
    for (int c0 = 0; c0 < 2;  ++c0) {
    for (int n1 = 0; n1 < N;  ++n1) {
    for (int c1 = 0; c1 < 2;  ++c1) {

      // compute number of differences of each type
      diff[0]=diff[1]=diff[2]=diff[3]=0;
      if(n1!=id){
	for (int u = 0; u < num_unknown; ++u) {	
	  r=pop[id].get_unknown_pos(u);
	  if((loci_type[r] == 'S') && (pop[id].n_missing(r)==0)){
	  from0 = pop[n0].get_haplotype(c0, r);
	  from1 = pop[n1].get_haplotype(c1, r);
	  targ0 = pop[id].get_allele(0, r);
	  targ1 = pop[id].get_allele(1, r);	  
	  ctype=1; // this is the type of copy that is occurring (0 - 3)
	  // see look-up table above (but actually more efficient to calculate
	  // on the fly.)
	  if((from0==from1) && (targ0==targ1)){
	    if(from0==targ0)
	      ctype=0;
	    else
	      ctype=2;
	  } else{
	    if ((from0!=from1) && (targ0 !=targ1))
	      ctype=3;
	  }	
	  diff[ctype]++;	
	  }	  
	}
      }

      for (int t0 = 0; t0 < SS; ++t0) {
	for (int t1 = 0; t1 < SS; ++t1) {  
    
	  vecff[ptr]=WEIGHTS[t0] * WEIGHTS[t1] * CC(id,0,n0,c0,t0) * CC(id,1,n1,c1,t1) *
	    DiploidDiffProb(Nindminus1,t0,t1,diff[0],0) *
	    DiploidDiffProb(Nindminus1,t0,t1,diff[1],1) *
	    DiploidDiffProb(Nindminus1,t0,t1,diff[2],2) *
	    DiploidDiffProb(Nindminus1,t0,t1,diff[3],3);
	    
	  for (int u = 0; u < num_unknown; ++u) {	
	    r=pop[id].get_unknown_pos(u);
     
	    // Number of missing alleles at r
	    switch ( pop[id].n_missing(r) ) {
	    case 0: 
	      from0 = pop[n0].get_haplotype(c0, r);
	      targ0 = pop[id].get_allele(0, r);
	      targ1 = pop[id].get_allele(1, r);
	      if(loci_type[r] != 'S'){
		if (n1 == id) {
		  vecff[ptr] *=
		    PrHitTarg(r, nchr, t0, from0, targ0, Q) *
		    PrHitTarg(r, nchrplus1, t1, targ0, targ1, Q) + 
		    PrHitTarg(r, nchr, t0, from0, targ1, Q) *
		    PrHitTarg(r, nchrplus1, t1, targ1, targ0, Q);	         
		} else {
		  from1 = pop[n1].get_haplotype(c1, r);
		  vecff[ptr] *=
		    PrHitTarg(r, nchr, t0, from0, targ0, Q) *
		    PrHitTarg(r, nchrplus1, t1, from1, targ1, Q) + 
		    PrHitTarg(r, nchr, t0, from0, targ1, Q) *
		    PrHitTarg(r, nchrplus1, t1, from1, targ0, Q);
		}
	      }
	      else{ // for SNP loci
		if ( n1 == id ) {
		  vecff[ptr] *=
		    PrHitTarg(r, nchr, t0, from0, targ0, Q) *
		    PrHitTarg(r, nchrplus1, t1, targ0, targ1, Q) + 
		    PrHitTarg(r, nchr, t0, from0, targ1, Q) *
		    PrHitTarg(r, nchrplus1, t1, targ1, targ0, Q);	         
		} 
	      }
	      break;
	      
	      
	    case 1: 
	      if ( n1 != id ) {
		/** One allele is missing
		 *  allele 0 is observed and allele 
		 *  1 is imputed so is ignored. 
		 */
		targ0 = pop[id].get_allele(0, r);
		from1 = pop[n1].get_haplotype(c1, r);
		from0 = pop[n0].get_haplotype(c0, r);
		vecff[ptr] *= 
		  PrHitTarg(r, nchr, t0, from0, targ0, Q) +
		  PrHitTarg(r, nchr, t1, from1, targ0, Q);
	      } else {
		vecff[ptr] = 0.0;
	      }
	      break;
	      
	    default:                                    
	      // Both alleles are missing
	      // Do nothing
	      break;
	    }
	  }	  
	  
	  total_sum += vecff[ptr++];
	  
	}
      }
    }
    }
    }
    }   
}




// {{{ Log

/*
 * $Log: arrayFF.cpp,v $
 * Revision 1.12  2003/06/14 00:24:04  stephens
 * Adding files, and committing a lot of changes
 * that have resulted in version 2.0 of PHASE
 *
 * Revision 1.11  2001/06/23 15:56:59  stephens
 * started to check missing data code. Corrected some bugs. Seems to give sensible results in a couple of small examples.
 *
 * Revision 1.10  2001/06/19 17:02:25  stephens
 * Changes to computation of arrayFF to make more efficient
 * Added facility to store "original phenotype" in indnode,
 * in preparation for allowing genotyping error.
 *
 * Revision 1.9  2001/05/31 16:26:13  stephens
 * Added DiploidDiffProb look-up table class (to be used to make
 * computation of arrayFF for SNPs more efficient)
 *
 * Revision 1.8  2001/05/30 06:02:18  stephens
 * Updated to be considerably more efficient, via introduction of
 * various new methods, including introduction of ArrayDiffProb
 * and ArrayDiffCount to improve computation al efficiency for SNP
 * data. Speedtest.inp now runs in about 7secs.
 * Also corrected several bugs. Output now looks more promising
 * and major bugs appear to have been eliminated. Convergence of chain
 * can now be monitored more easily by the output in temp.monitor, which
 * gives the pseudo-likelihood every Nthin repetitions.
 *
 * Revision 1.7  2001/05/23 02:16:55  stephens
 * Corrected some bugs in computation of CC. Reduced number of
 * loci we update by one when the total number of unknowns is <=5.
 * Confirmed program was giving vaguely sensible resutls, and is
 * competitive on the speed test example (36secs vs 31 for old version)
 *
 * Revision 1.6  2001/05/18 21:47:42  stephens
 * added facility to update just 5 at a time. Set it to do this
 * so we can test for speed with the current version.
 *
 * Revision 1.5  2001/05/08 16:58:23  stephens
 * Added class arrayCC, which computes quantities similar
 * to arrayFF, but only for sites whose phase is known.
 * This is then used to improve the efficiency of computation
 * for arrayFF by re-using the calculations for known sites.
 * Checked that the output was unchanged on an example.
 * Might slightly improve efficiency of computation of CC if
 * we do the loop through known sites more efficiently.
 *
 * Revision 1.4  2001/04/26 18:29:51  stephens
 * Fixed bug in computation of ArrayFF (total_sum now initialized
 * to 0 each time FF is computed)
 *
 * Revision 1.3  2001/04/24 19:37:07  nali
 * Minor changes.
 *
 * Revision 1.2  2001/04/20 00:32:45  nali
 * Put reading loci types into the contructor of ClassPop
 *
 * Revision 1.1  2001/04/19 19:43:02  nali
 * Class for FF array.
 *
 */

// }}}
