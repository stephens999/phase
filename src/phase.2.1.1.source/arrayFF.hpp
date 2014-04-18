/*
 * Define a class for array FF
 *
 * $Id: arrayFF.hpp,v 1.9 2003/06/14 00:24:04 stephens Exp $
 */

#ifndef ARRAY_FF_H
#define ARRAY_FF_H

#include "constants.hpp"
#include "arrayQ.hpp"
#include "arrayDiffProb.hpp"
#include "arrayDiploidDiffProb.hpp"
#include "indnode.hpp"
//#include "crand/ranlib.h"
#include "arrayCC.hpp"
#include "utility.hpp"

#include <cmath>
#include <vector>
using namespace::std;

class ArrayFF {
    int N;                        // total number of individuals
    int D6;  
    int D5;
    int D4;
    int D3;
    int D2;
    int D1;                      // total length of the vector
    std::vector<double> vecff;
    double total_sum;

    // Compute FF for id
    void compute   ( int id, 
                     const vector<ArrayQ *> &, const ArrayDiffProb &,
		     const ArrayDiploidDiffProb &,
                     const std::vector<CIndividual> &, const std::string & );

    // copy constructor disabled.
    // Can only declare a vector of pointers to FF
    ArrayFF ( const ArrayFF & ) {}
public:
    ArrayFF ( int );
    ~ArrayFF ();
    
    /** Direct access function to the elements of FF 
     */
  // double & operator() (int, int, int, int, int, int);

    /** Sample n0, c0, t0, n1, c1, t1 according to the probabilties
        proportional to
        Pr( observe id genotype | n0, c0, t0, n1, c1, t1 )
    */
    void sample_froms ( int id,
                        int *, int *, int *,
                        int *, int *, int *,
                        const vector<ArrayQ *> &, const ArrayDiffProb &,
			const ArrayDiploidDiffProb &,
                        const std::vector<CIndividual> &, const std::string & );
};


#endif // ARRAY_FF_H

// {{{ Log
/*
 * $Log: arrayFF.hpp,v $
 * Revision 1.9  2003/06/14 00:24:04  stephens
 * Adding files, and committing a lot of changes
 * that have resulted in version 2.0 of PHASE
 *
 * Revision 1.8  2001/10/12 23:49:25  stephens
 * Various updates, particularly cosmetic, but some bug fixes.
 * this version tested on MS and SNP data without missing alleles
 * gives very similar answers to the original phase.
 *
 * Revision 1.7  2001/06/19 17:02:25  stephens
 * Changes to computation of arrayFF to make more efficient
 * Added facility to store "original phenotype" in indnode,
 * in preparation for allowing genotyping error.
 *
 * Revision 1.6  2001/05/30 06:02:18  stephens
 * Updated to be considerably more efficient, via introduction of
 * various new methods, including introduction of ArrayDiffProb
 * and ArrayDiffCount to improve computational efficiency for SNP
 * data. Speedtest.inp now runs in about 7secs.
 * Also corrected several bugs. Output now looks more promising
 * and major bugs appear to have been eliminated. Convergence of chain
 * can now be monitored more easily by the output in temp.monitor, which
 * gives the pseudo-likelihood every Nthin repetitions.
 *
 * Revision 1.5  2001/05/18 21:47:42  stephens
 * added facility to update just 5 at a time. Set it to do this
 * so we can test for speed with the current version.
 *
 * Revision 1.4  2001/05/08 16:58:23  stephens
 * Added class arrayCC, which computes quantities similar
 * to arrayFF, but only for sites whose phase is known.
 * This is then used to improve the efficiency of computation
 * for arrayFF by re-using the calculations for known sites.
 * Checked that the output was unchanged on an example.
 * Might slightly improve efficiency of computation of CC if
 * we do the loop through known sites more efficiently.
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
