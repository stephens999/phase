/* Constants Definitions
 *
 * $Id: constants.hpp,v 1.17 2003/06/14 00:24:05 stephens Exp $
 */
  
#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

// Definition for constants


const int BIGINT    = 9999;     // max size for MS allele
const int NMAXD     = 4;        // used when outputting ids 
const int KMAX      = 50;       // Number of alleles
const int ITPRINT    = 10;      // how often to print iteration number

const int MAXLISTLENGTH = 11;

const double EPSILON = 1e-100; 
const double EMCONVERGE = 1e-4;
const double BIGNEGATIVE = -1.0e300;
const int RANDMAX = 2147483647;

const double DIRPRIOR = 1e-100; // parameter of dirichlet prior for Naive Gibbs
const double RHOMULTSIGMA = 1.15; // sd of prior on log(RhoMult)
const double MAXRHOMEAN = 1000;
const double MINRHOMEAN = 1E-8;

// Representation of missing alleles
const char MISSNP = '?';
const int  MISSMS = -1;

// Number of Ancestral Haps used in AncHap method
const int NANCHAP = 15;

const char POSITIONLINEINDICATOR = 'P';

// space holder when print out results for unambiguous sites
const char SPACEHOLDER = '='; 
const char MISSCHAR = '?'; // char to use to replace missing sites

const char UNKNOWNPHASECHAR = '*';

// Vector of points in quadrature; Laguerre polynomials, from 
// appendix B of Evans;


const int SS        = 2;        // Number of quadrature points
const int NPOWER    = 50;       // Number of powers to use in approximation
//  const double TIMES   [] = { 0.322547689619392312,
//                              1.74576110115834658,
//                              4.53662029692112798,
//                              9.39507091230113313 };
//  const double WEIGHTS [] = { 0.603154104341633602,
//                              0.357418692437799687,
//                              0.0388879085150053843,
//                              0.000539294705561327450 };

const double TIMES []={0.58578643,3.41421356}; //vector of points in quadrature; Laguerre polynomials, from Appendix B of Evans;
const double WEIGHTS []={0.85355339,0.14644661};


#endif 

// {{{ Log
/* 
   $Log: constants.hpp,v $
   Revision 1.17  2003/06/14 00:24:05  stephens
   Adding files, and committing a lot of changes
   that have resulted in version 2.0 of PHASE

   Revision 1.16  2002/02/27 18:56:45  stephens
   Commiting the current source, which is essentially that released as
   PHASE 1.0. The number of quadrature points is reduced to 2. Some code has
   been added to cope with recombination, but this is not included in the release
   and is still under development.

   Revision 1.15  2001/10/21 18:33:29  stephens
   fixed bug in GibbsUpdate

   Revision 1.13  2001/10/16 20:44:16  stephens
   Added way to input position of each marker, by preceding line
   of SSMMMSS with a line beginning with P followed by position of
   each marker.

   Revision 1.12  2001/10/16 18:12:06  stephens
   More minor bugs fixed. (including -m bug)

   Revision 1.11  2001/10/12 23:49:26  stephens
   Various updates, particularly cosmetic, but some bug fixes.
   this version tested on MS and SNP data without missing alleles
   gives very similar answers to the original phase.

   Revision 1.10  2001/10/02 18:02:59  stephens
   just some small additions. haven't committed these for a long time
   because cvs was not working propersly. (now know why: they remounted
   you from /user3 to /home)

   Revision 1.9  2001/04/24 19:37:07  nali
   Minor changes.

   Revision 1.8  2001/04/19 19:44:28  nali
   New constant SPACEHOLDER defined for output phases and haplotypes.

   Revision 1.7  2001/04/09 16:28:35  nali
   Most part of original ResolvePhase (without recombination) implemented.
   Now compiles and runs.

   Revision 1.6  2001/02/28 04:54:31  nali
   Make EM working and haplotype list right

   Revision 1.5  2001/02/20 16:15:29  nali
   Use assert ( )  for debugging.

   Revision 1.4  2001/02/16 17:13:17  nali
   New functions: InputRandom, InputHusdonData, etc..
   New member function: make_haplist.

   Revision 1.3  2001/02/16 02:41:48  nali

   Remove the character representation of phenotypes.

   Revision 1.2  2001/02/14 23:39:49  nali
   New constant BIGINT added.

   Revision 1.1  2001/02/13 06:02:59  nali
   Global constances

*/
// }}}
