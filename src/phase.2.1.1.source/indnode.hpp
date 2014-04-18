/*
 * CVS :  $Id: indnode.hpp,v 1.41 2003/06/14 00:24:05 stephens Exp $
 */
 
#ifndef CLASS_INDIVIDUAL_H
#define CLASS_INDIVIDUAL_H

#include "constants.hpp"
#include "arrayQ.hpp"
#include "Haplotype.hpp"
#include "Summary.hpp"


#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <functional>
using namespace::std;

class CIndividual {

    // Static members
    static int counts;                   // Number of individuals
  
    string id;                           // id

    vector<Haplotype> phenotype; // stores current guess of genotype
    vector<vector<int> > orig_phenotype; // stores original genotype 
    vector<vector<int> > saved_hap; // stores an old guess for the phenoty
 
    vector<vector<int> > Z; // stores which anc hap copied by this ind

    vector< vector< vector<float> > > AlleleCount;
    vector< vector<float> > PhaseCount;
  
    // phase[r] = {0, 1} the phase of the r th locus 
    std::vector<int> phase;


    /* missing[r] = {0, 1, 2} for r = 0 .. (nloci - 1)
       number of missing alleles at locus r.
       If there is one missing allele, then phenotype[r][0]
       is observed and phenotype[r][1] is missing  */
  std::vector<int> missing;
  std::vector<int> recom; // list of sites after which a recom occurs in transmission from this ind
  // (used in trios case)

  std::vector<int>  notmissing_list; //list of loci with no missing data
   
    /* unknown[i] = {0 .. (nloci - 1)}, index of the i th
       unknown phase locus */   
    std::vector<int>  unknown; // phase unknown loci 
    std::vector<int>  known; // phase known loci
    double overall_correct_prob; // overall prob that phase is correct (used by HapList2)


public:
  const CIndividual & operator= (const CIndividual & );

    CIndividual ( string loci_type="" );                  // Constructor
    CIndividual ( const CIndividual & );  // Copy constructor
    CIndividual ( const CIndividual &, int, int ); // individual, between 2 loci 
  CIndividual ( const CIndividual &, const CIndividual &); // concatenate 2 inds 
    ~CIndividual ();                      // Destructor    

    // ID accessors
    std::string get_id () const;
    void set_id ( std::istream & );
    void set_id ( const std::string & nid );
    void print_id (std::ostream &) const;
    void set_id ();                         

  // Z setter
  void set_Z( int chr, int locus, int value ); // inline below
  int get_Z( int chr, int locus ); // inline below

    // Allele
    int  get_allele ( int chr, int locus ) const;    
    int get_orig_allele (int chr, int locus ) const;
    int get_orig_nonmissing_allele (int locus) const;
    int missingchr (int locus) const;

    void set_allele ( int chr, int locus, int allele );   // sets current alleles
    void set_original_allele ( int chr, int locus, int allele );   // sets original alleles

    void input_orig_allele ( istream &, char locus_type,
                      int chr, int locus);
    // Return number of missing alleles at locus
    int n_missing ( int locus ) const;
    const int NMissingLoci () const; // inline

    // Get haplotype
    int              get_haplotype ( int c, int locus ) const;
  //std::vector<int> get_haplotype ( int c ) const;
    Haplotype get_haplotype (int c ) const;

  void update_haplotype( int chr, int locus, int allele ); // changes current allele for haplotype chr (leaves original unchanged)
  void update_haplotypes( const Haplotype & h0, const Haplotype & h1); 
  void update_haplotype( int chr, const Haplotype & h); // changes current
  //haplotype for chromosom chr (inline)

    // Get phase information
    
    const std::vector<int> & get_unknown_pos ( ) const; 
    const std::vector<int> & get_recom() const;
    const std::vector<int> & get_known_pos ( ) const;
   
   
  std::vector<int>  get_notmissing_list ();
  
  const vector<int> & get_nmissing() const;
   const int nmissing(int) const;
    
  char get_locus_type( int locus) const; // inline
    int              get_unknown_pos ( int count ) const;
    std::vector<int> get_phase ( ) const;
    int              get_phase ( int locus ) const;
    void        print_phase ( std::ostream &, bool = false ) const;
  void         print_phase_prob ( std::ostream &, bool = false ) const;
    void        set_phase ( int locus, int ph );

    void flip_phase();
    void flip_phase(int);
  
    int get_nloci () const;

    // Number of positions with unknown phase
    int numunknown () const;
    bool is_unknown (int) const;

  // number of possible phase reconstructions
    int numphase() const;


    // Phenotype IO
    int read_orig_phenotypes  ( std::istream &,
                           const std::string &,
                           int, int );
    
    void print_phenotypes ( std::ostream &,                
                            const std::string &,
                            const std::vector<int> * ) const;
    /** 
     * Print out both haplotypes, ignore loci with no ambiguity
     * (both alleles are observed and are the same).
     */
    void print_allele ( std::ostream &, int, int,             
                            const std::string &,
                            const std::vector<int> *,
                            bool=true, bool=true, bool=false, double=0.5, double=0.0) const;
  
    void print_haplotypes ( std::ostream &,                
                            const std::string &,
                            const std::vector<int> *,
                            bool, bool, bool, double, double) const;
  
   void print_haplotype ( int, std::ostream &,                
                            const std::string &,
                            const std::vector<int> *,
                            bool=true, bool=true, bool=false, double=0.5, double=0.0) const;
  
    // Initialize unknown phases and missing alleles
    // Both parameters can be omitted. 
    // However when startinfo = 1 or 2, the second parameter
    // must be valid input stream
    int initialize ( int knowninfo , std::istream &, int initmetho, std::istream &, const string & loci_type);
 
    double flipprob(); // computes the probability of the current phases
  // using the current phase counts

  void ResetCounts(); // resets the counts of phases and alleles (to 1)
    void UpdateCounts(); // updates the counts of phases and alleles 
  void TransferCounts(Summary & sum);
    
    double BestPhaseProb(int) const; // return probability of best phase at given pos
    int BestAllele(int c,int locus) const; // return best allele for chromosome chr at position locus, based on AlleleCount
    double BestAlleleProb(int c,int locus) const; // return prob of best allele, based on AlleleCount

  int BestPhase(int) const;
  int BestHaplotype(int,int) const;
  Haplotype BestHaplotype(int) const;

  void SaveCurrentState();
  void RestoreSavedState();
  void set_phenotypes();
  void set_overall_correct_prob(double p); // inline
  void output_overall_correct_prob(ostream & ostr); //inline
  double ObservedDataProbGivenParents(const Haplotype & h1, const Haplotype & h2, const Haplotype & h3, const Haplotype & h4, const vector<int> & recom1, const vector<int> & recom2);
  double PrOrigGenotypeData(int locus, int allele0, int allele1);
};

Haplotype GetCompHap(const Haplotype & h, const CIndividual & ind, bool & found, bool checkmissing = true); // returns Complementary haplotype in individual ind, if can be found (otherwise found is returned as false)

Haplotype GetCompHap(const Haplotype & h, const CIndividual & ind, bool & found,vector<int> & matches, bool checkmissing = true); // returns Complementary haplotype in individual ind, if can be found (otherwise found is returned as false)
//matches says which chrom in ind matches h at each unknown position

// return whether h can possibly be found in pop
bool CanBeFoundAtAll(const Haplotype & h, const vector<CIndividual> & pop);

int NDiff(const std::vector<CIndividual> & pop, int n0, int c0, int n1, int c1, const std::vector<int> & uselist); // finds number of diffs

double CCProb(const std::vector<CIndividual> & pop, int n0, int c0, int n1, int c1, int t, int nchr, const vector<ArrayQ *>  & Q, 
const std::vector<int> & uselist); 
// compute the probability of seeing n0,c0 
// given it copied n1,c1, separated by time t, 
// for loci in uselist


std::vector<int> subrange(const vector<int> & vec,int min, int max);

// Inline member functions

inline void CIndividual::set_overall_correct_prob(double p){
  overall_correct_prob  = p;
}

inline void CIndividual::output_overall_correct_prob(ostream & ostr){
  ostr << overall_correct_prob;
}

inline void CIndividual::SaveCurrentState(){
  for(int chr=0;chr<2;chr++){
    for(int r=0; r<phase.size();r++){
      saved_hap[chr][r] = get_haplotype(chr,r);
    }
  }
}

inline void CIndividual::RestoreSavedState(){
  for(int chr=0;chr<2;chr++){
    for(int r=0; r<phase.size();r++){
      update_haplotype(chr,r,saved_hap[chr][r]);
    }
  }
}

inline std::string CIndividual::get_id () const { return id; } 
inline void CIndividual::set_id ( istream & istr ) { istr >> id; }
inline void CIndividual::set_id ( const string & nid ) { id = nid; }

// Set Z value
inline void CIndividual::set_Z ( int chr, int locus, int a ) 
{ 
  Z[chr][locus] = a;
}

inline char CIndividual::get_locus_type( int locus) const
{
  return phenotype[0].get_locus_type()[locus];
}

inline int CIndividual::get_Z ( int chr, int locus) 
{ 
  return Z[chr][locus];
}


// Get allelic types
inline int CIndividual::get_allele ( int chr, int locus ) const { 
    return phenotype[chr].get_allele(locus); 
}
 
inline int CIndividual::get_orig_allele ( int chr, int locus ) const { 
    return orig_phenotype[chr][locus]; 
}
 
inline int CIndividual::get_orig_nonmissing_allele ( int locus ) const { 
  return (orig_phenotype[0][locus] != MISSMS) ? orig_phenotype[0][locus] : orig_phenotype[1][locus]; 
}
 
inline int CIndividual::missingchr ( int locus ) const { 
    return (orig_phenotype[0][locus] == MISSMS) ? 0 : 1; 
}
 
// Set allelic types
inline void CIndividual::set_allele ( int chr, int locus, int allele ) 
{ 
  phenotype[chr].set_allele(locus,allele);
  //orig_phenotype[chr][locus] = allele;  
}

inline void CIndividual::set_original_allele ( int chr, int locus, int allele ) 
{ 
  orig_phenotype[chr][locus] = allele;  
}



inline void CIndividual::update_haplotype ( int chr, int locus, int allele ) 
{ 
  //cout << "Setting: " << chr << "," << locus << "," << allele << endl;
  //cout << "Phase = " << phase[locus] << endl;
  if(chr==0)
    phenotype[phase[locus]].set_allele(locus,allele);
  else
    phenotype[1-phase[locus]].set_allele(locus,allele);
}

inline void CIndividual::update_haplotype ( int chr, const Haplotype & h) 
{ 
  int Nloci = phase.size(); 
  
  for(int locus = 0; locus < Nloci; locus++){
    if((phase[locus] != 0) && (phase[locus]!=1)){
      cerr << "Error: phase is not right at locus " << (locus+1) << endl;
      cout << phase[locus] << endl;
      cout << missing[locus] << endl;
      exit(1);
    }
    if(phase[locus] == 0)
      phenotype[chr].set_allele(locus,h.get_allele(locus));
    else
      phenotype[1-chr].set_allele(locus,h.get_allele(locus));
  }
}

inline void CIndividual::update_haplotypes ( const Haplotype & h0, const Haplotype & h1) 
{ 
  bool consistentwithknown = true;
  int Nloci = phase.size();
  for(vector<int>::iterator u = known.begin(); u!=known.end(); u++){
    if( ((h0.get_allele(*u) != orig_phenotype[0][*u]) && (orig_phenotype[0][*u]!=MISSMS)) || 
	((h1.get_allele(*u) != orig_phenotype[1][*u]) && (orig_phenotype[1][*u]!=MISSMS)) ){
      consistentwithknown = false;
      break;
    }
  }
  // cout << "Here!" << endl;
  // cout << "consistent with known = " << consistentwithknown << endl;
  if(consistentwithknown){ // set chr0 to be h0 and chr1 to be h1
    for(int locus = 0; locus < Nloci; locus++){
      // cout << "locus = " << locus << endl;
      // cout << "phase = " << phase[locus] << endl;
      // cout << "allele0 = " << h0.get_allele(locus) << endl;
      //	cout << "allelel = " << h1.get_allele(locus) << endl;
 
      phenotype[phase[locus]].set_allele(locus,h0.get_allele(locus));
      phenotype[1-phase[locus]].set_allele(locus,h1.get_allele(locus));
    }
  }
  else { // must switch to be consistent with known (chr0 = h1, chr1 = h0)
    for(int locus = 0; locus < Nloci; locus++){
      phenotype[phase[locus]].set_allele(locus,h1.get_allele(locus));
      phenotype[1-phase[locus]].set_allele(locus,h0.get_allele(locus));
    }
  }
}

// haplotypes
inline int CIndividual::get_haplotype ( int c, int locus ) const 
{
  if(c==0){
    return phenotype[phase[locus]].get_allele(locus);}
  else{
    return phenotype[1-phase[locus]].get_allele(locus);}
}

inline int CIndividual::get_nloci () const 
{
  return phase.size();
}

inline int CIndividual::n_missing ( int locus ) const 
{ 
    assert ( locus < missing.size() );
    assert ( locus >= 0 );
    return missing[locus]; 
}

inline const int CIndividual::NMissingLoci () const 
{ 
  return phase.size() - notmissing_list.size();
}



inline std::vector<int> CIndividual::get_notmissing_list ( )
{ 
  return notmissing_list; 
}



// Get list of unknown positions
inline const std::vector<int> & CIndividual::get_unknown_pos ( ) const 
{ 
    return unknown; 
}

inline const std::vector<int> & CIndividual::get_recom ( ) const 
{ 
    return recom; 
}

inline const std::vector<int> & CIndividual::get_nmissing ( ) const 
{ 
    return missing; 
}

inline const int CIndividual::nmissing (int locus ) const 
{ 
    return missing[locus]; 
}

inline const std::vector<int> & CIndividual::get_known_pos ( ) const 
{
    return known;
}


// return whether locus has unknown phase
inline bool CIndividual::is_unknown (int locus) const
{
  return std::binary_search ( unknown.begin(), unknown.end(), locus );
}


inline int CIndividual::get_unknown_pos ( int count ) const 
{
    assert ( count < unknown.size() );
    assert ( count >= 0 );
    return unknown[count];
}



inline std::vector<int> CIndividual::get_phase ( ) const 
{ 
    return phase; 
}
    
inline int CIndividual::get_phase ( int locus ) const 
{ 
    return phase[locus];
}

inline void CIndividual::set_phase (int locus, int ph) 
{
    phase[locus] = ph;
}

inline void CIndividual::flip_phase () 
{
  for(int locus=0; locus < phase.size(); locus++)
    phase[locus] = 1-phase[locus];
}

inline void CIndividual::flip_phase (int locus) 
{
  phase[locus] = 1-phase[locus];
}

inline int CIndividual::numunknown () const 
{ 
    return unknown.size(); 
}


inline int CIndividual::numphase () const
{
  int u=unknown.size();
  if(u>1)
    return 1<<(u-1); // number of phases for individual with 
  else
    return 1;
}

inline void CIndividual::set_phenotypes() 
{
  int Nloci = phase.size();
  for(int r = 0; r<Nloci; r++){
    phenotype[0].set_allele(r, orig_phenotype[0][r]);
    phenotype[1].set_allele(r, orig_phenotype[1][r]);
  }

}

inline double CIndividual::PrOrigGenotypeData(int locus, int allele0, int allele1)
{
   if(missing[locus] == 2) return 1;
   if(missing[locus] == 1) return 0.5* ((allele0==get_orig_nonmissing_allele(locus)) + (allele1==get_orig_nonmissing_allele(locus)));
   return ((allele0+allele1) == (orig_phenotype[0][locus] + orig_phenotype[1][locus]));


}
#endif // CLASS_INDIVIDUAL_H


// {{{ Log  
/* 
   $Log: indnode.hpp,v $
   Revision 1.41  2003/06/14 00:24:05  stephens
   Adding files, and committing a lot of changes
   that have resulted in version 2.0 of PHASE

   Revision 1.40  2002/02/27 18:56:45  stephens
   Commiting the current source, which is essentially that released as
   PHASE 1.0. The number of quadrature points is reduced to 2. Some code has
   been added to cope with recombination, but this is not included in the release
   and is still under development.

   Revision 1.39  2001/10/21 18:33:29  stephens
   fixed bug in GibbsUpdate

   Revision 1.37  2001/10/16 20:44:16  stephens
   Added way to input position of each marker, by preceding line
   of SSMMMSS with a line beginning with P followed by position of
   each marker.

   Revision 1.36  2001/10/16 18:12:06  stephens
   More minor bugs fixed. (including -m bug)

   Revision 1.34  2001/10/06 00:35:40  stephens
   Impute missing positions added- changed update_allele to update_haplotype.

   Revision 1.33  2001/10/02 18:02:59  stephens
   just some small additions. haven't committed these for a long time
   because cvs was not working propersly. (now know why: they remounted
   you from /user3 to /home)

   Revision 1.32  2001/06/23 15:57:00  stephens
   started to check missing data code. Corrected some bugs. Seems to give sensible results in a couple of small examples.

   Revision 1.31  2001/06/19 17:02:25  stephens
   Changes to computation of arrayFF to make more efficient
   Added facility to store "original phenotype" in indnode,
   in preparation for allowing genotyping error.

   Revision 1.30  2001/05/30 06:02:18  stephens
   Updated to be considerably more efficient, via introduction of
   various new methods, including introduction of ArrayDiffProb
   and ArrayDiffCount to improve computation al efficiency for SNP
   data. Speedtest.inp now runs in about 7secs.
   Also corrected several bugs. Output now looks more promising
   and major bugs appear to have been eliminated. Convergence of chain
   can now be monitored more easily by the output in temp.monitor, which
   gives the pseudo-likelihood every Nthin repetitions.

   Revision 1.29  2001/05/23 02:16:56  stephens
   Corrected some bugs in computation of CC. Reduced number of
   loci we update by one when the total number of unknowns is <=5.
   Confirmed program was giving vaguely sensible resutls, and is
   competitive on the speed test example (36secs vs 31 for old version)

   Revision 1.28  2001/05/21 20:16:34  nali

   Add another member vector<int> known for phase known loci, and
   corresponding accessor and initialization.

   Revision 1.27  2001/05/18 21:47:42  stephens
   added facility to update just 5 at a time. Set it to do this
   so we can test for speed with the current version.

   Revision 1.26  2001/05/08 16:58:23  stephens
   Added class arrayCC, which computes quantities similar
   to arrayFF, but only for sites whose phase is known.
   This is then used to improve the efficiency of computation
   for arrayFF by re-using the calculations for known sites.
   Checked that the output was unchanged on an example.
   Might slightly improve efficiency of computation of CC if
   we do the loop through known sites more efficiently.

   Revision 1.25  2001/04/26 18:29:51  stephens
   Fixed bug in computation of ArrayFF (total_sum now initialized
   to 0 each time FF is computed)

   Revision 1.24  2001/04/26 17:09:17  stephens
   Added function ClassPop::print_allele
   and made a few minor changes

   Revision 1.23  2001/04/24 22:00:08  stephens
   Added flag in ClassPop::output_hap and output_phase
   to indicate whether to output all positions, or only
   unknown positions.

   Revision 1.22  2001/04/24 19:34:58  nali
   Merge several print_haplotypes member functions into one with an option of
   printing all alleles or not (default not).
   Add default argument 0 to the constructor so default constructor can be called.

   Revision 1.21  2001/04/20 00:32:45  nali
   Put reading loci types into the contructor of ClassPop

   Revision 1.20  2001/04/19 19:47:35  nali
   Output SPACEHOLDER when there is no ambiguity in haplotype or phase.

   Revision 1.19  2001/04/17 22:08:44  nali

   Major revisement in overall structure. Created new class ClassPop and
   almost all global functions now became member functions of ClassPop, most
   of them private.

   "mult" removed in update_phase_NR. No other changes in terms of algorithm.
   Haven't check the results yet.

   proc_args() is moved to utility.cpp, which also defines a couple of other
   global functions.

   Revision 1.18  2001/04/12 17:13:10  stephens
   edited haplotype output routines, and added line to output
   haplotypes to output file. Output from test files indicates
   a bug as output haplotypes are inconsistent with input genotypes

   Revision 1.17  2001/04/12 01:53:26  stephens
   Added -D option for reading in several datasets
   Debugged -H option for inputting data from Hudson's simulations
   Reduced Nthin to 1 to reduce time for examples, and
   added a test example for the -H option and -D options

   Revision 1.16  2001/04/09 16:28:35  nali
   Most part of original ResolvePhase (without recombination) implemented.
   Now compiles and runs.

   Revision 1.15  2001/04/07 07:18:10  nali
   Revised to concentrate on implementing SD method.
   Discarded support for EM algorithm, etc..

   Revision 1.14  2001/04/04 06:26:48  nali
   OutFreq compiled.

   Revision 1.13  2001/02/28 04:54:32  nali
   Make EM working and haplotype list right

   Revision 1.12  2001/02/27 08:16:49  nali
   Make use of the new haplotype list class.

   Revision 1.11  2001/02/21 18:36:23  nali
   indnode.cpp

   Revision 1.10  2001/02/20 16:15:29  nali
   Use assert ( )  for debugging.

   Revision 1.9  2001/02/16 17:13:17  nali
   New functions: InputRandom, InputHusdonData, etc..
   New member function: make_haplist.

   Revision 1.8  2001/02/16 02:41:48  nali

   Remove the character representation of phenotypes.

   Revision 1.7  2001/02/14 23:40:40  nali

   Make some functions inline and no range checking.

   Revision 1.6  2001/02/14 02:00:21  nali

   Want get_phenotype() and get_phenotype_ch() to be public member functions
   eventually.

   Revision 1.5  2001/02/14 01:57:07  nali

   Add member function to print phenotype (mainly for debugging).

   Revision 1.4  2001/02/13 20:40:16  nali

   Two bugs fixed.
   1. Needs a copy constructor if a vector of objects is to be created.
      Since the vector is initialized by first creating a temporary object
      by the default constructor and copying it to every object in the
      vector. Especially when new operator is used in the constructor.
   2. For string stream, char str[10], where the length is necessary.
      Then ostrstream ostr(str, 10), where 10 is the length of char array.

   Revision 1.3  2001/02/12 05:28:42  nali
   Now it can read in and output phenotype data.
   No further processing yet.

   Revision 1.2  2001/02/10 01:13:52  nali
   New files added.

   Revision 1.1  2001/02/09 23:45:33  nali
   Reorganized

   Revision 1.1.1.1  2001/02/09 23:40:10  nali
   initial checkin
*/
// }}}

