/*
 * A class for vectors of individuals
 * $Id: classpop.hpp,v 1.30 2003/06/14 00:24:05 stephens Exp $
 */

#ifndef CLASS_POP_H
#define CLASS_POP_H

#include "indnode.hpp"
#include "arrayFF.hpp"
#include "arrayCC.hpp"
#include "arrayDiffCount.hpp"
#include "arrayDiffProb.hpp"
#include "arrayDiploidDiffProb.hpp"
#include "HapList2.hpp"
#include "Summary.hpp"

#include "tnt/tnt.h"
#include "tnt/cmat.h"

#include <map>
//#ifdef DMALLOC
//#include <dmalloc.h>
//#endif

using namespace::std;

class ClassPop {
  

  int Nloci;
  int Nind;
  int Nchild;

  vector<CIndividual> pop;    // Vector of individuals
  
  vector<int> childindex; // childindex[i] says which child is child of individual i

  vector<int> casecontrol; // vector of whether each ind is case or control
  vector<int> groupsize; // number of chromosomes in each group

  vector< vector<int> > buddy;

  string loci_type;           // Loci types
  
  vector<double> position;    // Map position

  vector<int> order; // order in which inds are considered

  vector<double> vecRho;
  vector<double> vecRhoDeriv;
  double RhoMean;
  vector<double> RhoMult; // vecRho[i] = RhoMean * RhoMult[i]

  vector<double> right; // parameters for simple hotspot model
  vector<double> lambda;
  vector<double> left;

  int RecomModel; // 0 = general model, 1 = simple hotspot, 2 = const recom
  int NindToUseForRho; // max number of inds to use in estimating Rho

  double CurrentLogProb;

  vector<ArrayQ *> Qptr;
  static vector<ArrayQ> Qvalues;
  static ArrayQ QSNP;

  vector<double> vecdiffprob; // prob that 2 haps differ by 0,1,2,... under SD approx (SNPs only)

  vector<double> vecTheta;

        /* For SNP sites, coding[0][r] are the integer value
       of the allele type (recoded as 0). Same for 1.
       For MS sites, coding[0][r] is the offset, 
       and coding[1][r] is range of repeats.  */
    vector<int> coding [2];  
  
    double BestLogProb; // best match so far

    int NSNP;
    int ALLSNP;
    int TREATSNPSASMS;
    vector<int> SNPlist;
    vector<int> nonSNPlist;


  ArrayCC CC; // 
  ArrayDiffCount DiffCount; //DiffCount stores the number of differences between each pair of haplotypes at SNP positions 

  // ArrayCC knownCC;               //as above, but based only on positions with known phase in id
  // ArrayDiffCount knownDiffCount; //

  // create lists of the haplotypes

  HapList haplist;

  
    // {{{ Private functions

  // functions to make a list of haplotypes in pop
  
    // compute pseudo-likelihood (=product of conditionals)
    double monitor_prob(const ArrayDiffProb &);

    // IO related functions
    void input_hudson_data    ( std::istream & );
    void input_random  ( std::istream &, int = 0 );

  void InitialiseSimpleHotspot(int nhotspot,  bool fixedpos = false);

    // Draw a random allele from the population at a locus
    int  draw_random_allele ( int ) const;
    // Initialization, called by input functions
    
    void normalize (int format = 0 );

    double calc_heterozygosity ( int locus ) const;
    void calc_theta ( );
  

    // Computations
     
  void GibbsUpdate ( int n1, double dirprior);

    void update_NR ( int n1, ArrayFF &, const ArrayDiffProb &, const ArrayDiploidDiffProb &, int fAncUpdate = 0, bool fNaiveGibbs=false);

   double update_R (vector<double> & rho, double theta, double delta, double temperature, int switchallowed);
  
  // simulated annealing version, updating just one strand at a time
  // (need to improve fast Forwards and fast Backwards algorithms
    double update_R_SA (vector<double> & rho, double theta,  double delta, double temperature, int switchallowed);

   //  void UpdateMissing_NR ( int ,
// 			    const ArrayDiffProb &);


    void ImputeMissingPositions(int n0, 
				int c0, 
				int n1, 
				int c1, 
				int t1, 
				int nchr, 
				const vector<int> & uselist);

    void diff_calc_phase_prob(int, 
			 std::vector<int>, 
			 std::vector<double>::iterator &,
			 const ArrayDiffProb &);

    void calc_phase_prob(int, 
			 std::vector<int>, 
			 std::vector<double>::iterator &, 
			 const ArrayDiffProb &);

  double DiploidForwardsAlg (std::vector<std::vector<double> > &,
                              int,
                              std::vector<double> &);

  void DiploidBackwardsAlg ( int n, double theta, double delta, vector<vector<double> > &, vector<double> & rho, double temperature, vector <vector<int> > & CopiedInd, vector <vector<int> > & CopiedChr, vector <vector<int> > & CopiedTime);


    double FastForwardsAlg ( int n, int c, double theta, double delta, vector<double> * FF, vector<double> & rho);
    void FastBackwardsAlg ( int n, int c, double theta, double delta, vector<double> * FF, vector<double> & rho, double temperature, vector <int> & CopiedInd, vector<int> & CopiedChr, vector<int> & CopiedTime);

  void ConditionalDiploidSim(int ind, double theta, double delta, double temperature,
vector<int> * CopiedInd, vector<int> * CopiedChr, vector<int> * CopiedTime);

    int impute_allele (int r, int nchr,
                       int t0, int from0) const;

    void update_phase_R  ( vector<double> & rho, double theta, double delta);    
    void update_phase_NR ( int, ArrayFF &, const ArrayDiffProb &, const ArrayDiploidDiffProb & );
   
    int update_phase_NR_fastestforsmallr ( int,  const ArrayDiffProb &, int );
  
    void UpdateCounts();
  void TransferCounts(vector<Summary> & sum);

    // }}}
 public:
  ClassPop ();  
  ClassPop ( const ClassPop & ); // copy constructor
  ClassPop ( const ClassPop &, int firstlocus, int lastlocus ); //extract sub class
  ClassPop ( const ClassPop &,  const ClassPop &, double minfreq = 0.0); //concatenate 2 classes
  ClassPop ( const ClassPop & CP1, const ClassPop & CP2, map<string,int> & cmds, int Niter, double theta, vector<double> & vecDelta, map<string, double> & d_cmds, double minfreq = 0.0);
    
    ~ClassPop ();
  const ClassPop & operator= ( const ClassPop & );
  void AugmentHapList(const ClassPop & c, double minfreq);

  double GetBestLogProb() const;
  void SaveCurrentState(); // save current state of all inds
  void SaveCurrentState(int); // save current state of ind i
  void RestoreSavedState(int); // restore current state 
  void RestoreSavedState(); // restore current state of all inds

  void Computevecdiffprob(); // compute vector of probs for fast SSD method

  // for making and outputting summary of best guesses
  void MakeHapList(bool);
  void MakeAllPossHapList();
 
  void ClearHapList(); // inline
  void OutputList(std::ostream &);
  void OutputHaplistSummary(std::ostream &);
  void OutputPositions(ostream &);
 
  void OutputPhaseProbs(std::ostream & ostr, bool=false);

  void OutputRho(ostream & ostr);
  void OutputHotspotParams(ostream & ostr);

  void OutputCurrentLogProb(ostream & ostr);

  void ResetCounts(); // reset counts in each individual
  
  // Access functions for private members
  // mostly coded inline below
  double get_position(int r) const;
    int get_allele (int ind, int chr, int locus ) const;  
    int get_phase (int ind, int locus ) const;  
    int get_haplotype( int , int, int) const;
    int get_nloci() const;
    int get_nind() const;
    std::vector<int> get_coding ( int ) const;
    char get_locus_type(int) const; // return type of rth locus
    std::string get_loci_type() const;
  double get_first_position() const;
  double get_last_position() const;
  double get_physical_length() const;

  // Manipulation functions; 
  // these are public because they are 
  // used by the class Summary in postprocessing.
  // Coded inline below
    void set_allele( int ind, int chr, int locus, int allele);
    void set_phase (int ind, int locus, int phase );
    void flip(int); 
    void flip_phase(int,int);
  
  // Recode to have same coding as another ClassPop
  // (used in postprocessing)
    void renormalize ( ClassPop & );

  // Read in data file     
    void read_data  ( std::istream &,
                      std::map<string, int> & cmds , std::map<string, double> & d_cmds, std::map<string,string> & filenames);
  
    // output routines
    void output_all_haps ( std::ostream & , bool=true, bool=true, bool=false, bool=true, double=0.5, double=0.0) const;
    void output_all_phases ( std::ostream & , bool=false) const;    
    void output ( std::ostream & , bool= false, bool= false) const;
  
    void print_id(std::ostream &, int) const;
    void print_allele (std::ostream &, int, int, int, bool=false) const;
     
    // Recode alleles and impute missing data
    void initialize ( std::istream &, std::istream &, int knowninfo, int initmethod, double theta, vector<double> & vecDelta, int format =0);

  // set up vecRho, RhoMean=rho, RhoMult = 1
  //void InitialiseRho(double rho);

  // Set the inds in pop to their best guesses according to Haplist
  // use = 'P' says to use the HapList probs; 'F' (default) to use the Freqs
  //void SetHaplistBestGuesses(char use = 'F');

  // Computations

  // EM-type algorithm
  double ListResolvePhase(char method, int Niter, vector<double> & vecDelta, double rho, bool randomise = false, bool initialise = true);
  
  // updated PHASE algorithm
  double MCMCListResolvePhase(map<string, int> & cmds, int Niter, int Nthin, int Nburn, vector<double> & vecDelta, map<string,double> & d_cmds, string filename, bool initialise = true, bool collectdata=false);

 // naive Gibbs algorithm
    void GibbsResolvePhase( int Niter, double dirprior);
  
  // original PHASE algorithm
    double resolve_phase_NR  (
                         int Nburn,
                         int Nthin,
                         int Niter,
                         std::ostream &,
                         //std::ostream &,
			 std::ostream &,
			 vector<double> &,
			 int verbose,
			 int fAncUpdate, int fNaiveGibbs);

    void resolve_phase_R  ( double theta, double totalrho,
                         int Nburn,
                         int Nthin,
                         int Niter,
                         std::ostream &,
                         //std::ostream &,
			 std::ostream &,
			 vector<double> &,
			 int verbose);

  void output_all_correct_probs( ostream & ostr); // inline
  void OutputHapList( ostream & ostr, double probcutoff = 0.1, bool printheader=true);

 
  double logProb(char method, vector<double> & rho, double Dprior, int group = -1); // call routine below, with casecontrols as groupvec
  double logProb(char method, vector<double> & rho, vector<int> & groupvec, double Dprior, int group = -1); // compute the log prob of current config, at specified value for rho 
  // (rho is ignored if method is not R)
  
  // these separate routines also allow computation of derivatives wrt rho, for updating rho
  double logFDLSProb(vector<double> & rho, vector<double> & rhoderiv, bool computederiv, int group = -1); // call routine below, with casecontrols as groupvec
  double logFDLSProb(vector<double> & rho, vector<double> & rhoderiv, bool computederiv, vector<int> & groupvec, int group = -1); // compute the log prob of current config, at specified value for rho (also compute derivatives and store in rhoderiv, if bool is true)

  void ComputeRho();
  void ComputeRho(vector<double> &,vector<double> &,vector<double> &);
 
  void ComputeRhoDerivAndCurrentLogProb();
  void TestComputeRhoDerivNaively();

  double EffectiveLength(vector<double> & , vector<double> & ,vector<double> &);
  bool AcceptOrRejectSimpleHotspot(vector<double> & newleft , vector<double> & newright, vector<double> & newlambda,bool fixedpos, double MAXLAMBDA, double Hotspotsd, double MinHotspotSize, double HOTSPOTRATE, double meanRhoMean, double sdRhoMean);
  
  bool updateRhoSimpleHotspot(bool fixedpos, map<string,double> & d_cmds);
  double logpriorprobRhoMean(double oldval, double newval, double mu, double sigma) const; // inline

  bool updateRhoMultRandomWalk(double sigma = 0.1);
  bool updateRhoMultLangevin(double sigma = 0.1);
  bool updateRhoMeanRandomWalk(double sigma, map<string,double> & d_cmds);
  bool updateRhoMeanLangevin(double sigma, map<string,double> & d_cmds);
  
  void RandomiseOrder();
  void MHUpdateOrder();


  void UpdateRho(double sigmaRhoMean, double sigmaRhoMult, int & , int &, map<string,double> & d_cmds);

  double InferRho(int Niter, double & sigmaRhoMean, double & sigmaRhoMult, int verbose, map<string,double> & d_cmds );
  
  void FastHapMapResolve(int Niter, int Nburn);
  void FastHapMapUpdate(int ind, bool burnin);
  double PriorProbFromCopyProb(int allele, int locus, vector< vector<double> > & CopyProb);

  double BuddyHapListMCMCResolvePhaseRemove(map<string, int> & cmds, int Niter, int Nthin, int Nburn, map<string,double> & d_cmds, string filename, bool collectdata);
  double HapListMCMCResolvePhaseRemove(map<string, int> & cmds, int Niter, int Nthin, int Nburn, map<string,double> & d_cmds, string filename, bool collectdata); 
  double FuzzyHapListMCMCResolvePhaseRemove(map<string, int> & cmds, int Niter, int Nthin, int Nburn, map<string,double> & d_cmds, string filename, bool collectdata);
  int GetGroupCount(const Haplotype & h, int g, int omit);
  void HapListImputeMissing(int ind);
};

double PrObserveGivenTruth(int observed0, int observed1, int true0, int true1, double delta);

double PrObserve(int observed0,int observed1,int imputedallele,int fromallele, int nchr, double theta, double delta, int time, vector<ArrayQ *> & Q, int locus);


// Inline member functions
// Get allelic types

inline double ClassPop::PriorProbFromCopyProb(int allele, int locus, vector<vector<double> > & CopyProb){
   return CopyProb[locus][0] * PrHitTarg(locus, 2*Nind - 2, 0, 0,allele,Qptr) + CopyProb[locus][1] * PrHitTarg(locus, 2*Nind - 2, 1, 0, allele, Qptr) + CopyProb[locus][2] * PrHitTarg(locus, 2*Nind - 2, 0, 1, allele, Qptr) + CopyProb[locus][3] * PrHitTarg(locus,2*Nind - 2, 1, 1, allele, Qptr);
}

inline void ClassPop::ClearHapList(){
  haplist.RemoveAll();
}

inline void ClassPop::OutputPositions(ostream & ostr){
  for(int r=0;r<Nloci;r++)
     ostr << position[r] << ' ';
  ostr << endl;
} 

inline void ClassPop::output_all_correct_probs ( ostream & ostr){
  for(int ind = 0; ind<Nind; ind++){
    pop[ind].output_overall_correct_prob(ostr);
    ostr << endl;
  }
}

inline double ClassPop::get_first_position ( ) const	{ 
  return position[0];
}

inline double ClassPop::get_last_position ( ) const	{ 
  return position[get_nloci()-1];
}

inline double ClassPop::get_position ( int r) const	{ 
  return position[r];
}


inline double ClassPop::logpriorprobRhoMean(double oldval, double newval, double mu, double sigma) const {
  double ov = log(oldval);
  double nv = log(newval);
  return (0.5/(sigma*sigma)) * ((ov - mu) * (ov -mu) - (nv - mu) * (nv - mu));
}

inline double ClassPop::get_physical_length ( ) const	{ 
  return position[get_nloci()-1] - position[0];
}


inline int ClassPop::get_allele (int ind, int chr, int locus ) const	{ 
    return pop[ind].get_allele(chr,locus); 
}

inline double ClassPop::GetBestLogProb () const	{ 
    return BestLogProb; 
}

inline int ClassPop::get_haplotype (int ind, int chr, int locus ) const	{ 
    return pop[ind].get_haplotype(chr,locus); 
}

inline int ClassPop::get_phase (int ind, int locus ) const	{ 
    return pop[ind].get_phase(locus); 
}

inline int ClassPop::get_nloci() const {
    return loci_type.size();
}

inline int ClassPop::get_nind() const {
    return pop.size();
}

inline char ClassPop::get_locus_type(int r) const {
    return loci_type[r];
}

inline std::string ClassPop::get_loci_type() const {
    return loci_type;
}

inline std::vector<int> ClassPop::get_coding(int c) const {
    return coding[c];
}

inline void ClassPop::set_allele (int ind, int chr, int locus, int allele )
{ 
    pop[ind].set_allele(chr,locus,allele); 
}

inline void ClassPop::set_phase (int ind, int locus, int phase )
{ 
    pop[ind].set_phase(locus,phase); 
}

inline void ClassPop::flip(int ind){
    pop[ind].flip_phase();
}

inline void ClassPop::flip_phase(int ind, int pos){
    pop[ind].flip_phase(pos);
}

inline void ClassPop::OutputList(std::ostream & ostr)
{ 
  haplist.Output(ostr,coding);
}

#endif // CLASS_POP_H

// {{{ Log
// 
// $Log: classpop.hpp,v $
// Revision 1.30  2003/06/14 00:24:05  stephens
// Adding files, and committing a lot of changes
// that have resulted in version 2.0 of PHASE
//
// Revision 1.29  2002/02/27 18:56:44  stephens
// Commiting the current source, which is essentially that released as
// PHASE 1.0. The number of quadrature points is reduced to 2. Some code has
// been added to cope with recombination, but this is not included in the release
// and is still under development.
//
// Revision 1.28  2001/10/21 18:33:29  stephens
// fixed bug in GibbsUpdate
//
// Revision 1.26  2001/10/16 20:44:16  stephens
// Added way to input position of each marker, by preceding line
// of SSMMMSS with a line beginning with P followed by position of
// each marker.
//
// Revision 1.25  2001/10/16 18:12:05  stephens
// More minor bugs fixed. (including -m bug)
//
// Revision 1.23  2001/10/06 00:35:40  stephens
// Impute missing positions added- changed update_allele to update_haplotype.
//
// Revision 1.22  2001/10/02 18:02:59  stephens
// just some small additions. haven't committed these for a long time
// because cvs was not working propersly. (now know why: they remounted
// you from /user3 to /home)
//
// Revision 1.21  2001/06/23 15:57:00  stephens
// started to check missing data code. Corrected some bugs. Seems to give sensible results in a couple of small examples.
//
// Revision 1.20  2001/06/19 17:02:25  stephens
// Changes to computation of arrayFF to make more efficient
// Added facility to store "original phenotype" in indnode,
// in preparation for allowing genotyping error.
//
// Revision 1.19  2001/05/31 16:26:14  stephens
// Added DiploidDiffProb look-up table class (to be used to make
// computation of arrayFF for SNPs more efficient)
//
// Revision 1.18  2001/05/30 06:02:18  stephens
// Updated to be considerably more efficient, via introduction of
// various new methods, including introduction of ArrayDiffProb
// and ArrayDiffCount to improve computation al efficiency for SNP
// data. Speedtest.inp now runs in about 7secs.
// Also corrected several bugs. Output now looks more promising
// and major bugs appear to have been eliminated. Convergence of chain
// can now be monitored more easily by the output in temp.monitor, which
// gives the pseudo-likelihood every Nthin repetitions.
//
// Revision 1.17  2001/05/23 02:16:55  stephens
// Corrected some bugs in computation of CC. Reduced number of
// loci we update by one when the total number of unknowns is <=5.
// Confirmed program was giving vaguely sensible resutls, and is
// competitive on the speed test example (36secs vs 31 for old version)
//
// Revision 1.16  2001/05/21 20:17:29  nali
// No real changes
//
// Revision 1.15  2001/05/18 21:47:42  stephens
// added facility to update just 5 at a time. Set it to do this
// so we can test for speed with the current version.
//
// Revision 1.14  2001/05/08 16:58:23  stephens
// Added class arrayCC, which computes quantities similar
// to arrayFF, but only for sites whose phase is known.
// This is then used to improve the efficiency of computation
// for arrayFF by re-using the calculations for known sites.
// Checked that the output was unchanged on an example.
// Might slightly improve efficiency of computation of CC if
// we do the loop through known sites more efficiently.
//
// Revision 1.13  2001/04/27 00:20:48  stephens
// First version of postprocess done - usage is
// postprocess -n file.in file.sum dummyfilename
//
// Revision 1.12  2001/04/26 18:29:51  stephens
// Fixed bug in computation of ArrayFF (total_sum now initialized
// to 0 each time FF is computed)
//
// Revision 1.11  2001/04/26 17:09:17  stephens
// Added function ClassPop::print_allele
// and made a few minor changes
//
// Revision 1.10  2001/04/24 22:00:08  stephens
// Added flag in ClassPop::output_hap and output_phase
// to indicate whether to output all positions, or only
// unknown positions.
//
// Revision 1.9  2001/04/24 20:25:08  stephens
// Minor changes to renormalize function
//
// Revision 1.8  2001/04/24 19:54:43  stephens
// Added access function get_coding0 to ClassPop
//
// Revision 1.7  2001/04/24 19:31:31  nali
// Move data input out of the constructor. Member functions read_data 
// and initialize have to be called explicitly. Put everything related
// to hudson data set into a single function so that it might go away 
// one day.
//
// Revision 1.6  2001/04/23 18:55:16  stephens
// Added member function ClassPop::renormalize(newcoding0)
// to renormalize a population; for use when postprocessing
// results to make sure all samples have been normalized the
// same way.
//
// Revision 1.5  2001/04/21 01:20:32  stephens
// added a couple of public functions to ClassPop: get_allele
// and get_nloci
//
// Revision 1.4  2001/04/20 00:32:45  nali
// Put reading loci types into the contructor of ClassPop
//
// Revision 1.3  2001/04/19 19:45:28  nali
// Add update missing data, use class ArrayFF.
//
// Revision 1.2  2001/04/17 22:53:27  stephens
// Moved calculation of FF in non-recom case to its own function
// CalcFF_NR
//
// Revision 1.1  2001/04/17 22:08:44  nali
//
// Major revisement in overall structure. Created new class ClassPop and
// almost all global functions now became member functions of ClassPop, most
// of them private.
//
// "mult" removed in update_phase_NR. No other changes in terms of algorithm.
// Haven't check the results yet.
//
// proc_args() is moved to utility.cpp, which also defines a couple of other
// global functions.
//
// }}}

