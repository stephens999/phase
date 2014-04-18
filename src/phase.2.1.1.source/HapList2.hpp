/*
 * A class for list of haplotypes
 */

#ifndef HAPLIST_H
#define HAPLIST_H

#include "Haplotype.hpp"
#include "indnode.hpp"
#include "arrayQ.hpp"
#include "utility.hpp"
#include "Summary.hpp"
#include <map>
#include <iostream>
using namespace::std;

struct HapRecord
{double Freq; double PseudoCount; double Prob; vector<double> GroupFreq; vector<double> GroupFreqSq; double SqPseudoCount;

};

typedef map<Haplotype,HapRecord> ListType;
typedef pair<Haplotype,HapRecord> PairType;

class HapList {

public:
  ListType  haplist; // list of haplotypes and frequencies
  vector<ListType::iterator> PositiveHaps; // pointers to haps with positive frequency
  
public:
  HapList ();
  HapList (const HapList &);
  HapList (istream & istr);
  HapList (const HapList &, int firstlocus, int lastlocus); // extract
  //list from first up to (but not including) last
  HapList (const HapList &, const HapList &, double minfreq=0.0); // concatenate two haplists
  HapList (const HapList &, const HapList &, vector<CIndividual> & pop, double minfreq=0.0); // concatenate two haplists, but only include haps that can exist in pop
  
 
  const HapList & operator = (const HapList &);
  
  ~HapList();

  // the next two functions respectively
  // add a record whose Freq, PseudoCounts, and Prob are 
  // the product, or minimum, of 2 other records
  // (used in concatenating two haplists)
  void AddMultiple(const Haplotype &, const HapRecord &, const HapRecord &);
  void AddMin(const Haplotype &, const HapRecord &, const HapRecord &);
  void Add (const HapList & h1, double minfreq);
  void Add(const Haplotype &, const HapRecord &);  
  void Add(const Haplotype &, double freq = 1.0);  
  void Add(CIndividual &,int,bool,double freq = 1.0);
  void Add(CIndividual &,double freq = 1.0,bool = false);
  void Add(vector<CIndividual> & pop,int, bool = false); 

  // Add routines that return pointers to the added haplotype
  ListType::iterator Add (CIndividual & ind, int chr, double freq, bool & isnewhap, bool usebestguess = false);
  ListType::iterator Add (const Haplotype & h, double freq, bool & isnewhap);

  ListType::iterator Find(const Haplotype & h, const vector<int> & uselist);
  
  void AddAllPossible(CIndividual &,vector<int> unknown_list);
  void AddAllPossible(CIndividual &); 
  void AddAllPossible(vector<CIndividual>, int);

  void HardRemove(Haplotype &); // remove haplotype from list, without ref to frequency

  void Remove(Haplotype &, double freq = 1.0);
  void Remove(CIndividual &, double freq = 1.0);
  
  // soft remove only reduces freq, not actually remove from list
  void SoftRemove(const Haplotype &, double freq = 1.0);
  void SoftRemove(CIndividual &, double freq = 1.0);

  // output list of all haps with freq >= minfreq
  void Output(ostream & ostr, const vector<int> * coding, double min = 0.0, bool printheader = true); 
  void OutputProbs(ostream & ostr, const vector<int> * coding, double min = 0.0); 

  // output/input haps only, without coding info
  
  double FullDataPseudoLogLikelihood(char method, vector<ArrayQ *> & Q, int nchr, vector<double> & rho, double DPRIOR);

 
  void MakePairsIndex(vector < vector< pair< ListType::iterator,ListType::iterator> > > & index, vector <CIndividual> & pop, bool RemoveDuds, bool CheckMissing = true );
  void PrunePairsIndex(vector < vector< pair< ListType::iterator,ListType::iterator> > > & index, vector < vector< double > > & phaseprobs,vector <CIndividual> & pop,
		       double minthreshold);

  void ExtendPairsIndex(ListType::iterator & newrecord, vector < vector< pair< ListType::iterator,ListType::iterator> > > & index, vector < vector< double > > & phaseprobs,vector <CIndividual> & pop, bool CheckMissing = true);
  void ExtendPairsIndex(ListType::iterator & newrecord, vector< pair< ListType::iterator,ListType::iterator> >  & index, vector< double >  & phaseprobs, CIndividual & ind, bool CheckMissing = true);

  void MakePositiveHaps();

  // resolving phase function; returns likelihood
  double ResolvePhase(char method, int Niter, vector <CIndividual> & pop, vector<ArrayQ *> & Q , vector<double> & rho, bool randomise = false);
  
  vector<Summary> ProduceSummary(vector< vector< pair< ListType::iterator,ListType::iterator> > > & index, vector< vector<double> > & phaseprobs, int startlocus, int endlocus, vector<CIndividual> & pop, bool allowsplit = true);

 // functions to make a list of haplotypes in pop

  void RemoveAll(); // Clear whole list - deletes all members (inline)
 
  void ClearPseudoCounts(); // set all PseudoCounts to 0 (inline)
  void ClearProbs(); // set all Probs to 0 (inline)
  void ClearFreqs(); // set all Freqs to 0 (inline)
  void SetupGroupFreqs(int); // Set all freqs to 0 (inline)

  void NormaliseFreqs(); // divide Freqs by constant to sum to 1 (inline)
  void NormaliseSqPseudoCounts(double norm); // divide SqPseudoCount by norm
  void NormaliseGroupFreqs(); // normalise GroupFreq and GroupFreqSq 

    void RandomiseFreqs(); // randomly assign frequencies according to uniform
  // on simplex

  int Find(CIndividual &,int,bool=false);

  // these functions not actually to be used (see the functions below)
  void SetBestGuesses(vector<CIndividual> & pop, vector < vector< pair< ListType::iterator,ListType::iterator> > > & index, char use);
  void SetBestGuesses(vector<CIndividual> & pop, char use);

  // these functions are used instead
  pair <ListType::iterator, ListType::iterator> FindBestPair(vector< pair< ListType::iterator,ListType::iterator> > & index, vector< double > & phaseprobs, double & bestphaseprob, int startlocus = -1, int endlocus = -1 );
  void SetBestGuesses(vector<CIndividual> & pop, vector < vector< pair< ListType::iterator,ListType::iterator> > > & index, vector < vector< double > > & phaseprobs) ;
  
  // compute probability of each hap in list, in light of list
  // and place in the Haplotype record (.Prob) 
  double CalcProb(const Haplotype & h, char method, vector<ArrayQ *> & Q, int nchr, vector<double> & rho, double DPRIOR, const vector<int> & nmissing = vector<int>(0), bool ALLSNP = false, const vector<double> & diffprobvec = vector<double>(0), bool fuzzy = false);
  
  void ComputeProbs(char method, vector<ArrayQ *> & Q, int nchr, vector<double> & rho); 
  void ComputeEMProbs();
  void ComputeSDProbs(vector<ArrayQ *> & Q, int nchr);
  void ComputeFDLSProbs(const vector<ArrayQ *> & Q, int nchr, vector<double> & rho);
		  
  double EMProb(const Haplotype & h, double DPRIOR); 
  double SDProb(const Haplotype & h, vector<ArrayQ *> & Q, int nchr); 
                                // computes prob of haplotype h, 
                                // given haplotypes in list (+ freqs)
                                // under SD approximation

  double SNPSDProb(const Haplotype & h, const vector<double> & diffprobs); 
  double FDLSProb(const Haplotype & h, const vector<ArrayQ *> & Q, int nchr, vector<double> & vecRho, vector<double> & vecRhoDerivs, bool computederivs, bool usequad, vector<int> nmissing = vector<int>(0), const vector<double> & vecTheta = vector<double>(0),int Nforcorrection = 0, bool fuzzy = false );
 
  double ForwardsAlgorithm(const Haplotype & h, const vector<ArrayQ *> & Q, int nchr, vector<double> & vecRho, vector< vector<double> > & Alpha, vector<double> & AlphaSum, bool usequad, const vector<int> & nmissing , bool gofromrhs = false, const vector<double> & vecTheta = vector<double>(0), int Nforcorrection = 0);
  double FuzzyForwardsAlgorithm(const Haplotype & h, const vector<ArrayQ *> & Q, int nchr, vector<double> & vecRho, vector< vector<double> > & Alpha, vector<double> & AlphaSum, bool usequad, const vector<int> & nmissing , bool gofromrhs = false, const vector<double> & vecTheta = vector<double>(0), int Nforcorrection = 0);


  void BackwardsAlgorithm(const Haplotype & h, int nchr, vector<double> & vecRho, vector< vector<double> > & Alpha, vector<double> & AlphaSum, vector<int> & copiedtime, vector<int> & copiedallele, vector<int> & copiedhap, bool usequad = true);
  
  
  void ComputeHiddenStateProbs(vector<vector<double> > & CopyProb, const Haplotype & h, const vector<ArrayQ *> & Q, int nchr,  vector<double> & vecRho, bool usequad, const vector<int> & ismissing, const vector<double> & vecTheta, int Nforcorrection = 0);

  void ComputeVectorOfNaiveGibbsProbs(CIndividual &, vector<double> &, double &, double); // computes prob of individual being made up of each in list

  //output pair of haplotypes for individual
  void OutputPair(ostream & ostr, pair< ListType::iterator,ListType::iterator > hpair,  const vector<int> * coding);
  Haplotype get_haplotype(int listpos); //inline

  int get_listlength(); // inline
  int get_positivelistlength(); // inline
  
  int get_nloci(); // inline
  int get_ngroups(); // inline
  void CopyPseudoCountsToFreqs(); // inline
  
};

inline int HapList::get_listlength(){ return haplist.size();}
inline int HapList::get_positivelistlength(){ return PositiveHaps.size();}

inline int HapList::get_nloci(){ return haplist.begin()->first.get_nloci();}
inline int HapList::get_ngroups(){ return haplist.begin()->second.GroupFreq.size();}

inline void HapList::RemoveAll(){ haplist.clear();}


inline void HapList::ClearPseudoCounts(){
  for(ListType::iterator h = haplist.begin(); h!=haplist.end(); h++){
    h->second.PseudoCount = 0;
    h->second.SqPseudoCount = 0;
  }
}

inline void HapList::ClearProbs(){
  for(ListType::iterator h = haplist.begin(); h!=haplist.end(); h++){
    h->second.Prob = 0;
  }
}

inline void HapList::ClearFreqs(){
  for(ListType::iterator h = haplist.begin(); h!=haplist.end(); h++){
    h->second.Freq = 0;
  }
}

inline void HapList::SetupGroupFreqs(int ngroups){
  for(ListType::iterator h = haplist.begin(); h!=haplist.end(); h++){
    h->second.GroupFreq = vector<double>(ngroups,0);
    h->second.GroupFreqSq = vector<double>(ngroups,0);
  }
}



inline void HapList::NormaliseFreqs(){
  double sum = 0;
  for(ListType::iterator h = haplist.begin(); h!=haplist.end(); h++){
    sum += h->second.Freq;
  }
  for(ListType::iterator h = haplist.begin(); h!=haplist.end(); h++){
    h->second.Freq /= sum;
  }  
}

inline void HapList::NormaliseSqPseudoCounts(double norm){
  for(ListType::iterator h = haplist.begin(); h!=haplist.end(); h++){
    h->second.SqPseudoCount /= norm;
  }
}

inline void HapList::NormaliseGroupFreqs(){
  
  for(int g=0; g< haplist.begin()->second.GroupFreq.size(); g++){
    double norm = 0;
    for(ListType::iterator h = haplist.begin(); h!=haplist.end(); h++)
      norm += h->second.GroupFreq[g];
    for(ListType::iterator h = haplist.begin(); h!=haplist.end(); h++){
      h->second.GroupFreq[g] /= norm;
      h->second.GroupFreqSq[g] /= norm;
    }
  }
    
}

inline void HapList::RandomiseFreqs(){
  for(ListType::iterator h = haplist.begin(); h!=haplist.end(); h++){
    h->second.Freq=ranf();
  }
  NormaliseFreqs();
}

inline void HapList::CopyPseudoCountsToFreqs(){
  for(ListType::iterator h1 = haplist.begin(); h1 !=haplist.end(); h1++){
    (*h1).second.Freq = (*h1).second.PseudoCount;
  }
}

inline Haplotype HapList::get_haplotype(int listpos){
  int pos = 0;
  ListType::const_iterator h = haplist.begin();
  while(pos<listpos){
    h++;
    pos++;
  }
  return (*h).first;
}

#endif
