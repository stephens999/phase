/*
 * A class for list of pairs of haplotypes
 * (used for summarising output particularly)
 */

#ifndef HAPPAIRLIST_H
#define HAPPAIRLIST_H

#include "Haplotype.hpp"
#include "indnode.hpp"
#include "arrayQ.hpp"
#include "utility.hpp"
#include "HapList2.hpp"
#include <map>
using namespace::std;

typedef pair<Haplotype,Haplotype> HPairType;
typedef map<HPairType, double> PairListType ;

class HapPairList {

private:
  PairListType hpairlist;  // list of pairs of haplotypes and probabilities
  const HapPairList & operator = (const HapPairList &);  // declared as private so can't use for the moment

public:
  HapPairList ();
  HapPairList (const HapPairList &);
  HapPairList (const HapPairList &, int firstlocus, int lastlocus); // extract
  //list from first up to (but not including) last
  
  HapPairList (vector< pair< ListType::iterator,ListType::iterator> > & index, vector<double> & phaseprobs, int startlocus=-1, int endlocus=-1);

  ~HapPairList();

  int get_nloci(); //inline

  //add record to list
  void Add (const HPairType & hpair, const double prob);
 
  // find pair with highest prob in list
  HPairType BestPair(double & bestprob);
  HPairType BestPair(); // inline
  double BestKLdivergence(const vector<int> & nmissing);
  double KLsplitdivergence(int splitlocus, const vector <int> & nmissing);
  double KLdivergence(const HPairType & guesspair, vector<double> & flipprob, vector< vector<double> > & errprob, vector< vector< vector<double> > > & alleleprob, const vector<int> & nmissing);
  void ComputeFlipProbErrProbAlleleProb(const HPairType & guesspair, vector<double> & flipprob, vector< vector<double> > & errprob, vector< vector< vector<double> > > & alleleprob, const vector<int> & nmissing);
  
  // return summary of HapPairList
  Summary Summarise(const vector<int> & nmissing, bool allowsplit = true);
  Summary Summarise(int startlocus, int endlocus, const vector<int> & nmissing, bool allowsplit = true);

};

inline HPairType HapPairList::BestPair(){ double bp = 0; return BestPair(bp);}

//find a vector saying whether or not to flip each site to get from truepair to matchpair, 
//between startlocus and endlocus (uses the minimum number of flips)
void ComputeBestFlipBestErrBestAllele(const HPairType  & truepair, const HPairType & guesspair, vector<int> & bestflip, vector< vector<int> > & besterr, vector< vector<int> > & newallele, const vector<int> & nmissing);

// compute the probability of getting matchpair, from truepair, when the flip probs are flipprob
double MatchProb ( const HPairType & truepair, const HPairType & guesspair, vector<double> & flipprob, vector< vector<double> > & errprob, vector< vector < vector<double> > > & alleleprob, const vector<int> & nmissing);
bool IsHeterozygous(const HPairType & hpair);


inline int HapPairList::get_nloci(){ return hpairlist.begin()->first.first.get_nloci(); }


#endif
