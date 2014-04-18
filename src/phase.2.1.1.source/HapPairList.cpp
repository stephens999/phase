#include "HapPairList.hpp"
#include "indnode.hpp"
#include "constants.hpp"

#include <iomanip>

using namespace std;

// constructor
HapPairList::HapPairList () :
  hpairlist ()
{
}

// copy constructor
HapPairList::HapPairList (const HapPairList & h2):
  hpairlist ( h2.hpairlist)
{
}

HapPairList::HapPairList (const HapPairList & h2, int firstlocus, int secondlocus):
  hpairlist ()
{
  for(PairListType::const_iterator i = h2.hpairlist.begin(); i != h2.hpairlist.end(); i++){
    Haplotype hap1((i->first).first,firstlocus,secondlocus);
    Haplotype hap2((i->first).second,firstlocus,secondlocus);
    pair<Haplotype,Haplotype> hpair(hap1,hap2);
    Add( hpair, i->second ); 
  }
}

HapPairList::HapPairList (vector< pair< ListType::iterator,ListType::iterator> > & index, vector<double> & phaseprobs, int firstlocus , int endlocus ):
  hpairlist ()
{
  // default behaviour - use whole haplotypes
  if(firstlocus<0){
    firstlocus =0;
    endlocus = phaseprobs.size();
  }

  vector<double>::iterator phasepointer = phaseprobs.begin();
  for(vector< pair< ListType::iterator,ListType::iterator> >:: iterator hpair =  index.begin(); hpair !=index.end(); hpair ++){
    Haplotype hap1(hpair->first->first,firstlocus,endlocus);
    Haplotype hap2(hpair->second->first,firstlocus,endlocus);
    HPairType temppair(hap1,hap2);
    Add(temppair,*phasepointer);
    phasepointer ++;
  }
    


}


HapPairList::~HapPairList()
{
}

const HapPairList & HapPairList::operator=(const HapPairList & rhs)
{
  if(this != &rhs){
    this->hpairlist = PairListType(rhs.hpairlist);
  }
  return *this;
}

void HapPairList::Add (const HPairType & hpair, const double prob)
{
  // impose hap1<=hap2
  Haplotype hap1 ((hpair.first));
  Haplotype hap2 ((hpair.second));
  if(hap2<hap1){
    Haplotype temp(hap1);
    hap1 = hap2;
    hap2 = temp;
  }
  PairListType::iterator i = hpairlist.lower_bound(hpair);
  if((i!=hpairlist.end()) && (*i).first == hpair){
    (*i).second += prob;
  }
  else
  {
    hpairlist.insert(i,pair< HPairType, double >( hpair, prob));
  }
}

HPairType HapPairList::BestPair(double & bestprob)
{
  bestprob = 0;
  HPairType bestpair = hpairlist.begin()->first;
  double sum = 0;
  for(PairListType::iterator i = hpairlist.begin(); i!=hpairlist.end(); i++){
    sum += (i->second);
    if((i->second) > bestprob){
      bestprob = i->second;
      bestpair = i->first;
    }
  }
  bestprob /=sum; // return prob of best pair in bestprob;

  return bestpair;

}
      

// return KL divergence (actually integral p log q)
// from the distribution p
// on haplotype pairs defined by the list
// to q = (bestpair, flipprob)
double HapPairList::BestKLdivergence(const vector<int> & nmissing)
{
  
  int nloci = nmissing.size();
  pair<Haplotype,Haplotype> bestpair = BestPair();
  vector< vector<double>  > errprob ( nloci, vector< double >(2,0) );
  vector< vector< vector<double> > > alleleprob( nloci, vector< vector<double> >(2, vector<double>(KMAX,0) ) );
  vector<double> flipprob(nloci,0);
  ComputeFlipProbErrProbAlleleProb(bestpair, flipprob, errprob, alleleprob, nmissing);
  return KLdivergence(bestpair,flipprob, errprob, alleleprob, nmissing);
}


// return KL divergence (actually integral p log q)
// from the distribution p
// on haplotype pairs defined by the list
// to q = (guesspair, flipprob, errprob, alleleprob)
double HapPairList::KLdivergence(const HPairType & guesspair, vector<double> & flipprob, vector< vector<double> > & errprob, vector< vector< vector<double> > > & alleleprob, const vector<int> & nmissing)
{
  double sum = 0;
  for(PairListType::const_iterator listpointer =  hpairlist.begin(); listpointer !=hpairlist.end(); listpointer ++){
    sum += (listpointer->second) * log(MatchProb(guesspair, listpointer->first, flipprob, errprob, alleleprob, nmissing));
  }
  return sum;
}


// compute the vector of flip probabilities for a particular guesspair
// based on the distribution defined by the list
void HapPairList::ComputeFlipProbErrProbAlleleProb(const HPairType & guesspair, vector<double> & flipprob, vector< vector<double> > & errprob, vector< vector< vector<double> > > & alleleprob, const vector<int> & nmissing) {
  
  vector<int> flipvec(get_nloci()); // stores whether to flip locus i in 
  //best way to get from guesspair to truepair
  vector< vector<int> > errvec( get_nloci(), vector<int>(2,0) ); // stores whether locus i, chr j, is in error 
  vector< vector<int> > allelevec( get_nloci(), vector<int>(2,0) );
  // stores allele that replaces the erroneous allele at error positions
  
  flipprob = vector<double> (get_nloci(),0);
  errprob = vector< vector<double> > (get_nloci(), vector<double>(2,0) );
  alleleprob = vector< vector< vector<double> > > (get_nloci(), vector< vector<double> > ( 2, vector<double>(KMAX,0)));

  for(PairListType::const_iterator listpointer =  hpairlist.begin(); listpointer !=hpairlist.end(); listpointer ++){
    ComputeBestFlipBestErrBestAllele(listpointer->first, guesspair, flipvec, errvec, allelevec, nmissing);
    for(int locus = 0; locus < get_nloci(); locus++){
      flipprob[locus] += listpointer->second * flipvec[locus];
      for(int chr = 0; chr<2; chr++){
	errprob[locus][chr] += listpointer->second * errvec[locus][chr];
	for(int allele = 0; allele < KMAX; allele ++)
	  alleleprob[locus][chr][allelevec[locus][chr] ]+= listpointer->second * errvec[locus][chr];
      }
    }
  }

// normalise alleleprob
  for(int locus = 0; locus < get_nloci(); locus++){
    for(int chr = 0; chr<2; chr++){
      double sum = 0;
      for(int allele = 0; allele < KMAX; allele ++)
	sum += alleleprob[locus][chr][allele];
      if(sum)
	for(int allele = 0; allele < KMAX; allele ++)
	  alleleprob[locus][chr][allele] /= sum;
    }
  }

}

double HapPairList::KLsplitdivergence(int splitlocus, const vector <int> & nmissing)
{
  int nloci = get_nloci();
  HapPairList lhs(*this, 0, splitlocus);
  HapPairList rhs(*this, splitlocus, get_nloci());
  
  vector< vector<double>  > errprob1 ( splitlocus, vector< double >(2,0) );
  vector< vector< vector<double> > > alleleprob1( splitlocus, vector< vector<double> >(2, vector<double>(KMAX,0) ) );
  vector<double> flipprob1(splitlocus,0);
  vector< vector<double>  > errprob2 (nloci - splitlocus , vector< double >(2,0) );
  vector< vector< vector<double> > > alleleprob2( nloci - splitlocus, vector< vector<double> >(2, vector<double>(KMAX,0) ) );
  vector<double> flipprob2(nloci - splitlocus,0);
  
  HPairType guesspair1=lhs.BestPair();
  HPairType guesspair2=rhs.BestPair();
  

  vector< int > nmissing1 = vector<int> (nmissing.begin(), nmissing.begin()+splitlocus);
  vector< int > nmissing2 = vector<int> (nmissing.begin()+splitlocus, nmissing.end());


  lhs.ComputeFlipProbErrProbAlleleProb(guesspair1, flipprob1, errprob1, alleleprob1, nmissing1);
  rhs.ComputeFlipProbErrProbAlleleProb(guesspair2, flipprob2, errprob2, alleleprob2, nmissing2); 
    
  // note that if both sides of the split are heterozygous then there's a factor of 0.5 in q, which
  // turns into a log(0.5) below. If either side of the split is homozygous at all positions then this
  // factor doesn't enter into q
  return ((log(0.5)* IsHeterozygous(guesspair1) * IsHeterozygous(guesspair2)) + lhs.KLdivergence(guesspair1, flipprob1, errprob1, alleleprob1, nmissing1) + rhs.KLdivergence(guesspair2, flipprob2, errprob2, alleleprob2, nmissing2));

}

Summary HapPairList::Summarise(const vector<int> & nmissing, bool allowsplit )
{
  int bestsplitlocus = 0;
  if(allowsplit){
    double kl = BestKLdivergence(nmissing);
    // compute KL divergence for splitting at each possible locus
    double bestkl = kl;    
    for(int locus = 1; locus < get_nloci(); locus++){
      //double klsplit = KLsplitdivergence(index[ind], phaseprobs[ind],locus);
      double klsplit = KLsplitdivergence(locus, nmissing);
      if(klsplit > bestkl){
	bestkl = klsplit;
	bestsplitlocus = locus;
      }
    }
  }
  
  if((bestsplitlocus>0) && allowsplit){ // split the list into 2, and summarise each side separately
    Summary s1 = Summarise(0, bestsplitlocus, nmissing);
    Summary s2 = Summarise(bestsplitlocus, get_nloci(), nmissing);
    return Summary(s1,s2);
  }
  else{ // just summarise the list by its bestguess and flipprob
    HPairType bestguess = BestPair();
    
    vector<double> flipprob;
    vector< vector <double> > errprob;
    vector< vector < vector <double> > > alleleprob;

    ComputeFlipProbErrProbAlleleProb(bestguess, flipprob, errprob, alleleprob, nmissing);
    return Summary(bestguess,flipprob,errprob, alleleprob);
  }

}

Summary HapPairList::Summarise(int startlocus, int endlocus, const vector<int> & nmissing, bool allowsplit)
{
  HapPairList partiallist(*this, startlocus, endlocus);
  vector<int> partialnmissing (nmissing.begin()+startlocus, nmissing.begin()+endlocus);
  return partiallist.Summarise(partialnmissing, allowsplit);
}



//returns a vector saying whether or not to flip each site to get from 
// guesspair to truepair, 
//between startlocus and endlocus (uses the minimum number of flips)
void ComputeBestFlipBestErrBestAllele(const HPairType  & truepair, const HPairType & guesspair, vector<int> & bestflip, vector< vector<int> > & besterr, vector< vector<int> > & bestallele, const vector<int> & nmissing){
  int nloci = truepair.first.get_nloci();
  vector<int> f0 (nloci,0); // flips required to match first to first
  vector<int> f1 (nloci,0); // flips required to match first in truepair to second in guesspair
  vector< vector<int> > e0 (nloci, vector<int>(2,0));
  vector< vector<int> > e1 (nloci, vector<int>(2,0));
  vector< vector<int> > newallele0 (nloci, vector<int> (2,0));
  vector< vector<int> > newallele1 (nloci, vector<int> (2,0));
  
  int sum0 = 0; // sum of f0+e0
  int sum1 = 0; // sum of f1+e1
  
  for(int i = 0; i < nloci; i++){
    int t1 = truepair.first.get_allele(i);
    int t2 = truepair.second.get_allele(i);
    int g1 = guesspair.first.get_allele(i);
    int g2 = guesspair.second.get_allele(i);

    if(nmissing[i]==0){
      if( (t1 != g1 ) && (t2 != g2) && (t1 == g2) && (t2 == g1) ){
	f0[i] = 1;
	f1[i] = 0;
	sum0++;
      } else if( (t1 != g2 ) && (t2 != g1) && (t1 == g1) && (t2 == g2) ){
	f0[i] = 0;
	f1[i] = 1;
	sum1++;
      }
    } else {
      if(t1 != g1) { e0[i][0]++; sum0++; newallele0[i][0] = t1; }
      if(t2 != g2) { e0[i][1]++; sum0++; newallele0[i][1] = t2; }
      if(t1 != g2) { e1[i][1]++; sum1++; newallele1[i][1] = t1; }
      if(t2 != g1) { e1[i][0]++; sum1++; newallele1[i][0] = t2; }
    }
  }
  if(sum0<sum1){
    bestflip = f0;
    besterr = e0;
    bestallele = newallele0;
  }
  else{
    bestflip = f1;
    besterr = e1;
    bestallele = newallele1;
  }
}

 

// compute the probability of getting truepair, from guesspair, when the flip probs are flipprob
double MatchProb ( const HPairType & truepair, const HPairType & guesspair, vector<double> & flipprob, vector< vector<double> > & errprob, vector< vector < vector<double> > > & alleleprob, const vector<int> & nmissing){
  int nloci = truepair.first.get_nloci();
  double p0 = 1; // prob of getting match first->first and second->second
  double p1 = 1; // prob of getting match first->second and second->first
  for(int i = 0; i < nloci; i++){
    int t1 = truepair.first.get_allele(i);
    int t2 = truepair.second.get_allele(i);
    int g1 = guesspair.first.get_allele(i);
    int g2 = guesspair.second.get_allele(i);
    if(nmissing[i]==0){
      if( (t1 != g1) && (t2 != g2) && (t1 == g2) && (t2 == g1)){
	p0 *= flipprob[i]; 
	p1 *= 1-flipprob[i];
      } else if ( (t1 != g2) && (t2 != g1) && (t1 == g1) && (t2 == g2) ){
	p0 *= 1 - flipprob[i];
	p1 *= flipprob[i];
      } 
    } else {
      p0 *= ( (1-errprob[i][0]) * (g1 == t1) + errprob[i][0] * alleleprob[i][0][t1] ) * ( (1-errprob[i][1]) * (g2 == t2) + errprob[i][1] * alleleprob[i][1][t2] );
      p1 *=  ( (1-errprob[i][0]) * (g1 == t2) + errprob[i][0] * alleleprob[i][0][t2] ) * ( (1-errprob[i][1]) * (g2 == t1) + errprob[i][1] * alleleprob[i][1][t1] );
    }
  }
  if(IsHeterozygous(guesspair))
    return p0+p1;
  else
    return p0;
}

bool IsHeterozygous(const HPairType & hpair)
{
  Haplotype hap1(hpair.first);
  Haplotype hap2(hpair.second);
  return(hap1!=hap2);
}
