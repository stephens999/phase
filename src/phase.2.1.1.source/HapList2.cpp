#include "HapList2.hpp"
#include "HapPairList.hpp"
#include "indnode.hpp"
#include "constants.hpp"
#include "Summary.hpp"
#include "utility.hpp"
#include "pacb.lookup.hpp"

#include <iomanip>


using namespace std;
extern int NHL;
// constructor
HapList::HapList () :
  haplist (),
  PositiveHaps ()
{
  NHL++;
}

// copy constructor
HapList::HapList (const HapList & h2):
  haplist ( h2.haplist),
  PositiveHaps ( h2.PositiveHaps)
{
  NHL++;
}

// read a haplist in from a file
HapList::HapList (istream & istr):
  haplist(),
  PositiveHaps ()
{
  cerr<< "Not implemented yet" << endl;
  exit(1);

  int Nhap;
  istr >> Nhap;
  for(int i = 0; i < Nhap; i++){
    Haplotype h;
    //h.read(istr);
    Add(h,0);
  }
}

HapList::HapList (const HapList & h2, int firstlocus, int secondlocus):
  haplist (),
  PositiveHaps ()
{
  NHL++;
  for(ListType::const_iterator i = h2.haplist.begin(); i != h2.haplist.end(); i++){
    Add( Haplotype(i->first,firstlocus,secondlocus), i->second ); 
  }
  MakePositiveHaps();
}


HapList::HapList (const HapList & h1, const HapList & h2, double minfreq ):
  haplist (),
  PositiveHaps()
{
  NHL++;
  if(minfreq>=0){
    for(ListType::const_iterator i = h1.haplist.begin(); i != h1.haplist.end(); i++){
      if((i->second.Freq)>minfreq){
	for(ListType::const_iterator j = h2.haplist.begin(); j!= h2.haplist.end(); j++){
	  if((j->second.Freq)>minfreq){
	    AddMin( Haplotype(i->first, j->first), i->second, j->second);
	  }
	}
      }
    }
  } else {
    // find 50th largest freq in each list
    vector<double> freq1;
    for(ListType::const_iterator i = h1.haplist.begin(); i != h1.haplist.end(); i++){
      freq1.push_back(i->second.Freq);
    }
    sort(freq1.begin(),freq1.end());

    vector<double> freq2;
    for(ListType::const_iterator i = h2.haplist.begin(); i != h2.haplist.end(); i++){
      freq2.push_back(i->second.Freq);
    }
    sort(freq2.begin(),freq2.end());

    for(ListType::const_iterator i = h1.haplist.begin(); i != h1.haplist.end(); i++){
      if((freq1.size()<50) || ((i->second.Freq)>freq1[freq1.size()-50])){
	for(ListType::const_iterator j = h2.haplist.begin(); j!= h2.haplist.end(); j++){
	  if((freq2.size()<50) || ((j->second.Freq)>freq2[freq2.size()-50])){
	    AddMin( Haplotype(i->first, j->first), i->second, j->second);
	  }
	}
      }
    }
    
  }
  MakePositiveHaps();
}


HapList::HapList (const HapList & h1, const HapList & h2, vector<CIndividual> & pop, double minfreq):
  haplist (),
  PositiveHaps()
{
  NHL++;
  for(ListType::const_iterator i = h1.haplist.begin(); i != h1.haplist.end(); i++){
    if((i->second.Freq)>minfreq){
      for(ListType::const_iterator j = h2.haplist.begin(); j!= h2.haplist.end(); j++){
	if((j->second.Freq)>minfreq){
	  Haplotype h(i->first, j->first);
	  if(CanBeFoundAtAll(h,pop)) // if h can be found in pop
	    AddMin( h, i->second, j->second);
	}
      }
    }
  } 
  MakePositiveHaps();
}




HapList::~HapList()
{
  NHL--;
}

const HapList & HapList::operator=(const HapList & rhs)
{
  if(this != &rhs){
    this->haplist = ListType(rhs.haplist);
    this->PositiveHaps = vector<ListType::iterator>(rhs.PositiveHaps);
  }
  return *this;
}


void HapList::AddMin (const Haplotype & h, const HapRecord & hrecord1, const HapRecord & hrecord2)
{
  ListType::iterator i = haplist.lower_bound(h);
  if((i!=haplist.end()) && (*i).first == h){
    if(hrecord1.Freq<hrecord2.Freq){
      (*i).second.Freq += hrecord1.Freq;
      (*i).second.PseudoCount += hrecord1.PseudoCount;
      (*i).second.Prob += hrecord1.Prob;
    }
    else{
      (*i).second.Freq += hrecord2.Freq;
      (*i).second.PseudoCount += hrecord2.PseudoCount;
      (*i).second.Prob += hrecord2.Prob;
    }
  }
  else
  {
    HapRecord emptyrecord = {0,0,0,vector<double>(0)};
    haplist.insert(i,pair<Haplotype,HapRecord>(h,emptyrecord));
    AddMin(h,hrecord1,hrecord2);
  }
}

void HapList::AddMultiple (const Haplotype & h, const HapRecord & hrecord1, const HapRecord & hrecord2)
{
  ListType::iterator i = haplist.lower_bound(h);
  if((i!=haplist.end()) && (*i).first == h){
    (*i).second.Freq += hrecord1.Freq * hrecord2.Freq;
    (*i).second.PseudoCount += hrecord1.PseudoCount * hrecord2.PseudoCount;
    (*i).second.Prob += hrecord1.Prob * hrecord2.Prob;
  }
  else
  {
    HapRecord emptyrecord = {0,0,0,vector<double>(0)};
    haplist.insert(i,pair<Haplotype,HapRecord>(h,emptyrecord));
    AddMultiple(h,hrecord1,hrecord2);
  }
}


// add the contents of h1 to the haplist
void HapList::Add (const HapList & h1,double minfreq)
{
  for(ListType::const_iterator i = h1.haplist.begin(); i!=h1.haplist.end(); i++){
    if(i->second.Freq>minfreq)
      Add(i->first,i->second.Freq);
  }
}

void HapList::Add (const Haplotype & h, const HapRecord & hrecord)
{
  ListType::iterator i = haplist.lower_bound(h);
  if((i!=haplist.end()) && (*i).first == h){
    (*i).second.Freq += hrecord.Freq;
    (*i).second.PseudoCount += hrecord.PseudoCount;
    (*i).second.SqPseudoCount += hrecord.SqPseudoCount;   
    (*i).second.Prob += hrecord.Prob; 
    int ngroups = haplist.begin()->second.GroupFreq.size();
    (*i).second.GroupFreq = vector<double>(ngroups,0);
    (*i).second.GroupFreqSq = vector<double>(ngroups,0);
    
  }
  else  {
    HapRecord emptyrecord = {0,0,0,vector<double>(0)};
    haplist.insert(i,pair<Haplotype,HapRecord>(h,emptyrecord));	    
    Add(h,hrecord);
  }
}

// adds haplotype, returning an iterator pointing to the Haplotype in
// the list. If the haplotype is new in the list, isnewhap returns true
ListType::iterator HapList::Add (const Haplotype & h, double freq, bool & isnewhap)
{
  ListType::iterator i = haplist.lower_bound(h);
  if((i!=haplist.end()) && (*i).first == h){
    isnewhap = false;
    (*i).second.Freq+=freq;   
    return i;
  }
  else
  {
    isnewhap = true;
    int ngroups = 0;
    if(haplist.begin()!=haplist.end())
      ngroups = haplist.begin()->second.GroupFreq.size();
    HapRecord emptyrecord = {0,0,0,vector<double>(ngroups,0),vector<double>(ngroups,0),0};
    haplist.insert(i,pair<Haplotype,HapRecord>(h,emptyrecord));	  
    return (Add(h,freq,isnewhap));   
  }
}

void HapList::Add (const Haplotype & h, double freq)
{
  ListType::iterator i = haplist.lower_bound(h);
  if((i!=haplist.end()) && (*i).first == h){
    (*i).second.Freq+=freq;
  }
  else
  {
    HapRecord emptyrecord = {0,0,0,vector<double>(0)};
    haplist.insert(i,pair<Haplotype,HapRecord>(h,emptyrecord));	  
    Add(h,freq);
  }
}

ListType::iterator HapList::Add (CIndividual & ind, int chr, double freq, bool & isnewhap, bool usebestguess )
{
  Haplotype hap(ind.get_haplotype(chr));
  if(usebestguess)
    hap = ind.BestHaplotype(chr);
  return( Add(hap,freq, isnewhap));
}

void HapList::Add (CIndividual & ind, int chr, bool usebestguess, double freq  )
{
  Haplotype hap(ind.get_haplotype(chr));
  if(usebestguess)
    hap = ind.BestHaplotype(chr);
  Add(hap,freq);
}

void HapList::Add (CIndividual & ind,  double freq , bool usebestguess)
{
  Haplotype hap0(ind.get_haplotype(0));
  Haplotype hap1(ind.get_haplotype(1));
  
  if(usebestguess){
    hap0 = ind.BestHaplotype(0);
    hap1 = ind.BestHaplotype(1);
  }

  Add(hap0,freq);
  Add(hap1,freq);
}

void HapList::Add (vector<CIndividual> & pop, int Nind, bool usebestguess )
{
  for(int ind = 0; ind < Nind; ind++)
    Add(pop[ind], 1.0, usebestguess);
}

void HapList::AddAllPossible(CIndividual & ind,vector<int> unknown_list){

  if(unknown_list.size()==0){
    Add(ind);
  }
  else{
    int lastpos = unknown_list.back(); 
    unknown_list.pop_back(); // remove end from update list
    int knownallele;

    // for microsats make list based on current imputed missing
    // alleles (list may be augmented later, so not crucial
    // that these imputed alleles may not be so great)
    if(ind.get_locus_type(lastpos)=='M'){
      AddAllPossible(ind,unknown_list);
      ind.flip_phase(lastpos);
      AddAllPossible(ind,unknown_list);
    }

    else{

      switch(ind.n_missing(lastpos)){ // for each possible set of values at lastpos, add all possible
      case 0:
	AddAllPossible(ind,unknown_list);
	ind.flip_phase(lastpos);
	AddAllPossible(ind,unknown_list);
        break;
      case 1: // 3 possibilites: homozygous for knownallele, or 2 possible het phases	
	knownallele = ind.get_orig_nonmissing_allele(lastpos);
	ind.update_haplotype(0, lastpos, knownallele);
	ind.update_haplotype(1, lastpos, knownallele);
	AddAllPossible(ind,unknown_list);
	ind.update_haplotype(1, lastpos, 1-knownallele);
	AddAllPossible(ind,unknown_list);
	ind.flip_phase(lastpos);
	AddAllPossible(ind,unknown_list);
	break;
      case 2: // 4 possibilites (for SNPs)
	
	ind.update_haplotype(0, lastpos, 0);
	ind.update_haplotype(1, lastpos, 1);      
	AddAllPossible(ind,unknown_list);
	ind.flip_phase(lastpos);
	AddAllPossible(ind,unknown_list);
	ind.update_haplotype(0, lastpos, 1);
	ind.update_haplotype(1, lastpos, 1); 
	AddAllPossible(ind,unknown_list);
	ind.update_haplotype(0, lastpos, 0);
	ind.update_haplotype(1, lastpos, 0); 
	AddAllPossible(ind,unknown_list);
	break;
      }
    }
  }
}

void HapList::AddAllPossible(CIndividual & ind){
  AddAllPossible(ind, ind.get_unknown_pos());
}

void HapList::AddAllPossible(vector<CIndividual> pop, int Nind){
  for(int i = 0; i< Nind; i++){
    AddAllPossible(pop[i]);
  }
}


// find a match for h in the list, but insist on a match only
// at the positions in uselist
ListType::iterator HapList::Find(const Haplotype & h, const vector<int> & uselist)
{ 
  ListType::iterator listpointer;
  ListType::iterator firsttry = haplist.find(h); // try to match exactly,
  if(firsttry != haplist.end())                  // ignore uselist
    listpointer = firsttry;
  else{
    listpointer = haplist.begin();
    while((!h.Matches(listpointer -> first, uselist)) 
	  && (listpointer != haplist.end()) )
      listpointer++;
  }
  return listpointer;
}

// remove haplotypes entirely from list
void HapList::HardRemove(Haplotype & h)
{
  haplist.erase(h);
}


void HapList::Remove (Haplotype & h, double freq)
{
  haplist[h].Freq-=freq;
  if(haplist[h].Freq<=0)
    haplist.erase(h);
}

void HapList::Remove (CIndividual & ind, double freq )
{
  Haplotype hap0(ind.get_haplotype(0));
  Haplotype hap1(ind.get_haplotype(1));
  Remove(hap0,freq);
  Remove(hap1,freq);
}

// remove by subtracting freq, but not remove from list
void HapList::SoftRemove (const Haplotype & h, double freq)
{
  haplist[h].Freq-=freq;
}

void HapList::SoftRemove (CIndividual & ind, double freq )
{
  Haplotype hap0(ind.get_haplotype(0));
  Haplotype hap1(ind.get_haplotype(1));
  SoftRemove(hap0,freq);
  SoftRemove(hap1,freq);
}

void HapList::Output (ostream & ostr, const vector<int> * coding,  double min ,bool printheader )
{
  if(haplist.size()>0){
  int count=1;
  int plen = haplist.begin()->first.get_printedlen();
  int ngroups = get_ngroups();

  if(printheader){
    ostr.setf(ios::right);
    ostr << setw(10) << "index";
    for(int i = 10; i< plen; i++)
      ostr << ' ';
    ostr << "  haplotype";
    ostr << setw(12) << " E(freq)";
    ostr << setw(12) << " S.E";
    if(ngroups >1 ){
      for(int i = 0; i < ngroups ; i++){
	ostr << setw(9) << " E[Freq(" << i << ")]";
	ostr << setw(10) << " S.E.(" << i << ")";
      }
    }
    ostr << endl;
  }
 
  for(ListType::const_iterator h = haplist.begin(); h != haplist.end(); h++){
    if((*h).second.Freq >= min){

      ostr.setf(ios::right);
      ostr << setw(10) << count++;

      for(int i = plen; i< 11; i++)
	ostr << ' ';
      ostr << ' ';
      (*h).first.print_haplotype(ostr,coding);           
      

      ostr.setf(ios::fixed);
      ostr.setf(ios::showpoint);
      ostr.precision(6); 
      ostr <<  setw(12) << (*h).second.Freq;
      if(printheader){
	double var = (h->second.SqPseudoCount) - (h->second.Freq * h->second.Freq);
	if(var<0) var=0;
	ostr <<  setw(12) << sqrt(var);
     
	if(ngroups>1){
	  for(int i = 0; i < ngroups; i++){
	    ostr.setf(ios::fixed);
	    ostr.setf(ios::showpoint);
	    ostr.precision(6); 
	    ostr <<  setw(12) << (*h).second.GroupFreq[i];
	    
	    double var = (h->second.GroupFreqSq[i]) - (h->second.GroupFreq[i] * h->second.GroupFreq[i]);
	    if(var<0) var=0;
	    ostr << setw(12) << sqrt(var);
	  }
	}
      }
      ostr << endl;
    }
  }}
}

void HapList::OutputProbs (ostream & ostr, const vector<int> * coding, double min )
{
 
  int count=1;
  for(ListType::const_iterator h = haplist.begin(); h != haplist.end(); h++){
    if((*h).second.Freq >= min){
      ostr << count++ << " : ";
      (*h).first.print_haplotype(ostr,coding);
      ostr << "(" << (*h).second.Prob << ")" << endl;
    }
  }
}


// Set up an index of the pairs in list that can make up each individual.
// index[ind] contains a vector of  pairs of pointers to haps in the list 
// that can make up each individual
//
// RemoveDuds: whether to remove haplotypes that cannot
// occur in at least one individual in pop
//
// CheckMissing: whether to require a match at the missing positions
// (if not, only unique representatives of each pair of haps at the
// non-missing positions are added)
//
// note: if missing data, standard procedure is simply to add one
// second hap for each first hap compatible with observed data (rather
// than all pairs compatible). When doing trios, want to change
// this to add all pairs compatible with the observed data.
void HapList::MakePairsIndex(vector < vector< pair< ListType::iterator,ListType::iterator> > > & index, vector <CIndividual> & pop, bool RemoveDuds,  bool CheckMissing )
{
 
  vector<ListType::iterator> removelist;
  bool found;
  int Nind = index.size();
 
  for(ListType::iterator h1 = haplist.begin(); h1 !=haplist.end(); h1++){
    bool foundatall = false; // check whether h1 can occur in pop
    for(int ind = 0; ind < Nind; ind++){
      Haplotype CompHap = GetCompHap((*h1).first, pop[ind], found, CheckMissing);
      if(found)
	foundatall = true;

       if((CompHap>=(*h1).first) || !CheckMissing ){ // impose ordering h2>=h1 to avoid duplicates
       
	ListType::iterator h2;
	
	if(found){ //find pointer to complementary haplotype in list 
	  //(returns length of list if not there)
	  if(!CheckMissing){ // add all pairs h1, h2 consistent with observed data to list
	    //h2 = Find(CompHap,pop[ind].get_notmissing_list());
	    
	    //cout << "Adding all pairs consistent with data for individual " << ind << endl;
	    
	    for(h2 = h1; h2 != haplist.end(); h2++){
	      if(CompHap.Matches(h2 -> first, pop[ind].get_notmissing_list())){
		pair<ListType::iterator,ListType::iterator> h(h1,h2);
		index[ind].push_back(h);
	      }
	    }
	  }
	
	  else{
	    h2 = haplist.find(CompHap);
	  } 
	} else
	  h2 = haplist.end();
	
	
	pair<ListType::iterator,ListType::iterator> h(h1,h2);
	
	if(h2!=haplist.end())
	  index[ind].push_back(h);
      }
    }
    if(!foundatall)
      removelist.push_back(h1);
  }
  // remove haplotypes in removelist
  if(RemoveDuds)
    for(vector<ListType::iterator>::iterator i = removelist.begin(); i != removelist.end(); i++)
      haplist.erase(*i);
}


// prune previously created index by removing those whose
// phaseprob is less than threshold
// this needn't be a member of HapList actually...
void HapList::PrunePairsIndex(vector < vector< pair< ListType::iterator,ListType::iterator> > > & index, vector < vector< double > > & phaseprobs,vector <CIndividual> & pop,
double minthreshold)
{
  bool found;
  int Nind = index.size();
  vector < vector< pair< ListType::iterator,ListType::iterator> > > oldindex(index);
  vector < vector< double > > oldphaseprobs(phaseprobs);
  index =  vector < vector< pair< ListType::iterator,ListType::iterator> > > (Nind);
  phaseprobs = vector < vector< double > > (Nind);
  
  for(int ind = 0; ind < Nind; ind++){   
    int i=0;
    for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair = oldindex[ind].begin(); hpair != oldindex[ind].end(); hpair++){
      if(oldphaseprobs[ind][i]>minthreshold){
	index[ind].push_back(*hpair);
	phaseprobs[ind].push_back(0.0); // reset phaseprobs to 0 for final set of iterations...
      }
      i++;
    }
    // if pruning gets rid of all possibilities (ie highly uncertain individual), reinstate them all
    if(phaseprobs[ind].size()==0){
      index[ind] =  vector< pair< ListType::iterator,ListType::iterator> >  (oldindex[ind]);
      phaseprobs[ind] = vector< double >  (oldphaseprobs[ind].size(),0.0);
      
    }
  }
}


// extend previously created index for new record added to list
// also extend phaseprobs (by adding 0s to end) as appropriate
void HapList::ExtendPairsIndex(ListType::iterator & newrecord, vector < vector< pair< ListType::iterator,ListType::iterator> > > & index, vector < vector< double > > & phaseprobs,vector <CIndividual> & pop, bool CheckMissing )
{
  bool found;
  int Nind = index.size();
  
  for(int ind = 0; ind < Nind; ind++){
    Haplotype CompHap = GetCompHap((*newrecord).first, pop[ind], found, CheckMissing);
    
    ListType::iterator h2;
    
    if(found){ //find pointer to complementary haplotype in list 
      //(returns length of list if not there)
      if(!CheckMissing)
	h2 = Find(CompHap,pop[ind].get_notmissing_list());
      else
	h2 = haplist.find(CompHap);
    } else
      h2 = haplist.end();
    
    pair<ListType::iterator,ListType::iterator> h(newrecord,h2);
    if(CompHap<=(*newrecord).first) // impose ordering h2>=newrecord
      h= pair<ListType::iterator,ListType::iterator>(h2,newrecord);
      
    if(h2!=haplist.end()){
      index[ind].push_back(h);
      phaseprobs[ind].push_back(0);
    }
  }
}


// extend previously created index for a particular individual
// for new record added to list
// also extend phaseprobs (by adding 0s to end) as appropriate
// differs from above in that it checks whether the record
// is already in the list
void HapList::ExtendPairsIndex(ListType::iterator & newrecord, vector< pair< ListType::iterator,ListType::iterator> > & index, vector< double > & phaseprobs,
CIndividual & ind, bool CheckMissing )
{
  bool found;
  
  Haplotype CompHap = GetCompHap((*newrecord).first, ind, found, CheckMissing);
    
  ListType::iterator h2;
  
  if(found){ //find pointer to complementary haplotype in list 
    //(returns length of list if not there)
    if(!CheckMissing)
      h2 = Find(CompHap,ind.get_notmissing_list());
    else
      h2 = haplist.find(CompHap);
  } else
    h2 = haplist.end();
  
  pair<ListType::iterator,ListType::iterator> h(newrecord,h2);
  if(CompHap<=(*newrecord).first) // impose ordering h2>=newrecord
    h= pair<ListType::iterator,ListType::iterator>(h2,newrecord);
    
  
  if(h2!=haplist.end() &&  // and h  not in index already
     ( find( index.begin(), index.end(), h) == index.end() ) ){
    index.push_back(h);
    phaseprobs.push_back(0);
  }
}

// use list method to resolve phase in pop
// Only haplotypes in the list are considered possible
// haplotypes for the reconstruction
double HapList::ResolvePhase(char method, int Niter, vector <CIndividual> & pop, vector<ArrayQ *> & Q, vector<double> & vecRho,bool randomise )
{
  if(randomise){
    cout << "Randomising List Frequencies " << endl;
    RandomiseFreqs();
  }

  int Nind = pop.size();

  // Set up an index of the pairs in list that can make up each individual.
  // index[ind] contains a vector of  pairs of pointers to haps in the list 
  // that can make up individual ind
  vector < vector< pair< ListType::iterator,ListType::iterator> > > index(Nind); 

  MakePairsIndex(index, pop, true); // true removes haps from list
  // that cannot be included in sample
  MakePositiveHaps(); // make list of haps with positive freqs
  //(probably all of them!)
  
// output list
//    vector<int> coding [2];
//    int nloci = (haplist.begin()->first).Nloci();
//    for(int i = 0; i<nloci; i++){
//      coding[0].push_back( (int) 'A');
//      coding[1].push_back( (int) 'G');
//    }
//    Output(cout,coding);

  int nchr = Nind+Nind;

  NormaliseFreqs(); // normalise the frequencies to sum to 1

  int iter=0;
  double oldloglik, loglik;
  loglik = 0;
  do{
    iter++;
    cout << "List method, iteration " << iter << endl;
    oldloglik = loglik;
    loglik = 1;
    ClearPseudoCounts(); 
    ComputeProbs(method,Q,nchr,vecRho);
    // make list of pseudo-counts by, for each individual,
    // computing the probability that it is made up of each possible
    // pair; normalising; and adding to Counts
    for(int ind = 0; ind < Nind; ind++){
      double sumprob = 0;
      vector<double> Prob;
      for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair = index[ind].begin(); hpair != index[ind].end(); hpair++){
	double p = ((*hpair).first->second).Prob * 
	  ((*hpair).second->second).Prob;
	sumprob +=p;
	Prob.push_back(p);	
      }
      loglik += log(sumprob);
      int i = 0;
      for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair = index[ind].begin(); hpair != index[ind].end(); hpair++){
	(((*hpair).first)->second).PseudoCount += Prob[i]/sumprob;
	(((*hpair).second)->second).PseudoCount += Prob[i]/sumprob;
	i++;
      }
    } 

    // copy PsedoCounts to Freqs (suitably normalised)
    CopyPseudoCountsToFreqs();
    NormaliseFreqs();

    cout << "Log Likelihood: " << loglik << endl;
    
  } while ( ( abs((loglik-oldloglik)) > EMCONVERGE ) &&  (iter!=Niter) );

  // set the individuals in pop to their best guesses
  SetBestGuesses(pop, index, 'P');

  return loglik;

}


// set each ind in pop to it's best guess, based on the
// haplist (with an index that has been computed already
// to say which pairs can make up each individual)
// if use = 'P' use the Probs in the haplist
// if use = 'F' use the frequencies
// Note: this did not work so well in preliminary tests,
// and is not really the correct thing to do. See the revised
// version that follows
void HapList::SetBestGuesses(vector<CIndividual> & pop, vector < vector< pair< ListType::iterator,ListType::iterator> > > & index, char use)
{
  int Nind = index.size();
  for(int ind = 0; ind < Nind; ind++){
    vector< pair<ListType::iterator,ListType::iterator> >::iterator bestpair;
    double bestprob = 0;
    
    for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair = index[ind].begin(); hpair != index[ind].end(); hpair++){

      double p;
      if(use == 'P')
	p = ((*hpair).first->second).Prob * ((*hpair).second->second).Prob;
      else if (use == 'F')
	p = ((*hpair).first->second).Freq * ((*hpair).second->second).Freq;
      else{
	cerr << "Error in value for use in SetBestGuesses" << endl;
	exit(1);
      }

      if(p>bestprob){
	bestpair = hpair;
	bestprob = p;
      }
    }   
    pop[ind].update_haplotype(0,(*bestpair).first->first);
    pop[ind].update_haplotype(1,(*bestpair).second->first);
    pop[ind].set_overall_correct_prob(bestprob);
  }
}

// as above, but where index has not been computed
void HapList::SetBestGuesses(vector<CIndividual> & pop, char use)
{
  int Nind = pop.size();
  vector < vector< pair< ListType::iterator,ListType::iterator> > > index(Nind); 
  MakePairsIndex(index, pop, false);
  SetBestGuesses(pop, index, use);
}
 

// return the pair in index that has the highest corresponding phaseprob
pair<ListType::iterator,ListType::iterator> HapList::FindBestPair(vector< pair< ListType::iterator,ListType::iterator> > &index, vector< double > & phaseprobs, double & bestphaseprob, int startlocus, int endlocus )
{
  // if startlocus < 0, then include all loci, so just find the pair
  // with the highest phase value
  bestphaseprob = 0;
  pair<ListType::iterator,ListType::iterator> bestpair;
  if(startlocus < 0){
    bestpair = index[0];
    vector< pair<ListType::iterator,ListType::iterator> >::iterator currentpair = index.begin();
    double sum = 0;
    for(vector<double>::iterator phasepointer = phaseprobs.begin();
	phasepointer != phaseprobs.end(); phasepointer++){
      sum += *phasepointer;
      if(*phasepointer > bestphaseprob){
	bestphaseprob = *phasepointer;
	bestpair = *currentpair;
      } 
      currentpair++;
    }
    bestphaseprob /= sum;
    return bestpair;
  }
  else {  
    HapPairList pairprobs; // a map of pairs of haplotypes 
    // (each hap of length endlocus-startlocus) to probabilities
    vector<double>::iterator phasepointer = phaseprobs.begin();
    for(vector< pair< ListType::iterator,ListType::iterator> >::iterator hpair=index.begin(); hpair !=index.end(); hpair++){
      Haplotype hap1 ((hpair->first)->first, startlocus,endlocus);
      Haplotype hap2 ((hpair->second)->first, startlocus,endlocus);
      pair<Haplotype,Haplotype> happair(hap1,hap2);
      pairprobs.Add(happair,*phasepointer);
      phasepointer++;
    }

    pair<Haplotype,Haplotype> besthappair = pairprobs.BestPair(bestphaseprob);

    for(vector< pair< ListType::iterator,ListType::iterator> >::iterator hpair=index.begin(); hpair !=index.end(); hpair++){
      Haplotype hap1((hpair->first)->first, startlocus,endlocus);
      Haplotype hap2((hpair->second)->first, startlocus,endlocus);
      if(((hap1 == besthappair.first) && (hap2 == besthappair.second)) || 
	 ((hap2 == besthappair.first) && (hap1 == besthappair.second))){
	bestpair = *hpair; hpair = index.end();
	hpair--;
      }
    }
    return bestpair;   
  }

}

// set each ind in pop to it's best guess, based on the
// phaseprobs (with an index that has been computed already
// to say which pairs can make up each individual)
void HapList::SetBestGuesses(vector<CIndividual> & pop, vector < vector< pair< ListType::iterator,ListType::iterator> > > & index, vector < vector< double > > & phaseprobs)
{
  int Nind = index.size();
  
  for(int ind = 0; ind < Nind; ind++){
    double bestphaseprob = 0;
    pair<ListType::iterator,ListType::iterator> bestpair = FindBestPair(index[ind],phaseprobs[ind],bestphaseprob);   
    //cout << "Bestphaseprob = " << bestphaseprob << endl;
    //cout << "ind = " << ind << endl;
    pop[ind].update_haplotypes(bestpair.first->first,bestpair.second->first);
    //pop[ind].update_haplotype(1,bestpair.second->first);
    pop[ind].set_overall_correct_prob(bestphaseprob);
  }
  
  
}

//
// compute the product over haplotypes of Pr(Hi | H-i)
//
double HapList::FullDataPseudoLogLikelihood(char method, vector<ArrayQ *> & Q, int nchr, vector<double> & vecRho, double DPRIOR)
{
  double removeamount = (1.0)/nchr;
  double addamount = (1.0)/(nchr-1);
  NormaliseFreqs();

  double loglik = 0;
  for(ListType::const_iterator h = haplist.begin(); h!=haplist.end(); h++){
    cout << h->second.Freq << endl;
    if(h->second.Freq>0){
      double freq = h->second.Freq;      
      SoftRemove(h->first,removeamount);
      NormaliseFreqs();
      loglik += freq * log(CalcProb(h->first,method,Q,nchr-1,vecRho,DPRIOR));
      Add(h->first,addamount);
      NormaliseFreqs();
    }
  }
  return loglik;
}


void HapList::MakePositiveHaps()
{
  PositiveHaps.clear();
  for( ListType::iterator h = haplist.begin(); h!=haplist.end(); h++)
    if(h->second.Freq>0)
      PositiveHaps.push_back(h);
}

// use MCMC-based list method to resolve phase in pop
// Only haplotypes in the list (plus those in the current guess)
//  are considered possible
// haplotypes for the reconstruction
// And at each stage each individual is removed; Hap Probs for that ind
// are computed, and the ind is replaced with randomly reconstructed
// haps according to the appropriate probabilities
// double HapList::MCMCResolvePhaseRemove(char method, int Niter, vector <CIndividual> & pop, vector<ArrayQ *> & Q, vector<double> & vecRho, double betastart, double betaend, bool collectdata)
// {
//   ClearFreqs(); //
//   Add(pop); // start with list of all haplotypes in current guess
//   NormaliseFreqs(); // normalise the frequencies to sum to 1
  
//   ClearPseudoCounts(); // PseudoCounts are used to store posterior means 
//   //of freqs of each hap


//   int Nind = pop.size();
//   vector<int> perm (Nind,0);
//   bool found;


//   // Set up an index of the pairs in list that can make up each individual.
//   // index[ind] contains a vector of  pairs of pointers to haps in the list 
//   // that can make up individual ind
//   vector < vector< pair< ListType::iterator,ListType::iterator> > > index(Nind); 
//   MakePairsIndex(index,pop,true); // true also removes those haps from list
//   // that cannot make up any individuals in pop

//   vector < vector< double > > phaseprobs(Nind); // set up a matrix for storing the prob of each phase for each ind
//   for(int ind = 0; ind< Nind; ind++){
//     phaseprobs[ind] = vector<double>(index[ind].size(),0.0);
//   }

//   int nchr = Nind+Nind;
//   int nchrminus2 = nchr-2;
//   double removeamount = (1.0)/nchr;
//   double addamount = (1.0)/nchrminus2;
  
//   double loglik;

//   // burn-in iterations
//   for(int iter = 0; iter < Niter; iter++){
    
//     //    cout << "MCMC List (Remove) method, burn-in iteration " << iter << endl;
//     loglik = 0;
    
//     //double DPRIOR = (betaend/nchr) + ((Niter-1-iter)*1.0/(Niter-1))*((betastart-betaend)/nchr);
//     double DPRIOR = (betaend) + ((Niter-1-iter)*1.0/(Niter-1))*((betastart-betaend));    

//     rperm(perm,Nind); // take random permutation of inds, to update in random order
//     for(int i = 0; i<Nind; i++){    
//       int ind = perm[i];     
//       //int ind = i;
//       //cout << "updating " << ind << endl;
//       SoftRemove(pop[ind],removeamount);
//       NormaliseFreqs();
//       MakePositiveHaps();
// #ifdef DEBUG
//       cout << "Removed ind " << ind << endl;
//       cout << "Freqs" << endl;
//       Output(cout, coding);
// #endif

//       double sumprob = 0;
//       vector<double> Prob;
//       int randomise = (ranf()<DPRIOR); // here DPRIOR is prob of randomising
//       // the choice of PHASE

//       if(randomise){
// 	for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair = index[ind].begin(); hpair != index[ind].end(); hpair++){
	
// 	// no prior annealing (but DPRIOR comes in later as prob of randomly
// 	// choosing haplotype)
// 	  double p = CalcProb(((*hpair).first->first),method,Q,nchrminus2,vecRho,0) + CalcProb(((*hpair).second->first),method,Q,nchrminus2,vecRho,0);
// 	  Prob.push_back(p);
// 	  sumprob+=p;
// 	} 
//       } else {
// 	for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair = index[ind].begin(); hpair != index[ind].end(); hpair++){
	  
// 	  // no prior annealing at this point (because of possible randomisation above)	
// 	  double p = CalcProb(((*hpair).first->first),method,Q,nchrminus2,vecRho,0) * CalcProb(((*hpair).second->first),method,Q,nchrminus2,vecRho,0);
	  
// 	  // original linear prior annealing (adding DPRIOR to p)
// 	  //double p = (CalcProb(((*hpair).first->first),method,Q,nchrminus2,vecRho,0)+DPRIOR)/(1+DPRIOR*get_listlength()) * 
// 	  //(CalcProb(((*hpair).second->first),method,Q,nchrminus2,vecRho,0)+DPRIOR)/(1+DPRIOR*get_listlength());
	  
// 	  // "powered" version of prior annealing, raising p to power DPRIOR
// 	  //double p = exp(DPRIOR * log(CalcProb(((*hpair).first->first),method,Q,nchrminus2,vecRho,0)*
// 	  //CalcProb(((*hpair).second->first),method,Q,nchrminus2,vecRho,0)));
	  
// #ifdef DEBUG       
// 	  ((*hpair).first->first).print_haplotype(cout,coding);
// 	  cout << ",";
// 	  ((*hpair).second->first).print_haplotype(cout,coding);
// 	  cout << ":" << p << endl;
// #endif
	  
// 	  sumprob +=p;
// 	  Prob.push_back(p);	
// 	}
//       }

//       // store phase probabilities in phaseprobs
//      //   vector<double>::iterator probpointer = Prob.begin();
// //        for (vector<double>::iterator phase = phaseprobs[ind].begin(); phase != phaseprobs[ind].end(); phase++){
// //  	*phase += (1.0/Niter)*(*probpointer/sumprob);
// //  	probpointer++;
// //        }

//       int choice = rint2(Prob,sumprob);
//       loglik += log(sumprob);
//       //loglik += log(Prob[choice]); 
//       vector< pair<ListType::iterator,ListType::iterator> >::iterator sampledpair = index[ind].begin();
//       for(int i = 0; i<choice; i++)
// 	sampledpair++;
      
// #ifdef DEBUG
//       cout << "Chosen Pair:" << endl;
//       OutputPair(cout,*sampledpair);
// #endif

//       pop[ind].update_haplotype(0,(*sampledpair).first->first);
//       pop[ind].update_haplotype(1,(*sampledpair).second->first);
//       Add(pop[ind],addamount);
//       NormaliseFreqs();

// #ifdef DEBUG
//       cout << "Adding ind:" << endl;
//       Output(cout,coding);
// #endif
      
//     }
//     cout << "Log Likelihood: " << loglik << endl;    
//   } 

//   // actual iterations
//   double meanloglik = 0;

//   for(int iter = 0; iter < Niter; iter++){
//     //    cout << "MCMC List (Remove) method, iteration " << iter << endl;

//     loglik = 0;
   
//     //double DPRIOR = betaend/nchr;
//     double DPRIOR = 0;
//     rperm(perm,Nind); // take random permutation of inds, to update in random order
//     for(int i = 0; i<Nind; i++){
//       int ind = perm[i];
//       //int ind = i;
//       SoftRemove(pop[ind],removeamount);
//       NormaliseFreqs();
//       MakePositiveHaps();
// #ifdef DEBUG
//       cout << "Removed ind " << ind << endl;
//       cout << "Freqs" << endl;
//       Output(cout, coding);
// #endif

//       //cout << "Individual " << ind << endl;
//       // make list of probabilities of each hap pair
//       double sumprob = 0;
//       vector<double> Prob;
//       for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair = index[ind].begin(); hpair != index[ind].end(); hpair++){
// 	double p = (CalcProb(((*hpair).first->first),method,Q,nchrminus2,vecRho,0)+DPRIOR)/(1+DPRIOR*get_listlength()) * 
// 	  (CalcProb(((*hpair).second->first),method,Q,nchrminus2,vecRho,0)+DPRIOR)/(1+DPRIOR*get_listlength());

// #ifdef DEBUG       
// 	((*hpair).first->first).print_haplotype(cout,coding);
// 	cout << ",";
// 	((*hpair).second->first).print_haplotype(cout,coding);
// 	cout << ":" << p << endl;
// #endif

// 	sumprob +=p;
// 	Prob.push_back(p);	
//       }

//       //      cout << "Number of possibilities = " << Prob.size() << endl;

//       // add the probabilities of each to PseudoCounts, thus storing
//       // the (Rao-Blackwellised) estimate of the posterior mean of  the freqs
//       int probpos = 0;
//       for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair = index[ind].begin(); hpair != index[ind].end(); hpair++){
// 	haplist[(*hpair).first->first].PseudoCount += Prob[probpos]/sumprob;
// 	haplist[(*hpair).second->first].PseudoCount += Prob[probpos++]/sumprob;
//       }

//       // store phase probabilities in phaseprobs
//       vector<double>::iterator probpointer = Prob.begin();
//       for (vector<double>::iterator phase = phaseprobs[ind].begin(); phase != phaseprobs[ind].end(); phase++){
// 	*phase += (1.0/Niter)*(*probpointer/sumprob);
// 	probpointer++;
//       }

//       int choice = rint2(Prob,sumprob);
//       loglik += log(sumprob);
//       //loglik += log(Prob[choice]); 
//       vector< pair<ListType::iterator,ListType::iterator> >::iterator sampledpair = index[ind].begin();
//       for(int i = 0; i<choice; i++)
// 	sampledpair++;
      
// #ifdef DEBUG
//       cout << "Chosen Pair:" << endl;
//       OutputPair(cout,*sampledpair);
// #endif

//       pop[ind].update_haplotype(0,(*sampledpair).first->first);
//       pop[ind].update_haplotype(1,(*sampledpair).second->first);
//       Add(pop[ind],addamount);
//       NormaliseFreqs();

// #ifdef DEBUG
//       cout << "Adding ind:" << endl;
//       Output(cout,coding);
// #endif
      
//       pop[ind].UpdateCounts();
//     }
//     cout << "Log Likelihood: " << loglik << endl; 
//     meanloglik +=loglik;
//   }

//   // copy PsedoCounts (posterior means of Freqs) to Freqs (suitably normalised)
//   CopyPseudoCountsToFreqs();
//   NormaliseFreqs();

//   SetBestGuesses(pop, index, phaseprobs);
  
  
//   vector<Summary> summary = ProduceSummary(index,phaseprobs,0,get_nloci(),pop);
//   vector<int> coding [2];
//     for(int i = 0; i < get_nloci(); i++){
//       coding[0].push_back( (int) '0');
//       coding[1].push_back( (int) '1');
//     }
    
//     for(vector<Summary>::iterator s = summary.begin(); s<summary.end(); s++)
//       s->Output(cout, coding);
 
//   return meanloglik;

// }



vector<Summary> HapList::ProduceSummary(vector< vector< pair< ListType::iterator,ListType::iterator> > > & index, vector< vector<double> > & phaseprobs, int startlocus, int endlocus, vector<CIndividual> & pop,bool allowsplit )
{
  
  int Nind = index.size();
  vector<Summary> finalanswer(Nind);

  // compute KL divergence for splitting at each possible locus
  vector<double> klsplit (endlocus-startlocus,0.0);
  

  // klsplit[0] has kl div for not splitting at all
  // klsplit[i] has kl div for splitting at locus i
  for(int ind = 0; ind < Nind; ind++){
    HapPairList PossiblePairs(index[ind], phaseprobs[ind], startlocus, endlocus);
    klsplit[0] += PossiblePairs.BestKLdivergence(pop[ind].get_nmissing());
  }

  int bestsplitlocus = 0;

  if(allowsplit){
    for(int ind = 0; ind < Nind; ind++){
      HapPairList PossiblePairs(index[ind], phaseprobs[ind],startlocus,endlocus);
      for(int locus = (startlocus+1); locus < endlocus; locus++){
	klsplit[locus-startlocus] += PossiblePairs.KLsplitdivergence(locus-startlocus, pop[ind].get_nmissing());
      }
    }
    
    int bestsplitlocus = 0;
    double bestkl = klsplit[0];
    for(int locus = (startlocus+1); locus < endlocus; locus++){
      if(klsplit[locus-startlocus]>bestkl){
	bestkl = klsplit[locus-startlocus];
	bestsplitlocus = locus;
      }
    }
  }

  if(bestsplitlocus>0){
    vector<Summary> lhssummary = ProduceSummary(index,phaseprobs,startlocus,bestsplitlocus, pop);
    vector<Summary> rhssummary = ProduceSummary(index,phaseprobs,bestsplitlocus,endlocus, pop);
    for(int ind = 0; ind < Nind; ind++){
      finalanswer[ind] = Summary(lhssummary[ind],rhssummary[ind]);
    }
    return finalanswer;
  }
  else {
    for(int ind = 0; ind < Nind; ind++){
      HapPairList PossiblePairs(index[ind], phaseprobs[ind], startlocus, endlocus);
      finalanswer[ind] = PossiblePairs.Summarise(pop[ind].get_nmissing(), false);
    }
    return finalanswer;
  }
}



// As above..
// use MCMC-based list method to resolve phase in pop
// Only haplotypes in the list are considered possible
// haplotypes for the reconstruction
// But.. at each stage don't update inds until end of iteration
// (quicker per iteration, but mixing worse per iteration)
// double HapList::MCMCResolvePhaseNoRemove(char method, int Niter, vector <CIndividual> & pop, vector<ArrayQ *> & Q, vector<double> & vecRho, bool collectdata)
// {
//   cerr << "Error: this routine not written!" << endl;
//   exit(1);
// }


// output pair of haplotypes
void HapList::OutputPair(ostream & ostr, pair< ListType::iterator,ListType::iterator > hpair, const vector<int> * coding)
{
  hpair.first->first.print_haplotype(ostr,coding);
  ostr << " , ";
  hpair.second->first.print_haplotype(ostr,coding);
}






// returns the prob of a given haplotype, according to method
double HapList::CalcProb(const Haplotype & h, char method, vector<ArrayQ *> & Q, int nchr, vector<double> & vecRho, double DPRIOR, const vector<int> & nmissing, bool ALLSNP, const vector<double> & vecdiffprob , bool fuzzy){
  
  vector<double> dummy; // used in FDLSProb below
 
  switch(method) {

  case 'E': //classic EM method 
    DPRIOR = 0.001;
    return EMProb(h, DPRIOR);
    break;
 
  case 'S': // SSD method
    if(ALLSNP){
      
      // if(abs(SNPSDProb(h, vecdiffprob) - SDProb(h, Q, nchr)) > 0.0001){
// 	cout << "WARNING!" << endl;
// 	cout << "SNP probs compare:" << SNPSDProb(h, vecdiffprob) << "," <<  SDProb(h, Q, nchr) << endl;
//       }
      
      return SNPSDProb(h, vecdiffprob);
    }
    else
      return SDProb(h, Q, nchr);
    break;

  case 'R': // Recom method, with Fearnhead-Donnelly-Li-Stephens conditional
    return FDLSProb(h, Q, nchr, vecRho, dummy, false, true, nmissing, vector<double>(0), 0, fuzzy);
    break;
				    
  default:
    cerr << "Error in method for computing Probs" << endl;
    exit(1);
  
  }

}

void HapList::ComputeProbs(char method, vector<ArrayQ *> & Q, int nchr, vector<double> & vecRho){
  
  switch(method) {

  case 'E': //classic EM method 
    ComputeEMProbs();
    break;
 
  case 'S': // SSD method
    ComputeSDProbs(Q, nchr);
    break;

  case 'R': // Recom method, with Fearnhead-Donnelly-Li-Stephens conditional
    ComputeFDLSProbs(Q, nchr, vecRho);
    break;
				    
  default:
    cerr << "Error in method for computing Probs" << endl;
    exit(1);
  
  }

}

void HapList::ComputeEMProbs(){
  static double DPRIOR = 1.0/20;
  for(ListType::iterator h = haplist.begin(); h!=haplist.end(); h++){
    (h->second).Prob = (h->second).Freq + DPRIOR;
  }
}

void HapList::ComputeSDProbs(vector<ArrayQ *> & Q, int nchr){
  for(ListType::iterator h = haplist.begin(); h!=haplist.end(); h++){
    (h->second).Prob = SDProb( h->first, Q, nchr);
  }
}

void HapList::ComputeFDLSProbs(const vector<ArrayQ *> & Q, int nchr, vector<double> & vecRho){
  vector<double> dummy; // dummy parameter for derivatives in FDLSProb
  for(ListType::iterator h = haplist.begin(); h!=haplist.end(); h++){
    (h->second).Prob = FDLSProb( h->first, Q, nchr, vecRho, dummy, false, true);
  }
}


double HapList::EMProb(const Haplotype & h, double DPRIOR){
  ListType::iterator hiter = haplist.find(h);
  if(hiter!=haplist.end())
    return (hiter->second).Freq + DPRIOR;
  else
    return DPRIOR;
}


// return the SD prob of haplotype h, based on list,
// when all sites are SNPs (more efficient method based on counting
// number of differences)
double HapList::SNPSDProb(const Haplotype & h, const vector<double> & diffprobs){
  //static int NSD=0;
  //cout << "Number of Evaluations = " << NSD++ << endl;
  if(PositiveHaps.begin() == PositiveHaps.end())
    return 1;
  else{
    double p = 0;
    for(vector<ListType::iterator>:: iterator hcopy = PositiveHaps.begin(); hcopy!=PositiveHaps.end(); hcopy++){
      p += ((*hcopy)->second).Freq * diffprobs[NDiff(h, (*hcopy)->first)];
    }
    return p;
  }
}

// return the SD prob of haplotype h, based on list
double HapList::SDProb(const Haplotype & h, vector<ArrayQ *> & Q, int nchr){
  //static int NSD=0;
  //cout << "Number of Evaluations = " << NSD++ << endl;
  if(PositiveHaps.begin() == PositiveHaps.end())
    return 1;
  else{
    double p = 0;
    for(vector<ListType::iterator>:: iterator hcopy = PositiveHaps.begin(); hcopy!=PositiveHaps.end(); hcopy++){
      double indprob = 0; // the prob of h given it copied hcopy

      
      for(int t=0; t<SS; t++){ // add up over possible times
	double tempprob = 1;
	for(int locus = 0; locus < h.Nloci(); locus++){
	  int from = ((*hcopy)->first).get_allele(locus);
	  int targ = h.get_allele(locus);
	  tempprob *= PrHitTarg(locus, nchr, t, from, targ, Q);
	}
	indprob += WEIGHTS[t] * tempprob;
      }

      p += ((*hcopy)->second).Freq * indprob ;
    }
    return p;
  }
}

// return the FDLS prob of haplotype h, based on list
// if computederivs, then also computes derivative of prob with respect
// to each rho_i
//
// at some point probably want to recode to prevent underflow
//
double HapList::FDLSProb(const Haplotype & h, const vector<ArrayQ *> & Q, int nchr, vector<double> & vecRho, vector<double> & vecRhoDeriv, bool computederivs, bool usequad,  vector<int> nmissing, const vector<double> & vecTheta, int Nforcorrection, bool fuzzy)
{
  if(!usequad && vecTheta.size()==0){
    cerr << "Error in call to FDLSProb: if not using quadrature, must specify vecTheta" << endl;
    exit(1);
  }
  
  if(nmissing == vector<int>(0))
    nmissing = vector<int> (h.Nloci(),0);

  //cout << "Number of Evaluations = " << NF++ << endl; 
  if(nchr==0){
    if(computederivs){
      for(vector<double>::iterator i = vecRhoDeriv.begin(); i!=vecRhoDeriv.end(); i++)
	*i = 0;
    }
    return 1;
  }
  else{
    int Nloci = h.Nloci();
    int listlength = PositiveHaps.size();

    vector< vector<double> > Alpha( vector< vector<double> > (Nloci, vector<double>(SS*listlength,0.0)) );
    vector<double> AlphaSum( Nloci );
    vector<double>::iterator AlphaPointer;

    double prob;
    if(!fuzzy)
      prob = ForwardsAlgorithm( h, Q, nchr, vecRho, Alpha, AlphaSum, usequad,  nmissing, false, vecTheta, Nforcorrection);
    else
      prob = FuzzyForwardsAlgorithm( h, Q, nchr, vecRho, Alpha, AlphaSum, usequad,  nmissing, false, vecTheta, Nforcorrection);

    //Alpha[r][j] = Pr(Y1...Yr and Xr = j), where Y1..Yr are observed
    // alleles at first r sites, and Xr is position in list of copied 
    // chromosome at locus r
    // final probability we return is sum over j of Alpha[last locus][j]
  
    if(computederivs){
      vector< vector<double> > Beta( vector< vector<double> > (Nloci, vector<double>(SS*listlength,0.0)) );
      vector<double> BetaSum (Nloci);
  
      //Beta[r][j] = Pr(Ylastlocus...Y(lastlocus -r) and X(lastlocus-r) = j), where Yr..Ylastlocus
      // are observed
      // alleles at loci r to end, and Xr is position in list of copied 
      // chromosome at locus r
      // ie Beta is version of Alpha, starting at Rhs
      if(!fuzzy)
	ForwardsAlgorithm( h, Q, nchr, vecRho, Beta, BetaSum, usequad, nmissing, true, vecTheta, Nforcorrection );
      else
	FuzzyForwardsAlgorithm( h, Q, nchr, vecRho, Beta, BetaSum, usequad, nmissing, true, vecTheta, Nforcorrection );

      // do work to compute derivatives
      vector<double> TransProb = vector<double>(vecRho.size());
     
      int Sforcorrection = Nloci;
      for(int locus=0; locus<TransProb.size(); locus++){
	TransProb[locus] = 1 - exp(-( correction(Nforcorrection, Sforcorrection, vecRho[locus]) * vecRho[locus]/nchr));
      }

      vector<double>::iterator BetaPointer;
      
      for(int locus = 0; locus< Nloci-1; locus++){
	BetaPointer = Beta[(locus+1)].begin();
	AlphaPointer = Alpha[locus].begin();
	double AlphaBetaSum=0; // find sum of Alpha * Beta /weights * freq
	for(vector<ListType::iterator>:: iterator hcopy = PositiveHaps.begin(); hcopy!=PositiveHaps.end(); hcopy++){
	  if(usequad){
	    for(int t=0; t<SS; t++){
	      AlphaBetaSum+= (*AlphaPointer) * (*BetaPointer)/ ( (((*hcopy)->second).Freq / nchr) * WEIGHTS[t]);
	      AlphaPointer++;
	      BetaPointer++;
	    }
	  }
	  else{
	    AlphaBetaSum+= (*AlphaPointer) * (*BetaPointer)/ ( ((*hcopy)->second).Freq / nchr);
	    AlphaPointer++;
	    BetaPointer++;
	  }
	}

	 //cout << "Complex sum is " << AlphaBetaSum *(1-TransProb[locus]) + TransProb[locus]*AlphaSum[locus] * BetaSum[locus+1] << endl;

	 vecRhoDeriv[locus] = derivcorrection(Nforcorrection, Sforcorrection, vecRho[locus]) * (1- TransProb[locus]) * (AlphaSum[locus] * BetaSum[(locus+1)] - AlphaBetaSum) /  nchr;
      } 

      //cout << "Difference in sums is" << BetaSum[0] - AlphaSum[lastlocus] << endl;
      //cout << "AlphaSum is " << AlphaSum[lastlocus] << endl;
      //cout << "BetaSum is " << BetaSum[0] << endl;
      
    }

    //cout << "Prob compare: " << sum << "," << SDProb( h, Q, nchr) << endl;
    return prob;
  }
}


//compute Alpha[r][j] = Pr(Y1...Yr and Xr = j), where Y1..Yr are observed
// alleles at first r sites, and Xr is position in list of copied 
// chromosome at locus r
// final probability we return is sum over j of Alpha[last locus][j]
// (note: actually not a probability, unless Freq is a set of counts (with totol nchr))

double HapList::ForwardsAlgorithm(const Haplotype & h, const vector<ArrayQ *> & Q, int nchr,  vector<double> & vecRho, vector< vector<double> > & Alpha, vector<double> & AlphaSum, bool usequad, const vector<int> & ismissing, bool gofromrhs , const vector<double> & vecTheta, int Nforcorrection )
{
 
  if(!usequad && vecTheta.size()==0){
    cerr << "Error in call to Forwards Algorithm: if not using quadrature, must specify vecTheta" << endl;
    exit(1);
  }

  int Nloci = h.Nloci();
  vector<double> TransProb = vector<double>(Nloci - 1);
  int Sforcorrection = Nloci;
  
  for(int locus=0; locus<TransProb.size(); locus++)
    TransProb[locus] = 1 - exp(-( correction(Nforcorrection, Sforcorrection, vecRho[locus]) * vecRho[locus]/nchr));
  
  int firstlocus = 0;
  if(gofromrhs)
    firstlocus = Nloci - 1;

  // compute Alpha[firstlocus][.]
  vector<double>::iterator AlphaPointer = Alpha[firstlocus].begin();
  vector<double>::iterator OldAlphaPointer;
  AlphaSum[firstlocus] = 0;

  for(vector<ListType::iterator>::iterator hcopy = PositiveHaps.begin(); hcopy!=PositiveHaps.end(); hcopy++){	
    int from = ((*hcopy)->first).get_allele(firstlocus);
    int targ = h.get_allele(firstlocus);
    if(usequad){
      for(int t = 0; t<SS; t++){ // loop over the SS quadrature points
	(*AlphaPointer) = ((*hcopy)->second).Freq/nchr * WEIGHTS[t];
	if(ismissing[firstlocus]==0)
	  (*AlphaPointer) *= PrHitTarg(firstlocus, nchr, t, from, targ, Q);
	AlphaSum[firstlocus] += (*AlphaPointer);
	AlphaPointer++;
      }
    } else {
      (*AlphaPointer) = ((*hcopy)->second).Freq/nchr;
      if(ismissing[firstlocus]==0)
	(*AlphaPointer) *= PrHitTarg(nchr, from, targ, vecTheta[firstlocus] );
      AlphaSum[firstlocus] += (*AlphaPointer);
      AlphaPointer++;
    }
  }
  
 
  int locus=0;
  int previouslocus = firstlocus;
  // now compute each of the other alphas
  for(int r = 1; r < Nloci; r++){
    locus = r;
    if(gofromrhs)
      locus = Nloci - 1 - r;

    AlphaSum[locus] = 0; // stores sum of Alpha[locus]
    
    double TProb = TransProb[previouslocus]; // prob of transition
    if(gofromrhs)
      TProb = TransProb[locus];

    AlphaPointer = Alpha[locus].begin();
    OldAlphaPointer = Alpha[previouslocus].begin();
    
    for(vector<ListType::iterator>:: iterator hcopy = PositiveHaps.begin(); hcopy!=PositiveHaps.end(); hcopy++){	
      int from = ((*hcopy)->first).get_allele(locus);
      int targ = h.get_allele(locus);
      if(usequad){
	for(int t=0; t<SS; t++){
	  *AlphaPointer = (*OldAlphaPointer * (1-TProb) + AlphaSum[previouslocus] * TProb * (((*hcopy)->second).Freq / nchr) * WEIGHTS[t]);
	  if (ismissing[locus]==0)
	    *AlphaPointer *= PrHitTarg(locus, nchr, t, from, targ, Q);
	  AlphaSum[locus] += *AlphaPointer;
	  OldAlphaPointer++;
	  AlphaPointer++;
	}
      }
      else{
	*AlphaPointer = (*OldAlphaPointer * (1-TProb) + AlphaSum[previouslocus] * TProb * (((*hcopy)->second).Freq / nchr));
	if (ismissing[locus]==0)
	  *AlphaPointer *= PrHitTarg(nchr, from, targ, vecTheta[locus]);
	AlphaSum[locus] += *AlphaPointer;
	OldAlphaPointer++;
	AlphaPointer++;
      }
    }
    previouslocus = locus;
  }
  
  return AlphaSum[previouslocus];

}


// simulate backwards along chromosome, the haplotype that was copied by Haplotype h
void HapList::BackwardsAlgorithm(const Haplotype & h, int nchr, vector<double> & vecRho, vector< vector<double> > & Alpha, vector<double> & AlphaSum, vector<int> & copiedtime, vector<int> & copiedallele, vector<int> & copiedhap, bool usequad )
{
  int Nloci = h.Nloci();
  vector<double> TransProb = vector<double>(Nloci-1);
  for(int locus = 0; locus < TransProb.size(); locus++){ 
    TransProb[locus] = 1 - exp(-(vecRho[locus]/nchr));
  }

  for(int locus = Nloci - 1; locus>0; locus--){    
    int copied = rint2(Alpha[locus]);   

#ifdef DEBUG       
    cout << "First copied =" << copied << endl; 
    cout << "Score at copied = " << Alpha[locus][copied] << endl;
    cout << "Score at previous locus at copied = " << Alpha[locus-1][copied] << endl;
#endif

    copiedhap[locus] = copied / SS;
    copiedallele[locus] = PositiveHaps[copied / SS]->first.get_allele(locus);
    copiedtime[locus] = copied % SS;
    
    vector<double>::iterator AlphaPointer = Alpha[locus-1].begin();
    double TProb = TransProb[locus-1]; // prob of transition
    AlphaSum[locus-1] = 0;

    int ptr = 0;
    for(vector<ListType::iterator>:: iterator hcopy = PositiveHaps.begin(); hcopy!=PositiveHaps.end(); hcopy++){	
      if(usequad){
	for(int t=0; t<SS; t++){
	  double transprob = TProb * ((*hcopy)->second).Freq/nchr * WEIGHTS[t];
	  if(ptr == copied)
	    transprob += (1-TProb);
	  //cout << "AlphaSum =" << AlphaSum[locus] << endl;
	  *AlphaPointer = (*AlphaPointer/AlphaSum[locus]) * transprob;
	  //cout << ptr << "," << *AlphaPointer << endl;
	  AlphaSum[locus-1] += *AlphaPointer;
	  AlphaPointer++;
	  ptr ++;
	}
      }
      else{
	double transprob = TProb * ((*hcopy)->second).Freq/nchr;
	if(ptr == copied)
	  transprob += (1-TProb);
	*AlphaPointer = (*AlphaPointer/AlphaSum[locus]) * transprob;
	AlphaSum[locus-1] += *AlphaPointer;
	AlphaPointer++;
	ptr ++;
      }	
    }
  }

  int copied = rint2(Alpha[0]);
  //cout << "copied =" << copied << endl;
  copiedallele[0] = PositiveHaps[copied / SS]->first.get_allele(0);
  copiedtime[0] = copied % SS;
}

// compute the probability that Haplotype h copied each
// allele at each time
void HapList::ComputeHiddenStateProbs(vector<vector<double> > & CopyProb, const Haplotype & h, const vector<ArrayQ *> & Q, int nchr,  vector<double> & vecRho, bool usequad, const vector<int> & isunknown , const vector<double> & vecTheta, int Nforcorrection)
{

  int Nloci = h.Nloci();
  int listlength = PositiveHaps.size();
  //cout << "listlength = " << listlength << endl;

  vector< vector<double> > Alpha( vector< vector<double> > (Nloci, vector<double>(SS*listlength,0.0)) );
  vector<double> AlphaSum( Nloci );
  vector<double>::iterator AlphaPointer;

  vector< vector<double> > Beta( vector< vector<double> > (Nloci, vector<double>(SS*listlength,0.0)) );
  vector<double> BetaSum( Nloci );
  vector<double>::iterator BetaPointer;

  double a=ForwardsAlgorithm(h, Q, nchr, vecRho, Alpha, AlphaSum, usequad, isunknown, false , vecTheta, Nforcorrection );
  double b=ForwardsAlgorithm(h, Q, nchr, vecRho, Beta, BetaSum, usequad, isunknown, true , vecTheta, Nforcorrection );

  // cout << "a=" << a << endl;
//   cout << "b=" << b << endl;
//   cout << Alpha[Nloci-1][0] << "," << Alpha[Nloci-1][1] << "," << Alpha[Nloci-1][2] << "," << Alpha[Nloci-1][3] << endl;
//   cout << Beta[0][0] << "," << Beta[0][1] << "," << Beta[0][2] << "," << Beta[0][3] << endl;

  
  for(int locus = 0; locus < Nloci; locus++){    

    //cout << "locus: " << locus << endl;

    AlphaPointer = Alpha[locus].begin();
    BetaPointer = Beta[locus].begin();
    AlphaSum[locus]=0;
    BetaSum[locus]=0;

    for(vector<ListType::iterator>:: iterator hcopy = PositiveHaps.begin(); hcopy!=PositiveHaps.end(); hcopy++){	
      int from = ((*hcopy)->first).get_allele(locus);
      int targ = h.get_allele(locus);
      //cout << "from = " << from << "; targ = " << targ << endl;
      if(usequad){
	for(int t=0; t<SS; t++){
	  double prob =1.0;
	  if (isunknown[locus]==0)
	    prob = PrHitTarg(locus, nchr, t, from, targ, Q);
	  //cout << *AlphaPointer << "," << *BetaPointer << endl;
	  
	  *AlphaPointer *= (*BetaPointer);
	  *AlphaPointer /= (prob*prob);
	  AlphaSum[locus] += *AlphaPointer;
	  BetaSum[locus] += *BetaPointer;
	  AlphaPointer++;
	  BetaPointer++;
	}
      }
      else{
	double prob = 1.0;
	if (isunknown[locus]==0)
	  prob = PrHitTarg(nchr, from, targ, vecTheta[locus]);
	*AlphaPointer *= (*BetaPointer);
	*AlphaPointer /= (prob*prob);
	AlphaSum[locus] += *AlphaPointer;
	BetaSum[locus] += *BetaPointer;
	AlphaPointer++;
	BetaPointer++;
      }

    }
    
  }

  // OUTPUT ALPHA AND ALPHASUM
  // for(int locus = 0; locus < Nloci; locus++){    
//     cout << "locus: " << locus << endl;
//     AlphaPointer = Alpha[locus].begin();
//     for(vector<ListType::iterator>:: iterator hcopy = PositiveHaps.begin(); hcopy!=PositiveHaps.end(); hcopy++){
//       if(usequad){
// 	for(int t=0; t<SS; t++){
//     	  cout << *AlphaPointer/AlphaSum[locus] << "+" ;
// 	  AlphaPointer++;
// 	}
//       }
//       else{
// 	cout << *AlphaPointer << "," ;
// 	AlphaPointer++;
//       }
//     }
//     cout << endl;
//   }
  
  for(int locus = 0; locus < Nloci; locus++){    
    AlphaPointer = Alpha[locus].begin();
    
    for(vector<ListType::iterator>:: iterator hcopy = PositiveHaps.begin(); hcopy!=PositiveHaps.end(); hcopy++){	
      int from = ((*hcopy)->first).get_allele(locus);
      //cout << "from: " << from << endl;
      if(usequad){
	for(int t = 0; t< SS; t++){
	  CopyProb[locus][from*SS + t] += *AlphaPointer/AlphaSum[locus];
	  AlphaPointer++;
	}
      } else {
	  CopyProb[locus][from] += *AlphaPointer/AlphaSum[locus];
	  AlphaPointer++;
      }
    }
  }

  // for(int locus = 0; locus < Nloci; locus++){ 
//     cout << "locus: " << locus << endl;
    
//     for(int allele = 0; allele < 2; allele++){
//       if(usequad){
// 	for(int t = 0; t< SS; t++)
// 	  cout << CopyProb[locus][allele * SS + t] << " ";
//       } else
// 	cout << CopyProb[locus][allele] << " ";
//     }
    
//     cout << endl;
//   }
  
      

}


//compute Alpha[r][j] = Pr(Y1...Yr and Xr = j), where Y1..Yr are observed
// alleles at first r sites, and Xr is position in list of copied 
// chromosome at locus r
// final probability we return is sum over j of Alpha[last locus][j]
// (note: actually not a probability, unless Freq is a set of counts (with totol nchr))

double HapList::FuzzyForwardsAlgorithm(const Haplotype & h, const vector<ArrayQ *> & Q, int nchr,  vector<double> & vecRho, vector< vector<double> > & Alpha, vector<double> & AlphaSum, bool usequad, const vector<int> & ismissing, bool gofromrhs , const vector<double> & vecTheta, int Nforcorrection )
{
 
  if(!usequad && vecTheta.size()==0){
    cerr << "Error in call to Forwards Algorithm: if not using quadrature, must specify vecTheta" << endl;
    exit(1);
  }

  int Nloci = h.Nloci();
  vector<double> TransProb = vector<double>(Nloci - 1);
  int Sforcorrection = Nloci;
  
  for(int locus=0; locus<TransProb.size(); locus++)
    TransProb[locus] = 1 - exp(-( correction(Nforcorrection, Sforcorrection, vecRho[locus]) * vecRho[locus]/nchr));
  
  int firstlocus = 0;
  if(gofromrhs)
    firstlocus = Nloci - 1;

  // compute Alpha[firstlocus][.]
  vector<double>::iterator AlphaPointer = Alpha[firstlocus].begin();
  vector<double>::iterator OldAlphaPointer;
  AlphaSum[firstlocus] = 0;

  for(vector<ListType::iterator>::iterator hcopy = PositiveHaps.begin(); hcopy!=PositiveHaps.end(); hcopy++){	
    float from = ((*hcopy)->first).get_fuzzyallele(firstlocus);
    int targ = h.get_allele(firstlocus);
    if(usequad){
      for(int t = 0; t<SS; t++){ // loop over the SS quadrature points
	(*AlphaPointer) = ((*hcopy)->second).Freq/nchr * WEIGHTS[t];
	if(ismissing[firstlocus]==0)
	  (*AlphaPointer) *= FuzzyPrHitTarg(firstlocus, nchr, t, from, targ, Q);
	AlphaSum[firstlocus] += (*AlphaPointer);
	AlphaPointer++;
      }
    } else {
      (*AlphaPointer) = ((*hcopy)->second).Freq/nchr;
      if(ismissing[firstlocus]==0)
	(*AlphaPointer) *= FuzzyPrHitTarg(nchr, from, targ, vecTheta[firstlocus] );
      AlphaSum[firstlocus] += (*AlphaPointer);
      AlphaPointer++;
    }
  }
  
 
  int locus=0;
  int previouslocus = firstlocus;
  // now compute each of the other alphas
  for(int r = 1; r < Nloci; r++){
    locus = r;
    if(gofromrhs)
      locus = Nloci - 1 - r;

    AlphaSum[locus] = 0; // stores sum of Alpha[locus]
    
    double TProb = TransProb[previouslocus]; // prob of transition
    if(gofromrhs)
      TProb = TransProb[locus];

    AlphaPointer = Alpha[locus].begin();
    OldAlphaPointer = Alpha[previouslocus].begin();
    
    for(vector<ListType::iterator>:: iterator hcopy = PositiveHaps.begin(); hcopy!=PositiveHaps.end(); hcopy++){	
      float from = ((*hcopy)->first).get_fuzzyallele(locus);
      int targ = h.get_allele(locus);
      if(usequad){
	for(int t=0; t<SS; t++){
	  *AlphaPointer = (*OldAlphaPointer * (1-TProb) + AlphaSum[previouslocus] * TProb * (((*hcopy)->second).Freq / nchr) * WEIGHTS[t]);
	  if (ismissing[locus]==0)
	    *AlphaPointer *= FuzzyPrHitTarg(locus, nchr, t, from, targ, Q);
	  AlphaSum[locus] += *AlphaPointer;
	  OldAlphaPointer++;
	  AlphaPointer++;
	}
      }
      else{
	*AlphaPointer = (*OldAlphaPointer * (1-TProb) + AlphaSum[previouslocus] * TProb * (((*hcopy)->second).Freq / nchr));
	if (ismissing[locus]==0)
	  *AlphaPointer *= FuzzyPrHitTarg(nchr, from, targ, vecTheta[locus]);
	AlphaSum[locus] += *AlphaPointer;
	OldAlphaPointer++;
	AlphaPointer++;
      }
    }
    previouslocus = locus;
  }
  
  return AlphaSum[previouslocus];

}


void HapList::ComputeVectorOfNaiveGibbsProbs(CIndividual & ind, vector<double> & TempProb, double & TempProbSum, double dirprior){ 
  
  //  int found0,found1;
  bool found=false; 
  vector<double> CompFreq (get_listlength(),0); // list of numbers of comp haps.
  
  int listpos = 0;

  
  for(ListType::const_iterator h = haplist.begin(); h != haplist.end(); h++){

    // attempt to extract ith haplotype in list from genotype [n1]
    
    //  found=1; 
//      const vector<int> unknown_list =  ind.get_unknown_pos();
//      for (vector<int>::const_iterator u = unknown_list.begin();
//  	 u != unknown_list.end(); ++u ) {
//        //cout << *u << endl;

//        int haplistallele = (*h).first.get_allele(*u);
      
//        if(ind.get_haplotype(0,*u)!=haplistallele){
//  	if(ind.get_haplotype(1,*u)==haplistallele)
//  	  ind.flip_phase(*u); // get hap 0 in n1 to be haplist[listpos] 
//  	else{
//  	  found=0;
//  	  found0=0;
//  	  found1=0;
//  	  goto JUMP1; // end loop
//  	}
//        }
//      }

//      if(found==1){ // check the "known" phases to see they fit in
//        found0=1; // found0 indicates whether the inds hap0 matches
//        found1=1; // found1 indicates whether the inds hap1 matches
//        const vector<int> known_list =  ind.get_known_pos();
//        for (vector<int>::const_iterator u = known_list.begin();
//  	   u != known_list.end(); ++u ) {
//  	int haplistallele = (*h).first.get_allele(*u);
	
//  	if(ind.get_haplotype(0,*u) != haplistallele)
//  	  found0=0;
//  	if(ind.get_haplotype(1,*u) != haplistallele)
//  	  found1=0;
//  	if((found0 == 0) && (found1 ==0))
//  	  goto JUMP1; // jump out of loop
//        }
//      }
    
//    JUMP1:
//      if(found0!=1){
//        if(found1==1){
//  	const vector<int> known_list =  ind.get_known_pos();
//  	for (vector<int>::const_iterator u = known_list.begin();
//  	     u != known_list.end(); ++u )
//  	  ind.flip_phase(*u);
//        }
//        else
//  	found=0;
//      }
    
    //  Haplotype CompHap = ind.get_haplotype(1);

    Haplotype CompHap = GetCompHap((*h).first, ind, found);

    if(found){ //find frequency of complementary haplotype in list 
      //(returns length of list if not there)
      ListType::iterator i = haplist.lower_bound(CompHap);
     //   if(i!=haplist.end()) 
      if((i!=haplist.end()) && ((i->first) == CompHap))
	CompFreq[listpos] = i->second.Freq;
      else
	CompFreq[listpos] = 0;	
    }
    else
      CompFreq[listpos]=-1;
    listpos++;
  }
      
  TempProbSum=0;
  
  listpos = 0;
  for(ListType::const_iterator h = haplist.begin(); h != haplist.end(); h++){
    // compute TempProb as product of freq of h * CompFreq
    if(CompFreq[listpos]>=0){ //if that haplotype was found
      TempProb[listpos]=(*h).second.Freq+dirprior;
      if(CompFreq[listpos]>0)
	TempProb[listpos] *= (CompFreq[listpos]+dirprior);
      else
	TempProb[listpos]*=(2*dirprior);
      TempProb[listpos]-=dirprior*dirprior;
    }
    else
      TempProb[listpos]=0;
    TempProbSum+=TempProb[listpos];
    listpos++;
  }
}


int HapList::Find(CIndividual & ind, int chr, bool usebestguess){
  Haplotype hap = ind.get_haplotype(chr);
  if(usebestguess)
    hap = ind.BestHaplotype(chr);

  int listpos = 0;
  for(ListType::const_iterator h = haplist.begin(); h != haplist.end(); h++){
    if((*h).first == hap){
      return listpos;
    } 
    listpos++;
  }
  return (-1);
}

    
