#include "Summary.hpp"
#include "Haplotype.hpp"

#include <iomanip>

using namespace std;

// constructor
Summary::Summary () :
  bestguess(),
  flipprob(),
  errorprob(),
  alleleprob()
{
}

// copy constructor
Summary::Summary (const Summary & s2):
  bestguess(s2.bestguess),
  flipprob(s2.flipprob),
  errorprob(s2.errorprob),
  alleleprob(s2.alleleprob)
{
}

// merge constructor
Summary::Summary (const Summary & s1, const Summary & s2):
  bestguess(s1.bestguess),
  flipprob(s1.flipprob),
  errorprob(s1.errorprob),
  alleleprob(s1.alleleprob)
{
  bestguess.insert(bestguess.end(),s2.bestguess.begin(),s2.bestguess.end());
  flipprob.insert(flipprob.end(),s2.flipprob.begin(),s2.flipprob.end());
  errorprob.insert(errorprob.end(),s2.errorprob.begin(),s2.errorprob.end());
  alleleprob.insert(alleleprob.end(),s2.alleleprob.begin(),s2.alleleprob.end());
  
}


Summary::Summary ( const HPairType & hpair, const vector<double> & fp,
const vector< vector<double> > & ep, const vector< vector< vector <double> > > & ap):
  bestguess(),
  flipprob(),
  errorprob(),
  alleleprob()
{
  bestguess.push_back(hpair);
  flipprob.push_back(fp);
  errorprob.push_back(ep);
  alleleprob.push_back(ap);
}


const Summary & Summary::operator=(const Summary & rhs)
{
  if(this != &rhs){
    this->bestguess = vector< HPairType >(rhs.bestguess);
    this->flipprob = vector< vector<double> >(rhs.flipprob);
    this->errorprob = vector < vector< vector<double> > > (rhs.errorprob);
    this->alleleprob = vector< vector < vector< vector<double> > > > (rhs.alleleprob);
  }
  return *this;
}

Summary::~Summary()
{
}

void Summary::Output( ostream & ostr, vector<int> * coding)
{
  vector<int> tempcoding [2];
  double startlocus = 0;
  vector<int>::iterator codingstart0 = coding[0].begin();
  vector<int>::iterator codingstart1 = coding[1].begin();
  
  for(vector<HPairType>::iterator h=bestguess.begin(); h!=bestguess.end(); h++){
    tempcoding[0] = vector<int>(codingstart0,coding[0].end());
    tempcoding[1] = vector<int>(codingstart1,coding[1].end());
    h->first.print_haplotype(ostr, tempcoding);
    codingstart0 += h->first.get_nloci();
    codingstart1 += h->first.get_nloci();
    ostr << " || ";
  }
  ostr << endl;

  codingstart0 = coding[0].begin();
  codingstart1 = coding[1].begin();
  for(vector<HPairType>::iterator h=bestguess.begin(); h!=bestguess.end(); h++){
    tempcoding[0] = vector<int>(codingstart0,coding[0].end());
    tempcoding[1] = vector<int>(codingstart1,coding[1].end());
    h->second.print_haplotype(ostr, tempcoding);
    codingstart0 += h->second.get_nloci();
    codingstart1 += h->second.get_nloci();
    ostr << " || ";
  }
  ostr << endl;
  for(vector< vector<double> >::iterator pvec = flipprob.begin(); pvec!=flipprob.end(); pvec++){
    for(vector< double >::iterator p = pvec->begin(); p!= pvec->end(); p++){
      ostr.setf(ios::fixed);
      ostr.setf(ios::showpoint);
      ostr.precision(2); 
      ostr << *p << ";" ;
    }
    ostr << " || ";
  }
  ostr << endl;

  for(vector< vector< vector<double> > >::iterator pvec = errorprob.begin(); pvec!=errorprob.end(); pvec++){

    for(vector< vector<double> >::iterator p = (*pvec).begin(); p!= (*pvec).end(); p++){
      ostr.setf(ios::fixed);
      ostr.setf(ios::showpoint);
      ostr.precision(2); 
      ostr << (*p)[0] << ":" ;
    }
    ostr << " || ";
  }
  ostr << endl;

  for(vector< vector< vector<double> > >::iterator pvec = errorprob.begin(); pvec!=errorprob.end(); pvec++){
    for(vector< vector<double> >::iterator p = (*pvec).begin(); p!= (*pvec).end(); p++){
      ostr.setf(ios::fixed);
      ostr.setf(ios::showpoint);
      ostr.precision(2); 
      ostr << (*p)[1] << ":" ;
    }
    ostr << " || ";
  }

  ostr << endl;



}
