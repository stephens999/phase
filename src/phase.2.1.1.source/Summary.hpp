/*
 * A class for summarising a distribution on haplotypes
 */

#ifndef SUMMARY_H
#define SUMMARY_H

#include "Haplotype.hpp"
#include <map>
using namespace::std;

typedef pair<Haplotype,Haplotype> HPairType;

// summarises an individual as independent sections
// each section has a bestguess (pair of haps), and
// a flipprob vector (estimates of individual phase error probs)

class Summary {

public:
  vector< HPairType > bestguess; 
  vector< vector<double> > flipprob; // indices are (segment, locus)
  vector< vector< vector<double> > > errorprob;// indices are (segment, locus, chrom)
  vector< vector < vector < vector < double > > > > alleleprob; // indices are (segment, locus, chrom, allele)
  
public:
  Summary ();
  Summary (const Summary &);
  Summary (const Summary & s1, const Summary & s2);
  Summary ( const HPairType & hpair, const vector<double> & fp, const vector< vector<double> > & ep, const vector< vector< vector <double> > > & ap);
  
  
  const Summary & operator = (const Summary &);  

  ~Summary();
 

  void Output(ostream & ostr, vector<int> * coding);
  
};



#endif
