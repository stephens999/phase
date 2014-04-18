// Implementation of class CIndividual 
//
// CVS : $Id: indnode.cpp,v 1.33 2003/06/14 00:24:05 stephens Exp $
//

#include "indnode.hpp"
#include "constants.hpp"
#include "utility.hpp"
#include "Haplotype.hpp"
#include "Summary.hpp"

#include <strstream>
#include <algorithm>
#include <iomanip>

using namespace std;

extern int NIND;

// static member
int CIndividual::counts = 0;

// Constructors and destructor 

// Default constructor
CIndividual::CIndividual ( string loci_type ) : 
    phase   ( loci_type.size() ),
    recom (0),
    missing ( loci_type.size() ),
    notmissing_list ( 0 ),
    unknown ( 0 ),
    known   ( 0 ),
    phenotype(2,Haplotype(loci_type)),
    orig_phenotype(2,vector<int>(loci_type.size())),
    saved_hap(2,vector<int>(loci_type.size())),
    Z(2,vector<int>(loci_type.size())),
    AlleleCount(2,vector< vector<float> > (loci_type.size(), vector<float> ())),
    PhaseCount(loci_type.size(), vector<float>(2, 1.0) )
{ 
  NIND ++;

  

#if DEBUG == 5
    cout << "Constructor called ... ";
    cout.flush();
#endif // DEBUG
    int nloci = loci_type.size();
    //cout << "Creating Ind, nloci = " << nloci << endl;
    for(int c=0; c < 2; c++){
      for(int locus = 0; locus < nloci; locus++){
	int n_allele = (loci_type[locus] == 'S' ) ? 2 : KMAX;
	AlleleCount[c][locus] = vector<float> (n_allele,0.0);
      }
    }
       
    for (int i = 0; i < nloci; ++i) {
      phenotype[0].set_allele(i,MISSMS);
      phenotype[1].set_allele(i,MISSMS);
      orig_phenotype[0][i] = MISSMS;
      orig_phenotype[1][i] = MISSMS;
      saved_hap[0][i] = MISSMS;
      saved_hap[1][i] = MISSMS;
      Z[0][i] = (int) (NANCHAP * ranf());
      Z[1][i] = (int) (NANCHAP * ranf());
    }
    

#if DEBUG == 5
    cout << "OK" << endl;
#endif // DEBUG
}

CIndividual::CIndividual( const CIndividual & C ) :
  id (C.id),
  phase   ( C.phase ),
  recom (C.recom ),
  missing ( C.missing ),
  notmissing_list ( C.notmissing_list ),
  unknown ( C.unknown ),
  known   ( C.known ),   
  phenotype ( C.phenotype ),
  orig_phenotype (C.orig_phenotype),
  saved_hap (C.saved_hap),
  Z (C.Z),
  PhaseCount (C.PhaseCount),
  AlleleCount (C.AlleleCount)
{ 
  NIND ++;
#if DEBUG == 5
    cout << "Copy constructor called ... ";
    cout.flush();
#endif // DEBUG
    int nloci = C.phase.size();

#if DEBUG == 5
    cout << "OK" << endl;
#endif // DEBUG
}


CIndividual::CIndividual( const CIndividual & C, int firstlocus, int lastlocus ) :
  id (C.id),
  phase   ( C.phase.begin()+firstlocus, C.phase.begin()+lastlocus ), 
  recom (0),
  missing (C.missing.begin()+firstlocus, C.missing.begin()+lastlocus),
  notmissing_list(0),
  unknown(0),
  known(0),
  phenotype(2),
  orig_phenotype(2,vector<int>(lastlocus - firstlocus)),
  saved_hap(2,vector<int>(lastlocus - firstlocus)),
  Z(2,vector<int>(lastlocus - firstlocus)),
  AlleleCount(2,vector< vector<float> > (lastlocus - firstlocus, vector<float> ())),
  PhaseCount(lastlocus - firstlocus, vector<float>(2, 1.0) )
  
{ 
  NIND ++;
  phenotype[0] = Haplotype(C.phenotype[0],firstlocus,lastlocus);
  phenotype[1] = Haplotype(C.phenotype[1],firstlocus,lastlocus);

 
  for(vector<int>::const_iterator i = C.notmissing_list.begin(); i!=C.notmissing_list.end(); i++){
    if( ((*i)>=firstlocus) && ((*i)<lastlocus)){
      notmissing_list.push_back((*i)-firstlocus);
    }
  }
  
  for(vector<int>::const_iterator i = C.unknown.begin(); i!=C.unknown.end(); i++){
    if( ((*i)>=firstlocus) && ((*i)<lastlocus)){
      unknown.push_back((*i)-firstlocus);
    }
  }
  for(vector<int>::const_iterator i = C.known.begin(); i!=C.known.end(); i++){
    if( ((*i)>=firstlocus) && ((*i)<lastlocus)){
      known.push_back((*i)-firstlocus);
    }
  }
  for(vector<int>::const_iterator i = C.recom.begin(); i!=C.recom.end(); i++){
    if( ((*i)>=firstlocus) && ((*i)<lastlocus)){
      recom.push_back((*i)-firstlocus);
    }
  }


  int nloci = phase.size(); 
  for ( int locus = 0; locus < nloci; ++locus) {
    int oldloc = locus + firstlocus;
    phenotype[0].set_allele(locus,C.phenotype[0].get_allele(oldloc));
    phenotype[1].set_allele(locus,C.phenotype[1].get_allele(oldloc));
    orig_phenotype[0][locus] = C.orig_phenotype[0][oldloc];
    orig_phenotype[1][locus] = C.orig_phenotype[1][oldloc];
    saved_hap[0][locus] = C.saved_hap[0][oldloc];
    saved_hap[1][locus] = C.saved_hap[1][oldloc];
    Z[0][locus] = C.Z[0][oldloc];
    Z[1][locus] = C.Z[1][oldloc];
  }
    
  
  for(int c=0; c < 2; c++){
      for(int locus = 0; locus < nloci; locus++){
	AlleleCount[c][locus] = vector<float>(C.AlleleCount[c][locus+firstlocus]);
        //for(int j=0; j< n_allele; j++){
	//  AlleleCount[c][locus][j]=C.AlleleCount[c][locus+firstlocus][j];
	//}	
      }
  }
  
  for(int locus=0; locus<nloci; locus++){     
    for(int j=0; j < 2; j++)
      PhaseCount[locus][j]=C.PhaseCount[locus+firstlocus][j];	
  }

}

// construct individual from concatenation of two inds
CIndividual::CIndividual( const CIndividual & C1,  const CIndividual & C2 ) :
  id (C1.id),
  phase   ( C1.phase ),
  missing ( C1.missing ),
  notmissing_list ( C1.notmissing_list ),
  unknown ( C1.unknown ),
  known   ( C1.known ),
  recom (C1.recom),
  phenotype ( C1.phenotype ),
  orig_phenotype(2,vector<int>(C1.phase.size()+C2.phase.size())),
  saved_hap(2,vector<int>(C1.phase.size()+C2.phase.size())),
  Z(2,vector<int>(C1.phase.size()+C2.phase.size())),
  AlleleCount(2,vector< vector<float> > (C1.phase.size()+C2.phase.size(), vector<float> ())),
  PhaseCount(C1.phase.size()+C2.phase.size(), vector<float>(2, 1.0) )
{ 
  NIND ++;
  phase.insert(phase.end(),C2.phase.begin(),C2.phase.end());
  missing.insert(missing.end(),C2.missing.begin(),C2.missing.end()); 
  
  phenotype[0] = Haplotype(C1.phenotype[0],C2.phenotype[0]);
  phenotype[1] = Haplotype(C1.phenotype[1],C2.phenotype[1]);

  int offset = C1.phase.size();
 
  for(vector<int>::const_iterator i = C2.notmissing_list.begin(); i!=C2.notmissing_list.end(); i++)
    notmissing_list.push_back((*i)+offset);

  for(vector<int>::const_iterator i = C2.unknown.begin(); i!=C2.unknown.end(); i++)
    unknown.push_back((*i)+offset);
  
  for(vector<int>::const_iterator i = C2.known.begin(); i!=C2.known.end(); i++)
    known.push_back((*i)+offset);

  for(vector<int>::const_iterator i = C2.recom.begin(); i!=C2.recom.end(); i++)
    recom.push_back((*i)+offset);


  int nloci = phase.size();
  //copy C1 elements, and then C2 elements
  for ( int locus = 0; locus < offset; ++locus) {
    orig_phenotype[0][locus] = C1.orig_phenotype[0][locus];
    orig_phenotype[1][locus] = C1.orig_phenotype[1][locus];
    saved_hap[0][locus] = C1.saved_hap[0][locus];
    saved_hap[1][locus] = C1.saved_hap[1][locus];
    Z[0][locus] = C1.Z[0][locus];
    Z[1][locus] = C1.Z[1][locus];
  }
  for ( int locus = offset; locus < nloci; ++locus) {
    orig_phenotype[0][locus] = C2.orig_phenotype[0][locus-offset];
    orig_phenotype[1][locus] = C2.orig_phenotype[1][locus-offset];
    saved_hap[0][locus] = C2.saved_hap[0][locus-offset];
    saved_hap[1][locus] = C2.saved_hap[1][locus-offset];
    Z[0][locus] = C2.Z[0][locus-offset];
    Z[1][locus] = C2.Z[1][locus-offset];
  }
  
  
  for(int c=0; c < 2; c++){
    for(int locus = 0; locus < offset; locus++){
      //int n_allele = (loci_type[locus] == 'S' ) ? 2 : KMAX;
      AlleleCount[c][locus] = vector<float>(C1.AlleleCount[c][locus]);
      //for(int j=0; j< n_allele; j++){
      //	AlleleCount[c][locus][j]=C1.AlleleCount[c][locus][j];
      //}	
    }
    for(int locus = offset; locus < nloci; locus++){
      //int n_allele = (loci_type[locus] == 'S' ) ? 2 : KMAX;
      AlleleCount[c][locus] = vector<float>(C2.AlleleCount[c][locus-offset]);
      //for(int j=0; j< n_allele; j++){
      //	AlleleCount[c][locus][j]=C2.AlleleCount[c][locus-offset][j];
      //}	
    }
  }
  

  for(int locus=0; locus<offset; locus++){
    for(int j=0; j < 2; j++)
      PhaseCount[locus][j]=C1.PhaseCount[locus][j];	
  }
  for(int locus=offset; locus<nloci; locus++){
    for(int j=0; j < 2; j++)
      PhaseCount[locus][j]=C2.PhaseCount[locus-offset][j];	
  }
}


CIndividual::~CIndividual()
{ 
  NIND --;
#if DEBUG == 5
    cout << "Destructor called ... ";
    cout.flush();
#endif // DEBUG
 
    int nloci = phase.size();
    
    //cout << "DELETING IND nloci = " << nloci << "; NIND = " << NIND << endl;

#if DEBUG == 5
    cout << "OK" << endl;
#endif // DEBUG
}

const CIndividual & CIndividual::operator= (const CIndividual & rhs)
{
  //NIND ++;
  if(this != &rhs){
    this->phase =  ( rhs.phase );
    this->missing = ( rhs.missing );
    this->notmissing_list = ( rhs.notmissing_list );
    this->recom = (rhs.recom);
    this->unknown =( rhs.unknown );
    this->known   = ( rhs.known );   
    this->phenotype = (rhs.phenotype);
    this->id = (rhs.id);

    int nloci = rhs.phase.size();
    this->orig_phenotype = ( rhs.orig_phenotype );
    this->saved_hap = ( rhs.saved_hap );
    this->Z = ( rhs.Z );
    for ( int i = 0; i < nloci; ++i) {
        this->phenotype[0].set_allele(i,rhs.phenotype[0].get_allele(i));
        this->phenotype[1].set_allele(i,rhs.phenotype[1].get_allele(i));
    }

    this->AlleleCount =  (rhs.AlleleCount);
    this->PhaseCount= (rhs.PhaseCount);
    
  }

  return *this;

}
// Public functions

// This function sets the id to be 001, 002, ...
// when id is not present in the data set

void CIndividual::set_id ()
{
    char idstring[NMAXD + 3];
    ostrstream ostr ( idstring, NMAXD + 3 );
    ostr.fill('0');
    ostr.width(NMAXD);
    ostr << CIndividual::counts << ends;
    id = idstring;
}

void CIndividual::print_id ( ostream & ostr) const
{
    ostr << id;
}

// Return haplotype with phase c = 0 or 1.
Haplotype CIndividual::get_haplotype( int c ) const
{
  Haplotype haplo (phenotype[0].get_locus_type());
  for (int locus = 0; locus < phase.size(); ++locus) {
    haplo.set_allele(locus, get_haplotype(c, locus));
  }
  return haplo;
}
  
//  vector<int> CIndividual::get_haplotype ( int c ) const {
//      vector<int> haplo ( phase );
//      for (int locus = 0; locus < phase.size(); ++locus) {
//          haplo[locus] = get_haplotype(c, locus);
//      }
//      return haplo;
//  }


// Readin phenotypes
int CIndividual::read_orig_phenotypes ( istream & istr,
                                   const string & loci_type,
                                   int idpresent, 
                                   int format )
{
  ++CIndividual::counts;
    if ( idpresent ) {
        set_id ( istr );        // Read in id from file.
    } else {        
        set_id ();              // Assign id automatically
    }

    switch ( format ) {
    case 0:                  // 2 multilocus phenotypes on two lines
        for ( int chr = 0; chr < 2; ++chr ) {
            for ( int locus = 0; locus < loci_type.size(); ++locus) {
                input_orig_allele ( istr, loci_type[locus], chr, locus );
            }
        }
        break;
    case 1:                  // locus by locus phenotypes
        for ( int locus = 0; locus < loci_type.size(); ++locus ) {
            for ( int chr = 0; chr < 2; ++chr ) {
                input_orig_allele ( istr, loci_type[locus], chr, locus );	
            }
        }
        break;
    case 2: // locus by locus phenotypes, with one char for each, and with hets indicated by H
      for ( int locus = 0; locus < loci_type.size(); ++locus ) {	
	if(loci_type[locus] == 'M'){
	  cerr << "Error: format 2 not valid for data containing multiallelic markers" << endl;
	  exit(1);
	}
	input_orig_allele ( istr, loci_type[locus], 0, locus ); 
	orig_phenotype[1][locus] = orig_phenotype[0][locus];
        // set chromsome 1 to be the same as 0
	// note 'H's are replace by 01 in classpop.initialize
      }
      break;
      
    default:                 // Unknown format
      cerr << "Error: Unrecognized format requested" << endl;
      return 1;
      break;
    }
    
    return 0;

    
}

void CIndividual::input_orig_allele ( istream & istr, 
                                      char locus_type,
                                      int chr, int locus)
{
    switch ( locus_type ) {
    case 'S':
        char clt;
        // SNP loci
        istr >> clt;
        if ( clt == MISSNP ) {                 // Missing data
	  orig_phenotype[chr][locus] = MISSMS;
        } else {
	  orig_phenotype[chr][locus] = int (clt);
        }

	if( clt == '-'){
	  cerr << "Warning: use of - in input file at SNP locus may" << endl;
	  cerr << "indicate incorrect specification of missing allele" << endl;
	  cerr << "(Missing alleles at SNP locus should be specified as ?)." << endl;
	  cerr << "This warning may also occur if you use - to specify an indel allele," << endl;
	  cerr << "in which case you can ignore it." << endl;
	}
	
        break;
    case 'M':
        // Microsatellite loci
        istr >> orig_phenotype[chr][locus];
        break;
    default:
        cerr << "Error: Unrecognized locus type "
             << locus_type
             << "must be 'S' or 'M'. \n";
        exit(0);
    }
}

void CIndividual::print_phenotypes ( ostream & ostr,
                                     const string & loci_type,
                                     const vector<int> * coding ) const
{    
    for (int chr = 0; chr < 2; ++chr) {
        for (int locus = 0; locus < loci_type.size(); ++locus ) {
            int allele = get_allele ( chr, locus );       
            if ( loci_type[locus] == 'S' ) {
                assert ( allele == 0 || allele == 1 );
                allele = coding[allele][locus];
                ostr << (char) allele << ' ';
            } else {
                allele -= coding[0][locus];
                ostr << allele << ' ';
            }
        }
        ostr << endl;
    }
}

int CIndividual::BestAllele(int chr,int locus) const
{
  int best = 0;
  for(int j=1; j<AlleleCount[chr][locus].size(); j++){
    if(AlleleCount[chr][locus][j]>AlleleCount[chr][locus][best])
      best = j;
  }
  return best;
}

double CIndividual::BestAlleleProb(int chr,int locus) const
{
  double sum=0;
  for(int j=0; j<AlleleCount[chr][locus].size(); j++){
    sum += AlleleCount[chr][locus][j];
  }
  return AlleleCount[chr][locus][BestAllele(chr,locus)] / sum;
  
}


double CIndividual::BestPhaseProb(int locus) const
{
 if(PhaseCount[locus][0]>PhaseCount[locus][1])
   return PhaseCount[locus][0]/(PhaseCount[locus][0]+PhaseCount[locus][1]);
 else
   return PhaseCount[locus][1]/(PhaseCount[locus][0]+PhaseCount[locus][1]);
}

int CIndividual::BestPhase(int locus) const
{
 if(PhaseCount[locus][0]>=PhaseCount[locus][1])
   return 0;
 else
   return 1;
}

int CIndividual::BestHaplotype(int chr, int locus) const
{
  int bestphase = BestPhase(locus);
  if(bestphase == 0)
    return BestAllele(chr, locus);
  else
    return BestAllele(1-chr, locus);
}

Haplotype CIndividual::BestHaplotype(int chr) const
{
  Haplotype haplo(phenotype[0].get_locus_type());
  for (int locus = 0; locus < phase.size(); ++locus) {
    int bestphase = BestPhase(locus);
    if(bestphase == 0)
      haplo.set_allele(locus, BestAllele(chr, locus));
    else
      haplo.set_allele(locus, BestAllele(1-chr, locus));
  }
  return haplo;
}

void CIndividual::print_allele ( ostream & ostr, int locus, int chr,
				 const string & loci_type,
				 const vector<int> * coding,
				 bool PrintKnownPhase, bool PrintMissing, bool PrintBestGuess, double PhaseThreshold, double AlleleThreshold) const
{    
  double allele_prob,alternate_allele_prob; // stores probability of best allele
  bool UnknownPhase = std::binary_search ( unknown.begin(), 
					   unknown.end(), locus);
  bool IsHeterozygote;

  if ( PrintKnownPhase || UnknownPhase || missing[locus]>0) {                
                // print out the allele
    if((missing[locus]==0) || PrintMissing){
      int allele = get_haplotype ( chr, locus ); 
      int alternate_allele = get_haplotype( 1-chr, locus);

      if(PrintBestGuess){
	if(PhaseCount[locus][0]>=PhaseCount[locus][1])
	  {
	    allele = BestAllele(chr,locus);
	    alternate_allele = BestAllele(1-chr,locus);
	    allele_prob = BestAlleleProb(chr,locus);
	    alternate_allele_prob = BestAlleleProb(1-chr,locus);
	  }
	else
	  {
	    allele = BestAllele(1-chr,locus);
	    alternate_allele = BestAllele(chr,locus);
	    allele_prob = BestAlleleProb(1-chr,locus);
	    alternate_allele_prob = BestAlleleProb(chr,locus);	
	  }

	IsHeterozygote = (allele!=alternate_allele);
	//cout << ":IsHet" << IsHeterozygote << endl;
	
	if(allele_prob < AlleleThreshold)
	  ostr << '[';
	else if(alternate_allele_prob < AlleleThreshold)
	  ostr << ' ';

	// test if phase unknown, and exceeds threshold, and heterozygote
	if((BestPhaseProb(locus)<PhaseThreshold) && UnknownPhase && (IsHeterozygote))
	  ostr << '(';
       
      }

      if ( loci_type[locus] == 'S' ) {
	assert ( allele == 0 || allele == 1 );
	allele = coding[allele][locus];
	ostr << (char) allele;
      } else {
	allele -= coding[0][locus];
	ostr << allele       ;
      }

      if(PrintBestGuess && (BestPhaseProb(locus)<PhaseThreshold) && UnknownPhase && IsHeterozygote )
	  ostr << ')';
      if(PrintBestGuess)
	{
	  if(allele_prob < AlleleThreshold)
	    ostr << ']';
	  else if (alternate_allele_prob < AlleleThreshold)
	    ostr << ' ';
	}
      ostr << ' ';
    }
    else {
      ostr << MISSCHAR << ' ';
    }
  }
  else {
    ostr << SPACEHOLDER << ' ';
  } 
  ostr.flush();
}



void CIndividual::print_haplotype ( int chr, ostream & ostr,
                                     const string & loci_type,
                                     const vector<int> * coding,
                                     bool PrintKnownPhase, bool PrintMissing, bool PrintBestGuess, double PhaseThreshold, double AlleleThreshold ) const
{    
  
  for (int locus = 0; locus < loci_type.size(); ++locus ) 
    print_allele(ostr, locus, chr, loci_type, coding, PrintKnownPhase, PrintMissing, PrintBestGuess, PhaseThreshold, AlleleThreshold);
}

void CIndividual::print_haplotypes ( ostream & ostr,
                                     const string & loci_type,
                                     const vector<int> * coding,
                                     bool PrintKnownPhase, bool PrintMissing, bool PrintBestGuess, double PhaseThreshold, double AlleleThreshold  ) const
{    
    for (int chr = 0; chr < 2; ++chr)
      {
	print_haplotype(chr, ostr, loci_type, coding, PrintKnownPhase, PrintMissing, PrintBestGuess, PhaseThreshold, AlleleThreshold);
	ostr << endl;
      }
}

void CIndividual::print_phase ( ostream & ostr,
                                bool bfPrintAll ) const 
{
    for (int r = 0; r < phase.size(); ++r) {
        if ( bfPrintAll ||
	     ( is_unknown(r) && n_missing(r) < 2 ) ) {
	  ostr << phase[r];
        } else {
	  ostr << SPACEHOLDER;
        }
    }
}

void CIndividual::print_phase_prob ( ostream & ostr,
                                bool bfPrintAll ) const 
{
    for (int r = 0; r < phase.size(); ++r) {
        if ( bfPrintAll ||
                ( is_unknown(r)
                  && n_missing(r) < 2 ) ) {
	  ostr << setprecision(2) << BestPhaseProb(r) << ' ';
        } else {
	  if(n_missing(r) == 0)
	    ostr << SPACEHOLDER << ' ';
	  else
	    ostr << "? ";
        }
    }
    ostr << endl;
}

int CIndividual::initialize ( int knowninfo, istream & istr_known, int initmethod, istream & istr_init, const string & loci_type)
{
    char tempchar;
    int temp;

    for (int locus = 0; locus < phase.size(); ++locus) {
        if ( orig_phenotype[0][locus] == MISSMS &&
             orig_phenotype[1][locus] == MISSMS ) {
            // Missing two alleles
            missing[locus] = 2;
            // Also needs update
            unknown.push_back(locus);
	} else if ( orig_phenotype[0][locus] == MISSMS || 
                    orig_phenotype[1][locus] == MISSMS ) {
            // Missing one alleles
	  cerr << "Warning: there is a genotype missing one allele at locus " << (locus+1) << endl;
	  missing[locus] = 1;
            // Swap so that chr 1 is missing
	    // Note: Need to be careful when phase is specified
	    //if ( phenotype[0].get_allele(locus) == MISSMS ) {
	    //phenotype[0].set_allele(locus,phenotype[1].get_allele(locus));
	    // phenotype[1].set_allele(locus,MISSMS);
	    // orig_phenotype[0][locus] = orig_phenotype[1][locus];
	    // orig_phenotype[1][locus] = MISSMS;
            //} 
            // Phase unknown
	  unknown.push_back(locus);
	    
        } else {
	  
	  notmissing_list.push_back(locus);
	  
	  if ( phenotype[0].get_allele(locus) != phenotype[1].get_allele(locus) ) {
            // Phase unknown
            unknown.push_back ( locus );
	  } else {
            // phase known loci
            known.push_back (locus );
	  }
	}
    }
 
    //if ( unknown.size() == 1 ) {
      // Take out individual with one unknown phase
    //  unknown.clear();
    //}

    // If known info is supplied, use it
    // Data in "knownfile", with UNKNOWNPHASECHAR to
    // indicate unknown phase, 0 to indicate
    // phase is as input and 1 to indicate opposite
                     
    if(knowninfo > 0 ){
      if((knowninfo == 1) & !istr_known ) {
        cerr << "Need to specify file of known phases!" << endl;
        exit(1);
      }
      if(knowninfo == 999)
	cerr << "Warning: assuming all phases as input in file" << endl;

      unknown.clear();
      known.clear();
      
      for (int locus = 0; locus < phase.size(); ++locus) {
	if(knowninfo == 999)
	  tempchar = '0';
	else
	  istr_known >> tempchar;
	double r = ranf();
	//if(missing[locus]>0)
	//  tempchar = UNKNOWNPHASECHAR;

	switch(tempchar){
	case '0':
	  phase[locus]=0;
	  known.push_back(locus);
	  break;
	case '1':
	  phase[locus]=0;
	  temp = orig_phenotype[0][locus];
	  orig_phenotype[0][locus] = orig_phenotype[1][locus];
	  orig_phenotype[1][locus] = temp;
	  phenotype[0].set_allele(locus, orig_phenotype[0][locus]);
	  phenotype[1].set_allele(locus, orig_phenotype[1][locus]);
	  known.push_back(locus);
	  break;
	case UNKNOWNPHASECHAR:	         
	  if((orig_phenotype[0][locus] != orig_phenotype[1][locus]) || (missing[locus]>0)){
	    unknown.push_back(locus);	      
	    if ( r < 0.5 ) {
	      phase[locus] = 0;
	    } else {
	      phase[locus] = 1;
	    }
	  }
	  else{
	    phase[locus]=0;
	    known.push_back(locus);
	  }
	  break;
	default:
	  cerr << "Error: illegal character ('" << tempchar << "') in file for start info" << endl;
	  exit(1);
	} 
           
      }
    }

    // this actually in just to debug: shouldn't happen!
    // cerr << "Known info: " << endl;
//     for (vector<int>::const_iterator i = known.begin();
// 	   i != known.end(); ++i) {
      
//       cerr << *i << ":" << missing[*i] << endl;

//       if(missing[*i]>0){
// 	cerr << "Error: missing an allele at a position indicated as KNOWN; individual" << id << " position " << *i << endl;
// 	exit(1);
//       }
//     }

    // initialise phase

    switch ( initmethod ) {
    case 0: // randomly initialise phase
      for (vector<int>::const_iterator i = unknown.begin();
	   i != unknown.end(); ++i) {
	double r = ranf();
	if ( r < 0.5 ) {
	  phase[*i] = 0;
	} else {
	  phase[*i] = 1;
	}
      }
      break;

    case 1: // Use phase in inputfile as starting point
      for (int locus = 0; locus < phase.size(); ++locus) {
	istr_init >> tempchar;
	phase[locus] = tempchar - (int) '0';
	if(phase[locus]!=0 && phase[locus]!=1){
	  cerr << "Error: in initial phase file" << endl;
	  cerr << "Phase must be 0 or 1, but read in as " << tempchar << endl;
	  exit(1);
	}
      }
      break;

    case 2: // set phase to be 0 as starting point
      for (vector<int>::const_iterator i = unknown.begin();
	   i != unknown.end(); ++i) {
	//cout << *i << phase[*i] << endl;
	phase[*i] = 0;
      }
      break;   
    }
   
    return 0;
}


double CIndividual::flipprob()
{
  double fp=0;
  double nfp=0;
  int nloci = phase.size();
  
  for(int r=0; r< nloci; r++){
    int allele0 = get_haplotype(0,r);
    int allele1 = get_haplotype(1,r);
    double PSEUDOCOUNT = 1.0;
    nfp += log( (AlleleCount[0][r][ allele0 ]+PSEUDOCOUNT) * (AlleleCount[1][r][ allele1 ] + PSEUDOCOUNT) * PhaseCount[r][0] + (AlleleCount[0][r][ allele1 ]  + PSEUDOCOUNT)* (AlleleCount[1][r][ allele0 ] + PSEUDOCOUNT) * PhaseCount[r][1]);
    fp += log((AlleleCount[0][r][ allele1 ] + PSEUDOCOUNT) * (AlleleCount[1][r][ allele0 ]  + PSEUDOCOUNT) * PhaseCount[r][0] + (AlleleCount[0][r][ allele0 ] + PSEUDOCOUNT) * (AlleleCount[1][r][ allele1 ] + PSEUDOCOUNT) * PhaseCount[r][1] );
  }

  return 1/(1+exp(nfp-fp));
}

void CIndividual::ResetCounts()
{
  int nloci = get_nloci();
  for(int c=0; c < 2; c++){
    for(int r=0; r<nloci; r++){
      for(int j=0; j< AlleleCount[c][r].size(); j++){
	AlleleCount[c][r][j]=0.0;
      }	
    }
  }
  
  for(int r=0; r<nloci; r++){
    for(int j=0; j < 2; j++)
      PhaseCount[r][j]=1.0;	
  }

}


void CIndividual::UpdateCounts()
{
  double fp = flipprob();

  //*NCount += 1;
  int nloci = get_nloci();
  for(int r=0; r < nloci ; r++){
    int allele0 = get_haplotype(0,r);
    int allele1 = get_haplotype(1,r);
  
    if(allele0 == allele1){ // for homozygotes, just update allele counts
      AlleleCount[0][r][allele0] += 1.0;
      AlleleCount[1][r][allele1] += 1.0;
    }

    else{ // for heterozygotes
      double PSEUDOCOUNT = 1.0;
      // following test if allele0 is identified with bestallele0 and 
      // allele1 is identified with bestallele1,
      if( (AlleleCount[0][r][allele0] + PSEUDOCOUNT) * (AlleleCount[1][r][allele1] + PSEUDOCOUNT)
	  >= (AlleleCount[1][r][allele0] + PSEUDOCOUNT) * (AlleleCount[0][r][allele1] + PSEUDOCOUNT) ){
	AlleleCount[0][r][allele0] += 1.0;
	AlleleCount[1][r][allele1] += 1.0;
	
	if(fp<=0.5)
	  PhaseCount[r][0] += 1.0;
	else
	  PhaseCount[r][1] += 1.0;
      }
      else{
	AlleleCount[0][r][allele1] += 1.0;
	AlleleCount[1][r][allele0] += 1.0;
	if(fp<=0.5)
	  PhaseCount[r][1] += 1.0;
	else
	  PhaseCount[r][0] += 1.0;
      }
    }   
  }
}

double CIndividual::ObservedDataProbGivenParents(const Haplotype & h1, const Haplotype & h2, const Haplotype & h3, const Haplotype & h4, const vector<int> & recom1, const vector<int> & recom2){
  int chr1 = 0;
  int chr2 = 0;
  int nloci = get_nloci();
  double prob = 1;

  for(int r = 0; r < nloci; r++){
    int allele0, allele1;
    if(chr1 == 0)
      allele0 = h1.get_allele(r);
    else
      allele0 = h2.get_allele(r);
    if(chr2 == 0)
      allele1 = h3.get_allele(r);
    else
      allele1 = h4.get_allele(r);

    if(missing[r] == 0){
      if(!( 
	   ( ( (allele0 == orig_phenotype[0][r]) || (orig_phenotype[0][r] == MISSMS) )  
	    && ( (allele1 == orig_phenotype[1][r]) || (orig_phenotype[1][r] == MISSMS) ) )
	   || 
	   ( ( (allele0 == orig_phenotype[1][r]) || (orig_phenotype[1][r] == MISSMS) ) 
	     && ( (allele1 == orig_phenotype[0][r]) || (orig_phenotype[0][r] == MISSMS) ))
	     )){
	prob = 0;
	break;
      }
    }
    if(std::binary_search(recom1.begin(),recom1.end(),r))
      chr1 = 1-chr1;
    if(std::binary_search(recom2.begin(),recom2.end(),r))
      chr2 = 1-chr2;
  }
  return prob;
}


// transfer counts (allelecounts, phasecounts) from summary to individual
void CIndividual::TransferCounts(Summary & sum)
{
  int locus = 0;
  for (vector< vector<double> >::iterator p1 = sum.flipprob.begin(); p1 !=sum.flipprob.end(); p1++){
    for(vector<double>::iterator p2 = (*p1).begin(); p2 != (*p1).end(); p2++){
      PhaseCount[locus][0] = 1-(*p2);
      PhaseCount[locus++][1] = (*p2);
    }
  }
  
  locus = 0;
  for (int segment = 0; segment < sum.alleleprob.size(); segment++){
    for (int pos = 0; pos < sum.alleleprob[segment].size(); pos++){
      for(int allele = 0; allele < AlleleCount[0][locus].size(); allele++){
	if(allele != sum.bestguess[segment].first.get_allele(pos))
	  AlleleCount[0][locus][allele] = sum.errorprob[segment][pos][0] * sum.alleleprob[segment][pos][0][allele];
	else
	  AlleleCount[0][locus][allele] = 1 - sum.errorprob[segment][pos][0];
	if(allele != sum.bestguess[segment].second.get_allele(pos))
	  AlleleCount[1][locus][allele] = sum.errorprob[segment][pos][1] * sum.alleleprob[segment][pos][1][allele];
	else
	  AlleleCount[1][locus][allele] = 1 - sum.errorprob[segment][pos][1];
      }
      locus++;
    }
  }
}
  
// returns whether it is possible for h to exist
// in the population pop
bool CanBeFoundAtAll(const Haplotype & h, const vector<CIndividual> & pop)
{
  bool found;
  for(int ind = 0; ind < pop.size(); ind++){
    Haplotype CompHap = GetCompHap(h, pop[ind], found, false);
    if(found)
      return true;
  }
  return false;
}

// calls routine below when
// you aren't bothered about returning which positions
// match
Haplotype GetCompHap(const Haplotype & h, const CIndividual & ind, bool & found,bool checkmissing ){
  vector<int> dummy;
  return GetCompHap(h,ind,found,dummy,checkmissing);
}





// returns complementary haplotype to h in ind
// (returns found = false if couldn't find h in ind)
//
// if checkmissing is false, h is looked for
// only at the positions where there are no missing alleles
// In this case. if there's a match at the non-missing positions,
// found is returned as true, and CompHap at the missing
// positions is as close as possible to the comp of h
// matches returns which chrom (0 or 1) ind matches h at each
// unknown position
Haplotype GetCompHap(const Haplotype & h, const CIndividual & ind, bool & found,vector<int> & matches, bool checkmissing )
{
  Haplotype Comp = ind.get_haplotype(1);
  found=false;  
  matches = vector<int>(ind.get_nloci(),0);

  const vector<int> unknown_list =  ind.get_unknown_pos();
  for (vector<int>::const_iterator u = unknown_list.begin();
       u != unknown_list.end(); ++u ) {   
    int hapallele = h.get_allele(*u);
    int allele0 = ind.get_haplotype(0,*u);
    if(allele0!=hapallele){
      if( (ind.get_haplotype(1,*u)==hapallele) // if h is in chrom 1 in ind,
	  || ((checkmissing==false) && (ind.n_missing(*u)>0)) ) {
	                                   //or this is a missing pos, and we
	                            // are not checking the missing positions
	Comp.set_allele(*u,allele0); //  set Comp to be chrom 0 in ind
	matches[*u] = 1;
      }
      else
	return Comp; // couldn't find it
    }
  }
  
  int found0=1; // found0 indicates whether the inds hap0 matches
  int found1=1; // found1 indicates whether the inds hap1 matches
  const vector<int> known_list =  ind.get_known_pos();
  for (vector<int>::const_iterator u = known_list.begin();
       u != known_list.end(); ++u ) {
    int hapallele = h.get_allele(*u);
   
    if(ind.get_haplotype(0,*u) != hapallele)
      found0=0;
    if(ind.get_haplotype(1,*u) != hapallele)
      found1=0;
    
    if((found0 == 0) && (found1 ==0))
      return Comp; // couldn't find it
  }
    
  
  if(found0!=1){
    if(found1==1){
      const vector<int> known_list =  ind.get_known_pos();
      for (vector<int>::const_iterator u = known_list.begin();
	   u != known_list.end(); ++u ){
	Comp.set_allele(*u,ind.get_haplotype(0,*u));
      }
    }
    else
      return Comp; // found = false
  }

  found = true;
  return Comp;
    
}


int NDiff(const std::vector<CIndividual> & pop, int n0, int c0, int n1, int c1, const std::vector<int> & uselist)
{
  int nd=0;
  for (vector<int>::const_iterator u = uselist.begin();
          u != uselist.end(); ++u ) {
        // compute DiffCount based only on loci in uselist	             
    nd += (pop[n0].get_haplotype(c0, *u) != pop[n1].get_haplotype(c1, *u));
  }
  return nd;

}


// compute the probability of seeing n0,c0 
// given it copied n1,c1, separated by time t, 
// for loci in uselist
double CCProb(const std::vector<CIndividual> & pop, int n0, int c0, int n1, int c1, int t, int nchr, const vector<ArrayQ *>  & Q, const std::vector<int> & uselist)
{

  int from,targ;
  double prod=1;

  for ( vector<int>::const_iterator u = uselist.begin();
	  u != uselist.end(); ++u ) {
      // compute based only on loci in uselist	
    from = pop[n0].get_haplotype (c0, *u);
    targ = pop[n1].get_haplotype (c1, *u);       
    prod *= PrHitTarg(*u, nchr , t, from, targ, Q);
  }  
  return prod;

}

// return subvector whose values fall within min and max-1
std::vector<int> subrange(const vector<int> & vec,int min, int max)
{ 
  return vector<int>(0);
  // return vector(start,finish);
}


//

// {{{ Log
/*
  $Log: indnode.cpp,v $
  Revision 1.33  2003/06/14 00:24:05  stephens
  Adding files, and committing a lot of changes
  that have resulted in version 2.0 of PHASE

  Revision 1.32  2002/02/27 18:56:45  stephens
  Commiting the current source, which is essentially that released as
  PHASE 1.0. The number of quadrature points is reduced to 2. Some code has
  been added to cope with recombination, but this is not included in the release
  and is still under development.

  Revision 1.31  2001/10/21 18:33:29  stephens
  fixed bug in GibbsUpdate

  Revision 1.29  2001/10/16 20:44:16  stephens
  Added way to input position of each marker, by preceding line
  of SSMMMSS with a line beginning with P followed by position of
  each marker.

  Revision 1.28  2001/10/16 18:12:06  stephens
  More minor bugs fixed. (including -m bug)

  Revision 1.26  2001/10/06 00:35:40  stephens
  Impute missing positions added- changed update_allele to update_haplotype.

  Revision 1.25  2001/10/02 18:02:59  stephens
  just some small additions. haven't committed these for a long time
  because cvs was not working propersly. (now know why: they remounted
  you from /user3 to /home)

  Revision 1.24  2001/06/19 17:02:25  stephens
  Changes to computation of arrayFF to make more efficient
  Added facility to store "original phenotype" in indnode,
  in preparation for allowing genotyping error.

  Revision 1.23  2001/05/21 20:16:34  nali

  Add another member vector<int> known for phase known loci, and
  corresponding accessor and initialization.

  Revision 1.22  2001/04/26 17:09:17  stephens
  Added function ClassPop::print_allele
  and made a few minor changes

  Revision 1.21  2001/04/24 22:00:08  stephens
  Added flag in ClassPop::output_hap and output_phase
  to indicate whether to output all positions, or only
  unknown positions.

  Revision 1.20  2001/04/24 19:34:58  nali
  Merge several print_haplotypes member functions into one with an option of
  printing all alleles or not (default not).
  Add default argument 0 to the constructor so default constructor can be called.

  Revision 1.19  2001/04/20 00:32:45  nali
  Put reading loci types into the contructor of ClassPop

  Revision 1.18  2001/04/19 19:47:35  nali
  Output SPACEHOLDER when there is no ambiguity in haplotype or phase.

  Revision 1.17  2001/04/17 22:08:44  nali

  Major revisement in overall structure. Created new class ClassPop and
  almost all global functions now became member functions of ClassPop, most
  of them private.

  "mult" removed in update_phase_NR. No other changes in terms of algorithm.
  Haven't check the results yet.

  proc_args() is moved to utility.cpp, which also defines a couple of other
  global functions.

  Revision 1.16  2001/04/12 17:57:41  stephens
  Fixed bug in coding for SNPs (azero and aone are now defined
  inside the locus loop, rather than out of it.)
  This seems to fix the bug in previous version (output haplotypes are now
  consistent with input genotypes)

  Revision 1.15  2001/04/12 17:13:10  stephens
  edited haplotype output routines, and added line to output
  haplotypes to output file. Output from test files indicates
  a bug as output haplotypes are inconsistent with input genotypes

  Revision 1.14  2001/04/12 01:53:25  stephens
  Added -D option for reading in several datasets
  Debugged -H option for inputting data from Hudson's simulations
  Reduced Nthin to 1 to reduce time for examples, and
  added a test example for the -H option and -D options

  Revision 1.13  2001/04/09 16:28:35  nali
  Most part of original ResolvePhase (without recombination) implemented.
  Now compiles and runs.

  Revision 1.12  2001/04/07 07:18:10  nali
  Revised to concentrate on implementing SD method.
  Discarded support for EM algorithm, etc..

  Revision 1.11  2001/02/28 04:54:31  nali
  Make EM working and haplotype list right

  Revision 1.10  2001/02/27 08:16:49  nali
  Make use of the new haplotype list class.

  Revision 1.9  2001/02/21 18:36:23  nali
  indnode.cpp

  Revision 1.8  2001/02/16 17:13:17  nali
  New functions: InputRandom, InputHusdonData, etc..
  New member function: make_haplist.

  Revision 1.7  2001/02/16 02:41:48  nali

  Remove the character representation of phenotypes.

  Revision 1.6  2001/02/14 23:40:40  nali

  Make some functions inline and no range checking.

  Revision 1.5  2001/02/14 02:00:21  nali

  Want get_phenotype() and get_phenotype_ch() to be public member functions
  eventually.

  Revision 1.4  2001/02/14 01:57:07  nali

  Add member function to print phenotype (mainly for debugging).

*/
// }}}
