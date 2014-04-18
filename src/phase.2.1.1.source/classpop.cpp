/*
 * Implementation of ClassPop
 * $Id: classpop.cpp,v 1.30 2003/06/14 00:24:04 stephens Exp $
 */

#include "classpop.hpp"
#include "constants.hpp"
#include "utility.hpp"
#include "arrayDiffProb.hpp"
#include "arrayDiffCount.hpp"
#include "arrayCC.hpp"
#include "indnode.hpp"
#include "HapList2.hpp"
#include "arrayQ.hpp"
#include "errcheck.hpp"

#include <iomanip>
#include <fstream>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <cmath>

using namespace std;

vector<ArrayQ> ClassPop::Qvalues;
ArrayQ ClassPop::QSNP;

// {{{ Constructor and destructors

ClassPop::ClassPop () :
    pop (0),
    Nind(0),
    Nchild(0),
    childindex (0),
    casecontrol (0),
    groupsize (0),
    buddy (0),
    BestLogProb (0),
    haplist (),
    loci_type (""),
    position (0),
    order (0),
    vecRho(0),
    vecRhoDeriv(0),
    RhoMean(1),
    RhoMult(1),
    right(0),
    left(0),
    lambda(0),
    RecomModel(0),
    NindToUseForRho(0),
    CurrentLogProb(0),
    vecTheta(0),
    NSNP(0),
    ALLSNP(0),
    TREATSNPSASMS(0),
    Qptr (0),
    vecdiffprob (0),
    CC (0),
    //knownCC (0),
    DiffCount(0)
    //knownDiffCount(0)
{
  //coding = new vector<int> [2];
    coding[0] = vector<int> (0);
    coding[1] = vector<int> (0);
    SNPlist = vector<int> (0);
    nonSNPlist = vector<int> (0);
}

ClassPop::ClassPop ( const ClassPop & cpop2 ) :
    casecontrol( cpop2.casecontrol),
    groupsize( cpop2.groupsize),
    Nloci( cpop2.Nloci),
    Nind ( cpop2.Nind),
    Nchild (cpop2.Nchild),
    loci_type ( cpop2.loci_type ),
    pop(cpop2.pop),
    childindex(cpop2.childindex),
    buddy(cpop2.buddy),
    //pop (cpop2.Nind, CIndividual ( cpop2.loci_type )),
    position ( cpop2.position ),
    order ( cpop2.order ),
    vecRho ( cpop2.vecRho ),
    vecRhoDeriv ( cpop2.vecRhoDeriv ),
    RhoMean ( cpop2.RhoMean),
    RhoMult ( cpop2.RhoMult),
    right (cpop2.right),
    left (cpop2.left),
    lambda (cpop2.lambda),
    RecomModel (cpop2.RecomModel),
    NindToUseForRho(0),
    CurrentLogProb (cpop2.CurrentLogProb),
    vecTheta ( cpop2.vecTheta ),
    BestLogProb ( cpop2.BestLogProb ),
    NSNP (cpop2.NSNP),
    ALLSNP (cpop2.ALLSNP),
    TREATSNPSASMS (cpop2.TREATSNPSASMS),
    SNPlist (cpop2.SNPlist),
    nonSNPlist (cpop2.nonSNPlist),
    Qptr (cpop2.Qptr),
    vecdiffprob (cpop2.vecdiffprob),
    CC (0), //cpop2.Nind),
    DiffCount (0), //(cpop2.Nind),
    haplist ( cpop2.haplist )    
{    
  for(int ind = 0; ind < cpop2.pop.size(); ind++)
    pop[ind] = CIndividual(cpop2.pop[ind]);
  //coding = new vector<int> [2];
  coding[0] = vector<int>(cpop2.coding[0]);
  coding[1] = vector<int>(cpop2.coding[1]);
}

ClassPop::ClassPop ( const ClassPop & cpop2, int firstlocus, int lastlocus ) : 
  casecontrol(cpop2.casecontrol),
  groupsize(cpop2.groupsize),
  Nloci(lastlocus-firstlocus),
  Nind(cpop2.Nind),
  Nchild(cpop2.Nchild),
  loci_type(cpop2.loci_type.begin()+firstlocus,cpop2.loci_type.begin()+lastlocus),
  //pop(cpop2.Nind, CIndividual ( cpop2.loci_type )),
  pop(cpop2.Nind+cpop2.Nchild),
  childindex(cpop2.childindex),
  buddy(cpop2.buddy),
  position(cpop2.position.begin()+firstlocus,cpop2.position.begin()+lastlocus),
  order(cpop2.order),
  vecRho(cpop2.vecRho.begin()+firstlocus,cpop2.vecRho.begin()+lastlocus),
  vecRhoDeriv(cpop2.vecRhoDeriv.begin()+firstlocus,cpop2.vecRhoDeriv.begin()+lastlocus),
  RhoMean(cpop2.RhoMean),
  RhoMult(cpop2.RhoMult.begin()+firstlocus,cpop2.RhoMult.begin()+lastlocus),
  right (cpop2.right),
  left (cpop2.left),
  lambda (cpop2.lambda),
  RecomModel (0), // (cpop2.RecomModel),
  NindToUseForRho (cpop2.NindToUseForRho),
  CurrentLogProb(0),
  vecTheta(cpop2.vecTheta.begin()+firstlocus,cpop2.vecTheta.begin()+lastlocus),
  BestLogProb ( cpop2.BestLogProb ),
  haplist ( cpop2.haplist, firstlocus, lastlocus ),
  TREATSNPSASMS (cpop2.TREATSNPSASMS),
  Qptr (cpop2.Qptr.begin()+firstlocus,cpop2.Qptr.begin()+lastlocus),
  vecdiffprob (lastlocus - firstlocus + 1),
  CC (0), //cpop2.pop.size()),
  DiffCount (0), //(cpop2.pop.size()),
  NSNP(0),
  SNPlist(0),
  nonSNPlist(0),
  ALLSNP(0)
{
  for(int ind = 0; ind < cpop2.pop.size(); ind++)
    pop[ind] = CIndividual(cpop2.pop[ind],firstlocus,lastlocus);

  NSNP = 0;
    //coding = new vector<int> [2];

  coding[0] = vector<int>(cpop2.coding[0].begin()+firstlocus, cpop2.coding[0].begin()+lastlocus);
  coding[1] = vector<int>(cpop2.coding[1].begin()+firstlocus, cpop2.coding[1].begin()+lastlocus); 
  
 //   if(TREATSNPSASMS)
//      NSNP = 0;
//    else
//      for( int r = 0; r < Nloci; r++){
//        if(loci_type[r]=='S'){
//  	NSNP++;
//        }
//      }    

//    SNPlist = vector<int> (NSNP);
//    nonSNPlist = vector<int> (Nloci - NSNP);

//    NSNP = 0;
//    int nonSNP =0;
//    for( int r = 0; r < Nloci; r++){
//      if(TREATSNPSASMS)
//        nonSNPlist[nonSNP++] = r;
//      else
//        if(loci_type[r]=='S'){
//  	SNPlist[NSNP++] = r;
//        }
//        else{
//  	nonSNPlist[nonSNP++] = r;
//        }
//    }    
  
  for( int r = 0; r < Nloci; r++){
    if(TREATSNPSASMS)
      nonSNPlist.push_back(r);
    else
      if(loci_type[r]=='S'){
	SNPlist.push_back(r);
	NSNP++;
      }
      else{
	nonSNPlist.push_back(r);
      }
  }    
 
  if(Nloci == NSNP)
    ALLSNP=1;
  else
    ALLSNP=0;

  ComputeRho();
  ComputeRhoDerivAndCurrentLogProb();
  
  Computevecdiffprob();

}

// concatenates two classpops into a single one
// (only adds haplotypes to the haplist if the min of their
// frequencies exceeds minfreq)
ClassPop::ClassPop ( const ClassPop & CP1, const ClassPop & CP2, double minfreq ) :
    casecontrol (CP1.casecontrol ),
    groupsize(CP1.groupsize ),
    Nloci( CP1.Nloci+CP2.Nloci),
    Nind ( CP1.Nind),
    Nchild (CP1.Nchild),
    loci_type ( CP1.loci_type + CP2.loci_type ),
    pop( CP1.Nind + CP1.Nchild),
    childindex( CP1.childindex),
    buddy( CP1.buddy ),
    position ( CP1.position ),
    order ( CP1.order ),
    vecRho ( CP1.vecRho ),
    vecRhoDeriv (CP1.vecRhoDeriv ),
    RhoMean ((CP1.RhoMean * CP1.get_physical_length() + CP2.RhoMean * CP2.get_physical_length())/( CP1.get_physical_length() + CP2.get_physical_length() ) ),
    RhoMult ( CP1.RhoMult ),
    right (CP1.right),
    left (CP1.left),
    lambda (CP1.lambda),
    RecomModel (CP1.RecomModel),
    NindToUseForRho (CP1.NindToUseForRho),
    CurrentLogProb ( 0 ),
    vecTheta ( CP1.vecTheta ),
    BestLogProb ( BIGNEGATIVE ),
    NSNP (CP1.NSNP+CP2.NSNP),
    ALLSNP (CP1.ALLSNP*CP2.ALLSNP),
    TREATSNPSASMS (CP1.TREATSNPSASMS * CP2.TREATSNPSASMS),
    SNPlist (CP1.SNPlist),
    nonSNPlist (CP1.nonSNPlist),
    Qptr (CP1.Qptr),
    vecdiffprob ( CP1.loci_type.size() + CP2.loci_type.size() + 1),
    CC (0),//CP1.Nind),
    DiffCount (0) //CP1.Nind)
    //    haplist ( CP1.haplist, CP2.haplist, minfreq )
{    
  if(CP1.Nind != CP2.Nind){
    cerr << "Error: trying to concatenate classpops of different size!" << endl;
    exit(1);
  }	 

  for(int ind = 0; ind < pop.size(); ind++)
    pop[ind] = CIndividual(CP1.pop[ind],CP2.pop[ind]);

  position.insert(position.end(), CP2.position.begin(), CP2.position.end());
  Qptr.insert(Qptr.end(), CP2.Qptr.begin(), CP2.Qptr.end() );
  
  vecRho.insert(vecRho.end(), CP2.vecRho.begin(), CP2.vecRho.end());
  vecRhoDeriv.insert(vecRhoDeriv.end(), CP2.vecRhoDeriv.begin(), CP2.vecRhoDeriv.end());
  RhoMult.insert(RhoMult.end(), CP2.RhoMult.begin(), CP2.RhoMult.end());

  if(RecomModel == 0){
    for(int i = 0; i<CP1.RhoMult.size(); i++){
      if(RhoMean > 0)
	RhoMult[i] = CP1.RhoMean * CP1.RhoMult[i] / RhoMean;
      else
	RhoMult[i] = 1.0;
      if(RhoMult[i]<0.01)
	RhoMult[i] = 0.01;
    }
    for(int i = 0; i<CP2.RhoMult.size(); i++){
      if(RhoMean >0)
	RhoMult[i+CP1.get_nloci()] = CP2.RhoMean * CP2.RhoMult[i] / RhoMean;
      else
	RhoMult[i+CP1.get_nloci()] = 1.0;
      if(RhoMult[i+CP1.get_nloci()]<0.01)
	RhoMult[i+CP1.get_nloci()] = 0.01;
    }
    RhoMult[CP1.get_nloci()-1] = 1.0; // set to 1 the lambda for interval between blocks
  }

    if(RecomModel == 1 || RecomModel == 2){
    InitialiseSimpleHotspot(lambda.size(),RecomModel==2);
  }
  
  vecTheta.insert(vecTheta.end(), CP2.vecTheta.begin(), CP2.vecTheta.end());
   
  haplist = HapList(CP1.haplist, CP2.haplist, pop, minfreq);

  //coding = new vector<int> [2];
  coding[0] = vector<int>(CP1.coding[0]);
  coding[1] = vector<int>(CP1.coding[1]);
  coding[0].insert(coding[0].end(),CP2.coding[0].begin(),CP2.coding[0].end());
  coding[1].insert(coding[1].end(),CP2.coding[1].begin(),CP2.coding[1].end());

  ComputeRho();
  ComputeRhoDerivAndCurrentLogProb();
 
  
  Computevecdiffprob();

}

// concatenates two classpops into a single one
// (takes current states of each pop, and does a few quick iterations to
// decide how to line them up)
ClassPop::ClassPop ( const ClassPop & CP1, const ClassPop & CP2, map<string, int> & cmds, int Niter, double theta, vector<double> & vecDelta, map<string, double> & d_cmds, double minfreq ) :
    casecontrol (CP1.casecontrol),
    groupsize (CP1.groupsize),
    Nloci( CP1.Nloci+CP2.Nloci),
    Nind ( CP1.Nind),
    Nchild (CP1.Nchild),
    loci_type ( CP1.loci_type + CP2.loci_type ),
    pop( CP1.Nind  + CP1.Nchild),
    childindex (CP1.childindex),
    buddy ( CP1.buddy ),
    position ( CP1.position ),
    order ( CP1.order ),
    vecRho ( CP1.vecRho ),
    vecRhoDeriv ( CP1.vecRhoDeriv ),
    RhoMean ((CP1.RhoMean * CP1.get_physical_length() + CP2.RhoMean * CP2.get_physical_length())/( CP1.get_physical_length() + CP2.get_physical_length() ) ),
    RhoMult ( CP1.RhoMult ),
    RecomModel (CP1.RecomModel ),
    NindToUseForRho (CP1.NindToUseForRho ),
    CurrentLogProb ( 0 ),
    vecTheta ( CP1.vecTheta ),
    BestLogProb ( BIGNEGATIVE ),
    NSNP (CP1.NSNP+CP2.NSNP),
    ALLSNP (CP1.ALLSNP*CP2.ALLSNP),
    TREATSNPSASMS (CP1.TREATSNPSASMS * CP2.TREATSNPSASMS),
    SNPlist (CP1.SNPlist),
    nonSNPlist (CP1.nonSNPlist),
    Qptr (CP1.Qptr),
    vecdiffprob ( CP1.loci_type.size() + CP2.loci_type.size() + 1 ),
    CC (0), //CP1.Nind),
    DiffCount (0),// CP1.Nind),
    haplist ()
{    
  if(CP1.Nind != CP2.Nind){
    cerr << "Error: trying to concatenate classpops of different size!" << endl;
    exit(1);
  }	 

  for(int ind = 0; ind < CP1.pop.size(); ind++)
    pop[ind] = CIndividual(CP1.pop[ind],CP2.pop[ind]);

  position.insert(position.end(), CP2.position.begin(), CP2.position.end());
  Qptr.insert(Qptr.end(), CP2.Qptr.begin(), CP2.Qptr.end() );
  
  vecRho.insert(vecRho.end(), CP2.vecRho.begin(), CP2.vecRho.end());
  vecRhoDeriv.insert(vecRhoDeriv.end(), CP2.vecRhoDeriv.begin(), CP2.vecRhoDeriv.end());
  RhoMult.insert(RhoMult.end(), CP2.RhoMult.begin(), CP2.RhoMult.end());

  if(RecomModel == 0){
    for(int i = 0; i<CP1.RhoMult.size(); i++){
      RhoMult[i] = CP1.RhoMean * CP1.RhoMult[i] / RhoMean; 
      if(RhoMult[i]<0.01)
	RhoMult[i] = 0.01;
    }
    for(int i = 0; i<CP2.RhoMult.size(); i++){
      RhoMult[i+CP1.get_nloci()] = CP2.RhoMean * CP2.RhoMult[i] / RhoMean;
      if(RhoMult[i+CP1.get_nloci()]<0.01)
	RhoMult[i+CP1.get_nloci()] = 0.01;
    }
  }

  if(RecomModel == 1 || RecomModel == 2)
    InitialiseSimpleHotspot(lambda.size(),RecomModel==2);

  vecTheta.insert(vecTheta.end(), CP2.vecTheta.begin(), CP2.vecTheta.end());
    
  //coding = new vector<int> [2];
  coding[0] = vector<int>(CP1.coding[0]);
  coding[1] = vector<int>(CP1.coding[1]);
  coding[0].insert(coding[0].end(),CP2.coding[0].begin(),CP2.coding[0].end());
  coding[1].insert(coding[1].end(),CP2.coding[1].begin(),CP2.coding[1].end());

  ComputeRho();
  ComputeRhoDerivAndCurrentLogProb();
 
  

  // add all four combinations for each individual
  for(int ind = 0; ind < Nind; ind++){
    haplist.Add(Haplotype(CP1.pop[ind].get_haplotype(0),CP2.pop[ind].get_haplotype(0)));
    haplist.Add(Haplotype(CP1.pop[ind].get_haplotype(0),CP2.pop[ind].get_haplotype(1)));
    haplist.Add(Haplotype(CP1.pop[ind].get_haplotype(1),CP2.pop[ind].get_haplotype(0)));
    haplist.Add(Haplotype(CP1.pop[ind].get_haplotype(1),CP2.pop[ind].get_haplotype(1)));
  }
  
  Computevecdiffprob();

  MCMCListResolvePhase(cmds, Niter, 1, Niter, vecDelta, d_cmds, "", false, false);

  haplist = HapList( CP1.haplist, CP2.haplist, pop, minfreq );

}


ClassPop::~ClassPop ()
{
}

// }}}

// {{{ Public functions

// add to haplist the haplist in c
void ClassPop::AugmentHapList(const ClassPop & c,double minfreq)
{
  //haplist.Output(cout,coding);
  haplist.Add(c.haplist,minfreq);
  //haplist.Output(cout, coding);

}

const ClassPop & ClassPop::operator=(const ClassPop & rhs)
{
  if(this != &rhs){
    this->Nloci = rhs.Nloci;
    this->Nind = rhs.Nind;
    this->Nchild = rhs.Nchild;
    this->loci_type = string(rhs.loci_type);
    this->pop = vector<CIndividual>(rhs.pop);
    this->childindex = vector<int>(rhs.childindex);
    this->buddy = vector< vector<int> > (rhs.buddy);
    this->casecontrol = vector<int>(rhs.casecontrol);
    this->groupsize = vector<int>(rhs.groupsize);
    this->position = vector<double>(rhs.position);
    this->order = vector<int>(rhs.order);
    this->vecRho = vector<double>(rhs.vecRho);
    this->vecRhoDeriv = vector<double>(rhs.vecRhoDeriv);   
    this->RhoMean = rhs.RhoMean;
    this->RhoMult = vector<double>(rhs.RhoMult);
    this->right = rhs.right;
    this->left = rhs.left;
    this->lambda = rhs.lambda;
    this->RecomModel = rhs.RecomModel;
    this->CurrentLogProb = rhs.CurrentLogProb;
    this->vecTheta = vector<double>(rhs.vecTheta);
    this->BestLogProb = rhs.BestLogProb;
    this->NSNP = rhs.NSNP;
    this->ALLSNP = rhs.ALLSNP;
    this->TREATSNPSASMS = rhs.TREATSNPSASMS;
    this->SNPlist = vector<int>(rhs.SNPlist);
    this->nonSNPlist = vector<int>(rhs.nonSNPlist);
    this->Qptr = vector<ArrayQ *>(rhs.Qptr);
    this->vecdiffprob = vector<double> (rhs.vecdiffprob);
    this->CC = ArrayCC(rhs.CC);
    this->DiffCount = ArrayDiffCount(rhs.DiffCount);
    this->haplist = HapList( rhs.haplist );
    //this->coding = new vector<int> [2];
    this->coding[0] = vector<int>(rhs.coding[0]);
    this->coding[1] = vector<int>(rhs.coding[1]);
  }
  return *this;
  //  cerr << "Error: ClassPop::operator= called" << endl;
//    exit(1);

  //  for(int ind = 0; ind < Nind; ind++)
//      pop[ind] = CIndividual(cpop2.pop[ind]);
 
}

void ClassPop::InitialiseSimpleHotspot(int nhotspot, bool fixedpos)
{
  
  if(!fixedpos){
    right = vector<double>(nhotspot,position[Nloci-1]);
    left = vector<double>(nhotspot,position[0]);
  }
  lambda = vector<double>(nhotspot,0);
  ComputeRho(right, lambda, left);
}


void ClassPop::read_data ( istream & input,
                           map<string, int> & cmds, map<string, double> & d_cmds, map<string,string> & filenames )
{

    char lt = 'S';

    if ( cmds["usehudson"] == 1 ) {
        input_hudson_data ( input );
    } else {
      input >> Nind;
      NindToUseForRho = (Nind>cmds["NindToUseForRho"]) ? cmds["NindToUseForRho"] : Nind; // limit the number of individuals to use when inferring rho

      buddy = vector<vector<int> > (Nind, vector<int>());

      if(cmds["useped"]){
	childindex = vector<int>(Nind);
	Nchild = Nind / 2;
	for(int i = 0; i<Nchild; i++){
	  childindex[i*2] = i+Nind;
	  childindex[i*2+1] = i+Nind;
	  buddy[i*2].push_back(i*2+1);
	  buddy[i*2+1].push_back(i*2);
	}
	
	// ifstream pedinput (filenames["pedigreefile"].c_str());
// 	assure ( pedinput, filenames["pedigreefile"] );
// 	pedinput >> Nchild;
// 	childindex = vector<int>(Nind);
// 	for(int i = 0; i<Nchild; i++){
// 	  int parent1;
// 	  int parent2;
// 	  int ch;
// 	  pedinput >> ch;
// 	  pedinput >> parent1;
// 	  pedinput >> parent2;
//        buddy[parent1].push_back(parent2);
//        buddy[parent2].push_back(parent1);
// 	  childindex[parent1] = ch;
// 	  childindex[parent2] = ch;
// 	}
      }

      // Number of individuals

        
	
// Number of loci
        input >> Nloci;

        loci_type = "";
	NSNP = 0;
        // read in loci types
	input >> lt;
	if(lt == POSITIONLINEINDICATOR){ // read in positions
	  cerr << "Reading Positions of loci" << endl;
	  double p;
	  double lastp = 0;
	  while( position.size() < Nloci){
	    input >> p;
	    //cout << p << endl;
	    position.push_back(p);
	    if(p<lastp){
	      cerr << "Error: loci must be in correct order along chromosome" << endl;
	      cerr << "locus at position " << lastp << " should be after locus at " << p << endl;
	      exit(1);

	    }
	    lastp = p;
	  }
	}	
	else if(lt !='S' && lt != 'M'){
	  cerr << "Error in input file, character " << lt << " should be S, M or " << POSITIONLINEINDICATOR << endl;
	  exit(1);
	}
	else {
	    position = vector<double> ( Nloci, 0.0 );
	    if(cmds["method"] == 'R')
	      cerr << "Warning: no positions specified in input file; assuming loci equally spaced" << endl;

	    for(int r=0; r<Nloci; r++){
	      position[r] = r*1000; // assume the loci are 1kb apart
	    }
	    loci_type +=lt;
	}


	while ( loci_type.size() < Nloci ) {
	  input >> lt;
	  if ( lt != 'S' && lt != 'M' ) {
	    cerr << "Invalid locus type " << lt << endl
		 << "Must be either 'S' or 'M'" << endl;
	    exit (1);
	  }
	  loci_type += lt;
	}

	//cerr << "Nind = " << Nind << endl;
	//cerr << "Nchild = " << Nchild << endl;

	pop = vector<CIndividual> ( Nind + Nchild, CIndividual ( loci_type ) );
	casecontrol = vector<int> (Nind + Nchild,0);

	// if(cmds["testing"]){ // for testing that p value is uniform under H0
// 	  for(int i=0; i< Nind; i++)
// 	    casecontrol[i] = (ranf()<0.5);
// 	  cmds["nperm"] = 100;
// 	}


	order = vector<int>(Nind, 0);
	for(int i=0; i<Nind;i++){
	  order[i] = i;
	}
	RandomiseOrder();

	vecTheta = vector<double> (Nloci,0);
	
	vecRho = vector<double> (Nloci,0);
	vecRhoDeriv = vector<double> (Nloci,0);
	RhoMult = vector<double> (Nloci,1);
	RhoMean = d_cmds["rhostart"];

	lambda = vector<double>(cmds["nhotspot"],0);
	left = vector<double>(cmds["nhotspot"]);
	right = vector<double>(cmds["nhotspot"]);
 
	if(cmds["nhotspot"]>0){
	  left[0] = d_cmds["left"]; // set left and right of hotspot
	  right[0] = d_cmds["right"];
	}
	if(cmds["nhotspot"]>1){
	  left[1] = d_cmds["left2"]; // set left and right of hotspot
	  right[1] = d_cmds["right2"];
	}

	if(cmds["RecomModel"] == 1 || cmds["RecomModel"] == 2){
	  InitialiseSimpleHotspot(cmds["nhotspot"],cmds["RecomModel"]==2);	  
	}
	else
	  ComputeRho();
	
        SNPlist = vector<int> ();
	nonSNPlist = vector<int> ();
	for( int r = 0; r< Nloci; r++)
	  {
#ifdef TREAT_SNPS_AS_MS
	    nonSNPlist.push_back(r);
#endif

#ifndef TREAT_SNPS_AS_MS
	    if(loci_type[r]=='S'){
	      SNPlist.push_back(r);
	      NSNP++;
	    }
	    else{
	      nonSNPlist.push_back(r);
	    }
#endif
	  }

#ifdef TREAT_SNPS_AS_MS
	TREATSNPSASMS=1;
#endif
 
	if(Nloci == NSNP)
	  ALLSNP=1;
	else
	  ALLSNP=0;

	 
        // Read in the genotypes
        //  if ( cmds["randomise"] == 1 ) {
//              input_random( input, cmds["idpresent"] );
//          } else {
	int j=0;
	for (vector<CIndividual>::iterator i = pop.begin();
	     i != pop.end(); ++i) {
	  cerr << "Reading individual " << setw(6) << (j+1) << "\033[A" << endl;
	  if(cmds["casecontrol"]){
	    input >> casecontrol[j];
	  }
	  
	  j++;
	  i->read_orig_phenotypes( input, loci_type,
				   cmds["idpresent"], cmds["format"] );
	}
	
        // Initialize some other members
       
	BestLogProb = BIGNEGATIVE;
        coding[0] = vector<int> ( Nloci, 0);
        coding[1] = vector<int> ( Nloci, 0);

	//CC.resize(Nind);
	//DiffCount.resize(Nind);

	//knownDiffCount.resize(Nind);
	//	}

	//knownCC.resize(Nind);
	//for(int n=0;n<Nind;n++){
	//  CC[n].resize(Nind);
	//  knownCC[n].resize(Nind);
	//	}
       
	//knownDiffCount.resize(Nind);
	//for(int n=0; n<Nind; n++){
	//  DiffCount[n].resize(Nind);
	//  knownDiffCount.resize(Nind);
	//	}
    }
}

void ClassPop::initialize ( istream & istr_known, istream & istr_init, int knowninfo, int initmethod, double theta, vector<double> & vecDelta , int format )
{
   // Recoding the alleles
  normalize (format);
  
  // Impute missing alleles by random draws from the population
  // also initialize phase
  int count = 0;
  for (vector<CIndividual>::iterator i = pop.begin();
       i != pop.end(); ++i) {
    //if(!istr_start){ cout << "WRONG!" << endl;}
    //else{ cout << "RIGHT" << endl;}
    //cerr << "Initialising Individual " << ++count << endl;

    i->initialize( knowninfo, istr_known, initmethod, istr_init, loci_type );
    for (int locus = 0; locus < loci_type.size(); ++locus) {
      if ( i->n_missing(locus) == 1 ) {
	// Missing one allele
	i->set_allele ( i->missingchr(locus), locus, draw_random_allele (locus) );
      }
      else if ( i->n_missing(locus) == 2 ) {
	// Missing both allele
	i->set_allele ( 1, locus, draw_random_allele (locus) );
	i->set_allele ( 0, locus, draw_random_allele (locus) );
	
      }
    }
  }

  vecTheta= vector<double> ( loci_type.size(), theta );
   
  // Calculate theta if necessary
  if ( theta <= 0.0 ) calc_theta ( ); 
  
  // cout << "Theta Values" << endl;
  //for(int r=0; r<Nloci; r++)
  //  cout << vecTheta[r] << endl;
  //cout << endl;

  cerr << "Computing matrix Q, please wait" << endl;
  ClassPop::Qvalues = vector<ArrayQ>(loci_type.size());
  ClassPop::QSNP = ArrayQ('S', pop.size()+1, vecTheta[0], 0);
  
  Qptr = vector<ArrayQ *>(loci_type.size());
  vecdiffprob = vector<double> (loci_type.size()+1);

  for(int r=0; r<loci_type.size(); r++){
    cerr << "Locus " << setw(5) << (r+1) << "\033[A" << endl;
    //cer << "; Theta = " << vecTheta[r] << endl;
    if(loci_type[r] == 'M'){
      Qvalues[r] = ArrayQ(loci_type[r], pop.size()+1, vecTheta[r], vecDelta[r] );
      Qptr[r] = &Qvalues[r];
    } else
      Qptr[r] = &QSNP;
  }
  cerr << "Done computing Q" << endl;
   
  Computevecdiffprob();

  ComputeRhoDerivAndCurrentLogProb();

#ifdef DEBUG
  // Print out the genotype with fomat 1
  cout << Nind << endl;
  cout << Nloci << endl;
  cout << loci_type << endl;
  for (vector<CIndividual>::const_iterator i = pop.begin();
       i != pop.end(); ++i) {
    cout << i->get_id() << endl;
    i->print_phenotypes( cout, loci_type, coding );
  }
#endif // DEBUG
}


void ClassPop::Computevecdiffprob()
{
  for(int ndiff = 0; ndiff <= get_nloci(); ndiff++){
    vecdiffprob[ndiff] = 0;
    for(int t=0; t<SS; t++){ // add up over possible times
      double tempprob = 1;
      int nsame = get_nloci() - ndiff;
      for(int i = 0; i < ndiff; i++)
	tempprob *= PrHitTarg(0, 2*Nind - 2, t, 0, 1, Qptr);
      for(int i = 0; i < nsame; i++)
	tempprob *= PrHitTarg(0, 2*Nind - 2, t, 0, 0, Qptr);
      vecdiffprob[ndiff] += WEIGHTS[t] * tempprob;
    }
  }
}

//
// make a list of the haplotypes present in the sample
// if usebestguess is true, then use best guess to make list
// otherwise use current haps to make list
void ClassPop::MakeHapList(bool usebestguess) 
{
  haplist.RemoveAll();
  for(int ind =0; ind<Nind; ind++){
    haplist.Add(pop[ind],1.0,usebestguess);
  }
}

//
// make a list of all possible haplotypes present in the sample

void ClassPop::MakeAllPossHapList() 
{
  cerr << "Making List of all possible haplotypes" << endl;
  haplist.RemoveAll();
  haplist.AddAllPossible(pop,Nind);

}

void ClassPop::OutputHaplistSummary(ostream & ostr)
{
  for(int n=0; n<Nind; n++)
    {
      pop[n].print_id(ostr);
      ostr << ": (" << (haplist.Find(pop[n],0,true)+1) << "," << (haplist.Find(pop[n],1,true)+1) << ")" << endl;
    }
}

void ClassPop::OutputPhaseProbs(ostream & ostr, bool PrintAll)
{
  for(int n=0; n<Nind; n++)
    {
      pop[n].print_phase_prob(ostr, PrintAll);      
    }
  
}


double ClassPop::resolve_phase_NR  (
                               int Nburn,
                               int Nthin,
                               int Niter,
                               ostream & ostrHaplotypes,
                               //ostream & ostrPhases,
			       ostream & ostrMonitor,
			       vector<double> & vecDelta, 
			       int verbose, int fAncUpdate, int fNaiveGibbs)
{

#ifdef BIGDATASETS // if compiled for big data sets, 
                   //do not assign memory for FF
    ArrayFF FF ( 1 ); 
    
    ArrayDiffProb DiffProb (loci_type, 0, Qptr);
    ArrayDiploidDiffProb DiploidDiffProb(loci_type, 0, Qptr);
    DiffCount.resize(0);
    CC.resize(0);
#else
    ArrayFF FF ( pop.size() );
    cerr << "Allocating Memory for Arrays" << endl;

    ArrayDiffProb DiffProb (loci_type, pop.size(), Qptr);
    ArrayDiploidDiffProb DiploidDiffProb(loci_type, pop.size(), Qptr);
    DiffCount.resize(Nind);
    CC.resize(Nind);
    cerr << "Computing Initial Values for arrays for each individual" << endl;
    
    DiffCount.compute(pop, SNPlist);
    for(int n=0;n<pop.size(); ++n){
      cerr<< "Individual " << setw(6) << (n+1) << "\033[A";
      CC.compute(n,Qptr,pop,nonSNPlist,loci_type,DiffProb);
      cerr << " done" << endl;
    } 
#endif

   

    
    
    // Probablity of choosing individuals to update
    vector<double> ChooseProb ( pop.size() );
    // Sum of ChooseProb
    double ChooseSum = 0.0;
    for (int n = 0; n < pop.size(); ++n) {
      ChooseProb[n] = pop[n].numunknown() > 0 ? 1.0 : 0.0;
      ChooseSum += ChooseProb[n];
    }
    if ( ChooseSum == 0 ) {
      cerr << "No one is ambiguous!" << endl;
      return 0;
    }
    
    // The person to update at each iteration
    int n1 = 0;
    // {{{ Burn-in
  
    ChooseSum = 0;
    for (int n = 0; n < pop.size(); ++n) {
      ChooseProb[n] = pop[n].numunknown() > 1 ? 1.0 : 0.0;
      ChooseSum += ChooseProb[n];
    }

     
    if(ChooseSum>0){
      for (int iter = 0; iter < Nburn; ++iter) {
	if ( iter % ITPRINT == 0 ) {
	  cerr << "Burn-in " << iter << endl;
	}
	for (int thin = 0; thin < Nthin; ++thin) {	  
	  n1 = rint2 ( ChooseProb, ChooseSum );
	  update_NR ( n1, FF, DiffProb, DiploidDiffProb, fAncUpdate, fNaiveGibbs);
	}
	
	for(int n=0;n<pop.size(); ++n) //recompute CC to avoid rounding errors
	  CC.compute(n,Qptr,pop,nonSNPlist,loci_type,DiffProb);	    
	ostrMonitor << monitor_prob(DiffProb  ) << endl;
      }
    }

    
    // }}}
    // {{{ Main loop
    
    cerr << "Performing Main iterations" << endl;
    ChooseSum = 0;
    for (int n = 0; n < pop.size(); ++n) {
      ChooseProb[n] = pop[n].numunknown() > 1 ? 1.0 : 0.0;
      ChooseSum += ChooseProb[n];
    }
    
    double averagemonitorprob =0;
    for (int iter = 0; iter < Niter; ++iter) {
      
      if ( iter % ITPRINT == 0 ) {
	cerr << "Iteration " << iter << endl;
      }
      for (int thin = 0; thin < Nthin; ++thin) {	
	n1 = rint2 ( ChooseProb, ChooseSum );
	update_NR ( n1, FF, DiffProb, DiploidDiffProb, fAncUpdate, fNaiveGibbs);
      }
      
      //UpdateCounts();
      for(int n=0;n<pop.size(); ++n) //recompute CC to avoid rounding errors
	CC.compute(n,Qptr,pop,nonSNPlist,loci_type,DiffProb);

      double mp = monitor_prob(DiffProb);
      averagemonitorprob += mp;
      ostrMonitor << mp << endl;

      //output_all_phases ( ostrPhases );
      //output( ostrHaplotypes, true, true );
    }

    //  HapList ListOfBestGuesses;
//      ListOfBestGuesses.Add(pop,true);
//      ListOfBestGuesses.MakePositiveHaps();
//      int nchr = Nind+Nind;
//      ListOfBestGuesses.Output(cout,coding);
//      double loglik = ListOfBestGuesses.FullDataPseudoLogLikelihood('S',Qptr,nchr,vecRho,0);
//      cout << "Full Data LogLik: " << loglik << endl;

    averagemonitorprob/=Niter;

    return averagemonitorprob;
    // return loglik;
}

double ClassPop::ListResolvePhase(char method, int Niter, vector<double> & vecDelta, double rho , bool randomise, bool initialise )
{
  vecRho = vector<double> (position.size(),0);

  
  if(initialise){
    MakeAllPossHapList();
  }
  haplist.Output(cout, coding);

  double likelihood = haplist.ResolvePhase(method, Niter, pop, Qptr, vecRho, randomise);
  haplist.Output(cout, coding, 0.001);
  UpdateCounts();

  return likelihood;

}

void ClassPop::OutputHapList(ostream & ostr, double probcutoff , bool printheader )
{
  haplist.Output(ostr, coding, probcutoff*(1.0/(2*Nind)), printheader);
}

// void ClassPop::InitialiseRho(double rho)
// {
//   vecRho = vector<double>(position.size());
//   vecRhoDeriv = vector<double>(position.size());
//   RhoMult = vector<double>(position.size(),1);
//   RhoMean = rho;
  
//   ComputeRho();
//   ComputeRhoDerivAndCurrentLogProb();
// }

double ClassPop::MCMCListResolvePhase(map<string,int> & cmds, int Niter, int Nthin, int Nburn, vector<double> & vecDelta, map<string,double> & d_cmds, string filename, bool initialise , bool collectdata)
{
    
  if(initialise && (cmds["knowninfo"]!=999)){
    haplist.RemoveAll();
    MakeAllPossHapList();
  }

  
  if(cmds["verbose"]){
    cout << "Haplist being used:" << endl;
    haplist.Output(cout, coding);
  }

  
  double loglik;

  if(cmds["useped"])
    loglik = BuddyHapListMCMCResolvePhaseRemove(cmds, Niter, Nthin, Nburn, d_cmds, filename, collectdata );

  else if(cmds["usefuzzy"])
    loglik = FuzzyHapListMCMCResolvePhaseRemove(cmds, Niter, Nthin, Nburn, d_cmds, filename, collectdata );

  else
    loglik = HapListMCMCResolvePhaseRemove(cmds, Niter, Nthin, Nburn, d_cmds, filename, collectdata );

  if(cmds["verbose"]){
    cout << "Final HapList:" << endl;
    haplist.Output(cout, coding);
  }

  //  haplist.Output(cout, coding, 0.001);
  //UpdateCounts();
  
  //cout << "Original LogLik: " << loglik << endl;
 
  // Now compute loglik as the FullDataPseudoLogLikelihood for
  // the best guesses (which have been put into pop by the MCMC method)
  // (actually found this performed slightly worse)
  //  HapList ListOfBestGuesses;
//    ListOfBestGuesses.MakePositiveHaps();
//    ListOfBestGuesses.Add(pop);
//    int nchr = Nind+Nind;
//    loglik = ListOfBestGuesses.FullDataPseudoLogLikelihood(method,Qptr,nchr,vecRho,betaend);
//    cout << "Full Data LogLik: " << loglik << endl;

  return loglik;

}

// following sets up the best guesses for each individual
// using the frequencies in haplist
//void ClassPop::SetHaplistBestGuesses(char use )
//{
//  haplist.SetBestGuesses(pop, use);
//}


void ClassPop::GibbsResolvePhase  (int Niter, double dirprior)
{
     
    // }}}
    // {{{ Main loop
   
    MakeHapList(false); // make list of haplotypes

    for (int iter = 0; iter < Niter; ++iter) {

        if ( iter % ITPRINT == 0 ) {
            cerr << "Iteration " << iter << endl;
        }
	for(int n1=0; n1<Nind; n1++){
	  GibbsUpdate ( n1, dirprior );
        }

    }	
    // }}}
}

// update individual n1
void ClassPop::GibbsUpdate( int n1, double dirprior)
{
  // see if individual n1 can be made up from one or two haps in
  // the current list.
  
  int c,r,r0,newhappos;
 
  // contains -1 if that hap cannot be extracted
  
  vector<double> TempProb (haplist.get_listlength(),0);
  double TempProbSum=0;


  //  cout << "before removing" << n1 << endl;
  //OutputList(cout);

  haplist.Remove(pop[n1]);

  //haplist[h].Output(cout,coding);

  //cout << "after removing" << endl;
  //OutputList(cout);

// attempt to extract haplotypes in list from genotype [n1]
   
  haplist.ComputeVectorOfNaiveGibbsProbs(pop[n1],TempProb,TempProbSum,dirprior);

  // for(int listpos=0; listpos < listlength; listpos++){
  //   cout << listpos << ":" << CompList[listpos] << "," << TempProb[listpos] << endl;
  // }
  
  //cout << TempProbSum << endl;

  double randprob=exp(2*(log(dirprior)-Nloci*log(2.0)) + (pop[n1].numunknown()-1)*log(2.0));
  //cout << randprob << endl;

  // prob of randomising

  

  //if(ranf()>randprob/(randprob+TempProbSum)) {
 
  //Note: the above line is how things really ought to be done for
  //a Dirichlet prior. However, we found the following worked better
  //in practice (which does make intuitive sense:
  //why randomise when we can do better?)

  //Essentially we are doing the case dirprior tends to zero!
  //

  if(TempProbSum>0){
    newhappos=rint2(TempProb,TempProbSum);
    Haplotype newhap = haplist.get_haplotype(newhappos);
    const vector<int> unknown_list =  pop[n1].get_unknown_pos();
    for (int locus = 0; locus < Nloci; locus ++){	
      if(pop[n1].get_haplotype(0,locus) != newhap.get_allele(locus))
	pop[n1].flip_phase(locus);
    }
  }
  else{ 
    // randomise haplotype at unknown positions
    const vector<int> unknown_list =  pop[n1].get_unknown_pos();
    for (vector<int>::const_iterator u = unknown_list.begin();
	 u != unknown_list.end(); ++u ) {     
      pop[n1].set_phase( *u , (int) floor(2*ranf()) );
    }
  }

  haplist.Add(pop[n1]);

}


// }}}

// {{{ Private functions

void ClassPop::input_hudson_data ( istream & istr  )
{
    // Read in Hudson data
    char lt = '0';
    do {
        istr >> lt;
    } while ( lt != ':');

    // Number of loci
    istr >> Nloci;

    loci_type = string (Nloci, 'S');
    position = vector<double> ( Nloci, 0.0 );

    // Read in map positions
    istr >> position[0];
    if ( position[0] < 0 ) {
        cerr << "Error: All positions must be between 0 and 1 "
             << "inclusive" << endl;
        exit (1);
    }
    for (int r = 1; r < Nloci; ++r) {
        istr >> position[r];
        if ( position[r] < position[r-1] ) {
            cerr << "Error: in current implementation loci "
                 << "must be in order along sequence." << endl;
            exit (1);
        }
    }
    if ( position[Nloci - 1] > 1.0 ) {
        cerr << "Error: All positions must be between 0 and 1 "
             << "inclusive" << endl;
        exit (1);
    }
    //
    // Now read the phenotype data into a vector of strings,
    // without knowing how many individuals there are.
    vector<string> instrings;
    string tempbuf;
    // Pass the end of last line
    getline(istr, tempbuf);
    while (1) {
        getline(istr, tempbuf);
        // An empty line signals the end of one data set
        if ( tempbuf.size() == 0 || tempbuf[0] == ' ' ) break;
        instrings.push_back ( tempbuf );
    }

#ifdef DEBUG
    for (vector<string>::const_iterator i = instrings.begin();
         i != instrings.end(); ++i) {
        cout << *i << endl;
    }
#endif // DEBUG

    // Now we have the number of individuals
    pop = vector<CIndividual> ( instrings.size() / 2,
                                CIndividual ( loci_type ) );

    // Randomly pair haplotypes
    cerr << "Warning: Randomly pairing haplotypes!" << endl;
    vector<int> perm ( instrings.size() );
    for (int n = 0; n < perm.size(); ++n) {
        perm[n] = n;
    }
    random_shuffle ( perm.begin(), perm.end() );
    for (int n = 0; n < perm.size(); ++n) {
        pop[n/2].set_id ();
        int chr = perm[n] % 2;
        int n1  = perm[n] / 2;
        for ( int locus = 0; locus < loci_type.size(); ++locus ) {
            pop[n1].set_allele (chr, locus, instrings[n][locus]);
	    pop[n1].set_original_allele(chr, locus, instrings[n][locus]);
        }
    }
    // Initialize some other members
   
    coding[0] = vector<int> ( Nloci, 0);
    coding[1] = vector<int> ( Nloci, 0);
}

void ClassPop::input_random ( istream & istr, int idpresent )
{
    cerr << "Warning: Randomly pairing haplotypes!" << endl;
    vector<int> perm ( 2 * pop.size() );
    for (int n = 0; n < perm.size(); ++n) {
        perm[n] = n;
    }
    random_shuffle ( perm.begin(), perm.end() );

    // Read in haplotypes
    for (int n = 0; n < perm.size(); ++n) {
        if ( n % 2 == 0 ) {         // new individual
            if ( idpresent == 1 ) {
                pop[n/2].set_id ( istr );
            } else {
                pop[n/2].set_id ();
            }
        }
        int chr = perm[n] % 2;
        int n1  = perm[n] / 2;

        for ( int locus = 0; locus < loci_type.size(); ++locus ) {
            pop[n1].input_orig_allele (istr, loci_type[locus], chr, locus);
        }
    }
}

//
// following function normalizes to given coding
// (used when postprocessing to normalize all samples
// the same way)
//
void ClassPop::renormalize ( ClassPop & cp)
{
    int allele = 0;
    vector <int> newcoding0 = cp.get_coding(0);

    for (int locus = 0; locus < loci_type.size(); ++locus) {
        if ( loci_type[locus] == 'M' ) {
            // Center the number of repeats
            for (vector<CIndividual>::iterator i = pop.begin();
                 i != pop.end(); ++i) {
                for ( int chr = 0; chr < 2; ++chr ) {
                    allele = i->get_allele(chr, locus);
                    if ( allele != MISSMS ) {
                        i->set_allele (chr, locus, allele -
                            coding[0][locus] + newcoding0[locus]);
			i->set_original_allele (chr, locus, allele -
                            coding[0][locus] + newcoding0[locus]);
			
                    }
                }
            }
        } else {
            // For SNP sites, recode so that the two alleles are 0 and 1
            for (vector<CIndividual>::iterator i = pop.begin();
                 i != pop.end(); ++i) {
                for ( int chr = 0; chr < 2; ++chr ) {
                    if ( coding[0][locus] != newcoding0[locus] ) {
                        // Reverse coding 0 vs 1
                        allele = i->get_allele(chr, locus);
                        i->set_allele ( chr, locus, 1-allele );
			i->set_original_allele ( chr, locus, 1-allele );
			
                    }
                }
            }
        }
        coding[0][locus] = newcoding0[locus];
    }
}


void ClassPop::normalize ( int format )
{
    int max_repeats   = 0;
    int min_repeats   = 0;
    int range_repeats = 0;
    int offset = 0;
    int allele = 0;

    for (int locus = 0; locus < loci_type.size(); ++locus) {
        if ( loci_type[locus] == 'M' ) {
            // Find out the number of alleles for each MS site
            max_repeats = -BIGINT;
            min_repeats = BIGINT;
            for (vector<CIndividual>::const_iterator i = pop.begin();
                 i != pop.end(); ++i) {
                for ( int chr = 0; chr < 2; ++chr ) {
                    allele = i->get_orig_allele(chr, locus);
                    if ( allele  == MISSMS ) {
                        continue;
                    } else if ( allele < min_repeats ) {
                        min_repeats = allele;
                    } else if ( allele > max_repeats ) {
                        max_repeats = allele;
                    }
                }
            }

            range_repeats = max_repeats - min_repeats;

            if ( KMAX <= range_repeats ) {
                cerr << "Error: number of alleles too large ( K = "
                     << range_repeats + 1 << ", KMAX = "
                     << KMAX << ") at locus " << locus << endl;
                cerr << "Try increasing KMAX in constants.h and "
                     << "recompiling." << endl;
                exit(1);
            } else {
                // Center the number of repeats
                offset = (KMAX - max_repeats - min_repeats) / 2;
                for (vector<CIndividual>::iterator i = pop.begin();
                     i != pop.end(); ++i) {
                    for ( int chr = 0; chr < 2; ++chr ) {
                        allele = i->get_orig_allele(chr, locus);
                        if ( allele != MISSMS ) {
                            i->set_allele (chr, locus, allele + offset);
			    i->set_original_allele (chr, locus, allele + offset);

                        }
                    }
                } 
            }

            coding[0][locus] = offset;
            coding[1][locus] = range_repeats;

        } else {

	  
            // For SNP sites, recode so that the two alleles are 0 and
            // 1. Store the decode theme.
	  int firstnonmissing=-1;
	  int azero;
	  do{ // find first non missing character (and non 'H' if format =2)
	    firstnonmissing++;
	    azero = pop[firstnonmissing].get_orig_nonmissing_allele(locus);
	  }while( (azero==MISSMS)
		  || ((format == 2) && (azero == (int) 'H')) );
	  
	  int aone = -1;
	  for (vector<CIndividual>::iterator i = pop.begin();
	       i != pop.end(); ++i) {
	    for ( int chr = 0; chr < 2; ++chr ) {
	      allele = i->get_orig_allele(chr, locus);
	      if ( (allele  == MISSMS) || 
		   ( (allele == (int) 'H') && (format == 2) ) ) {
		continue;
	      } else if ( allele == azero ) {
		i->set_allele ( chr, locus, 0 );
		i->set_original_allele ( chr, locus, 0 );
	      } else {
		if ( aone == -1 ) {
		  aone = allele;
		} else if(aone != allele){
		  cerr << "Error in input file: more than 2 alleles at SNP locus." << endl;
		  cerr << "Individual = ";
		  i->print_id(cerr);
		  cerr << "; Locus = " << (locus+1) << endl; 
		  cerr << "Alleles are: " << (char) azero << "," << (char) aone << ", and " << (char) allele << endl;
		  cerr << "(Note that this error may indicate other problems with the formatting of" << endl;
		  cerr << "the input file, and not actually the use of multiple alleles at a SNP locus)" << endl;
		  exit(1);
		}
		i->set_allele ( chr, locus, 1 );
		i->set_original_allele ( chr, locus, 1 );
		
	      }
	    }
	  }
	  
	  if(aone == -1){ // if only been one allele set other allele
	    aone = azero+1; // to be that allele + 1 (used for format =2
	    // in case where no minor allele hets observed).
	  }

	  if(format == 2){ // go through sites, replacing H with 0 and 1
	    for (vector<CIndividual>::iterator i = pop.begin();
	       i != pop.end(); ++i) {
	      if(i->get_orig_allele(0,locus) == (int) 'H'){
		i -> set_allele( 0, locus, 0);
		i -> set_allele( 1, locus, 1);
		i -> set_original_allele( 0, locus, 0);
		i -> set_original_allele( 1, locus, 1);
	      }
	    }
	  }
	  coding[0][locus] = azero;
	  coding[1][locus] = aone;
        }
    }
}



int ClassPop::draw_random_allele ( int locus ) const
{
    // Draw a random number between 0 and 2 * Nind
    int r = 0;
    int allele = MISSMS;
    double ndraw = 0;
    while ( (allele == MISSMS) && (++ndraw < 1.0E06)) {         // nonmissing allele
        // Draw a random number between 0 and 2 * pop.size()
        r = (int) ( 2.0 * pop.size() * ranf() );
        allele = pop[r/2].get_allele( r % 2, locus );

    }
    if (allele == MISSMS){
      cerr << "Error: failed to find a non-missing allele at locus " << (locus+1) << endl;
      exit(1);
    }
    return allele;
    
}


void ClassPop::output_all_phases ( ostream & ostr , bool printknownphase) const
{
    for (vector<CIndividual>::const_iterator i = pop.begin();
         i != pop.end(); ++i) {
        i->print_phase( ostr, printknownphase );
        ostr << endl;
    }
    ostr << endl;
}


void ClassPop::output_all_haps    ( ostream & ostr  , bool printknownphase , bool printnames, bool printbestguess, bool printmissing, double PhaseThreshold, double AlleleThreshold) const
{
  int j = 0;
    for (vector<CIndividual>::const_iterator i = pop.begin();
         i != pop.end(); ++i) {
      if(printnames){
	ostr << casecontrol[j++] << " ";
	i->print_id( ostr );
	ostr << endl;
      }
      i->print_haplotypes( ostr, loci_type, coding, printknownphase, printmissing, printbestguess, PhaseThreshold, AlleleThreshold );
    }
    //ostr << endl;
}

void ClassPop::output ( ostream & ostr  , bool printknownphase, bool printnames) const
{
  ostr << pop.size() << endl;
  ostr << loci_type.size() << endl;
  ostr << loci_type << endl;
  output_all_haps(ostr,printknownphase,printnames);
}

void ClassPop::print_id ( ostream & outfile  , int n) const
{
  pop[n].print_id(outfile);
}

// print allele of nth ind, chromosome chr, locus locus
void ClassPop::print_allele( ostream & outfile , int n, int chr, int locus, bool printknownphase) const
{
  pop[n].print_allele(outfile,chr,locus,loci_type,coding,printknownphase);
}

double ClassPop::calc_heterozygosity ( int locus ) const
{
    // Heterozygosity = Pr(2 randomly chosen are different)
  double H = 0.0;
  vector<double> freq( KMAX, 0.0 );
  double temp = 1.0 / (2.0 * pop.size());
  for (vector<CIndividual>::const_iterator i = pop.begin(); i != pop.end(); ++i){
    for (int c = 0; c < 2; ++c) {
      freq[i->get_allele(c, locus)] += temp; 
      //cout << "allele: " << ( i->get_allele(c,locus) ) << endl;
    }
  }
  for (int k1 = 0; k1 < KMAX; ++k1) {
    for (int k2 = 0; k2 < KMAX; ++k2) {
      if ( k1 != k2) {
	H += freq[k1] * freq[k2];
      }
    }
  }
  return H;
}

void ClassPop::calc_theta ()
{
    for (int r = 0; r < vecTheta.size(); ++r) {
        if ( loci_type[r] == 'S' ) {
            vecTheta[r] = 1.0 / log (2.0 * pop.size());
        } else {
            vecTheta[r] = 1.0 / ( 1 - calc_heterozygosity ( r ) );
            vecTheta[r] = ( vecTheta[r] * vecTheta[r] - 1.0) * 0.5;
        }	
    }
}
//
// return pseudo-likelihood for current configuration
//
double ClassPop::monitor_prob(const ArrayDiffProb & DiffProb)
{
  
  //  int firstSNP= 0; //first SNP locus;
//    while((loci_type[firstSNP]!='S') && (firstSNP<(Nloci-1)))
//      firstSNP++;

//    double logprob = 10*Nind*NSNP*
//      (log(Qptr(firstSNP,2*Nind-1,0,0,0))+log(Qptr(firstSNP,2*Nind-2,0,0,0))); // a constant of prop

  double logprob = 0;
  for(int n1=0; n1 < Nind; n1++){
    double temp=log(DiffCount.CombineProb(n1,NSNP,DiffProb,CC));
    if(temp > -1000 && temp < 1000){ //prevent numerical problems
      logprob += temp;
    }
    
  }
  if(logprob > BestLogProb){
    SaveCurrentState();
    BestLogProb = logprob;
  }
  return logprob;
}


// void ClassPop::UpdateMissing_NR ( int id, const ArrayDiffProb & DiffProb)
// {
//   vector<double> tempprob(Nind*2*SS);
//   double total_sum=0;
//   int n0,c0,t0;

//   vector<int> missingSNP = pop[id].get_missingSNP_list();
//   vector<int> missingMS = pop[id].get_missingMS_list();
 
//   if(missingSNP.size()+missingMS.size()>0){

//     vector<int> notmissingSNP = pop[id].get_notmissingSNP_list();
//     vector<int> notmissingMS = pop[id].get_notmissingMS_list();
    
    
//     int num_notmissingSNP = notmissingSNP.size();
    
//     // update the missing data positions for chromosome 0 of id
    
//     int nchr = Nind + Nind -2;
//     int ptr=0;
//     total_sum=0;
//     for(n0 =0; n0 < Nind; n0++){
//       if(n0 != id){
// 	for(c0 = 0; c0 < 2; c0++){
// 	  int nd=NDiff(pop,id,0,n0,c0,notmissingSNP);
// 	  for(t0 = 0; t0 < SS; t0++){
// 	    tempprob[ptr] = DiffProb(nchr, t0, nd, 1) 
// 	      * DiffProb(nchr, t0, num_notmissingSNP - nd, 0) 
// 	      * CCProb(pop,id,0,n0,c0,t0,nchr,Qptr,notmissingMS);
// 	    total_sum += tempprob[ptr++];
// 	  }
// 	}
//       }
//     }
    
//     ptr = 0;
//     // sample n0,c0,t0 from tempprob

//     int temp = rint2 ( tempprob, total_sum );


//     for (n0 = 0; n0 < Nind; ++n0) {
//       if(n0 != id){
// 	for (c0 = 0; c0 < 2; ++c0) {
// 	  for (t0 = 0; t0 < SS; ++t0) {
// 	    if ( ptr == temp ) {	  
// 	      goto ENDLOOP1;
// 	    } else {
// 	      ++ptr;
// 	    }
// 	  }
// 	}
//       }
//     }
//   ENDLOOP1: 
    
//     //cout << pop[id].get_missingSNP_list().size() << endl;
//     //cout << pop[id].get_missingMS_list().size() << endl;
    
//     //cout << "HERE1" << id << " " << n0 << " " <<  c0 << " " << t0 << " " << nchr << " " << endl;
//     ImputeMissingPositions(id, 0, n0, c0, t0, nchr, missingSNP);
//     //cout << "HERE2" << id << " " << n0 << " " <<  c0 << " " << t0 << " " << nchr << " " << endl;
//     ImputeMissingPositions(id, 0, n0, c0, t0, nchr, missingMS);
    
//   // update missing data positions for chromosome 1
//     nchr++;

//     // this gives a more elegant way of sampling from tempprob, 
//     //but untested; and is it faster?
//  //     ptr=0;
// //      total_sum=0;
// //      for(n0 =0; n0 < Nind; n0++){
// //        for(c0 = 0; c0 < 2; c0++){
// //  	int nd=NDiff(pop,id,1,n0,c0,notmissingSNP);	
// //  	for(t0 = 0; t0 < SS; t0++){
// //  	  if((n0!=id) || (c0==0))
// //  	    tempprob[ptr] = DiffProb(nchr, t0, nd, 1) 
// //  	      * DiffProb(nchr, t0, num_notmissingSNP - nd, 0) 
// //  	      * CCProb(pop,id,1,n0,c0,t0,nchr,Qptr,notmissingMS);
// //  	  else
// //  	    tempprob[ptr] = 0;
// //  	  total_sum += tempprob[ptr++];
// //  	}
// //        }	
// //      }
    
    
// //      ptr = 0;
    
// //      // sample n0,c0,t0 from tempprob
// //      temp = rint2 ( tempprob, total_sum );
// //      t0 = temp % SS;
// //      temp = temp / SS;
// //      c0 = temp % 2;
// //      temp = temp / 2;
// //      n0 = temp;


//     ptr=0;
//     total_sum=0;
//     for(n0 =0; n0 < Nind; n0++){
//       for(c0 = 0; c0 < 2; c0++){
// 	if((n0!=id) || (c0==0)){
// 	  int nd=NDiff(pop,id,1,n0,c0,notmissingSNP);	
// 	  for(t0 = 0; t0 < SS; t0++){
// 	    tempprob[ptr] = DiffProb(nchr, t0, nd, 1) 
// 	      * DiffProb(nchr, t0, num_notmissingSNP - nd, 0) 
// 	      * CCProb(pop,id,1,n0,c0,t0,nchr,Qptr,notmissingMS);
// 	    total_sum += tempprob[ptr++];
// 	  }
// 	}
//       }
//     }
    
//     ptr = 0;
    
//     // sample n0,c0,t0 from tempprob
//     temp = rint2 ( tempprob, total_sum );
    

//     for (n0 = 0; n0 < Nind; ++n0) {
//       for (c0 = 0; c0 < 2; ++c0) { 
// 	if((n0!=id) || (c0==0)){
// 	  for (t0 = 0; t0 < SS; ++t0) {
// 	    if ( ptr == temp ) {	  
// 	      goto ENDLOOP2;
// 	    } else {
// 	      ++ptr;
// 	    }
// 	  }
// 	}
//       }
//     }
//   ENDLOOP2: 
       
//     //cout << pop[id].get_missingSNP_list().size() << endl;
//     //cout << pop[id].get_missingMS_list().size() << endl;
//     //cout << "HERE3" << id << " " << n0 << " " <<  c0 << " " << t0 << " " << nchr << " "  << endl;
//     ImputeMissingPositions(id, 1, n0, c0, t0, nchr, missingSNP);
//     //cout << "HERE4" << id << " " << n0 << " " <<  c0 << " " << t0 << " " << nchr << " "  << endl;
//     ImputeMissingPositions(id, 1, n0, c0, t0, nchr, missingMS);
    
//     Resymmetrize(DiffCount,id);
//     if(ALLSNP==0)
//       Resymmetrize(CC,id);
//   }
// }

// update the positions in uselist
// for individual n0,c0, given it copied
// n1,c1 at time t
void ClassPop::ImputeMissingPositions(int n0, int c0, int n1, int c1, int t, int nchr, const vector<int> & uselist)
{
  //  if(uselist.size()>0){
//      cout << "Imputing " << n0 << "," << c0  << " from " << n1 << "," << c1 << endl;
//      pop[n0].print_haplotype(c0, cout,loci_type,coding,true);
//      pop[n1].print_haplotype(c1, cout,loci_type,coding,true);
//    }

  for (vector<int>::const_iterator u = uselist.begin();
          u != uselist.end(); ++u ) {  
    //cout << "Locus " << *u << endl;
    
    int from = pop[n1].get_haplotype(c1,*u); 
    int new_allele = impute_allele (*u, nchr, t, from);
    //cout << "from " << from << "to " << new_allele << endl; 
    int oldtarg0 = pop[n0].get_haplotype(0,*u);
    int oldtarg1 = pop[n0].get_haplotype(1,*u);
    pop[n0].update_haplotype (c0, *u, new_allele);

    if(loci_type[*u]!='S' || TREATSNPSASMS)
      CC.Update(n0,Qptr,pop,*u,oldtarg0,oldtarg1);
    else
      DiffCount.Update(n0,pop,*u,oldtarg0,oldtarg1);
  }

  //  if(uselist.size()>0){
//      cout << "After: " << n0 << "," << c0  << " from " << n1 << "," << c1 << endl;
//      pop[n0].print_haplotype(c0,cout,loci_type,coding,true);
//      pop[n1].print_haplotype(c1,cout,loci_type,coding,true);
//      cout << endl;
//    }

}

void ClassPop::SaveCurrentState()
{
  for(int n=0; n<Nind; n++){
    pop[n].SaveCurrentState();
  }
}

void ClassPop::SaveCurrentState(int n)
{
  pop[n].SaveCurrentState();
}

void ClassPop::RestoreSavedState(int n)
{
  pop[n].RestoreSavedState();
}

void ClassPop::RestoreSavedState()
{
  for(int n=0; n<Nind; n++){
    pop[n].RestoreSavedState();
  }
}

void ClassPop::update_NR ( int n1,
                        ArrayFF & FF,
			const ArrayDiffProb & DiffProb,
			const ArrayDiploidDiffProb & DiploidDiffProb,
			int fAncUpdate,
                        bool fNaiveGibbs)
{     
  if (fNaiveGibbs ){
    GibbsUpdate(n1,DIRPRIOR);
  } else {
    //#ifdef DEBUG
    //CC[n1].compute(n1, Qptr, pop, SNPlist, loci_type, DiffProb);
    //cout << "CC probability " << CC[n1].calc_prob() << endl;
    //DiffCount[n1].compute(n1, pop, SNPlist);
    //cout << "Diffcount probability " << DiffCount[n1].calc_prob(n1,NSNP,DiffProb) << endl;
    //#endif
    
    
    //NOTE: could make more efficient by ignoring missing positions
    // in update_phase_NR_fastestforsmallr
    
#ifdef BIGDATASETS     
    update_phase_NR_fastestforsmallr ( n1, DiffProb, 5);  	
#else
    if ( pop[n1].numunknown() < 7)
      update_phase_NR_fastestforsmallr ( n1, DiffProb, 7);    
    else      
      update_phase_NR ( n1, FF, DiffProb, DiploidDiffProb);
#endif
    //UpdateMissing_NR ( n1, DiffProb);         
  }  
}

//
// compute probs of phase at SNP positions in updatelist, based on
// the DiffCount matrix. Put it in
// an array starting at the position pointed to by prob_pointer,
// and return the pointer so it points at the next position
// (implemented recursively)
//
void ClassPop::diff_calc_phase_prob(int n, vector<int> updatelist, vector<double>::iterator & prob_pointer, const ArrayDiffProb & DiffProb )
{
#ifdef DEBUG
  cout << "diff_calc_phase_prob called with n=" << n << " and updatelist length " << updatelist.size() << endl;
#endif
  if(updatelist.size()==1){   
    *prob_pointer++ = DiffCount.calc_prob(n,NSNP, DiffProb);
    pop[n].flip_phase(updatelist[0]);
    DiffCount.Update(n,pop,updatelist[0]);
    *prob_pointer++ = DiffCount.calc_prob(n,NSNP, DiffProb);
  }
  else {
    int lastpos = updatelist.back(); 
    updatelist.pop_back(); // remove end from update list
    diff_calc_phase_prob(n,updatelist,prob_pointer,DiffProb);
    pop[n].flip_phase(lastpos);	
    DiffCount.Update(n,pop,lastpos);
    diff_calc_phase_prob(n,updatelist,prob_pointer,DiffProb);
  }
}

//
// compute probs of phase at positions in updatelist, put it in
// an array starting at the position pointed to by prob_pointer,
// and return the pointer so it points at the next position
// (implemented recursively)
//
void ClassPop::calc_phase_prob(int n, vector<int> updatelist, vector<double>::iterator & prob_pointer, const ArrayDiffProb & DiffProb)
{
#ifdef DEBUG
  cout << "calc_phase_prob called with n=" << n << " and updatelist length " << updatelist.size() << endl;
#endif

  
  if(updatelist.size()==1){
    *prob_pointer++ = DiffCount.CombineProb(n,NSNP,DiffProb,CC);
    
    pop[n].flip_phase(updatelist[0]);
    
    if(loci_type[updatelist[0]]!='S' || TREATSNPSASMS)
      CC.Update(n,Qptr,pop,updatelist[0]);
    else
      DiffCount.Update(n,pop,updatelist[0]);	
    
    *prob_pointer++ = DiffCount.CombineProb(n,NSNP,DiffProb,CC);
    
  }
  else {
    int lastpos = updatelist.back(); 
    updatelist.pop_back(); // remove end from update list
    calc_phase_prob(n, updatelist, prob_pointer, DiffProb);
    pop[n].flip_phase(lastpos);	
    
    if(loci_type[lastpos]!='S' || TREATSNPSASMS )
      CC.Update(n,Qptr,pop,lastpos);
    else
      DiffCount.Update(n,pop,lastpos);
    
    calc_phase_prob(n, updatelist, prob_pointer, DiffProb);
  }
    
}





//
// this routine converts the phase number corresponding to position phasenum
// in the prob vector created by the recursive calls to calc_phase_prob and diff_calc_phase_prob
// into a phase call at position pos (pos=0 is at the start of the updatelist).
//
int convert_phase(int phasenum, int pos){
  return (((phasenum+(1 << pos))/ (1 << (pos+1))) % 2) ;
}

void ClassPop::update_phase_NR ( int n,
                                 ArrayFF & FF, const ArrayDiffProb & DiffProb, 
				 const ArrayDiploidDiffProb & DiploidDiffProb)
{

#ifdef BIGDATASETS
  cout << "Error: this version is compiled for big data sets, but is calling";

cout << " update_phase_NR, which is not allowed for big data sets" << endl;
  exit(1);
#endif

    int n0 = 0;
    int c0 = 0;
    int t0 = 0;
    int n1 = 0;
    int c1 = 0;
    int t1 = 0;

    // Now choose n0 -- t1 at random, according to FF
    // FF is computed first before sampling
    FF.sample_froms ( n, &n0, &c0, &t0, &n1, &c1, &t1, Qptr, DiffProb, DiploidDiffProb, pop, loci_type );

    //cout << n0 << "," << c0 << "," << t0 << "," << n1 << "," << c1 << "," << t1 << endl;

    vector<double> tempprob(2, 0.0);

    vector<double> errorprobs(4,0.0); // There are 4 types of error: no error, 
    //error in allele 0, error in allele 1, and  error in both!

    int r = 0;
    int targ0 = 0;
    int targ1 = 1;
    int from0 = 0;
    int from1 = 0;
    int oldtarg0=0;
    int oldtarg1=0;

    int nchr=Nind+Nind-2;
    int nchrplus1 = nchr+1;

    int newphase;
    // Now choose the phase at random at each position
    // conditional on n0 -- t1
    const vector<int> update_list =  pop[n].get_unknown_pos();
    for (vector<int>::const_iterator u = update_list.begin();
          u != update_list.end(); ++u ) {
      int r = *u;
      from0 = pop[n0].get_haplotype(c0, r);
      from1 = pop[n1].get_haplotype(c1, r);
      targ0 = pop[n].get_allele(0, r);
      targ1 = pop[n].get_allele(1, r);
      oldtarg0 = pop[n].get_haplotype(0, r);
      oldtarg1 = pop[n].get_haplotype(1, r);
	
      switch ( pop[n].n_missing(r) ) {
      case 0:            
	if (n1 == n) {
	  tempprob[0] = PrHitTarg(r, nchr, t0, from0, targ0, Qptr) *
	    PrHitTarg(r, nchrplus1, t1, targ0, targ1, Qptr);
	  tempprob[1] = PrHitTarg(r, nchr, t0, from0, targ1, Qptr) *
	    PrHitTarg(r, nchrplus1, t1, targ1, targ0, Qptr);
	} else {
	  tempprob[0] = PrHitTarg(r, nchr, t0, from0, targ0, Qptr) *
	    PrHitTarg(r, nchrplus1, t1, from1, targ1, Qptr);
	  tempprob[1] = PrHitTarg(r, nchr, t0, from0, targ1, Qptr) *
	    PrHitTarg(r, nchrplus1, t1, from1, targ0, Qptr);
	}
	// Update phase
	newphase=  rint2(tempprob,tempprob[0] + tempprob[1]);
	
	pop[n].set_phase (r, newphase);	      
	break;
      case 1:
	// allele 1 was imputed
	if (n1 == n) {
	  tempprob[0]=tempprob[1]=1;
	}
	else{
	  tempprob[0] = PrHitTarg(r, nchr, t0, from0, targ0, Qptr);
	  tempprob[1] = PrHitTarg(r, nchr, t1, from1, targ0, Qptr);
	}
	
	// Update phase
	newphase=  rint2(tempprob,tempprob[0] + tempprob[1]);
	pop[n].set_phase (r, newphase);    
	
	// Impute allele 1
	if(newphase==0)
	  pop[n].update_haplotype (1, r,
				impute_allele (r, nchrplus1, t1, from1));
	else
	  pop[n].update_haplotype (1, r,
				impute_allele (r, nchrplus1, t0, from0));
	
	break;
      case 2:
	// Impute both alleles
	//cout << "Updating alleles in ind " << n << ", position" << r << endl;
	//pop[n].set_phase (r, 0);
	pop[n].update_haplotype (0, r,
			      impute_allele (r, nchr, t0, from0));
	pop[n].update_haplotype (1, r,
			      impute_allele (r, nchrplus1, t1, from1));
	
	
      }
      
      if(loci_type[r]!='S' || TREATSNPSASMS)
	CC.Update(n,Qptr,pop,r,oldtarg0,oldtarg1);
      else
	DiffCount.Update(n,pop,r,oldtarg0,oldtarg1);
      
    }
    Resymmetrize(DiffCount,n);
    if(ALLSNP==0)
      Resymmetrize(CC,n);
    
}

int ClassPop::update_phase_NR_fastestforsmallr ( int n,
						  const ArrayDiffProb & DiffProb,
						  int maxnpos)
{
    vector<int> updatelist ( pop[n].get_unknown_pos() );
    vector<int> noupdatelist ( pop[n].get_known_pos() );

    //noupdatelist.push_back ( updatelist.back() ); // remove one locus
    //updatelist.pop_back(); // from update list (arbitrarily)
    // The above is OK if all noupdate are homozygotes, but
    // if some heterozygote sites are known then it is not
    // allowed. I've therefore removed this for safety right
    // now, though could reinstate it with the appropriate
    // test in the future.

    if( updatelist.size() < 1) // don't update empty list!
      return 0;

    if( updatelist.size() > maxnpos) {
        random_shuffle ( updatelist.begin(), updatelist.end() );

        // take the loci at the end of the update list and add them to the end of
        // the no-update list
        for(int i = updatelist.size(); i > maxnpos; --i){
            noupdatelist.push_back ( updatelist.back() );
            updatelist.pop_back();
        }
    }

    vector <int> start_phase(updatelist.size(),0);
    
    for(int u=0; u< updatelist.size(); u++){
      start_phase[u]=pop[n].get_phase(updatelist[u]);
    }

    int Nphase = ( 1 << updatelist.size() ); // 2^ number of unknown positions
    vector<double> tempprob(Nphase, 0.0);
    vector<double>::iterator p= tempprob.begin();
    calc_phase_prob(n, updatelist, p, DiffProb );
        
    int newphase = rint2(tempprob);

#ifdef DEBUG
    int j=0;
    for(vector<double>::iterator i = tempprob.begin(); i != tempprob.end(); i++){
      cout << *i << " ";
      for(int u=0; u < updatelist.size(); u++)
	cout << convert_phase(j,u);
      j++;
      cout << endl;
    }
    cout << "New Phase = " << newphase << endl;
    cout << convert_phase(newphase,0) << convert_phase(newphase,1) << convert_phase(newphase,2) << convert_phase(newphase,3) << endl;
    cout << endl;
#endif
       
    // set phase to be newphase
    int np;
    for ( int u = 0; u < updatelist.size(); ++u ) {
      int r  = updatelist[u];      
      np=((convert_phase(newphase,u) +  start_phase[u]) % 2);
      if(np != pop[n].get_phase(r)){
	pop[n].flip_phase(r);
	if(loci_type[r]!='S' || TREATSNPSASMS){
	  CC.Update(n,Qptr,pop,r);
	}
	else
	  DiffCount.Update(n,pop,r);	  
      }
    }
    Resymmetrize(DiffCount,n);
    if(ALLSNP==0)
      Resymmetrize(CC,n);
}



// chooses an allele at locus r, at
//random, given that it copied "from" at time t
int ClassPop::impute_allele (int r, int nchr,
                             int t, int from) const
{
    int n_allele = ( loci_type[r] == 'S' ) ? 2 : KMAX;
    vector<double> tempprob ( n_allele, 0.0 );
    double sum = 0.0;
    for (int a = 0; a < n_allele; ++a) {
        tempprob[a] = PrHitTarg(r, nchr, t, from, a, Qptr);
        sum += tempprob[a];
    }
    return rint2(tempprob, sum);
}

//update counts that are used to determine best phase and
//best alleles
void ClassPop::UpdateCounts()
{
  for(int n=0; n < Nind; n++)
    pop[n].UpdateCounts();
}

//reset counts that are used to determine best phase and
//best alleles
void ClassPop::ResetCounts()
{
  for(int n=0; n < Nind; n++)
    pop[n].ResetCounts();
}

void ClassPop::TransferCounts(vector<Summary> & sum)
{
  for(int n=0; n < Nind; n++)
    pop[n].TransferCounts(sum[n]);
}


// compute the log of the PAC likelihood of all individuals whose casecontrol label is "group" (all individuals if group < 0)
double ClassPop::logProb(char method, vector<double> & rho, double Dprior, int group )
{
  return logProb(method, rho, casecontrol, Dprior, group );
}

// compute the log of the PAC likelihood for individuals whose "groupvec" label is 
// "group" (all individuals if group < 0)
// (use the order in order[])
double ClassPop::logProb(char method, vector<double> & rho, vector<int> & groupvec, double Dprior, int group )
{
  HapList templist; // templist contains the list of haplotypes,
  // built up from an empty list until it contains them all
  // adding the conditional probabilities at each step
  double logprob = 0;
  int nchr = 0;

  templist.RemoveAll();

  for(int n=0; n<Nind; n++){
    int ind = order[n];
    for(int chr=0; chr<2; chr++){
      if((group < 0) || (groupvec[ind] == group)){
	double condprob = templist.CalcProb(pop[ind].get_haplotype(chr), method, Qptr, nchr, rho, Dprior, vector<int>(0), ALLSNP, vecdiffprob);

	logprob += log(condprob);
	//cout << "Templist before:" << endl;
	//templist.Output(cout, coding);
	templist.Add(pop[ind],chr,false,1.0);
	templist.MakePositiveHaps();  
	
	//cout << "Templist after:" << endl;
	//templist.Output(cout, coding);
	nchr++;
      }
    }
  }
  return logprob;
}

// compute the log of the Prob of all individuals whose casecontrol label is "group" (all individuals if group < 0)
double ClassPop::logFDLSProb(vector<double> & rho, vector<double> & rhoderiv, bool computederiv, int group )
{
  return logFDLSProb(rho, rhoderiv, computederiv, casecontrol, group);
}

// compute the log of the Prob of the individuals whose "groupvec" label is 
// "group" (all individuals if group < 0)
// uses non-quadrature version (for inference of rho)
double ClassPop::logFDLSProb(vector<double> & rho, vector<double> & rhoderiv, bool computederiv, vector<int> & groupvec, int group )
{
  HapList templist; // templist contains the list of haplotypes,
  // built up from an empty list until it contains them all
  // adding the conditional probabilities at each step

  templist.RemoveAll();
  
  double logprob = 0;
  int nchr = 0;

  vector<double> deriv(rhoderiv.size(),0);
  for(vector<double>::iterator i=rhoderiv.begin(); i!=rhoderiv.end(); i++)
    *i=0;

  for(int n=0; n<NindToUseForRho; n++){
    
    int ind = order[n];
    
    for(int chr=0; chr<2; chr++){
      
      //cout << "templist= " << endl;
      //templist.Output(cout, coding);
      
      if((group < 0) || (groupvec[ind] == group)){
	double condprob = templist.FDLSProb(pop[ind].get_haplotype(chr), Qptr, nchr, rho, deriv, computederiv, false, vector<int>(0), vecTheta, 2*NindToUseForRho);
	// false means do not use quadrature; vector<int>(0) says no missing data (ie use complete data)

	// test computation of rho derivs
	//  cout << "Testing computation of conditional derivatives" << endl;
	//        vector<double> tempRho = rho;
	//        vector<double> tempRhoDeriv = deriv;
	//        vector<double> dummy; 
	//        for(int i=0; i<(Nloci-1); i++){
	//  	double epsilon = 1e-10;
	//  	tempRho = rho;
	//  	tempRho[i] += epsilon;
	//  	double newcondprob = templist.FDLSProb(pop[ind].get_haplotype(chr), Qptr, nchr,tempRho,dummy,false);
	//  	double approxderiv = (newcondprob - condprob)/epsilon;
	//  	cout << i << ": " << approxderiv << " " << tempRhoDeriv[i] << endl;
	//        }	
	// end of test

	//cout << "condprob = " << condprob << endl;
	if(condprob == 0)
	  cerr << "Warning: underflow problem in computation of logFDLSProb" << endl;
	  
	logprob += log(condprob);
	for(int i = 0; i< rhoderiv.size(); i++)
	  rhoderiv[i] += deriv[i]/condprob;
	templist.Add(pop[ind],chr,false,1.0);
	templist.MakePositiveHaps();      
	nchr++;
      }
    }
  }
  return logprob;
}


// propose change to RhoMean, and accept or reject it
// currently uses MH random walk on log(rhomean).
// parameterization of rho is 
// log(rho(i)) = log(RhoMean) + log(rhomult(i))
// or rho(i) = RhoMean * rhomult(i)
// sigma is standard deviation of MH proposal
bool ClassPop::updateRhoMeanRandomWalk(double sigma ,map<string,double> & d_cmds )
{
  double rwfactor = exp( rnorm(0,sigma) );

  if(RhoMean == 0)
    cerr << "Warning: estimate of recom rate reached zero" << endl;
 
  bool accept;

  if(((RhoMean*rwfactor)>MAXRHOMEAN) || ((RhoMean*rwfactor)<MINRHOMEAN) ){
    accept = false;
  } else {
    // set up newRho to be rwfactor * old rho values
    vector<double> newRho = vecRho;
    for(vector<double>::iterator r = newRho.begin(); r!=newRho.end(); r++)
      (*r) *= rwfactor;
    
    vector<double> dummy;
    double newlogprob = logFDLSProb(newRho, dummy, false);
    
    double lpriorprob = logpriorprobRhoMean(RhoMean, rwfactor*RhoMean, d_cmds["meanRhoMean"], d_cmds["sdRhoMean"]);
   
    accept = (ranf()<exp(newlogprob - CurrentLogProb + lpriorprob)) ;
    if(accept){
      vecRho = newRho;
      RhoMean *= rwfactor;
      CurrentLogProb = newlogprob;
    } 
  } 
    //cout << RhoMean << "," << CurrentLogProb << endl;
  return accept;
}


// propose change to RhoMean, and accept or reject it
// using MH Langevin algorithm on log(rhomean).
// parameterization of rho is 
// log(rho(i)) = log(RhoMean) + log(rhomult(i))
// or rho(i) = RhoMean * rhomult(i)
//
// sigma is sd of proposal
//
bool ClassPop::updateRhoMeanLangevin(double sigma, map<string,double> & d_cmds)
{ 
  if(RhoMean == 0)
    cerr << "Warning: estimate of recom rate reached zero" << endl;
  
  // start with derivative of log(prior) for OmegaMean = log(RhoMean)
  double deriv = (d_cmds["meanRhoMean"]- log(RhoMean)) / (d_cmds["sdRhoMean"] * d_cmds["sdRhoMean"]);

  // compute derivative of log likelihood with respect to OmegaMean = log(RhoMean)
  for(int k=0; k<(Nloci-1); k++)
    deriv += vecRho[k] * vecRhoDeriv[k];

  //cout << "deriv = " << deriv << endl;

  // set up newRho to be langevin factor * rwfactor * old rho values
  vector<double> newRho = vecRho;
  
  double epsilon = rnorm(0,sigma);
  double factor = exp(0.5*sigma*sigma*deriv + epsilon );

  bool accept;

  //cout << "Current value of RhoMean = " << RhoMean << endl;
  // cout << "Proposed value of RhoMean = " << RhoMean*factor << endl;
  

  if( ((RhoMean*factor)>MAXRHOMEAN) || ((RhoMean*factor)<MINRHOMEAN) ){
    accept = false;
  } else {

    for(vector<double>::iterator r = newRho.begin(); r!=newRho.end(); r++)
      (*r) *= factor;
 
    vector<double> newRhoDeriv = vecRhoDeriv;
    double newlogprob = logFDLSProb(newRho,newRhoDeriv,true);
    // compute new derivative of likelihood with respect to OmegaMean = log(RhoMean)
    double newderiv = (d_cmds["meanRhoMean"]- log(RhoMean*factor)) / (d_cmds["sdRhoMean"] * d_cmds["sdRhoMean"]);

    for(int k=0; k<(Nloci-1); k++)
      newderiv += newRho[k] * newRhoDeriv[k];
    
    // prob of going from logRhoMean to NewlogRhoMean, and back
    double logforwardsprob = logdnorm(epsilon,0,sigma);
    double backwardsepsilon = -log(factor) - 0.5*sigma*sigma*newderiv;
    double logbackwardsprob = logdnorm(backwardsepsilon,0,sigma);
    double logHastingsRatio = logbackwardsprob - logforwardsprob;
    
    // cout << "HastingsRatio:" << exp(logHastingsRatio) << endl;
    
    //cout << "NewLogProb " << newlogprob << endl;
    //cout << "CurrentLogProb " << CurrentLogProb << endl;
    //cout << "Likelihood Ratio:" << exp(newlogprob - CurrentLogProb) << endl;
    double lpriorprob = logpriorprobRhoMean(RhoMean, factor*RhoMean, d_cmds["meanRhoMean"], d_cmds["sdRhoMean"]);
   
    accept = (ranf()<(exp(logHastingsRatio+newlogprob - CurrentLogProb + lpriorprob))) ;
    if(accept){
      vecRho = newRho;
      vecRhoDeriv = newRhoDeriv;
      RhoMean *= factor;
      CurrentLogProb = newlogprob;
    } 
  }

  //cout << "Final value of RhoMean = " << RhoMean << endl;
  //cout << "Final logprob = " << CurrentLogProb << endl;
  return accept;
}

void ClassPop::ComputeRho(vector<double> & newright,vector<double> & newlambda,vector<double> & newleft)
{
  for(int i=0; i < (Nloci-1);i++){
    RhoMult[i] = 1;
    for(int h = 0; h< newlambda.size(); h++){
      if(newlambda[h] > 0){
	if((newleft[h] <= position[i]) && (newright[h] >= position[i+1])){ //whole of interval is in hotspot
	  RhoMult[i] += exp(newlambda[h]) - 1;
	}
	else if((newleft[h] >= position[i]) && (newright[h] <= position[i+1])){ // special case where r and l within a single gap
	  RhoMult[i] += (exp(newlambda[h]) - 1) * (newright[h] - newleft[h]) / (position[i+1] - position[i]);
	}
	else if((newleft[h] <= position[i]) && (newright[h] >= position[i])){ //  left part of interval is in hotspot
	  RhoMult[i] += (exp(newlambda[h]) - 1) * (newright[h] - position[i]) / (position[i+1] - position[i]);
	}
	else if((newleft[h] <= position[i+1]) && (newright[h] >= position[i+1])){ // right part of interval is in hotspot
	  RhoMult[i] +=  (exp(newlambda[h]) - 1) * (position[i+1] - newleft[h]) / (position[i+1] - position[i]);
	}
      }
    }
  }

  ComputeRho();
  //OutputRho(cout);
}

double ClassPop::EffectiveLength(vector<double> & l,vector<double> & r,vector<double> & lambda)
{
  double len = get_physical_length();
  for(int h = 0; h < lambda.size(); h++){    
    if(lambda[h]>0)
      len += (r[h] - l[h]) * (exp(lambda[h]) - 1);
  }
  //cout << "effective len:" << len << endl;

  return len;

}

bool ClassPop::AcceptOrRejectSimpleHotspot(vector<double> & newleft, vector<double> & newright, vector<double> & newlambda, bool fixedpos, double MAXLAMBDA, double Hotspotsd, double MinHotspotSize, double HOTSPOTRATE, double meanRhoMean, double sdRhoMean)
{
  bool accept = true;
  double lpriorratio = 0;

  double MINLAMBDA = -1;
  double PriorProbHotspot;
  
  if(fixedpos)
    PriorProbHotspot = 0.5;
  else
    PriorProbHotspot = 1 - exp(-get_physical_length() * HOTSPOTRATE);

  double newRhoMean = RhoMean * ( EffectiveLength(left,right,lambda) / EffectiveLength(newleft,newright,newlambda));
  double oldRhoMean = RhoMean;

  if((newRhoMean>MAXRHOMEAN) || (newRhoMean<MINRHOMEAN) )
    accept = false;
  else{
    for(int h = 0; h< lambda.size(); h++){
      double newcenter = 0.5*(newright[h]+newleft[h]);
      if(((newright[h] - newleft[h]) < MinHotspotSize) || (newlambda[h]>MAXLAMBDA) || (newlambda[h]<MINLAMBDA) || (newcenter < get_first_position()) || (newcenter > get_last_position()) ) // (newleft[h] < get_first_position()) || (newright[h] > get_last_position())
	accept = false;
    }
  }

  if(accept){
    RhoMean = newRhoMean; 
    ComputeRho(newright,newlambda,newleft);
   
    // take prior on RhoMean into account
    lpriorratio = logpriorprobRhoMean(oldRhoMean, RhoMean, meanRhoMean, sdRhoMean);
  
    for(int h = 0; h< lambda.size(); h++){
      if(newlambda[h] < 0)
	lpriorratio += log((MAXLAMBDA/(-MINLAMBDA)) * (1-PriorProbHotspot)); 
      else
	lpriorratio += log(PriorProbHotspot);
    
      if(lambda[h] < 0)
	lpriorratio -= log((MAXLAMBDA/(-MINLAMBDA)) * (1-PriorProbHotspot));
      else
	lpriorratio -= log(PriorProbHotspot);

    
      lpriorratio += -(0.5/(Hotspotsd* Hotspotsd)) * ((newright[h] - newleft[h])*(newright[h]- newleft[h]) - (right[h]-left[h])*(right[h]-left[h]));
    
    }
        
    vector<double> dummy;
    double newlogprob = logFDLSProb(vecRho, dummy, false);

    //cout << "Accept = " << exp(newlogprob - CurrentLogProb + lpriorratio) << endl;
    accept = (ranf()<exp(newlogprob - CurrentLogProb + lpriorratio)) ;
    
    if(accept){
      right = newright;
      lambda = newlambda;
      left = newleft;
      CurrentLogProb = newlogprob;
    } else {
      RhoMean = oldRhoMean;
      ComputeRho(right,lambda,left);
    }
  }
  return accept;
}

bool ClassPop::updateRhoSimpleHotspot(bool fixedpos, map<string,double> & d_cmds)
{
  vector<double> newleft(left);
  vector<double> newlambda(lambda);
  vector<double> newright(right);

  bool accept;
  double lpriorratio;
  
  double MAXLAMBDA = d_cmds["maxlambda"]; //log(100.0);
  double Hotspotsd = d_cmds["Hotspotsd"]; //2000; // standard deviation of hotspot width
  double MinHotspotSize = d_cmds["MinHotspotSize"]; //200; // minimum size of hotspot
  double HOTSPOTRATE = d_cmds["Hotspotrate"]; //1.0/50000; // one hotspot per 50kb.
  double meanRhoMean = d_cmds["meanRhoMean"];
  double sdRhoMean = d_cmds["sdRhoMean"];
  
  for(int h = 0; h<lambda.size(); h++){

    if(!fixedpos){ // update position of hotspot (left and right positions)
      if(lambda[h] < 0){ // sample position of left and right from prior
	do{
	  double center = get_first_position() + ranf() * (get_last_position()-get_first_position());
	  double width = abs(rnorm(0,Hotspotsd));
	  left[h] = center - 0.5*width;
	  right[h] = center + 0.5*width;
	} while( ((right[h] - left[h]) < MinHotspotSize) ); // || (left[h] < get_first_position()) || (right[h] > get_last_position()) );
      } else { // propose via MH random walk
	
	newleft[h] = left[h] + 1000 * (ranf()-0.5);
	newlambda[h] = lambda[h];
	newright[h] = right[h];
	
	AcceptOrRejectSimpleHotspot(newleft,newright,newlambda,fixedpos,MAXLAMBDA,Hotspotsd,MinHotspotSize,HOTSPOTRATE,meanRhoMean,sdRhoMean);
	
	
	newleft[h] = left[h];
	newlambda[h]= lambda[h];
	newright[h] = right[h] + 1000 * (ranf()-.5);
	
	AcceptOrRejectSimpleHotspot(newleft,newright,newlambda,fixedpos,MAXLAMBDA,Hotspotsd,MinHotspotSize,HOTSPOTRATE,meanRhoMean,sdRhoMean);
      }
    }
  

    // propose new value for lambda

    newleft[h] = left[h];
    newlambda[h] = lambda[h] + 2 * (ranf()-.5);
    newright[h] = right[h];
    
    AcceptOrRejectSimpleHotspot(newleft,newright,newlambda,fixedpos,MAXLAMBDA,Hotspotsd,MinHotspotSize,HOTSPOTRATE,meanRhoMean,sdRhoMean);
  
  }
  return 0;

}


// propose change to RhoMult, and accept or reject it
// using MH Langevin algorithm on log(rhomult).
// parameterization of rho is 
// log(rho(i)) = log(RhoMean) + log(rhomult(i))
// or rho(i) = RhoMean * rhomult(i)
//
// sigma is sd of proposal
//
bool ClassPop::updateRhoMultLangevin(double sigma)
{ 
// set up newRho to be langevin factor * rwfactor * old rho values
  vector<double> newRho = vecRho;
  vector<double> newRhoDeriv = vecRhoDeriv;
  vector<double> factor(Nloci-1);
 
// prob of going from logRhoMult to NewlogRhoMult, and back 
  double logforwardsprob = 0;
  double logbackwardsprob = 0;
  double logpriorratio =0;
  bool accept;

  for(int locus = 0; locus<(Nloci-1); locus++){
      // compute derivative of loglikelihood+log(prior) with respect to 
      //Omegai = log(RhoMult[i])  
    double deriv = vecRho[locus] * vecRhoDeriv[locus] - log(RhoMult[locus])/(RHOMULTSIGMA*RHOMULTSIGMA);
      
    //cout << "deriv = " << deriv << endl;

    double epsilon = rnorm(0,sigma);
    logforwardsprob += logdnorm(epsilon,0,sigma);

    factor[locus] = exp(0.5*sigma*sigma*deriv + epsilon );
    logpriorratio += logdnorm( log(factor[locus])+log(RhoMult[locus]), 0, RHOMULTSIGMA) - logdnorm(log(RhoMult[locus]), 0, RHOMULTSIGMA);
    
    newRho[locus] *= factor[locus];
  }

  double newlogprob = logFDLSProb(newRho,newRhoDeriv,true);

  // compute backwards probabilities
  for(int locus = 0; locus<(Nloci-1); locus++){
    // compute new derivative of loglikelihood +logprior 
    // with respect to OmegaMean = log(RhoMean)
    double newderiv = newRho[locus] * newRhoDeriv[locus] - log(RhoMult[locus]*factor[locus])/(RHOMULTSIGMA*RHOMULTSIGMA);
    double backwardsepsilon = -log(factor[locus]) - 0.5*sigma*sigma*newderiv;
    logbackwardsprob += logdnorm(backwardsepsilon,0,sigma);
  }

  double logHastingsRatio = logbackwardsprob - logforwardsprob;

  //cout << "HastingsRatio:" << exp(logHastingsRatio) << endl;

  //cout << "NewLogProb " << newlogprob << endl;
  //cout << "CurrentLogProb " << CurrentLogProb << endl;
  //cout << "Likelihood Ratio:" << exp(newlogprob - CurrentLogProb) << endl;

  accept = (ranf()<(exp(logpriorratio + logHastingsRatio+newlogprob - CurrentLogProb))) ;
  if(accept){
    vecRho = newRho;
    vecRhoDeriv = newRhoDeriv;
    for(int locus = 0; locus<(Nloci-1); locus++){
      RhoMult[locus] *= factor[locus];
    }
    CurrentLogProb = newlogprob;
  }
 
  //  for(int locus = 0; locus<(Nloci-1); locus++){
//      cout << RhoMult[locus] << " ";
//    }
//    cout << CurrentLogProb << endl;
  return accept;
}

void ClassPop::ComputeRho()
{
  for(int i = 0; i<(position.size()-1); i++){
    vecRho[i] = RhoMean * RhoMult[i] * (position[i+1]-position[i]);
  }
}

void ClassPop::ComputeRhoDerivAndCurrentLogProb()
{
  //TestComputeRhoDerivNaively();
  CurrentLogProb = logFDLSProb(vecRho,vecRhoDeriv,true);
  //cout << "Current Log Prob is " << CurrentLogProb << endl;
}

void ClassPop::TestComputeRhoDerivNaively()
{
  vector<double> tempRho = vecRho;
  vector<double> tempRhoDeriv = vecRhoDeriv;
  vector<double> dummy;

  double logprob = logFDLSProb(tempRho,tempRhoDeriv,true);

  for(int i=0; i<position.size(); i++){
    double epsilon1 = 1e-10;
    tempRho = vecRho;
    tempRho[i] += epsilon1;
    double newlogprob = logFDLSProb(tempRho,dummy,false);
    double approxderiv1 = (newlogprob - logprob)/epsilon1;

    double epsilon2 = 1e-5;
    tempRho = vecRho;
    tempRho[i] += epsilon2;
    newlogprob = logFDLSProb(tempRho,dummy,false);
    double approxderiv2 = (newlogprob - logprob)/epsilon2;



    cout << i << ": " << approxderiv1 << " " << approxderiv2 << " " << tempRhoDeriv[i] << endl;
  }

  // compute derivative of likelihood with respect to OmegaMean = log(RhoMean)
  double deriv = 0;
  double epsilon = 1e-5;
  tempRho = vecRho;
  for(int i=0; i<(position.size()-1); i++){
    deriv += vecRho[i] * tempRhoDeriv[i];
    tempRho[i] *= exp(epsilon);
  }
  double newlogprob = logFDLSProb(tempRho,dummy,false);
  double approxderiv = (newlogprob-logprob)/epsilon;

  cout << "Mean: " << approxderiv << " " << deriv << endl;



}


// sigma is sd of proposal
bool ClassPop::updateRhoMultRandomWalk(double sigma)
{

  // set up newRho to be rwfactor * old rho values
  vector<double> newRho = vecRho;
  vector<double> rwfactor(vecRho.size());

  double logpriorratio = 0; //value of prior ratio for rhoMult
  // currently assuming a normal(0,RHOMULTSIGMA) prior on log(rhoMult)

  bool accept;

  int locus = 0;
  
  for(vector<double>::iterator r = newRho.begin(); r!=newRho.end(); r++)
  {
    rwfactor[locus] = exp( rnorm(0,sigma) );
    (*r) *= rwfactor[locus];
    logpriorratio += logdnorm( log(rwfactor[locus])+log(RhoMult[locus]), 0, RHOMULTSIGMA) - logdnorm(log(RhoMult[locus]), 0, RHOMULTSIGMA);
    locus++;
  }

  vector<double> dummy;
  double newlogprob = logFDLSProb(newRho,dummy,false);
  
  accept = (ranf()<exp(logpriorratio + newlogprob - CurrentLogProb)) ;
  if(accept){
    vecRho = newRho;
    for(int locus = 0; locus<RhoMult.size(); locus++)
      RhoMult[locus] *= rwfactor[locus];
    CurrentLogProb = newlogprob;
  } 

 //   cout << "RhoMult: ";
//    for(it locus = 0; locus<RhoMult.size(); locus++)
//      cout << RhoMult[locus] << ",";
//    cout << endl;

  return accept;
}

void ClassPop::RandomiseOrder()
{
  random_shuffle(order.begin(),order.end());
}

// propose new order; accept or reject based on probability
void ClassPop::MHUpdateOrder()
{
  // propose switching two chosen at random
  int i = (int) (NindToUseForRho * ranf());
  int j = (int) (NindToUseForRho * ranf());
  
  // for(int k=0; k< order.size(); k++){
//     cout << order[k] << ",";
//   }
//   cout << endl;
 
  vector<int> oldorder(order);
  order[i] = oldorder[j];
  order[j] = oldorder[i];
  vector<double> newRhoDeriv = vecRhoDeriv;
  double NewLogProb;
  //if(RecomModel ==0)
  NewLogProb = logFDLSProb(vecRho,newRhoDeriv,true);
    //else
    // NewLogProb = logFDLSProb(vecRho,newRhoDeriv,false);

  double A = exp(NewLogProb - CurrentLogProb);
  if(ranf() > A) // if reject
    order = oldorder;
  else{
    CurrentLogProb = NewLogProb;
    vecRhoDeriv = newRhoDeriv;
  }

}

void ClassPop::OutputHotspotParams(ostream & ostr)
{
  ostr << RhoMean << " ";
  for(int h = 0; h< lambda.size(); h++){
    ostr << left[h] << " " << right[h] << " ";
    if(lambda[h]>0)
      ostr << exp(lambda[h]);
    else
      ostr << 1.0;
    ostr << " ";
  }
  ostr << endl;
}

void ClassPop::OutputRho(ostream & ostr)
{
  ostr << RhoMean << " ";
  for(int locus = 0; locus<(Nloci-1); locus++)
    ostr << RhoMult[locus] << " ";
  ostr << endl;
}

void ClassPop::OutputCurrentLogProb(ostream & ostr)
{
  ostr << CurrentLogProb << endl;
}

void ClassPop::UpdateRho(double sigmaRhoMean, double sigmaRhoMult, int & nacceptRhoMean, int & nacceptRhoMult,map<string,double> & d_cmds)
{
  if(RhoMean>0){

    switch(RecomModel){
      
    case 0: 
      nacceptRhoMean += updateRhoMeanLangevin( sigmaRhoMean, d_cmds );
      nacceptRhoMult += updateRhoMultLangevin(sigmaRhoMult);
      break;
      
    case 1:  
      nacceptRhoMean += updateRhoMeanRandomWalk( sigmaRhoMean, d_cmds );
      nacceptRhoMult += updateRhoSimpleHotspot ( false, d_cmds );
      break;
      
    case 2:  
      nacceptRhoMean += updateRhoMeanRandomWalk( sigmaRhoMean, d_cmds );
      nacceptRhoMult += updateRhoSimpleHotspot ( true, d_cmds ); // fixed hotspot position
      break;
      
    case 3:  // uniform recom rate
      nacceptRhoMean += updateRhoMeanRandomWalk( sigmaRhoMean, d_cmds );
      break;
      
    default: break;
    }
  }
  
}

// tune sigmaRhoMean and sigmaRhoMult so that acceptance probs not too small
double ClassPop::InferRho(int Niter, double & sigmaRhoMean, double & sigmaRhoMult, int verbose,map<string,double> & d_cmds)
{

  int nacceptRhoMean=0;
  int nacceptRhoMult=0;
  int Nrep = 10;  
  int ntry = 0;

  while((nacceptRhoMean*1.0/Nrep > 0.7 || 
	nacceptRhoMean*1.0/Nrep < 0.3 || 
	nacceptRhoMult*1.0/Nrep > 0.7 || 
	nacceptRhoMult*1.0/Nrep < 0.3) && (ntry++ < 100) ){

  
    nacceptRhoMean=0;
    nacceptRhoMult=0;

    for(int count =0; count < Nrep; count++){
      //TestComputeRhoDerivNaively();
      UpdateRho(sigmaRhoMean,sigmaRhoMult,nacceptRhoMean, nacceptRhoMult, d_cmds);
      if(verbose){
	OutputRho(cout);
	OutputCurrentLogProb(cout);
      }
    }
 
    if(verbose){
      cout << "Acceptance Rate for RhoMean: " << nacceptRhoMean*1.0/Nrep << endl;
      cout << "Acceptance Rate for RhoMult: " << nacceptRhoMult*1.0/Nrep << endl; 
      cout << "SigmaMean = " << sigmaRhoMean << endl;
      cout << "SigmaMult = " << sigmaRhoMult << endl;
      cout << "RhoMean = " << RhoMean << endl;
    }	

    if(nacceptRhoMean*1.0/Nrep < 0.3 ) // tune sigmaRhoMean
      sigmaRhoMean /=(1+ranf());
    if(nacceptRhoMean*1.0/Nrep > 0.7 )
      sigmaRhoMean *=(1+ranf());
    
    if(nacceptRhoMult*1.0/Nrep < 0.3 ) // tune sigmaRhoMult
      sigmaRhoMult /=(1+ranf());
    if(nacceptRhoMult*1.0/Nrep > 0.7 )
      sigmaRhoMult *=(1+ranf()); 
    
   
   
  }

  if(ntry>100)
    cerr << "Warning: failed to find decent estimate of recombination parameters" << endl;


  // for(int iter=0; iter<(2*Niter); iter++){
   
//     MHUpdateOrder(); //RandomiseOrder(); ComputeRhoDerivAndCurrentLogProb();
    
//     for(int count =0; count < Nrep; count++){
//       //TestComputeRhoDerivNaively();
//       UpdateRho(sigmaRhoMean,sigmaRhoMult,nacceptRhoMean, nacceptRhoMult, d_cmds);
//       if(verbose){
// 	OutputRho(cout);
// 	OutputCurrentLogProb(cout);
//       }
//     }
 
//     if(iter < Niter){

//       if(verbose){
// 	cout << "nacceptRhoMean = " << nacceptRhoMean << endl; 
// 	cout << "sigmaRhoMean = " << sigmaRhoMean << endl;
//       }

//       if(nacceptRhoMean*1.0/Nrep < 0.3 ) // tune sigmaRhoMean
// 	sigmaRhoMean /=(1+ranf());
//       if(nacceptRhoMean*1.0/Nrep > 0.7 )
// 	sigmaRhoMean *=(1+ranf());
      
//       if(nacceptRhoMult*1.0/Nrep < 0.3 ) // tune sigmaRhoMult
// 	sigmaRhoMult /=(1+ranf());
//       if(nacceptRhoMult*1.0/Nrep > 0.7 )
// 	sigmaRhoMult *=(1+ranf());

//       if(verbose){
// 	cout << "New sigmaRhoMean = " << sigmaRhoMean << endl;
//    }


  //nacceptRhoMean = 0;
  //nacceptRhoMult = 0;
  //}
  //}
  
  if(verbose){
    cout << "Acceptance Rate for RhoMean: " << nacceptRhoMean*1.0/Nrep << endl;
    cout << "Acceptance Rate for RhoMult: " << nacceptRhoMult*1.0/Nrep << endl; 
    cout << "SigmaMean = " << sigmaRhoMean << endl;
    cout << "SigmaMult = " << sigmaRhoMult << endl;
    cout << "RhoMean = " << RhoMean << endl;
  }
  return RhoMean;
    
}
    
// do the same as HapListMCMCResolvePhaseRemove,
// but remove buddies (= partners of parents in trios) at same time, 
// and use the buddies haplotypes when updating individual.
 
double ClassPop::BuddyHapListMCMCResolvePhaseRemove(map<string,int> & cmds, int Niter, int Nthin, int Nburn, map<string,double> & d_cmds, string filename, bool collectdata)
{

  ///Note on how missing data are dealt with:
  // First, when list is initialised (before calling this fn)
  // all possibilities for SNPs are added. (for Multiallelics,
  // all possibilities based on current imputed alleles are added).
  // Second, when list of possible pairs is made up, all pairs 
  // consistent with observed alleles are added.

  static int firsttime = true;

  int verbose = cmds["verbose"];
  int outputsample = cmds["outputsample"];
  int method = cmds["method"];
  int nperm = cmds["nperm"];
  
  RecomModel = cmds["RecomModel"];
 
  int testcasecontrol = (nperm>0);
  
  double rho = d_cmds["rhostart"];
  double betastart = d_cmds["betastart"];
  double betaend = d_cmds["betaend"];
  double minfreqpairoutput = d_cmds["output_minfreq"]; //min prob that a particular
  // pair has to have before being output as a possible pair
  double minfreqpair = d_cmds["minfreq"];
  
  if(collectdata){
    Niter *= cmds["finalrepmult"];
    Nthin *= cmds["finalrepmult"];
    Nburn *= cmds["finalrepmult"];
  }

  if(method == 'Q') // for Q method, use no recom, except on the final run
    if(collectdata)
      method = 'R';
    else
      method = 'S';
  
  cout << "Method = " << (char) method << endl;
  
  if(collectdata)
    cout << "Performing Final Set of Iterations... nearly there!" << endl;

  ofstream recomfile;
  ofstream monitorfile; 
  ofstream samplefile; // sample from posterior
  ofstream pairsfile; // list of high probability haplotype pairs 
                      //for each individual
  // ofstream finalfile; // output final haplist used
  ofstream signiffile; // output final haplist used
  ofstream hotfile; // hotspot info
  
  if(collectdata){
    string recomfilename = filename+"_recom";
    string monitorfilename = filename+"_monitor";
    string  hotfilename = filename +"_hotspot";
      

       // this to make sure results for multiple datasets are sent to a single file
    if(cmds["NDatasets"]>1 && !firsttime){
      recomfile.open(recomfilename.c_str(),ios::app);
      assure (recomfile, recomfilename );
      monitorfile.open(monitorfilename.c_str(),ios::app);
      assure ( monitorfile, monitorfilename );
      if(RecomModel ==1 || RecomModel ==2){
	hotfile.open(hotfilename.c_str(),ios::app);
	assure (hotfile, hotfilename);
      }
    } else {
      recomfile.open(recomfilename.c_str());
      assure (recomfile, recomfilename );
      monitorfile.open(monitorfilename.c_str());
      assure ( monitorfile, monitorfilename );
      if(RecomModel ==1 || RecomModel ==2){
	hotfile.open(hotfilename.c_str());
	assure (hotfile, hotfilename);
      }
      firsttime = false;
    }
    
    string pairsfilename = filename+"_pairs";
    pairsfile.open(pairsfilename.c_str());
    assure (pairsfile, pairsfilename);

//     string finalfilename = filename + "_final";
//     finalfile.open(finalfilename.c_str());
//     assure (finalfile, finalfilename);
    
    if(testcasecontrol){
      string signiffilename = filename + "_signif";
      //if(cmds["testing"])
      //signiffile.open(signiffilename.c_str(),ios::app);
      //else
      signiffile.open(signiffilename.c_str());
      assure (signiffile, signiffilename);
    }
    
  }

  if(outputsample && collectdata){
    string samplefilename = filename+"_sample";
    samplefile.open(samplefilename.c_str()); 
    assure ( samplefile, samplefilename);
  }

  haplist.ClearFreqs(); //
  haplist.Add(pop,Nind); // start with list of all haplotypes in current guess
    
  haplist.ClearPseudoCounts(); // PseudoCounts are used to store posterior mean
  //of freqs of each hap

  bool found;

  // Set up an index of the pairs in list that can make up each individual.
  // index[ind] contains a vector of  pairs of pointers to haps in the list 
  // that can make up individual ind
  vector < vector< pair< ListType::iterator,ListType::iterator> > > index(Nind); 
  haplist.MakePairsIndex(index,pop,true, false); 
  // first bool (true) also removes those haps from list
  // that cannot make up any individuals in pop
  // second bool (false) means that missing positions are not
  // checked for match

  if(verbose){
    cout << "Possible pairs " << endl;
    for(int ind = 0; ind< Nind; ind++){
      cout << "Individual:" << ind << endl;
      for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair = index[ind].begin(); hpair != index[ind].end(); hpair++){
	(*hpair).first->first.print_haplotype(cout,coding);    
	cout << ",";
	(*hpair).second->first.print_haplotype(cout,coding); 
	cout << endl;
      }
    }
  }
 
  vector < vector< double > > phaseprobs(Nind); // set up a matrix for storing the prob of each phase for each ind
  for(int ind = 0; ind< Nind; ind++){
    phaseprobs[ind] = vector<double>(index[ind].size(),0.0);
  }


  int nchr = Nind+Nind;
  int nchrminus2 = nchr-2;
     
  double removeamount = 1.0;
  double addamount = 1.0;

  double loglik;

  double SigmaMean=1;
  double SigmaMult=1;
  ComputeRho();

  // burn-in iterations
  cerr << "Performing Burn-in iterations" << endl;
  cerr << setw(4) << 0 << "% done\033[A" << endl;
  int nacceptRhoMean = 0;
  int nacceptRhoMult = 0;



  vector< vector<double> > childprob(Nind);
  // set up a probability of child's data for each individual pair
  vector<double> totalp(Nind,0); 
  for(int ind = 0; ind<Nind;ind++){
    for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair = index[ind].begin(); hpair != index[ind].end(); hpair++){
      for(vector<int>::iterator b=buddy[ind].begin(); b!=buddy[ind].end(); b++){
	for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair2 = index[*b].begin(); hpair2 != index[*b].end(); hpair2++){
	  double pc = pop[childindex[ind]].ObservedDataProbGivenParents((*hpair).first->first, (*hpair).second->first, (*hpair2).first->first, (*hpair2).second->first, pop[ind].get_recom(), pop[*b].get_recom())
	    + pop[childindex[ind]].ObservedDataProbGivenParents( (*hpair).second->first, (*hpair).first->first, (*hpair2).first->first, (*hpair2).second->first, pop[ind].get_recom(), pop[*b].get_recom())
	    + pop[childindex[ind]].ObservedDataProbGivenParents((*hpair).first->first, (*hpair).second->first,  (*hpair2).second->first, (*hpair2).first->first, pop[ind].get_recom(), pop[*b].get_recom())
	    + pop[childindex[ind]].ObservedDataProbGivenParents((*hpair).second->first,(*hpair).first->first,  (*hpair2).second->first, (*hpair2).first->first, pop[ind].get_recom(), pop[*b].get_recom());
	  childprob[ind].push_back(pc);
	  totalp[ind] += (pc>0);
	
	  if(verbose){
	    cout << "Individual:" << ind << endl;
	    (*hpair).first->first.print_haplotype(cout,coding);    
	    cout << endl;
	    (*hpair).second->first.print_haplotype(cout,coding); 
	    cout << endl;
	    (*hpair2).first->first.print_haplotype(cout,coding);     
	    cout << endl;
	    (*hpair2).second->first.print_haplotype(cout,coding);  
	    cout << endl;
	    cout << "Childprob:" <<  pc << endl;
      	  }

	  
	}
      }
    }
    //    cout << "Total number of possibiliites =" << totalp[ind] << endl;
    // cout << "Sizes = " << index[ind].size() << "," << index[buddy[ind][0]].size() << endl;

    if(totalp[ind]==0){
      cerr << "Error: Program has reached state where child is incompatible with all parental possibilities; possibly due to Mendelian error in input file?" << endl;
      cerr << "Individual " << pop[ind].get_id() << endl;
      cerr << "partner is : " << pop[buddy[ind][0]].get_id() << endl;
      cerr << "child is: " << pop[childindex[ind]].get_id() << endl;
      exit(1);
    }
  }

  for(int iter = 0; iter < Nburn; iter++){
    
    if(method == 'R' && (RhoMean>0)){
      ComputeRhoDerivAndCurrentLogProb();
      MHUpdateOrder();
      
      if((iter == Nburn/2)){
	cerr << endl;
	cerr << "Estimating recom rates" << endl;
	if(verbose)
	  cout << "Estimating recom rates" << endl;
	InferRho(Niter/10,SigmaMean,SigmaMult,verbose, d_cmds); 
	cerr << "Continuing Burn-in" << endl;
      }
              
      if(iter > Nburn/2){
	UpdateRho(SigmaMean,SigmaMult,nacceptRhoMean, nacceptRhoMult, d_cmds);
	
	if(verbose){
	  OutputRho(cout);
	  OutputCurrentLogProb(cout);
	}
      }
 
    } else {
      RandomiseOrder();
    }
    
    loglik = 0;

    //double DPRIOR = (betaend/nchr) + ((Nburn-1-iter)*1.0/(Nburn-1))*((betastart-betaend)/nchr);
    double DPRIOR = (betaend) + ((Nburn-1-iter)*1.0/(Nburn-1))*((betastart-betaend));    
    
    for(int i = 0; i<Nind; i++){    
      int ind = order[i];     
      haplist.SoftRemove(pop[ind],removeamount);
      for(vector<int>::iterator b=buddy[ind].begin(); b!=buddy[ind].end(); b++){
	haplist.SoftRemove(pop[*b],removeamount);
      }
      haplist.MakePositiveHaps();

#ifdef DEBUG
      cout << "Removed ind " << ind << endl;
      cout << "Freqs" << endl;
      haplist.Output(cout, coding);
#endif

      double sumprob = 0;
      vector<double> Prob;
      bool addprobs = (ranf()<DPRIOR); // here DPRIOR is prob of adding,
      // rather than multiplying probs (adding gives more weight to
      // configurations where only one of the 2 haps has high prob
      double prob1, prob2;
     
      //if(index[ind].size()>1){
	vector<double>::iterator  pc = childprob[ind].begin();
	for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair = index[ind].begin(); hpair != index[ind].end(); hpair++){
	  double p;
	  vector<double> buddyprob1, buddyprob2, childlike;
	  double sumchildlike = 0;
      	  for(vector<int>::iterator b=buddy[ind].begin(); b!=buddy[ind].end(); b++){
	    for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair2 = index[*b].begin(); hpair2 != index[*b].end(); hpair2++){
	      if(*pc>0){ // if this pair is consistent with child's data
		buddyprob1.push_back(haplist.CalcProb(((*hpair2).first->first),method,Qptr,nchrminus2,vecRho,0, vector<int>(0), ALLSNP, vecdiffprob));
		buddyprob2.push_back(haplist.CalcProb(((*hpair2).second->first),method,Qptr,nchrminus2,vecRho,0, vector<int>(0), ALLSNP, vecdiffprob)); // pop[ind].get_nmissing() as last param would just check missing positions
		childlike.push_back(*pc);
	      } else{
		buddyprob1.push_back(0);
		buddyprob2.push_back(0);
		childlike.push_back(0);
	      }
	      sumchildlike += *pc;
	      pc++;
	    }
	  }
	
	  prob1 = 0;
	  prob2 = 0;
	  if(sumchildlike>0){
	    prob1 = haplist.CalcProb(((*hpair).first->first),method,Qptr,nchrminus2,vecRho,0, vector<int>(0), ALLSNP, vecdiffprob); // pop[ind].get_nmissing() as last param would just check missing positions
	    prob2 = haplist.CalcProb(((*hpair).second->first),method,Qptr,nchrminus2,vecRho,0, vector<int>(0), ALLSNP, vecdiffprob); // 
	  }
	  
	  for(int i = 0; i< buddyprob1.size(); i++){
	    if(addprobs)
	      p = childlike[i] * (prob1 + prob2 + buddyprob1[i] + buddyprob2[i]);
	    else
	      p = childlike[i] * prob1 * prob2 * buddyprob1[i] * buddyprob2[i];
	    Prob.push_back(p);
	    sumprob += p;
	  }
	
	  
#ifdef DEBUG       
	  ((*hpair).first->first).print_haplotype(cout,coding);
	cout << ",";
	((*hpair).second->first).print_haplotype(cout,coding);
	cout << ":" << p << endl;
#endif
	}
	//}	
      // else{
// 	Prob.push_back(1);
// 	sumprob = 1;
//       }


	//SOME OTHER PRIOR ANNEALING SCHEMES WE TRIED
      // original linear prior annealing (adding DPRIOR to p)
      //double p = (CalcProb(((*hpair).first->first),method,Qptr,nchrminus2,vecRho,0)+DPRIOR)/(1+DPRIOR*haplist.get_listlength()) * 
      //(CalcProb(((*hpair).second->first),method,Qptr,nchrminus2,vecRho,0)+DPRIOR)/(1+DPRIOR*haplist.get_listlength());
	
      // "powered" version of prior annealing, raising p to power DPRIOR
      //double p = exp(DPRIOR * log(CalcProb(((*hpair).first->first),method,Qptr,nchrminus2,vecRho,0)*
      //CalcProb(((*hpair).second->first),method,Qptr,nchrminus2,vecRho,0)));
      
      // store phase probabilities in phaseprobs -only necessary in order to do the later pruning...
      //if(index[ind].size()>1){
	vector<double>::iterator probpointer = Prob.begin();
	for (vector<double>::iterator phase1 = phaseprobs[ind].begin(); phase1 != phaseprobs[ind].end(); phase1++){
	  for(vector<int>::iterator b=buddy[ind].begin(); b!=buddy[ind].end(); b++){
	    for (vector<double>::iterator phase2 = phaseprobs[*b].begin(); phase2 != phaseprobs[*b].end(); phase2++){
	      *phase1 += 0.5*(1.0/Nburn)*(*probpointer/sumprob);
	      *phase2 += 0.5*(1.0/Nburn)*(*probpointer/sumprob);	    
	      probpointer++;
	    }
	  }
	}
      // }
//       else{
// 	phaseprobs[ind][0] = 1.0;
//       }

      int choice = rint2(Prob,sumprob);
      loglik += log(sumprob);
      //loglik += log(Prob[choice]); // a different version of the pseudologlik
      vector< pair<ListType::iterator,ListType::iterator> >::iterator sampledpair1 = index[ind].begin()+ choice / index[buddy[ind][0]].size();	
      vector< pair<ListType::iterator,ListType::iterator> >::iterator sampledpair2 = index[buddy[ind][0]].begin() + choice % index[buddy[ind][0]].size();	
      
      
      pop[buddy[ind][0]].update_haplotypes((*sampledpair2).first->first,(*sampledpair2).second->first);
      pop[ind].update_haplotypes((*sampledpair1).first->first,(*sampledpair1).second->first);

      //HapListImputeMissing(ind);  
      //HapListImputeMissing(buddy[ind][0]);

      for(int chr = 0; chr <2 ; chr++){
	bool isnewhap;
	ListType::iterator newrecord = haplist.Add(pop[ind], chr,addamount, isnewhap);
	// extend index just for the individual ind
	// if it has missing loci, as may have added haplotype pair
	// that wasn't on its list of allowed pairs
	if( pop[ind].NMissingLoci() > 0)
	  haplist.ExtendPairsIndex(newrecord, index[ind], phaseprobs[ind], pop[ind], true);

	if(isnewhap)
	  haplist.ExtendPairsIndex(newrecord, index, phaseprobs, pop, true);   
	// true param says to check missing positions
      }
      
      // update buddys
      for(vector<int>::iterator b=buddy[ind].begin(); b!=buddy[ind].end(); b++){
	for(int chr = 0; chr <2 ; chr++){
	  bool isnewhap;
	  ListType::iterator newrecord = haplist.Add(pop[*b], chr,addamount, isnewhap);
	   // extend index just for the individual ind
	   // if it has missing loci, as may have added haplotype pair
	   // that wasn't on its list of allowed pairs
	   if( pop[*b].NMissingLoci() > 0)
	     haplist.ExtendPairsIndex(newrecord, index[*b], phaseprobs[*b], pop[*b], true);
	   
	   if(isnewhap)
	     haplist.ExtendPairsIndex(newrecord, index, phaseprobs, pop, true);   
	   // true param says to check missing positions
	 }
      }
      

    
#ifdef DEBUG
      cout << "Adding ind:" << endl;
      haplist.Output(cout,coding);
#endif
      
    }

    cerr << setw(4) << (int) ((iter+1)*100.0/Nburn) << "% done\033[A" << endl;
    //cout << "Log Likelihood: " << loglik << endl;    
  }
  
  if(verbose){
      cout << "Acceptance Rate for RhoMean, burnin: " << nacceptRhoMean*2.0/(Nburn) << endl;
      cout << "Acceptance Rate for RhoMult, burnin:: " << nacceptRhoMult*2.0/(Nburn) << endl; 
      cout << "SigmaMean = " << SigmaMean << endl;
      cout << "SigmaMult = " << SigmaMult << endl;
  }

  
  // after burn-in, if necessary get better estimate of rho and appropriate sigmas
  if(nacceptRhoMean *2.0/Nburn < 0.1 || nacceptRhoMean *2.0/Nburn >0.9 || nacceptRhoMult *2.0/Nburn < 0.1 || nacceptRhoMult *2.0/Nburn >0.9){
    if(method == 'R' && (RhoMean >0)){
      cerr << "Re-Estimating recom rates" << endl;
      if(verbose)
	cout << "Estimating recom rates" << endl;
      InferRho(Niter/10,SigmaMean,SigmaMult,verbose, d_cmds); 
    }
  }

  nacceptRhoMean=nacceptRhoMult=0;

  // actual iterations
  double meanloglik = 0;

  // set up vecpermcc as vector of permutations of casecontrol labels    
  vector<double> meanpermloglik(nperm , 0); // mean log lik for controls

  vector< vector<int> > vecpermcc;

  if(testcasecontrol && collectdata){
    vecpermcc = vector<vector<int> > (nperm, casecontrol);
    for(int i = 1; i<nperm; i++){
      random_shuffle(vecpermcc[i].begin(),vecpermcc[i].end());
    }
  }

  
  // Here we prune index to remove very unlikely hap pairs
  // (also sets phaseprobs to be 0)
  haplist.PrunePairsIndex(index, phaseprobs,pop,minfreqpair*0.5);
  childprob = vector< vector<double> >(Nind);

  for(int ind = 0; ind<Nind;ind++){
    totalp[ind] = 0;
    for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair = index[ind].begin(); hpair != index[ind].end(); hpair++){
      for(vector<int>::iterator b=buddy[ind].begin(); b!=buddy[ind].end(); b++){
	for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair2 = index[*b].begin(); hpair2 != index[*b].end(); hpair2++){
	  double pc = pop[childindex[ind]].ObservedDataProbGivenParents((*hpair).first->first, (*hpair).second->first, (*hpair2).first->first, (*hpair2).second->first, pop[ind].get_recom(), pop[*b].get_recom())
	    + pop[childindex[ind]].ObservedDataProbGivenParents( (*hpair).second->first, (*hpair).first->first, (*hpair2).first->first, (*hpair2).second->first, pop[ind].get_recom(), pop[*b].get_recom())
	    + pop[childindex[ind]].ObservedDataProbGivenParents((*hpair).first->first, (*hpair).second->first,  (*hpair2).second->first, (*hpair2).first->first, pop[ind].get_recom(), pop[*b].get_recom())
	    + pop[childindex[ind]].ObservedDataProbGivenParents((*hpair).second->first,(*hpair).first->first,  (*hpair2).second->first, (*hpair2).first->first, pop[ind].get_recom(), pop[*b].get_recom());
	  childprob[ind].push_back(pc);
	  totalp[ind] += (pc>0);
	

	  if(verbose){
	    cout << "After pruning" << endl;
	    cout << "Individual:" << ind << endl;
	    cout << "Buddy: " << buddy[ind][0] << endl;
	    (*hpair).first->first.print_haplotype(cout,coding);    
	    cout << endl;
	    (*hpair).second->first.print_haplotype(cout,coding); 
	    cout << endl;
	    (*hpair2).first->first.print_haplotype(cout,coding);     
	    cout << endl;
	    (*hpair2).second->first.print_haplotype(cout,coding);  
	    cout << endl;
	    cout << "Childprob:" <<  pc << endl;
	  }
	}
      }
    }

    if(totalp[ind]==0){
      cerr << "Error: Program has reached state where child is incompatible " << endl;
      cerr << "with all parental possibilities" << endl;
      cerr << "Individual " << pop[ind].get_id() << endl;
      cerr << "partner is : " << pop[buddy[ind][0]].get_id() << endl;
      cerr << "child is: " << pop[childindex[ind]].get_id() << endl;
      exit(1);
    }
  }

  // set up group frequency vectors
  int ngroups = 0;
  for(int ind = 0; ind<Nind; ind++)
    if(casecontrol[ind]>ngroups)
      ngroups = casecontrol[ind];
  ngroups++;    
  haplist.SetupGroupFreqs(ngroups);
  
  // set up group size vector
  groupsize = vector<int>(ngroups,0);
  for(int ind = 0; ind<Nind; ind++)
    groupsize[casecontrol[ind]] += 2;


  cerr << "Performing Main iterations" << endl;
  cerr << setw(4) << 0 << "% done\033[A" << endl; 
   
  //int NumberOfSqPseudoCountUpdates = 0;
  
  haplist.ClearPseudoCounts(); // PseudoCounts are used to store posterior mean
  //of freqs of each hap

  if(method =='R' && collectdata){
   OutputPositions(recomfile);
  }

  for(int iter = 0; iter < Niter; iter++){
    loglik = 0;
      
    if(method == 'R' && (RhoMean>0)){
      ComputeRhoDerivAndCurrentLogProb();
      MHUpdateOrder();
      
      
      UpdateRho(SigmaMean,SigmaMult, nacceptRhoMean, nacceptRhoMult, d_cmds);
      if(verbose){
	OutputRho(cout);
	OutputCurrentLogProb(cout);
      }
      if(collectdata){
	OutputRho(recomfile);
	if(RecomModel == 1 || RecomModel ==2)
	  OutputHotspotParams(hotfile);
      }
    } else {
      RandomiseOrder();
    }

    if(testcasecontrol && collectdata){
      for(int i = 0; i< nperm; i++){
	meanpermloglik[i] += 
	  logProb( method, vecRho, vecpermcc[i], 1.0, 0) + 
	  logProb( method, vecRho, vecpermcc[i], 1.0, 1);
	// add likelihoods for group 0 and group 1
      }
    }
   
    for(int i = 0; i<Nind; i++){
      int ind = order[i];
      //cout << "ind = " << ind << endl;

      haplist.SoftRemove(pop[ind],removeamount);
      for(vector<int>::iterator b=buddy[ind].begin(); b!=buddy[ind].end(); b++){
	haplist.SoftRemove(pop[*b],removeamount);
      }
      haplist.MakePositiveHaps();

#ifdef DEBUG
      cout << "Removed ind " << ind << endl;
      cout << "Freqs" << endl;
      haplist.Output(cout, coding);
#endif

      // make list of probabilities of each hap pair
      double sumprob = 0;
      vector<double> Prob;
      double prob1, prob2;


      //if(index[ind].size()>1){
	vector<double>::iterator pc = childprob[ind].begin();
	for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair = index[ind].begin(); hpair != index[ind].end(); hpair++){
	  double p;
	  vector<double> buddyprob1, buddyprob2, childlike;
	  double sumchildlike = 0;
      	  for(vector<int>::iterator b=buddy[ind].begin(); b!=buddy[ind].end(); b++){
	    for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair2 = index[*b].begin(); hpair2 != index[*b].end(); hpair2++){
	      if(*pc>0){ // if this pair is consistent with child's data
		buddyprob1.push_back(haplist.CalcProb(((*hpair2).first->first),method,Qptr,nchrminus2,vecRho,0, vector<int>(0), ALLSNP, vecdiffprob));
		buddyprob2.push_back(haplist.CalcProb(((*hpair2).second->first),method,Qptr,nchrminus2,vecRho,0, vector<int>(0), ALLSNP, vecdiffprob)); // pop[ind].get_nmissing() as last param would just check missing positions
		childlike.push_back(*pc);
		double p = *pc * ((2-((*hpair2).first == (*hpair2).second)));
		for(int locus = 0; locus< Nloci; locus++){
		  if(pop[*b].n_missing(locus) == 1)
		    if(((*hpair2).first->first).get_allele(locus)!=((*hpair2).second->first).get_allele(locus))
		      p *= 0.5;
		}  
	      } else{
		buddyprob1.push_back(0);
		buddyprob2.push_back(0);
		childlike.push_back(0);
	      }
	      sumchildlike += *pc;
	      pc++;
	    }
	  }
	
	  prob1 = 0;
	  prob2 = 0;
	  if(sumchildlike>0){
	    prob1 = haplist.CalcProb(((*hpair).first->first),method,Qptr,nchrminus2,vecRho,0, vector<int>(0), ALLSNP, vecdiffprob); // pop[ind].get_nmissing() as last param would just check missing positions
	    prob2 = haplist.CalcProb(((*hpair).second->first),method,Qptr,nchrminus2,vecRho,0, vector<int>(0), ALLSNP, vecdiffprob); // 
	  }
	  
	  for(int i = 0; i< buddyprob1.size(); i++){
	    //cout << "i= " << i << endl;
	    double factor = (2-((*hpair).first == (*hpair).second)); // take acount of factor of 2 for case where there are missing positions
	    double p = factor * childlike[i] * prob1 * prob2 * buddyprob1[i] * buddyprob2[i];
	    // now account for Pr(gen | hap) at those positions where
	    // only one missing
	    for(int locus = 0; locus< Nloci; locus++){
	      if(pop[ind].n_missing(locus) == 1)
		if(((*hpair).first->first).get_allele(locus)!=((*hpair).second->first).get_allele(locus))
		  p *= 0.5;
	    }
	    Prob.push_back(p);
	    sumprob += p;
#ifdef DEBUG       
	    ((*hpair).first->first).print_haplotype(cout,coding);
	    cout << ",";
	    ((*hpair).second->first).print_haplotype(cout,coding);
	  
	    //cout << ",";
	    //((*hpair2).first->first).print_haplotype(cout,coding);
	    //cout << ",";
	    //((*hpair2).second->first).print_haplotype(cout,coding);
	    
	    cout << ":" << prob1 << "," << prob2 << "," << childlike[i] << "," << buddyprob1[i] << "," << buddyprob2[i]  << endl;
	    
#endif
	  
	  }	  
	}
      // }	
//       else{
// 	Prob.push_back(1);
// 	sumprob = 1;
//       }

      int probpos = 0;
   
      // store phase probabilities in phaseprobs 
      //if(index[ind].size()>1){
	vector<double>::iterator probpointer = Prob.begin();
	for (vector<double>::iterator phase1 = phaseprobs[ind].begin(); phase1 != phaseprobs[ind].end(); phase1++){
	  for(vector<int>::iterator b=buddy[ind].begin(); b!=buddy[ind].end(); b++){
	    for (vector<double>::iterator phase2 = phaseprobs[*b].begin(); phase2 != phaseprobs[*b].end(); phase2++){
	      *phase1 += (1.0/Niter)*(*probpointer/sumprob);
	      //*phase2 += (1.0/Niter)*(*probpointer/sumprob);	    
	      probpointer++;
	    }
	  }
	}
      // }
//       else{
// 	phaseprobs[ind][0] = 1;
//       }
      
      int choice = rint2(Prob,sumprob);
      loglik += log(sumprob);
      //loglik += log(Prob[choice]); // a different version of the pseudologlik
      vector< pair<ListType::iterator,ListType::iterator> >::iterator sampledpair1 = index[ind].begin()+ choice / index[buddy[ind][0]].size();	
      vector< pair<ListType::iterator,ListType::iterator> >::iterator sampledpair2 = index[buddy[ind][0]].begin() + choice % index[buddy[ind][0]].size();	
      
      
      // update buddy's haplotypes

      pop[buddy[ind][0]].update_haplotypes((*sampledpair2).first->first,(*sampledpair2).second->first);
      
      for(vector<int>::iterator b=buddy[ind].begin(); b!=buddy[ind].end(); b++){
	for(int chr = 0; chr <2 ; chr++){
	  bool isnewhap;
	  ListType::iterator newrecord = haplist.Add(pop[*b], chr,addamount, isnewhap);
	   // extend index just for the individual ind
	   // if it has missing loci, as may have added haplotype pair
	   // that wasn't on its list of allowed pairs
	   if( pop[*b].NMissingLoci() > 0)
	     haplist.ExtendPairsIndex(newrecord, index[*b], phaseprobs[*b], pop[*b], true);
	   
	   if(isnewhap)
	     haplist.ExtendPairsIndex(newrecord, index, phaseprobs, pop, true);   
	   // true param says to check missing positions
	}
      }


      // update estimated allele freqs

      for(ListType::iterator h = haplist.haplist.begin(); h != haplist.haplist.end(); h++){
	double hapcount = h->second.Freq;
	double grouphapcount = nchr-2;
	if(ngroups>1)
	  grouphapcount = (double) GetGroupCount(h->first,casecontrol[ind],ind);
	int groupn = groupsize[casecontrol[ind]];
	
	//cout << "nchr = " << nchr << endl;
	//cout << "hapcount = " << hapcount << endl;
	
	h->second.PseudoCount += (hapcount /nchr); // pseudocount stores posterior means of overall frequencies (bit of a misnomer for historical reasons! should change...)	 
	
	h->second.GroupFreq[casecontrol[ind]] += (grouphapcount / groupn);
	
	if(collectdata){
	  h->second.SqPseudoCount += (hapcount * hapcount)/(nchr*nchr);
	  h->second.GroupFreqSq[casecontrol[ind]] += (grouphapcount * grouphapcount)/(groupn*groupn) ;
	}
	
	probpos = 0;
	for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair = index[ind].begin(); hpair != index[ind].end(); hpair++){
	  for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair2 = index[buddy[ind][0]].begin(); hpair2 != index[buddy[ind][0]].end(); hpair2++){
	 
	    int indicator1 = (*hpair).first->first == (h->first);
	    int indicator2 = (*hpair).second->first == (h->first);
	    
	    //cout << "Prob=" << Prob[probpos]/sumprob << endl;
	    //cout << "tot = " << indicator1 + indicator2 << endl;
	    h->second.PseudoCount += (Prob[probpos]/sumprob) * (1.0/nchr) * (indicator1 + indicator2);
	    h->second.GroupFreq[casecontrol[ind]] += (Prob[probpos]/sumprob) * (1.0/groupn) * (indicator1 + indicator2);
	    
	    if(collectdata){
	      h->second.SqPseudoCount += (Prob[probpos]/sumprob) * (1.0/(nchr*nchr)) * (indicator1 + indicator2) *
		(indicator1 + indicator2 + 2*hapcount);
	      h->second.GroupFreqSq[casecontrol[ind]] += (Prob[probpos]/sumprob) * (1.0/(groupn*groupn)) * (indicator1 + indicator2) *
		(indicator1 + indicator2 + 2*grouphapcount);	    
	    }
	    probpos++;
	  } 
	}
      }
      
      pop[ind].update_haplotypes((*sampledpair1).first->first,(*sampledpair1).second->first);
     
      // this bit allows missing positions to be reconstructed
      // without reference to the HapList;
      // 
      // is not quite correct for case where one allele may be missing
      // 
      //HapListImputeMissing(ind);
//        // add new hap to list; extend Pairs index if necessary
      for(int chr = 0; chr <2 ; chr++){
	bool isnewhap;
	ListType::iterator newrecord = haplist.Add(pop[ind], chr,addamount, isnewhap);
	// extend index just for the individual ind
	// if it has missing loci, as may have added haplotype pair
	// that wasn't on its list of allowed pairs
	if(pop[ind].NMissingLoci()>0)
	   haplist.ExtendPairsIndex(newrecord, index[ind], phaseprobs[ind], pop[ind], true);
	if(isnewhap)
	  haplist.ExtendPairsIndex(newrecord, index, phaseprobs, pop, true);
      }
      
     
      

#ifdef DEBUG
      cout << "Adding ind:" << endl;
      haplist.Output(cout,coding);
#endif
      
      //pop[ind].UpdateCounts(); 

    }
   
    if(collectdata){
      monitorfile << loglik << " " << CurrentLogProb <<  endl; 
      if(outputsample)
	output(samplefile, true, true);
    }
    
    cerr << setw(4) << (int) ((iter+1)*100.0/(Niter)) << "% done\033[A" << endl; 
    meanloglik +=loglik;
  }
  
  if(verbose){
      cout << "Acceptance Rate for RhoMean, main: " << nacceptRhoMean*1.0/(Niter) << endl;
      cout << "Acceptance Rate for RhoMult, main: " << nacceptRhoMult*1.0/(Niter) << endl; 
      cout << "SigmaMean = " << SigmaMean << endl;
      cout << "SigmaMult = " << SigmaMult << endl;
  }
  
  vecpermcc = vector<vector<int> > (0);

  if(testcasecontrol && collectdata){
    double nbigger = 0;
    //cout << "meanpermloglik[0] = " << meanpermloglik[0] << endl;
    for(int i = 1; i<nperm; i++){
      nbigger += (meanpermloglik[i]>=meanpermloglik[0]);
      //if(meanpermloglik[i]==meanpermloglik[0])
      //nbigger += 0.5;
      //cout << "meanpermloglik[i] = " << meanpermloglik[i] << endl;
    }
    signiffile << "P-value for testing H0: Cases ~ Controls = " << (1.0*(nbigger+1))/(nperm) << endl;
  }

// copy PsedoCounts (posterior means of Freqs) to Freqs (suitably normalised)
  haplist.CopyPseudoCountsToFreqs();
  haplist.NormaliseFreqs();
  
  if(collectdata){
    haplist.NormaliseSqPseudoCounts(Nind * Niter); 

    haplist.NormaliseGroupFreqs();

    cerr << "Writing output to files " << endl;
    // Output list of plausible pairs for each individual 
    for(int ind = 0; ind<Nind; ind++){
      pairsfile << "IND: " << pop[ind].get_id() << endl;
      for(int j = 0; j< index[ind].size(); j++){
	if(phaseprobs[ind][j] > minfreqpairoutput){
	  haplist.OutputPair(pairsfile,index[ind][j],coding);
	  pairsfile << " , ";
	  pairsfile.setf(ios::fixed);
	  pairsfile.setf(ios::showpoint);
	  pairsfile.precision(3); 
	  pairsfile << phaseprobs[ind][j] << endl;
	}
      }
    }
    
    // Output final haplist used
    
    // haplist.Dump(finalfile);

    haplist.SetBestGuesses(pop,index,phaseprobs); // set all individuals in pop to their best guesses
    
    cerr << "Producing Summary, please wait " << endl;
    // each element of this vector is a summary for a single individual
    // the false at the end forces it to produce a single summary with no splits
    vector<Summary> summary = haplist.ProduceSummary(index,phaseprobs,0,get_nloci(),pop,false);
    
    // for(vector<Summary>::iterator s = summary.begin(); s<summary.end(); s++)
    //  s->Output(cout, coding);
    
    TransferCounts(summary);

    // produce separate lists of freqs for cases and controls
    // int ngroups = 0;
//     for(int ind = 0; ind<Nind; ind++)
//       if(casecontrol[ind]>ngroups)
// 	ngroups = casecontrol[ind];
//     ngroups++;
    
//     haplist.SetupGroupFreqs(ngroups);
//     for(int ind = 0; ind<Nind; ind++){
//       int probpos = 0;
//       for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair = index[ind].begin(); hpair != index[ind].end(); hpair++){  
// 	haplist.haplist[(*hpair).first->first].GroupFreq[casecontrol[ind]] += phaseprobs[ind][probpos]/(2*Nind);
// 	haplist.haplist[(*hpair).second->first].GroupFreq[casecontrol[ind]] += phaseprobs[ind][probpos]/(2*Nind);	  
// 	probpos++;
//       }
//     }
  }


  return meanloglik;

}

double ClassPop::HapListMCMCResolvePhaseRemove(map<string,int> & cmds, int Niter, int Nthin, int Nburn, map<string,double> & d_cmds, string filename, bool collectdata)
{

  ///Note on how missing data are dealt with:
  // First, when list is initialised (before calling this fn)
  // all possibilities for SNPs are added. (for Multiallelics,
  // all possibilities based on current imputed alleles are added).
  // Second, when list of possible pairs is made up, all pairs 
  // consistent with current imputed alleles are added.

  static int firsttime = true;
  
  int verbose = cmds["verbose"];
  int outputsample = cmds["outputsample"];
  int method = cmds["method"];
  int nperm = cmds["nperm"];
  
  RecomModel = cmds["RecomModel"];
 
  int testcasecontrol = (nperm>0);
  
  double rho = d_cmds["rhostart"];
  double betastart = d_cmds["betastart"];
  double betaend = d_cmds["betaend"];
  double minfreqpairoutput = d_cmds["output_minfreq"]; //min prob that a particular
  // pair has to have before being output as a possible pair
  double minfreqpair = d_cmds["minfreq"];
  
  if(collectdata){
    Niter *= cmds["finalrepmult"];
    Nthin *= cmds["finalrepmult"];
    Nburn *= cmds["finalrepmult"];
  }

  if(method == 'Q') // for Q method, use no recom, except on the final run
    if(collectdata)
      method = 'R';
    else
      method = 'S';
  
  cout << "Method = " << (char) method << endl;
  
  if(collectdata)
    cout << "Performing Final Set of Iterations... nearly there!" << endl;

  ofstream recomfile;
  ofstream monitorfile; 
  ofstream samplefile; // sample from posterior
  ofstream pairsfile; // list of high probability haplotype pairs 
                      //for each individual
  // ofstream finalfile; // output final haplist used
  ofstream signiffile; // output final haplist used
  ofstream hotfile; // hotspot info
  
  if(collectdata){
    string recomfilename = filename+"_recom";
    string monitorfilename = filename+"_monitor";
    string  hotfilename = filename +"_hotspot";
      

       // this to make sure results for multiple datasets are sent to a single file
    if(cmds["NDatasets"]>1 && !firsttime){
      recomfile.open(recomfilename.c_str(),ios::app);
      assure (recomfile, recomfilename );
      monitorfile.open(monitorfilename.c_str(),ios::app);
      assure ( monitorfile, monitorfilename );
      if(RecomModel ==1 || RecomModel ==2){
	hotfile.open(hotfilename.c_str(),ios::app);
	assure (hotfile, hotfilename);
      }
    } else {
      recomfile.open(recomfilename.c_str());
      assure (recomfile, recomfilename );
      monitorfile.open(monitorfilename.c_str());
      assure ( monitorfile, monitorfilename );
      if(RecomModel ==1 || RecomModel ==2){
	hotfile.open(hotfilename.c_str());
	assure (hotfile, hotfilename);
      }
      firsttime = false;
    }
    
    string pairsfilename = filename+"_pairs";
    pairsfile.open(pairsfilename.c_str());
    assure (pairsfile, pairsfilename);

//     string finalfilename = filename + "_final";
//     finalfile.open(finalfilename.c_str());
//     assure (finalfile, finalfilename);
    
    if(testcasecontrol){
      string signiffilename = filename + "_signif";
      //if(cmds["testing"])
      //signiffile.open(signiffilename.c_str(),ios::app);
      //else
      signiffile.open(signiffilename.c_str());
      assure (signiffile, signiffilename);
    }
    
  }

  
  if(outputsample && collectdata){
    string samplefilename = filename+"_sample";
    samplefile.open(samplefilename.c_str()); 
    assure ( samplefile, samplefilename);
  }

  haplist.ClearFreqs(); //
  haplist.Add(pop, Nind); // start with list of all haplotypes in current guess
    
  haplist.ClearPseudoCounts(); // PseudoCounts are used to store posterior mean
  //of freqs of each hap

  bool found;

  // Set up an index of the pairs in list that can make up each individual.
  // index[ind] contains a vector of  pairs of pointers to haps in the list 
  // that can make up individual ind
  vector < vector< pair< ListType::iterator,ListType::iterator> > > index(Nind); 
  haplist.MakePairsIndex(index,pop,true, true); 
  // first bool (true) also removes those haps from list
  // that cannot make up any individuals in pop
  // second bool (true) means that missing positions are
  // checked for match
 if(verbose){
    cout << "Possible pairs " << endl;
    for(int ind = 0; ind< Nind; ind++){
      cout << "Individual:" << ind << endl;
      for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair = index[ind].begin(); hpair != index[ind].end(); hpair++){
	(*hpair).first->first.print_haplotype(cout,coding);    
	cout << ",";
	(*hpair).second->first.print_haplotype(cout,coding); 
	cout << endl;
      }
    }
  }
 
  vector < vector< double > > phaseprobs(Nind); // set up a matrix for storing the prob of each phase for each ind
  for(int ind = 0; ind< Nind; ind++){
    phaseprobs[ind] = vector<double>(index[ind].size(),0.0);
  }


  int nchr = Nind+Nind;
  int nchrminus2 = nchr-2;
     
  double removeamount = 1.0;
  double addamount = 1.0;

  double loglik;

  double SigmaMean=1;
  double SigmaMult=1;  
  ComputeRho();


 
  
  // burn-in iterations
  cerr << "Performing Burn-in iterations" << endl;
  cerr << setw(4) << 0 << "% done\033[A" << endl;
  int nacceptRhoMean = 0;
  int nacceptRhoMult = 0;

  for(int iter = 0; iter < Nburn; iter++){

    if(method == 'R' && (RhoMean>0)){
      ComputeRhoDerivAndCurrentLogProb();
      MHUpdateOrder();
      
      if((iter == Nburn/2)){
	cerr << endl;
	cerr << "Estimating recom rates" << endl;
	if(verbose)
	  cout << "Estimating recom rates" << endl;
	InferRho(Niter/10,SigmaMean,SigmaMult,verbose, d_cmds); 
	cerr << "Continuing Burn-in" << endl;
      }
    
      if(iter > Nburn/2){
	UpdateRho(SigmaMean,SigmaMult,nacceptRhoMean, nacceptRhoMult, d_cmds);
	
	if(verbose){
	  OutputRho(cout);
	  OutputCurrentLogProb(cout);
	}
      }
 
    } else {
      RandomiseOrder();
    }
    
    loglik = 0;

    //double DPRIOR = (betaend/nchr) + ((Nburn-1-iter)*1.0/(Nburn-1))*((betastart-betaend)/nchr);
    double DPRIOR = (betaend) + ((Nburn-1-iter)*1.0/(Nburn-1))*((betastart-betaend));    

    for(int i = 0; i<Nind; i++){    
      int ind = order[i];     
      haplist.SoftRemove(pop[ind],removeamount);
      haplist.MakePositiveHaps();

#ifdef DEBUG
      cout << "Removed ind " << ind << endl;
      cout << "Freqs" << endl;
      haplist.Output(cout, coding);
#endif

      double sumprob = 0;
      vector<double> Prob;
      bool addprobs = (ranf()<DPRIOR); // here DPRIOR is prob of adding,
      // rather than multiplying probs (adding gives more weight to
      // configurations where only one of the 2 haps has high prob

      for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair = index[ind].begin(); hpair != index[ind].end(); hpair++){
		
	double prob1 = haplist.CalcProb(((*hpair).first->first),method,Qptr,nchrminus2,vecRho,0, vector<int>(0), ALLSNP, vecdiffprob); // pop[ind].get_nmissing() as last param would just check missing positions
	double prob2 = haplist.CalcProb(((*hpair).second->first),method,Qptr,nchrminus2,vecRho,0, vector<int>(0), ALLSNP, vecdiffprob); // 
		
	double p = prob1;
	if(addprobs)
	  p += prob2;
	else
	  p *= prob2;
	Prob.push_back(p);
	sumprob+=p;

#ifdef DEBUG       
	((*hpair).first->first).print_haplotype(cout,coding);
	cout << ",";
	((*hpair).second->first).print_haplotype(cout,coding);
	cout << ":" << p << endl;
#endif
	
      }

      //SOME OTHER PRIOR ANNEALING SCHEMES WE TRIED
      // original linear prior annealing (adding DPRIOR to p)
      //double p = (CalcProb(((*hpair).first->first),method,Qptr,nchrminus2,vecRho,0)+DPRIOR)/(1+DPRIOR*haplist.get_listlength()) * 
      //(CalcProb(((*hpair).second->first),method,Qptr,nchrminus2,vecRho,0)+DPRIOR)/(1+DPRIOR*haplist.get_listlength());
      
      // "powered" version of prior annealing, raising p to power DPRIOR
      //double p = exp(DPRIOR * log(CalcProb(((*hpair).first->first),method,Qptr,nchrminus2,vecRho,0)*
      //CalcProb(((*hpair).second->first),method,Qptr,nchrminus2,vecRho,0)));
      
      // store phase probabilities in phaseprobs -only necessary in order to do the later pruning...
      vector<double>::iterator probpointer = Prob.begin();
      for (vector<double>::iterator phase = phaseprobs[ind].begin(); phase != phaseprobs[ind].end(); phase++){
	*phase += (1.0/Nburn)*(*probpointer/sumprob);
	probpointer++;
      }

      int choice = rint2(Prob,sumprob);
      loglik += log(sumprob);
      //loglik += log(Prob[choice]); // a different version of the pseudologlik
      vector< pair<ListType::iterator,ListType::iterator> >::iterator sampledpair = index[ind].begin();
      for(int i = 0; i<choice; i++)
	sampledpair++;
      
#ifdef DEBUG
      cout << "Chosen Pair:" << endl;
      haplist.OutputPair(cout,*sampledpair, coding);
#endif
#ifdef DEBUG
      cout << "ind = " << ind << endl;
      cout << "Status before update:" << endl;
      pop[ind].print_haplotypes(cout,loci_type,coding,true,true,false,0,0);
#endif

      pop[ind].update_haplotypes((*sampledpair).first->first,(*sampledpair).second->first);

#ifdef DEBUG
      cout << "ind = " << ind << endl;
      cout << "Status after update:" << endl;
      pop[ind].print_haplotypes(cout,loci_type,coding,true,true,false,0,0);
#endif

      HapListImputeMissing(ind);  
 
#ifdef DEBUG
      cout << "ind = " << ind << endl;
      cout << "Status after Imputing missing:" << endl;
      pop[ind].print_haplotypes(cout,loci_type,coding,true,true,false,0,0);
#endif
     
      for(int chr = 0; chr <2 ; chr++){
	bool isnewhap;
	ListType::iterator newrecord = haplist.Add(pop[ind], chr,addamount, isnewhap);

	// extend index just for the individual ind
	// if it has missing loci, as may have added haplotype pair
	// that wasn't on its list of allowed pairs
	if( pop[ind].NMissingLoci() > 0)
	  haplist.ExtendPairsIndex(newrecord, index[ind], phaseprobs[ind], pop[ind], true);



	if(isnewhap)
	  haplist.ExtendPairsIndex(newrecord, index, phaseprobs, pop, true);   
	// true param says to check missing positions
      }

   if(verbose){
    cout << "Possible pairs " << endl;
    for(int ind = 0; ind< Nind; ind++){
      cout << "Individual:" << ind << endl;
      for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair = index[ind].begin(); hpair != index[ind].end(); hpair++){
        (*hpair).first->first.print_haplotype(cout,coding);
        cout << ",";
        (*hpair).second->first.print_haplotype(cout,coding);
        cout << endl;
      }
    }
  }

 
#ifdef DEBUG
      cout << "Adding ind:" << endl;
      haplist.Output(cout,coding);
#endif
      
    }

    cerr << setw(4) << (int) ((iter+1)*100.0/Nburn) << "% done\033[A" << endl;
    //cout << "Log Likelihood: " << loglik << endl;    
  }
  
  if(verbose){
      cout << "Acceptance Rate for RhoMean, burnin: " << nacceptRhoMean*1.0/(Nburn) << endl;
      cout << "Acceptance Rate for RhoMult, burnin:: " << nacceptRhoMult*1.0/(Nburn) << endl; 
      cout << "SigmaMean = " << SigmaMean << endl;
      cout << "SigmaMult = " << SigmaMult << endl;
  }
  

  // after burn-in, if necessary get better estimate of rho and appropriate sigmas
  if(nacceptRhoMean *2.0/Nburn < 0.1 || nacceptRhoMean *2.0/Nburn >0.9 || nacceptRhoMult *2.0/Nburn < 0.1 || nacceptRhoMult *2.0/Nburn >0.9){
    if(method == 'R' && (RhoMean >0)){
      cerr << "Re-Estimating recom rates" << endl;
      if(verbose)
	cout << "Estimating recom rates" << endl;
      InferRho(Niter/10,SigmaMean,SigmaMult,verbose, d_cmds); 
    }
  }

  nacceptRhoMean=nacceptRhoMult=0;

  // actual iterations
  double meanloglik = 0;

  // set up vecpermcc as vector of permutations of casecontrol labels    
  vector<double> meanpermloglik(nperm , 0); // mean log lik for controls

  vector< vector<int> > vecpermcc;

  if(testcasecontrol && collectdata){
    vecpermcc = vector<vector<int> > (nperm, casecontrol);
    for(int i = 1; i<nperm; i++){
      random_shuffle(vecpermcc[i].begin(),vecpermcc[i].end());
    }
  }

  
  // Here we prune index to remove very unlikely hap pairs
  // (also sets phaseprobs to be 0)
  haplist.PrunePairsIndex(index, phaseprobs,pop,minfreqpair*0.5);

  // set up group frequency vectors
  int ngroups = 0;
  for(int ind = 0; ind<Nind; ind++)
    if(casecontrol[ind]>ngroups)
      ngroups = casecontrol[ind];
  ngroups++;    
  haplist.SetupGroupFreqs(ngroups);
  
  // set up group size vector
  groupsize = vector<int>(ngroups,0);
  for(int ind = 0; ind<Nind; ind++)
    groupsize[casecontrol[ind]] += 2;


  cerr << "Performing Main iterations" << endl;
  cerr << setw(4) << 0 << "% done\033[A" << endl; 
   
  //int NumberOfSqPseudoCountUpdates = 0;
  
  haplist.ClearPseudoCounts(); // PseudoCounts are used to store posterior mean
  //of freqs of each hap

  if(method =='R' && collectdata){
   OutputPositions(recomfile);
  }

  for(int iter = 0; iter < Niter; iter++){
    loglik = 0;
      
    if(method == 'R' && (RhoMean>0)){
      ComputeRhoDerivAndCurrentLogProb();
      MHUpdateOrder();
      
      UpdateRho(SigmaMean,SigmaMult, nacceptRhoMean, nacceptRhoMult, d_cmds);
      if(verbose){
	OutputRho(cout);
	OutputCurrentLogProb(cout);
      }
      if(collectdata){
	OutputRho(recomfile);
	if(RecomModel == 1 || RecomModel ==2)
	  OutputHotspotParams(hotfile);
      }
    } else {
      RandomiseOrder();
    }

    if(testcasecontrol && collectdata){
      for(int i = 0; i< nperm; i++){
	meanpermloglik[i] += 
	  logProb( method, vecRho, vecpermcc[i], 1.0, 0) + 
	  logProb( method, vecRho, vecpermcc[i], 1.0, 1);
	// add likelihoods for group 0 and group 1
      }
    }
   
    for(int i = 0; i<Nind; i++){
      int ind = order[i];
    
      haplist.SoftRemove(pop[ind],removeamount);
      haplist.MakePositiveHaps();

#ifdef DEBUG
      cout << "Removed ind " << ind << endl;
      cout << "Freqs" << endl;
      haplist.Output(cout, coding);
#endif

      // make list of probabilities of each hap pair
      double sumprob = 0;
      vector<double> Prob;
      for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair = index[ind].begin(); hpair != index[ind].end(); hpair++){
	double prob1 = haplist.CalcProb(((*hpair).first->first),method,Qptr,nchrminus2,vecRho,0, vector<int>(0), ALLSNP, vecdiffprob);
	double prob2 = haplist.CalcProb(((*hpair).second->first),method,Qptr,nchrminus2,vecRho,0, vector<int>(0), ALLSNP, vecdiffprob);

	double factor = (2-((*hpair).first == (*hpair).second)); // take acount of factor of 2 for case where there are missing positions
	double p = factor * prob1 * prob2;

	// now account for Pr(gen | hap) at those positions where
	// only one missing
	for(int locus = 0; locus< Nloci; locus++){
	  if((pop[ind].n_missing(locus) == 1) && (pop[ind].is_unknown(locus)))
	    if(((*hpair).first->first).get_allele(locus)!=((*hpair).second->first).get_allele(locus))
	      p *= 0.5;
	}
	
	Prob.push_back(p);
	sumprob +=p;
	
	
#ifdef DEBUG  
	cout << "ind:" << ind << endl;     
	((*hpair).first->first).print_haplotype(cout,coding);
	cout << ",";
	((*hpair).second->first).print_haplotype(cout,coding);
	cout << ":" << p << endl;
#endif	
      }
    
    // Store data necessary to compute (Rao-Blackwellised) estimates of posterior mean of freqs, and freqs squared
      // (this is the old version - new version appears below)
      // for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair = index[ind].begin(); hpair != index[ind].end(); hpair++){
// 	haplist.haplist[(*hpair).first->first].PseudoCount += Prob[probpos]/sumprob;
// 	haplist.haplist[(*hpair).second->first].PseudoCount += Prob[probpos]/sumprob;	
// 	probpos++;
//       }

      int probpos = 0;
      for(ListType::iterator h = haplist.haplist.begin(); h != haplist.haplist.end(); h++){
	  double hapcount = h->second.Freq;
	  double grouphapcount = nchr-2;
	  if(ngroups>1)
	    grouphapcount = (double) GetGroupCount(h->first,casecontrol[ind],ind);
	  int groupn = groupsize[casecontrol[ind]];

	  h->second.PseudoCount += (hapcount /nchr); // pseudocount stores posterior means of overall frequencies (bit of a misnomer for historical reasons! should change...)	 
	  
	  h->second.GroupFreq[casecontrol[ind]] += (grouphapcount / groupn);

	  if(collectdata){
	    h->second.SqPseudoCount += (hapcount * hapcount)/(nchr*nchr);
	    h->second.GroupFreqSq[casecontrol[ind]] += (grouphapcount * grouphapcount)/(groupn*groupn) ;
	  }

	  probpos = 0;
	  for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair = index[ind].begin(); hpair != index[ind].end(); hpair++){
	    int indicator1 = (*hpair).first->first == (h->first);
	    int indicator2 = (*hpair).second->first == (h->first);
 
	    h->second.PseudoCount += (Prob[probpos]/sumprob) * (1.0/nchr) * (indicator1 + indicator2);
	    h->second.GroupFreq[casecontrol[ind]] += (Prob[probpos]/sumprob) * (1.0/groupn) * (indicator1 + indicator2);
	    
	    if(collectdata){
	      h->second.SqPseudoCount += (Prob[probpos]/sumprob) * (1.0/(nchr*nchr)) * (indicator1 + indicator2) *
		(indicator1 + indicator2 + 2*hapcount);
	      h->second.GroupFreqSq[casecontrol[ind]] += (Prob[probpos]/sumprob) * (1.0/(groupn*groupn)) * (indicator1 + indicator2) *
		(indicator1 + indicator2 + 2*grouphapcount);	    
	    }
	    probpos++;
	  } 
      }
      
      //store phase probabilities in phaseprobs
      vector<double>::iterator probpointer = Prob.begin();
      for (vector<double>::iterator phase = phaseprobs[ind].begin(); phase != phaseprobs[ind].end(); phase++){
	*phase += (1.0/Niter)*(*probpointer/sumprob);
	probpointer++;
      }

      int choice = rint2(Prob,sumprob);
      loglik += log(sumprob);      //other possibility: loglik += log(Prob[choice]); 

      vector< pair<ListType::iterator,ListType::iterator> >::iterator sampledpair = index[ind].begin();
      for(int i = 0; i<choice; i++)
	sampledpair++;
      
#ifdef DEBUG
      cout << "Chosen Pair:" << endl;
      haplist.OutputPair(cout,*sampledpair,coding);
#endif

      //pop[ind].update_haplotype(0,(*sampledpair).first->first);
      //pop[ind].update_haplotype(1,(*sampledpair).second->first);
      pop[ind].update_haplotypes((*sampledpair).first->first,(*sampledpair).second->first);
      
      // this bit allows missing positions to be reconstructed
      // without reference to the HapList;
      // 
      // is not quite correct for case where one allele may be missing
      // 
      HapListImputeMissing(ind);
//        // add new hap to list; extend Pairs index if necessary
      for(int chr = 0; chr <2 ; chr++){
	bool isnewhap;
	ListType::iterator newrecord = haplist.Add(pop[ind], chr,addamount, isnewhap);
	// extend index just for the individual ind
	// if it has missing loci, as may have added haplotype pair
	// that wasn't on its list of allowed pairs
	if(pop[ind].NMissingLoci()>0)
	   haplist.ExtendPairsIndex(newrecord, index[ind], phaseprobs[ind], pop[ind], true);
	if(isnewhap)
	  haplist.ExtendPairsIndex(newrecord, index, phaseprobs, pop, true);
      }
      
#ifdef DEBUG
      cout << "Adding ind:" << endl;
      haplist.Output(cout,coding);
#endif
      
      //pop[ind].UpdateCounts(); 

      if(verbose){
	OutputRho(cout);
	OutputCurrentLogProb(cout);
      }
    }
   
    if(collectdata){
      monitorfile << loglik << " " << CurrentLogProb <<  endl; 
      if(outputsample)
	output(samplefile, true, true);
    }
    
    cerr << setw(4) << (int) ((iter+1)*100.0/(Niter)) << "% done\033[A" << endl; 
    meanloglik +=loglik;
  }
  
  if(verbose){
      cout << "Acceptance Rate for RhoMean: " << nacceptRhoMean*1.0/(Niter) << endl;
      cout << "Acceptance Rate for RhoMult: " << nacceptRhoMult*1.0/(Niter) << endl; 
      cout << "SigmaMean = " << SigmaMean << endl;
      cout << "SigmaMult = " << SigmaMult << endl;
  }
  
  vecpermcc = vector<vector<int> > (0);

  if(testcasecontrol && collectdata){
    double nbigger = 0;
    //cout << "meanpermloglik[0] = " << meanpermloglik[0] << endl;
    for(int i = 1; i<nperm; i++){
      nbigger += (meanpermloglik[i]>=meanpermloglik[0]);
      //if(meanpermloglik[i]==meanpermloglik[0])
      //nbigger += 0.5;
      //cout << "meanpermloglik[i] = " << meanpermloglik[i] << endl;
    }
    signiffile << "P-value for testing H0: Cases ~ Controls = " << (1.0*(nbigger+1))/(nperm) << endl;
  }

// copy PsedoCounts (posterior means of Freqs) to Freqs (suitably normalised)
  haplist.CopyPseudoCountsToFreqs();
  haplist.NormaliseFreqs();
  
  if(collectdata){
    haplist.NormaliseSqPseudoCounts(Nind * Niter); 

    haplist.NormaliseGroupFreqs();

    cerr << "Writing output to files " << endl;
    // Output list of plausible pairs for each individual 
    for(int ind = 0; ind<Nind; ind++){
      pairsfile << "IND: " << pop[ind].get_id() << endl;
      for(int j = 0; j< index[ind].size(); j++){
	if(phaseprobs[ind][j] > minfreqpairoutput){
	  haplist.OutputPair(pairsfile,index[ind][j],coding);
	  pairsfile << " , ";
	  pairsfile.setf(ios::fixed);
	  pairsfile.setf(ios::showpoint);
	  pairsfile.precision(3); 
	  pairsfile << phaseprobs[ind][j] << endl;
	}
      }
    }
    
    // Output final haplist used
    
    // haplist.Dump(finalfile);

    haplist.SetBestGuesses(pop,index,phaseprobs); // set all individuals in pop to their best guesses
   
    //    cout << "TestOut0:" << pop[0].get_allele(0,1) << endl;
    //cout << "TestOut1:" << pop[0].get_allele(1,1) << endl;
 
    cerr << "Producing Summary, please wait " << endl;
    // each element of this vector is a summary for a single individual
    // the false at the end forces it to produce a single summary with no splits
    vector<Summary> summary = haplist.ProduceSummary(index,phaseprobs,0,get_nloci(),pop,false);
    
    // for(vector<Summary>::iterator s = summary.begin(); s<summary.end(); s++)
    //  s->Output(cout, coding);
    
    TransferCounts(summary);

    // produce separate lists of freqs for cases and controls
    // int ngroups = 0;
//     for(int ind = 0; ind<Nind; ind++)
//       if(casecontrol[ind]>ngroups)
// 	ngroups = casecontrol[ind];
//     ngroups++;
    
//     haplist.SetupGroupFreqs(ngroups);
//     for(int ind = 0; ind<Nind; ind++){
//       int probpos = 0;
//       for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair = index[ind].begin(); hpair != index[ind].end(); hpair++){  
// 	haplist.haplist[(*hpair).first->first].GroupFreq[casecontrol[ind]] += phaseprobs[ind][probpos]/(2*Nind);
// 	haplist.haplist[(*hpair).second->first].GroupFreq[casecontrol[ind]] += phaseprobs[ind][probpos]/(2*Nind);	  
// 	probpos++;
//       }
//     }
  }


  return meanloglik;

}

double ClassPop::FuzzyHapListMCMCResolvePhaseRemove(map<string,int> & cmds, int Niter, int Nthin, int Nburn, map<string,double> & d_cmds, string filename, bool collectdata)
{

  ///Note on how missing data are dealt with:
  // First, when list is initialised (before calling this fn)
  // all possibilities for SNPs are added. (for Multiallelics,
  // all possibilities based on current imputed alleles are added).
  // Second, when list of possible pairs is made up, all pairs 
  // consistent with observed alleles are added.


  HapList fuzzylist; // fuzzylist contains haplotypes to be used in computing
  // conditional probs
  fuzzylist.RemoveAll();
  fuzzylist.Add(pop[0]);
  fuzzylist.MakePositiveHaps();

  static int firsttime = true;

  int verbose = cmds["verbose"];
  int outputsample = cmds["outputsample"];
  int method = cmds["method"];
  int nperm = cmds["nperm"];
  
  RecomModel = cmds["RecomModel"];
 
  int testcasecontrol = (nperm>0);
  
  double rho = d_cmds["rhostart"];
  double betastart = d_cmds["betastart"];
  double betaend = d_cmds["betaend"];
  double minfreqpairoutput = d_cmds["output_minfreq"]; //min prob that a particular
  // pair has to have before being output as a possible pair
  double minfreqpair = d_cmds["minfreq"];
  
  if(collectdata){
    Niter *= cmds["finalrepmult"];
    Nthin *= cmds["finalrepmult"];
    Nburn *= cmds["finalrepmult"];
  }

  if(method == 'Q') // for Q method, use no recom, except on the final run
    if(collectdata)
      method = 'R';
    else
      method = 'S';
  
  cout << "Method = " << (char) method << endl;
  
  if(collectdata)
    cout << "Performing Final Set of Iterations... nearly there!" << endl;

  ofstream recomfile;
  ofstream monitorfile; 
  ofstream samplefile; // sample from posterior
  ofstream pairsfile; // list of high probability haplotype pairs 
                      //for each individual
  // ofstream finalfile; // output final haplist used
  ofstream signiffile; // output final haplist used
  ofstream hotfile; // hotspot info
  
  if(collectdata){
    string recomfilename = filename+"_recom";
    string monitorfilename = filename+"_monitor";
    string  hotfilename = filename +"_hotspot";
      

       // this to make sure results for multiple datasets are sent to a single file
    if(cmds["NDatasets"]>1 && !firsttime){
      recomfile.open(recomfilename.c_str(),ios::app);
      assure (recomfile, recomfilename );
      monitorfile.open(monitorfilename.c_str(),ios::app);
      assure ( monitorfile, monitorfilename );
      if(RecomModel ==1 || RecomModel ==2){
	hotfile.open(hotfilename.c_str(),ios::app);
	assure (hotfile, hotfilename);
      }
    } else {
      recomfile.open(recomfilename.c_str());
      assure (recomfile, recomfilename );
      monitorfile.open(monitorfilename.c_str());
      assure ( monitorfile, monitorfilename );
      if(RecomModel ==1 || RecomModel ==2){
	hotfile.open(hotfilename.c_str());
	assure (hotfile, hotfilename);
      }
      firsttime = false;
    }
    
    string pairsfilename = filename+"_pairs";
    pairsfile.open(pairsfilename.c_str());
    assure (pairsfile, pairsfilename);

//     string finalfilename = filename + "_final";
//     finalfile.open(finalfilename.c_str());
//     assure (finalfile, finalfilename);
    
    if(testcasecontrol){
      string signiffilename = filename + "_signif";
      //if(cmds["testing"])
      //signiffile.open(signiffilename.c_str(),ios::app);
      //else
      signiffile.open(signiffilename.c_str());
      assure (signiffile, signiffilename);
    }
    
  }

  if(outputsample && collectdata){
    string samplefilename = filename+"_sample";
    samplefile.open(samplefilename.c_str()); 
    assure ( samplefile, samplefilename);
  }

  haplist.ClearFreqs(); //
  haplist.Add(pop, Nind); // start with list of all haplotypes in current guess
    
  haplist.ClearPseudoCounts(); // PseudoCounts are used to store posterior mean
  //of freqs of each hap

  bool found;

  // Set up an index of the pairs in list that can make up each individual.
  // index[ind] contains a vector of  pairs of pointers to haps in the list 
  // that can make up individual ind
  vector < vector< pair< ListType::iterator,ListType::iterator> > > index(Nind); 
  haplist.MakePairsIndex(index,pop,true, true); 
  // first bool (true) also removes those haps from list
  // that cannot make up any individuals in pop
  // second bool (true) means that missing positions are
  // checked for match
 if(verbose){
    cout << "Possible pairs " << endl;
    for(int ind = 0; ind< Nind; ind++){
      cout << "Individual:" << ind << endl;
      for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair = index[ind].begin(); hpair != index[ind].end(); hpair++){
	(*hpair).first->first.print_haplotype(cout,coding);    
	cout << ",";
	(*hpair).second->first.print_haplotype(cout,coding); 
	cout << endl;
      }
    }
  }
 
  vector < vector< double > > phaseprobs(Nind); // set up a matrix for storing the prob of each phase for each ind
  for(int ind = 0; ind< Nind; ind++){
    phaseprobs[ind] = vector<double>(index[ind].size(),0.0);
  }


  int nchr = Nind+Nind;
  int nchrminus2 = nchr-2;
     
  double removeamount = 1.0;
  double addamount = 1.0;

  double loglik;

  double SigmaMean=1;
  double SigmaMult=1;
  if(method == 'R' && (RhoMean >0)){
    cerr << "Estimating recom rates" << endl;
    if(verbose)
      cout << "Estimating recom rates" << endl;
    InferRho(Niter/10,SigmaMean,SigmaMult,verbose, d_cmds); 
  }
  
  // burn-in iterations
  cerr << "Performing Burn-in iterations" << endl;
  cerr << setw(4) << 0 << "% done\033[A" << endl;
  int nacceptRhoMean = 0;
  int nacceptRhoMult = 0;

  for(int iter = 0; iter < Nburn; iter++){
    loglik = 0;
    
    if(method == 'R' && (RhoMean>0)){
      ComputeRhoDerivAndCurrentLogProb();
      MHUpdateOrder();
      
      UpdateRho(SigmaMean,SigmaMult,nacceptRhoMean, nacceptRhoMult, d_cmds);
      if(verbose){
	OutputRho(cout);
	OutputCurrentLogProb(cout);
      }
    } else {
      RandomiseOrder();
    }

    //double DPRIOR = (betaend/nchr) + ((Nburn-1-iter)*1.0/(Nburn-1))*((betastart-betaend)/nchr);
    double DPRIOR = (betaend) + ((Nburn-1-iter)*1.0/(Nburn-1))*((betastart-betaend));    

    for(int i = 0; i<Nind; i++){    
      int ind = order[i];     
      haplist.SoftRemove(pop[ind],removeamount);
      haplist.MakePositiveHaps();
      
      double sumprob = 0;
      vector<double> Prob;
      bool addprobs = (ranf()<DPRIOR); // here DPRIOR is prob of adding,
      // rather than multiplying probs (adding gives more weight to
      // configurations where only one of the 2 haps has high prob

      for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair = index[ind].begin(); hpair != index[ind].end(); hpair++){
		
	double prob1 = fuzzylist.CalcProb(((*hpair).first->first),method,Qptr,nchrminus2,vecRho,0, vector<int>(0), ALLSNP, vecdiffprob, true); // pop[ind].get_nmissing() as last param would just check missing positions
	double prob2 = fuzzylist.CalcProb(((*hpair).second->first),method,Qptr,nchrminus2,vecRho,0, vector<int>(0), ALLSNP, vecdiffprob, true); // 
		
	double p = prob1;
	if(addprobs)
	  p += prob2;
	else
	  p *= prob2;
	Prob.push_back(p);
	sumprob+=p;

#ifdef DEBUG       
	((*hpair).first->first).print_haplotype(cout,coding);
	cout << ",";
	((*hpair).second->first).print_haplotype(cout,coding);
	cout << ":" << p << endl;
#endif
	
      }

      // store phase probabilities in phaseprobs -only necessary in order to do the later pruning...
      vector<double>::iterator probpointer = Prob.begin();
      for (vector<double>::iterator phase = phaseprobs[ind].begin(); phase != phaseprobs[ind].end(); phase++){
	*phase += (1.0/Nburn)*(*probpointer/sumprob);
	probpointer++;
      }

      int choice = rint2(Prob,sumprob);
      loglik += log(sumprob);
      //loglik += log(Prob[choice]); // a different version of the pseudologlik
      vector< pair<ListType::iterator,ListType::iterator> >::iterator sampledpair = index[ind].begin();
      for(int i = 0; i<choice; i++)
	sampledpair++;
      
#ifdef DEBUG
      cout << "Chosen Pair:" << endl;
      haplist.OutputPair(cout,*sampledpair, coding);
#endif

      pop[ind].update_haplotypes((*sampledpair).first->first,(*sampledpair).second->first);
      
      HapListImputeMissing(ind);      
      for(int chr = 0; chr <2 ; chr++){
	bool isnewhap;
	ListType::iterator newrecord = haplist.Add(pop[ind], chr,addamount, isnewhap);
	// extend index just for the individual ind
	// if it has missing loci, as may have added haplotype pair
	// that wasn't on its list of allowed pairs
	if( pop[ind].NMissingLoci() > 0)
	  haplist.ExtendPairsIndex(newrecord, index[ind], phaseprobs[ind], pop[ind], true);

	if(isnewhap)
	  haplist.ExtendPairsIndex(newrecord, index, phaseprobs, pop, true);   
	// true param says to check missing positions
      }

    
#ifdef DEBUG
      cout << "Adding ind:" << endl;
      haplist.Output(cout,coding);
#endif
      
    }

    cerr << setw(4) << (int) ((iter+1)*100.0/Nburn) << "% done\033[A" << endl;
    //cout << "Log Likelihood: " << loglik << endl;    
  }
  
  if(verbose){
      cout << "Acceptance Rate for RhoMean, burnin: " << nacceptRhoMean*1.0/(Nburn) << endl;
      cout << "Acceptance Rate for RhoMult, burnin:: " << nacceptRhoMult*1.0/(Nburn) << endl; 
      cout << "SigmaMean = " << SigmaMean << endl;
      cout << "SigmaMult = " << SigmaMult << endl;
  }
  
  nacceptRhoMean=nacceptRhoMult=0;

  // actual iterations
  double meanloglik = 0;

  // set up vecpermcc as vector of permutations of casecontrol labels    
  vector<double> meanpermloglik(nperm , 0); // mean log lik for controls

  vector< vector<int> > vecpermcc;

  if(testcasecontrol && collectdata){
    vecpermcc = vector<vector<int> > (nperm, casecontrol);
    for(int i = 1; i<nperm; i++){
      random_shuffle(vecpermcc[i].begin(),vecpermcc[i].end());
    }
  }

  
  // Here we prune index to remove very unlikely hap pairs
  // (also sets phaseprobs to be 0)
  haplist.PrunePairsIndex(index, phaseprobs,pop,minfreqpair*0.5);

  // set up group frequency vectors
  int ngroups = 0;
  for(int ind = 0; ind<Nind; ind++)
    if(casecontrol[ind]>ngroups)
      ngroups = casecontrol[ind];
  ngroups++;    
  haplist.SetupGroupFreqs(ngroups);
  
  // set up group size vector
  groupsize = vector<int>(ngroups,0);
  for(int ind = 0; ind<Nind; ind++)
    groupsize[casecontrol[ind]] += 2;


  cerr << "Performing Main iterations" << endl;
  cerr << setw(4) << 0 << "% done\033[A" << endl; 
   
  //int NumberOfSqPseudoCountUpdates = 0;
  
  haplist.ClearPseudoCounts(); // PseudoCounts are used to store posterior mean
  //of freqs of each hap

  if(method =='R' && collectdata){
   OutputPositions(recomfile);
  }

  for(int iter = 0; iter < Niter; iter++){
    loglik = 0;
      
    if(method == 'R' && (RhoMean>0)){
      ComputeRhoDerivAndCurrentLogProb();
      MHUpdateOrder();
      
      UpdateRho(SigmaMean,SigmaMult, nacceptRhoMean, nacceptRhoMult, d_cmds);
      if(verbose){
	OutputRho(cout);
	OutputCurrentLogProb(cout);
      }
      if(collectdata){
	OutputRho(recomfile);
	if(RecomModel == 1 || RecomModel ==2)
	  OutputHotspotParams(hotfile);
      }
    } else {
      RandomiseOrder();
    }

    if(testcasecontrol && collectdata){
      for(int i = 0; i< nperm; i++){
	meanpermloglik[i] += 
	  logProb( method, vecRho, vecpermcc[i], 1.0, 0) + 
	  logProb( method, vecRho, vecpermcc[i], 1.0, 1);
	// add likelihoods for group 0 and group 1
      }
    }
   
    for(int i = 0; i<Nind; i++){
      int ind = order[i];
      haplist.SoftRemove(pop[ind],removeamount);
      haplist.MakePositiveHaps();

      // make list of probabilities of each hap pair
      double sumprob = 0;
      vector<double> Prob;
      for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair = index[ind].begin(); hpair != index[ind].end(); hpair++){
	double prob1 = fuzzylist.CalcProb(((*hpair).first->first),method,Qptr,nchrminus2,vecRho,0, vector<int>(0), ALLSNP, vecdiffprob, true);
	double prob2 = fuzzylist.CalcProb(((*hpair).second->first),method,Qptr,nchrminus2,vecRho,0, vector<int>(0), ALLSNP, vecdiffprob, true);

	double factor = (2-((*hpair).first == (*hpair).second)); // take acount of factor of 2 for case where there are missing positions
	double p = factor * prob1 * prob2;

	// now account for Pr(gen | hap) at those positions where
	// only one missing
	for(int locus = 0; locus< Nloci; locus++){
	  if(pop[ind].n_missing(locus) == 1)
	    if(((*hpair).first->first).get_allele(locus)!=((*hpair).second->first).get_allele(locus))
	      p *= 0.5;
	}
	
	Prob.push_back(p);
	sumprob +=p;
	
	
#ifdef DEBUG  
	cout << "ind:" << ind << endl;     
	((*hpair).first->first).print_haplotype(cout,coding);
	cout << ",";
	((*hpair).second->first).print_haplotype(cout,coding);
	cout << ":" << p << endl;
#endif	
      }
    

      int probpos = 0;
      for(ListType::iterator h = haplist.haplist.begin(); h != haplist.haplist.end(); h++){
	  double hapcount = h->second.Freq;
	  double grouphapcount = nchr-2;
	  if(ngroups>1)
	    grouphapcount = (double) GetGroupCount(h->first,casecontrol[ind],ind);
	  int groupn = groupsize[casecontrol[ind]];

	  h->second.PseudoCount += (hapcount /nchr); // pseudocount stores posterior means of overall frequencies (bit of a misnomer for historical reasons! should change...)	 
	  
	  h->second.GroupFreq[casecontrol[ind]] += (grouphapcount / groupn);

	  if(collectdata){
	    h->second.SqPseudoCount += (hapcount * hapcount)/(nchr*nchr);
	    h->second.GroupFreqSq[casecontrol[ind]] += (grouphapcount * grouphapcount)/(groupn*groupn) ;
	  }

	  probpos = 0;
	  for(vector< pair<ListType::iterator,ListType::iterator> >::iterator hpair = index[ind].begin(); hpair != index[ind].end(); hpair++){
	    int indicator1 = (*hpair).first->first == (h->first);
	    int indicator2 = (*hpair).second->first == (h->first);
 
	    h->second.PseudoCount += (Prob[probpos]/sumprob) * (1.0/nchr) * (indicator1 + indicator2);
	    h->second.GroupFreq[casecontrol[ind]] += (Prob[probpos]/sumprob) * (1.0/groupn) * (indicator1 + indicator2);
	    
	    if(collectdata){
	      h->second.SqPseudoCount += (Prob[probpos]/sumprob) * (1.0/(nchr*nchr)) * (indicator1 + indicator2) *
		(indicator1 + indicator2 + 2*hapcount);
	      h->second.GroupFreqSq[casecontrol[ind]] += (Prob[probpos]/sumprob) * (1.0/(groupn*groupn)) * (indicator1 + indicator2) *
		(indicator1 + indicator2 + 2*grouphapcount);	    
	    }
	    probpos++;
	  } 
      }
      
      //store phase probabilities in phaseprobs
      vector<double>::iterator probpointer = Prob.begin();
      for (vector<double>::iterator phase = phaseprobs[ind].begin(); phase != phaseprobs[ind].end(); phase++){
	*phase += (1.0/Niter)*(*probpointer/sumprob);
	probpointer++;
      }

int choice = rint2(Prob,sumprob);
loglik += log(sumprob);      //other possibility: loglik += log(Prob[choice]); 

vector< pair<ListType::iterator,ListType::iterator> >::iterator sampledpair = index[ind].begin();
for(int i = 0; i<choice; i++)
sampledpair++;

#ifdef DEBUG
cout << "Chosen Pair:" << endl;
haplist.OutputPair(cout,*sampledpair,coding);
#endif

//pop[ind].update_haplotype(0,(*sampledpair).first->first);
//pop[ind].update_haplotype(1,(*sampledpair).second->first);
pop[ind].update_haplotypes((*sampledpair).first->first,(*sampledpair).second->first);

// this bit allows missing positions to be reconstructed
// without reference to the HapList;
// 
// is not quite correct for case where one allele may be missing
// 
HapListImputeMissing(ind);
//        // add new hap to list; extend Pairs index if necessary
for(int chr = 0; chr <2 ; chr++){
bool isnewhap;
ListType::iterator newrecord = haplist.Add(pop[ind], chr,addamount, isnewhap);
// extend index just for the individual ind
// if it has missing loci, as may have added haplotype pair
// that wasn't on its list of allowed pairs
if(pop[ind].NMissingLoci()>0)
   haplist.ExtendPairsIndex(newrecord, index[ind], phaseprobs[ind], pop[ind], true);
if(isnewhap)
  haplist.ExtendPairsIndex(newrecord, index, phaseprobs, pop, true);
}

#ifdef DEBUG
cout << "Adding ind:" << endl;
haplist.Output(cout,coding);
#endif

//pop[ind].UpdateCounts(); 

if(verbose){
OutputRho(cout);
OutputCurrentLogProb(cout);
}
}

if(collectdata){
monitorfile << loglik << " " << CurrentLogProb <<  endl; 
if(outputsample)
output(samplefile, true, true);
}

cerr << setw(4) << (int) ((iter+1)*100.0/(Niter)) << "% done\033[A" << endl; 
meanloglik +=loglik;
}

if(verbose){
cout << "Acceptance Rate for RhoMean: " << nacceptRhoMean*1.0/(Niter) << endl;
cout << "Acceptance Rate for RhoMult: " << nacceptRhoMult*1.0/(Niter) << endl; 
cout << "SigmaMean = " << SigmaMean << endl;
cout << "SigmaMult = " << SigmaMult << endl;
}

vecpermcc = vector<vector<int> > (0);

if(testcasecontrol && collectdata){
double nbigger = 0;
//cout << "meanpermloglik[0] = " << meanpermloglik[0] << endl;
for(int i = 1; i<nperm; i++){
nbigger += (meanpermloglik[i]>=meanpermloglik[0]);
//if(meanpermloglik[i]==meanpermloglik[0])
//nbigger += 0.5;
//cout << "meanpermloglik[i] = " << meanpermloglik[i] << endl;
}
signiffile << "P-value for testing H0: Cases ~ Controls = " << (1.0*(nbigger+1))/(nperm) << endl;
}

// copy PsedoCounts (posterior means of Freqs) to Freqs (suitably normalised)
haplist.CopyPseudoCountsToFreqs();
haplist.NormaliseFreqs();

if(collectdata){
haplist.NormaliseSqPseudoCounts(Nind * Niter); 

haplist.NormaliseGroupFreqs();

cerr << "Writing output to files " << endl;
// Output list of plausible pairs for each individual 
for(int ind = 0; ind<Nind; ind++){
pairsfile << "IND: " << pop[ind].get_id() << endl;
for(int j = 0; j< index[ind].size(); j++){
if(phaseprobs[ind][j] > minfreqpairoutput){
  haplist.OutputPair(pairsfile,index[ind][j],coding);
  pairsfile << " , ";
  pairsfile.setf(ios::fixed);
  pairsfile.setf(ios::showpoint);
  pairsfile.precision(3); 
  pairsfile << phaseprobs[ind][j] << endl;
}
}
}

// Output final haplist used

// haplist.Dump(finalfile);

haplist.SetBestGuesses(pop,index,phaseprobs); // set all individuals in pop to their best guesses

cerr << "Producing Summary, please wait " << endl;
// each element of this vector is a summary for a single individual
// the false at the end forces it to produce a single summary with no splits
vector<Summary> summary = haplist.ProduceSummary(index,phaseprobs,0,get_nloci(),pop,false);

// for(vector<Summary>::iterator s = summary.begin(); s<summary.end(); s++)
//  s->Output(cout, coding);

TransferCounts(summary);
}


return meanloglik;

}

// count how many haplotypes of type h are in group g, 
// omitting individual "omit"
int ClassPop::GetGroupCount(const Haplotype & h, int g, int omit)
{
int count = 0;
for(int ind = 0; ind < Nind; ind++){
if((ind != omit) && (casecontrol[ind]==g)){
count += pop[ind].get_haplotype(0) == h;
count += pop[ind].get_haplotype(1) == h;
}
}
return count;
}



void ClassPop::FastHapMapResolve(int Niter, int Nburn)
{

  RandomiseOrder();
  
  cerr << endl;
  cerr << "Performing Burn-in iterations" << endl;
  cerr << setw(4) << 0 << "% done\033[A" << endl;
  for(int iter = 0; iter < Nburn; iter++){ 
    for(int i = 0; i < Nind; i++){
      int ind = order[i];
      FastHapMapUpdate(ind,true); 
    }
     cerr << setw(4) << (int) ((iter+1)*100.0/Nburn) << "% done\033[A" << endl;
  }
  cerr << endl;
  
  //cout << "Rho =" << RhoMean << endl; 
  // FastHapMapUpdate(0);
//   FastHapMapUpdate(1);
//   exit(1);
  cerr << "Performing Main iterations" << endl;
  cerr << setw(4) << 0 << "% done\033[A" << endl;
  for(int iter = 0; iter < Niter; iter++){
    for(int i = 0; i < Nind; i++){
      int ind = order[i];
      FastHapMapUpdate(ind,false);
    }
    UpdateCounts();
    cerr << setw(4) << (int) ((iter+1)*100.0/Niter) << "% done\033[A" << endl;
  }
  
  cerr << endl;
  

}


void ClassPop::FastHapMapUpdate(int ind, bool burnin)
{ 

  
  MakeHapList(false);
  //haplist.Output(cout,coding);

  int nchr = Nind+Nind;
  int nchrminus2 = nchr-2;
  vector<int>::iterator b = buddy[ind].begin();
  double removeamount = 1.0;
  double addamount = 1.0;

 //  for(int i = 0; i < Nind; i++){
//     cout << "ind: " << i << endl;
//     cout << "buddy: " << buddy[i][0] << endl;
//     cout << "child: " << childindex[i] << endl;
//   }
  
  haplist.SoftRemove(pop[ind],removeamount);
  haplist.SoftRemove(pop[*b],removeamount); 
  haplist.MakePositiveHaps();
  //haplist.Output(cout,coding);

  
  vector<vector<double> > CopyProb1( vector< vector<double> > (Nloci, vector<double>(SS*2,0.0)) );
  vector<vector<double> > CopyProb2( vector< vector<double> > (Nloci, vector<double>(SS*2,0.0)) );
  vector<vector<double> > CopyProb3( vector< vector<double> > (Nloci, vector<double>(SS*2,0.0)) );
  vector<vector<double> > CopyProb4( vector< vector<double> > (Nloci, vector<double>(SS*2,0.0)) );

  vector<int> ignore(Nloci,0);
  for(int i = 0; i< Nloci; i++)
    ignore[i] = 0; //burnin ? pop[ind].is_unknown(i) : 0; // condition on known positions if burnin is true; condition on everything if burnin is false
  haplist.ComputeHiddenStateProbs(CopyProb1,pop[ind].get_haplotype(0), Qptr, nchrminus2, vecRho, true, ignore, vecTheta); 
  haplist.ComputeHiddenStateProbs(CopyProb2,pop[ind].get_haplotype(1), Qptr, nchrminus2, vecRho, true, ignore, vecTheta); 
  
  for(int i = 0; i< Nloci; i++)
    ignore[i] = 0; //burnin ? pop[*b].is_unknown(i) : 0; 
  haplist.ComputeHiddenStateProbs(CopyProb3,pop[*b].get_haplotype(0), Qptr, nchrminus2, vecRho, true, ignore, vecTheta); 
  haplist.ComputeHiddenStateProbs(CopyProb4,pop[*b].get_haplotype(1), Qptr, nchrminus2, vecRho, true, ignore, vecTheta); 
  
 
  for(int locus = 0; locus < Nloci; locus++){
    if(pop[ind].is_unknown(locus)){
      vector<double> Prob(16,0.0);
      int i=0;
      double sum = 0;
      for(int p00 = 0; p00 < 2; p00++){ // pij = allele in true hap j of parent i
        for(int p01 = 0; p01 < 2; p01++){
          for(int p10 = 0; p10 < 2; p10++){
            for(int p11 = 0; p11 < 2; p11++){
              Prob[i] = PriorProbFromCopyProb(p00,locus,CopyProb1) 
                      * PriorProbFromCopyProb(p01,locus,CopyProb2) 
                      * PriorProbFromCopyProb(p10,locus,CopyProb3) 
                      * PriorProbFromCopyProb(p11,locus,CopyProb4) 
                      * pop[ind].PrOrigGenotypeData(locus,p00,p01) 
                      * pop[*b].PrOrigGenotypeData(locus,p10,p11) 
                      * pop[childindex[ind]].PrOrigGenotypeData(locus,p00,p10);
// note this last assumes no recom, as assumes parents transmit p00 and p10 to child
              sum += Prob[i++];
	    }
          }
        }
      }
      if(sum == 0){ 	
         cerr << "Error: Mendelian inconsistency, in ind " << (ind+1) << ", locus " << (locus+1) << endl;
         exit(1);
      }

      //output results
      i=0;
      for(int p00 = 0; p00 < 2; p00++){ // pij = allele in true hap j of parent i
        for(int p01 = 0; p01 < 2; p01++){
          for(int p10 = 0; p10 < 2; p10++){
            for(int p11 = 0; p11 < 2; p11++){
	      cout << ind << "," << locus << ","  << p00 << "," << p01 << "," << p10 << "," << p11 << ":" << Prob[i++]/sum << endl;
  	    }
          }
        }
      }



      int choice = rint2(Prob, sum);
      //cout << "choice:" << choice << endl;

      i=0;
      int p00, p01, p10, p11;
      for(p00 = 0; p00 < 2; p00++)
         for(p01 = 0; p01 < 2; p01++)
	   for(p10 = 0; p10 < 2; p10++)
	     for(p11 = 0; p11 < 2 ; p11++){
	       if(i==choice) goto GOTCHOICE;
	       else i++;
	     }
      
    GOTCHOICE:
      
      pop[ind].update_haplotype(0,locus,p00);
      pop[ind].update_haplotype(1,locus,p01);
      pop[*b].update_haplotype(0,locus,p10);
      pop[*b].update_haplotype(1,locus,p11);
      pop[childindex[ind]].update_haplotype(0,locus,p00);
      pop[childindex[ind]].update_haplotype(1,locus,p10);
      
    }

  }

  haplist.Add(pop[ind],addamount);
  haplist.Add(pop[*b],addamount);
  

}


void ClassPop::HapListImputeMissing(int ind)
{
  ///if(method != 'R'){
  //  cerr << "Error: missing data not implemented for non-R method" << endl;
  //  exit(1);
  //}
  if(pop[ind].NMissingLoci()>0){
    int nchrminus2 = Nind + Nind -2;
    int listlength = haplist.get_positivelistlength();
    vector< vector<double> > Alpha( vector< vector<double> > (Nloci, vector<double>(SS*listlength,0.0)) );
    vector<double> AlphaSum( Nloci );
    vector<int> copiedtime ( Nloci ); // stores time at which each locus copied
    vector<int> copiedallele ( Nloci ); // stores allele each locus copied
    vector<int> copiedhap (Nloci); // stores hap in haplist each locus copied
    //vector<int> missing = pop[ind].get_missing_list(); 
    
    vector< vector<int> > nmissing(2,pop[ind].get_nmissing());
 
    vector< vector<int> > missing(2); // list of positions considered missing in this chromosome

 // set up nmissing as a vector for that particular chromosome, to say whether or not
      // each allele is missing or not (to deal with case where one allele is missing)
      
    for(int locus = 0; locus < Nloci; locus++){
      if(pop[ind].nmissing(locus)==2){
	nmissing[0][locus] = 1; // this bit unnecessary, but makes it more obvious what is going on
	nmissing[1][locus] = 1; // this bit unnecessary, but makes it more obvious what is going on
      }

      if(pop[ind].nmissing(locus)==1){ // if only one missing, if only
	// one of the two alleles matches the known allele, the other 
	// must be the "missing" allele; otherwise choose one
	// at random to be missing

	if(!pop[ind].is_unknown(locus)){ // deal with phase known case: which is missing is unambiguous
	  nmissing[0][locus] = 0;
	  nmissing[1][locus] = 0;
	  nmissing[pop[ind].missingchr(locus)][locus] = 1;
	}
	else if(pop[ind].get_haplotype(0, locus) != pop[ind].get_orig_nonmissing_allele(locus)){
	  nmissing[1][locus] = 0;
	  nmissing[0][locus] = 1;
	}
	else if(pop[ind].get_haplotype(1, locus) != pop[ind].get_orig_nonmissing_allele(locus)){
	  nmissing[0][locus] = 0;
	  nmissing[1][locus] = 1;
	}
	else{
	  // choose one at random to be the missing one
	  int c = (ranf()<0.5);
	  nmissing[c][locus] = 0;
	  nmissing[1-c][locus] = 1;
	}
	if((pop[ind].get_haplotype(0,locus) != pop[ind].get_orig_nonmissing_allele(locus)) && (pop[ind].get_haplotype(1,locus) != pop[ind].get_orig_nonmissing_allele(locus)) ){
	  cout << "Error: haplotypes in individual " << ind << ", locus " << locus << " are not consistent with inputs" << endl;
	  exit(1);
	}
      }
      
      if(nmissing[0][locus]>0)
	missing[0].push_back(locus);
      if(nmissing[1][locus]>0)
	missing[1].push_back(locus);
      
    }

    
    
    // THIS PART TO BE CHANGED TO ALLOW FOR GENOTYPING ERROR?
    for(int chr = 0; chr < 2; chr++){
      haplist.ForwardsAlgorithm(pop[ind].get_haplotype(chr), Qptr, nchrminus2, vecRho, Alpha, AlphaSum, true, nmissing[chr], false);  
#ifdef DEBUG        
        cout << "Alpha values: " << endl;
            for(int locus = 0; locus<Nloci; locus++){
              for(int j=0; j<Alpha[locus].size(); j++)
        	cout << locus << "," << j << " : " << Alpha[locus][j] << endl;
              cout << "Sum =" << AlphaSum[locus] << endl;
            }
#endif
      
      haplist.BackwardsAlgorithm(pop[ind].get_haplotype(chr), nchrminus2, vecRho, Alpha, AlphaSum, copiedtime, copiedallele, copiedhap);

#ifdef DEBUG      
      cout << "Individual: " << ind << ", chromosome " << chr <<  endl;
      for(int i = 0; i<copiedallele.size(); i++)
        cout << copiedallele[i] << " , " << copiedtime[i] << endl;
#endif 
     
      for(vector<int>::iterator u = missing[chr].begin(); u!=missing[chr].end(); u++){
	//cout << "Calling impute_allele with " << *u << "," << nchrminus2 << "," << copiedtime[*u] << "," << copiedallele[*u] << endl;
	
	int newallele = impute_allele (*u, nchrminus2, copiedtime[*u], copiedallele[*u]);
	pop[ind].update_haplotype (chr, *u, newallele);
	
      }
    }
  }
}


// }}}

// {{{ Log
//
// $Log: classpop.cpp,v $
// Revision 1.30  2003/06/14 00:24:04  stephens
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
// Revision 1.25  2001/10/16 18:12:05  stephens
// More minor bugs fixed. (including -m bug)
//
// Revision 1.23  2001/10/09 23:03:25  stephens
// Modified to deal with very large files - put in checks for underflow
// and overflow, and changed computation of Q, DiffProb etc to only compute those elements that are necessary.
//
// Revision 1.22  2001/10/06 00:35:40  stephens
// Impute missing positions added- changed update_allele to update_haplotype.
//
// Revision 1.21  2001/10/02 18:02:58  stephens
// just some small additions. haven't committed these for a long time
// because cvs was not working propersly. (now know why: they remounted
// you from /user3 to /home)
//
// Revision 1.20  2001/06/23 15:57:00  stephens
// started to check missing data code. Corrected some bugs. Seems to give sensible results in a couple of small examples.
//
// Revision 1.19  2001/06/19 17:02:25  stephens
// Changes to computation of arrayFF to make more efficient
// Added facility to store "original phenotype" in indnode,
// in preparation for allowing genotyping error.
//
// Revision 1.18  2001/05/31 21:14:11  nali
//
// Implenting ArrayCC as double (*veccc)[2][SS][2] rather than double ****veccc.
//
// Revision 1.17  2001/05/31 16:26:14  stephens
// Added DiploidDiffProb look-up table class (to be used to make
// computation of arrayFF for SNPs more efficient)
//
// Revision 1.16  2001/05/30 06:02:18  stephens
// Updated to be considerably more efficient, via introduction of
// various new methods, including introduction of ArrayDiffProb
// and ArrayDiffCount to improve computation al efficiency for SNP
// data. Speedtest.inp now runs in about 7secs.
// Also corrected several bugs. Output now looks more promising
// and major bugs appear to have been eliminated. Convergence of chain
// can now be monitored more easily by the output in temp.monitor, which
// gives the pseudo-likelihood every Nthin repetitions.
//
// Revision 1.15  2001/05/23 02:16:55  stephens
// Corrected some bugs in computation of CC. Reduced number of
// loci we update by one when the total number of unknowns is <=5.
// Confirmed program was giving vaguely sensible resutls, and is
// competitive on the speed test example (36secs vs 31 for old version)
//
// Revision 1.14  2001/05/21 20:17:15  nali
// No real changes
//
// Revision 1.13  2001/05/18 21:47:42  stephens
// added facility to update just 5 at a time. Set it to do this
// so we can test for speed with the current version.
//
// Revision 1.12  2001/05/08 16:58:23  stephens
// Added class arrayCC, which computes quantities similar
// to arrayFF, but only for sites whose phase is known.
// This is then used to improve the efficiency of computation
// for arrayFF by re-using the calculations for known sites.
// Checked that the output was unchanged on an example.
// Might slightly improve efficiency of computation of CC if
// we do the loop through known sites more efficiently.
//
// Revision 1.11  2001/04/27 00:20:48  stephens
// First version of postprocess done - usage is
// postprocess -n file.in file.sum dummyfilename
//
// Revision 1.10  2001/04/26 18:29:51  stephens
// Fixed bug in computation of ArrayFF (total_sum now initialized
// to 0 each time FF is computed)
//
// Revision 1.9  2001/04/26 17:09:17  stephens
// Added function ClassPop::print_allele
// and made a few minor changes
//
// Revision 1.8  2001/04/24 22:00:08  stephens
// Added flag in ClassPop::output_hap and output_phase
// to indicate whether to output all positions, or only
// unknown positions.
//
// Revision 1.7  2001/04/24 20:25:07  stephens
// Minor changes to renormalize function
//
// Revision 1.6  2001/04/24 19:31:31  nali
// Move data input out of the constructor. Member functions read_data and initialize
// have to be called explicitly. Put everything related to hudson data set into a single
// function so that it might go away one day.
//
// Revision 1.5  2001/04/23 18:55:16  stephens
// Added member function ClassPop::renormalize(newcoding0)
// to renormalize a population; for use when postprocessing
// results to make sure all samples have been normalized the
// same way.
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
// "mult" removed in update_phase_NR No other changes in terms of algorithm.
// Haven't check the results yet.
//
// proc_args() is moved to utility.cpp, which also defines a couple of other
// global functions.
//
// }}}
