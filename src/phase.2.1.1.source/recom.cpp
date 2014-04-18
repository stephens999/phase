#include "classpop.hpp"
#include "constants.hpp"
#include "utility.hpp"
#include "arrayDiffProb.hpp"
#include "arrayDiffCount.hpp"
#include "arrayCC.hpp"
#include "indnode.hpp"

#include <iomanip>
#include <algorithm>
#include <cmath>

using namespace std;

void ClassPop::resolve_phase_R ( double thetainit, double totalrho,
                               int Nburn,
                               int Nthin,
                               int Niter,
                               ostream & ostrHaplotypes,
                               //ostream & ostrPhases,
			       ostream & ostrMonitor,
			       vector<double> & vecDelta, 
			       int verbose)
{

    // Mutation rate at each locus
    vecTheta = vector<double> ( loci_type.size(), thetainit );
   // Calculate theta if necessary
    if ( thetainit <= 0.0 ) calc_theta ();
        
    vector <double> rho;
    int Nloci = loci_type.size();

    //cout << "Computing rho" << endl;
    for(int r=0; r<(Nloci-1); r++){
      rho.push_back(totalrho * (position[r+1]-position[r]));
    }

   
 // {{{ Burn-in
    
    cout << "Beginning Burn-in" << endl;
 
    double thetastart = vecTheta[0];
    double thetafinish = vecTheta[0];
    double thetalambda = log(thetastart/thetafinish);
    double theta;

    int nchr = 2*Nind -1;

    // these params were working well before I put in diploid conditional sim
//      double deltastart = 10 * 0.99*0.5*thetastart/(nchr-1.0+thetastart);
//      double deltafinish = 0.1 * 0.99*0.5*thetastart/(nchr-1.0+thetastart);
//      double deltalambda = log(deltastart/deltafinish);
//      double delta;

//      double temperaturestart = 1;
//      double temperaturefinish = 0.05; 
//      double temperaturelambda = log(temperaturestart/temperaturefinish);
//      double temperature;
 
    double deltastart = 10 * 0.99*0.5*thetastart/(nchr-1.0+thetastart);
    double deltafinish = 0.1* 0.99*0.5*thetastart/(nchr-1.0+thetastart);
    double deltalambda = log(deltastart/deltafinish);
    double delta;

    double temperaturestart = 1; //10; 
    double temperaturefinish = 1; //0.01; 
    double temperaturelambda = log(temperaturestart/temperaturefinish);
    double temperature;
 
    double currentlogprob;
    double currentbest=-1e100;
    double lastlogprob=-1e100;

    SaveCurrentState();
    double oldlogprob = -1e100;

    for (int iter = 0; iter < Nburn; ++iter) {
      if ( iter % ITPRINT == 0 ) {
	cout << "Burn-in " << iter << endl;
      }
      delta = deltastart * exp(-deltalambda * iter / (Nburn));
      theta = thetastart * exp(-thetalambda * iter / (Nburn));
      temperature = temperaturestart * 	exp(-temperaturelambda * iter / (Nburn));

      //cout << "delta = " << delta << "; theta = " << theta << "; temperature = " << temperature << endl;

      for (int thin = 0; thin < Nthin; ++thin) {	
  	//currentlogprob = update_R_SA  ( rho, theta, delta, temperature, 1, Q );
	currentlogprob = update_R ( rho, theta, delta, temperature, 1);
	
	//cout << currentlogprob << endl; // << "old: " << lastlogprob << endl;
	ostrMonitor << currentlogprob << endl;
	//output( ostrHaplotypes, true, true );

	//double u = ranf();
	//if(u<exp((currentlogprob-lastlogprob+1))){
	//  if(currentlogprob>lastlogprob){
//  	  lastlogprob = currentlogprob;
//  	  SaveCurrentState();
//  	  oldlogprob = currentlogprob;
//  	}
//  	else{ 
//  	  RestoreSavedState();
//  	  currentlogprob = oldlogprob;
//  	}
	
      //      while(currentlogprob<currentbest && currentlogprob!=lastlogprob){
      //	lastlogprob = currentlogprob;
      //	currentlogprob = update_R  ( rho, theta, delta );
      //	cout << currentlogprob << endl;
      //}
      //currentbest = currentlogprob;
	//}
      }
      currentlogprob = update_R ( rho, theta, delta, temperature, 1);
      ostrMonitor << currentlogprob << endl;

    }
 
    temperaturestart = 1; 
    temperaturefinish = 0.01; //0.01; 
    temperaturelambda = log(temperaturestart/temperaturefinish);
 
    for (int iter = 0; iter < Nburn; ++iter) {
      if ( iter % ITPRINT == 0 ) {
	cout << "Burn-in 2 " << iter << endl;
      }
      //delta = deltastart * exp(-deltalambda * iter / (Nburn));
      //theta = thetastart * exp(-thetalambda * iter / (Nburn));
      temperature = temperaturestart * 	exp(-temperaturelambda * iter / (Nburn));

      //cout << "delta = " << delta << "; theta = " << theta << "; temperature = " << temperature << endl;

      for (int thin = 0; thin < Nthin; ++thin) {	
  	//currentlogprob = update_R_SA ( rho, theta, delta, temperature, 1, Q );
	currentlogprob = update_R ( rho, theta, delta, temperature, 1);
	//cout << currentlogprob << endl; // << "old: " << lastlogprob << endl;
	ostrMonitor << currentlogprob << endl;

      }
      currentlogprob = update_R ( rho, theta, delta, temperature, 1);
      ostrMonitor << currentlogprob << endl;
    }

    // }}}
    // {{{ Main loop

    //ostrHaplotypes << Niter << endl;

    delta = deltafinish;
    temperature = temperaturefinish;
    theta = thetafinish;

    cout << "Beginning main iterations" << endl;

    for (int iter = 0; iter < Niter; ++iter) {
      if ( iter % ITPRINT == 0 ) {
	cout << "Iteration " << iter << endl;
      }
      for (int thin = 0; thin < Nthin; ++thin) {	
	//currentlogprob = update_R_SA ( rho, theta, delta, temperature, 1, Q );
	currentlogprob = update_R ( rho, theta, delta, temperature, 1);
	//cout << currentlogprob << endl;
      }
      
      UpdateCounts();
      ostrMonitor << currentlogprob << endl;
      
      //output_all_phases ( ostrPhases );
      //output( ostrHaplotypes, true, true );
    }
}

double ClassPop::update_R_SA (vector<double> & rho, double theta, double delta, double temperature, int switchallowed)
{
  
  static double time =0;
  
  static vector <double> * FF;
  if(time==0){
    FF = new vector <double> [400]; //[loci_type.size()];
    for(int r=0; r<400; r++){
      FF[r] = vector<double> (2*SS*Nind,0.0);
    }
  }

  static vector <int> CopiedInd [2]; // holds which ind was copied by each hap
  static vector <int> CopiedChr [2]; // holds which chr was copied by each hap
  static vector <int> CopiedTime [2]; // holds time each chr was copied by each hap

  if(time==0){
    CopiedInd[0] = vector <int> (400,0);
    CopiedInd[1] = vector <int> (400,0);
    CopiedChr[0] = vector <int> (400,0);
    CopiedChr[1] = vector <int> (400,0);
    CopiedTime[0] = vector <int> (400,0);
    CopiedTime[1] = vector <int> (400,0);
  }

  time ++;
  //cout << "time:" << time << endl;

//    double logprob = 0;
//    for(int ind =0; ind< Nind; ind++){
//      for(int chr=0; chr<2; chr++){
//        logprob+=FastForwardsAlg ( ind, chr, theta, delta, FF, rho );
//        FastBackwardsAlg ( ind, chr, theta, delta, FF, rho, temperature, CopiedInd[chr], CopiedChr[chr] );
//      }

  double logprob = 0;
  for(int ind =0; ind< Nind; ind++){
    double oldlogprob = 0;
    double newlogprob = 0;
    
    for(int chr=0; chr<2; chr++){
      oldlogprob+=FastForwardsAlg ( ind, chr, theta, delta, FF, rho);
      FastBackwardsAlg ( ind, chr, theta, delta, FF, rho, temperature, CopiedInd[chr], CopiedChr[chr], CopiedTime[chr]);
    } 
    if(switchallowed)
      ConditionalDiploidSim(ind,theta,delta,temperature,CopiedInd,CopiedChr,CopiedTime);
  
    for(int chr=0; chr<2; chr++){
      newlogprob+=FastForwardsAlg ( ind, chr, theta, delta, FF, rho);     
    }
    //cout << "old: " << oldlogprob << ", new: " << newlogprob << endl;
    if(ranf() < exp( ((1/temperature) -1) * (newlogprob - oldlogprob - 1e-2) )){
      logprob += newlogprob;
      SaveCurrentState(ind);
    }
    else {
      logprob += oldlogprob;
      RestoreSavedState(ind);
    }

   //   cout << "Ind:" << ind << endl;
   //        for(int r=0;r<loci_type.size();r++)
   //          cout << pop[CopiedInd[0][r]].get_haplotype(CopiedChr[0][r],r) << " ";
//        cout << endl;
//        for(int r=0;r<loci_type.size(); r++)
//          cout << pop[CopiedInd[1][r]].get_haplotype(CopiedChr[1][r],r) << " ";
//        cout << endl;
 
  
  }

  return logprob;

  //cout << n << "," << c << endl;
  // for(int i=0; i<10;i++){
  //   for(int j=0;j<10;j++){
  //     cout << setw(2) << FF[i][j] << " ";
  //   }
  //  cout << endl;
  //}

  //output_all_haps    ( cout  , true , true, false, true, 0.5, 0.5);
}


double ClassPop::update_R (vector<double> & rho, double theta, double delta, double temperature, int switchallowed)
{
  
  static double time =0;
  
//    static vector <double> * FF;
//    if(time==0){
//      FF = new vector <double> [400]; //[loci_type.size()];
//      for(int r=0; r<400; r++){
//        FF[r] = vector<double> (2*SS*Nind,0.0);
//      }
//    }

//    static vector <int> CopiedInd [2]; // holds which ind was copied by each hap
//    static vector <int> CopiedChr [2]; // holds which chr was copied by each hap
//    static vector <int> CopiedTime [2]; // holds time each chr was copied by each hap

//    if(time==0){
//      CopiedInd[0] = vector <int> (400,0);
//      CopiedInd[1] = vector <int> (400,0);
//      CopiedChr[0] = vector <int> (400,0);
//      CopiedChr[1] = vector <int> (400,0);
//      CopiedTime[0] = vector <int> (400,0);
//      CopiedTime[1] = vector <int> (400,0);
//    }

//    time ++;


  static vector <vector<double> > FF (400);
  
  static vector <vector<int> > CopiedInd (2); // holds which ind was copied by each hap
  static vector <vector<int> > CopiedChr (2); // holds which chr was copied by each hap
  static vector <vector<int> > CopiedTime (2); // holds time each chr was copied by each hap

  
  if(time == 0){
    for(int r=0; r<400; r++){
      FF[r] = vector<double> (2*SS*Nind*2*SS*Nind,0.0);
    }
    CopiedInd[0] = vector <int> (400,0);
    CopiedInd[1] = vector <int> (400,0);
    CopiedChr[0] = vector <int> (400,0);
    CopiedChr[1] = vector <int> (400,0);
    CopiedTime[0] = vector <int> (400,0);
    CopiedTime[1] = vector <int> (400,0);
    time ++;
  }


   
  double logprob = 0;
  for(int ind =0; ind < Nind; ind++){
    cout << "ind:" << ind << endl;
    if(pop[ind].numunknown()>0){
      logprob += DiploidForwardsAlg ( FF, ind, rho );
      DiploidBackwardsAlg ( ind, theta, delta, FF, rho, temperature, CopiedInd, CopiedChr, CopiedTime);    
    }

    // print out copied alleles
    //  for(int r=0;r<loci_type.size();r++)
//        cout << CopiedInd[0][r] << "," << CopiedChr[0][r] << " : ";
//      cout << endl;
//      for(int r=0;r<loci_type.size();r++)
//        cout << CopiedInd[1][r] << "," << CopiedChr[1][r] << " : ";
//      cout << endl;

//      for(int r=0;r<loci_type.size();r++)
//        cout << pop[CopiedInd[0][r]].get_haplotype(CopiedChr[0][r],r) << " ";
//      cout << endl;    
//      for(int r=0;r<loci_type.size(); r++)
//        cout << pop[CopiedInd[1][r]].get_haplotype(CopiedChr[1][r],r) << " ";
//      cout << endl;    
  }

  return  logprob;

  //cout << n << "," << c << endl;
  // for(int i=0; i<10;i++){
  //   for(int j=0;j<10;j++){
  //     cout << setw(2) << FF[i][j] << " ";
  //   }
  //  cout << endl;
  //}

  //output_all_haps    ( cout  , true , true, false, true, 0.5, 0.5);
}




//
// update alleles in individual n, conditional on having copied
// chromosomes given in CopiedInd and CopiedChr at each position
//
void ClassPop::ConditionalDiploidSim(int n, double theta, double delta, double temperature,vector<int> * CopiedInd, vector<int> * CopiedChr, vector<int> * CopiedTime)
{
  for(int r=0; r < loci_type.size(); r++){
    int observedallele0 = pop[n].get_orig_allele(0,r);
    int observedallele1 = pop[n].get_orig_allele(1,r);
   
    int fromallele0 = pop[CopiedInd[0][r]].get_haplotype(CopiedChr[0][r],r);
    int fromallele1 = pop[CopiedInd[1][r]].get_haplotype(CopiedChr[1][r],r);
    int time0 = CopiedTime[0][r];
    int time1 = CopiedTime[1][r];

    vector<double> tempprob = vector<double>(4,0.0);
    
    int nchr = 2*Nind - 1;
    double pnomut0 = PrHitTarg(r,nchr,time0,fromallele0,fromallele0,Qptr);
    double pmut0 = PrHitTarg(r,nchr,time0,fromallele0,1-fromallele0,Qptr);
    double pnomut1 = PrHitTarg(r,nchr,time1,fromallele1,fromallele1,Qptr);
    double pmut1 = PrHitTarg(r,nchr,time1,fromallele1,1-fromallele1,Qptr);

    tempprob[0] = pnomut0*pnomut1*PrObserveGivenTruth(observedallele0, observedallele1, fromallele0, fromallele1, delta);
    tempprob[1] = pnomut0*pmut1*PrObserveGivenTruth(observedallele0, observedallele1, fromallele0, 1 - fromallele1, delta);
    tempprob[2] = pmut0*pnomut1*PrObserveGivenTruth(observedallele0, observedallele1, 1-fromallele0, fromallele1, delta);
    tempprob[3] = pmut0*pmut1*PrObserveGivenTruth(observedallele0, observedallele1, 1-fromallele0, 1-fromallele1, delta);

    double sum = tempprob[0]+tempprob[1]+tempprob[2]+tempprob[3];
    tempprob[0] /=sum;
    tempprob[1] /=sum;
    tempprob[2] /=sum;
    tempprob[3] /=sum;
    
   //   if(temperature!=1){
//        tempprob[0] = exp(log(tempprob[0])/temperature);
//        tempprob[1] = exp(log(tempprob[1])/temperature);
//        tempprob[2] = exp(log(tempprob[2])/temperature);
//        tempprob[3] = exp(log(tempprob[3])/temperature);
//      } 

   //   sum = tempprob[0]+tempprob[1]+tempprob[2]+tempprob[3];
//      tempprob[0] /=sum;
//      tempprob[1] /=sum;
//      tempprob[2] /=sum;
//      tempprob[3] /=sum;
    
    //cout << tempprob[0] << ":" << tempprob[1] << ":" << tempprob[2] << ":" << tempprob[3] << ":";

    int choice = rint2(tempprob);
    //cout << "prob" << tempprob[choice] << endl;

    //cout << choice << endl;

    if(choice==0){
      pop[n].update_haplotype(0,r,fromallele0);
      pop[n].update_haplotype(1,r,fromallele1);          
    }
    if(choice==1){
      pop[n].update_haplotype(0,r,fromallele0);
      pop[n].update_haplotype(1,r,1-fromallele1);          
    }
    if(choice==2){
      pop[n].update_haplotype(0,r,1-fromallele0);
      pop[n].update_haplotype(1,r,fromallele1);          
    }    
    if(choice==3){
      pop[n].update_haplotype(0,r,1-fromallele0);
      pop[n].update_haplotype(1,r,1-fromallele1);          
    }

  }
  
}

//
// version of the Forwards-backwards alg based on SLFD conditional,
// which keeps one chrom in individual fixed while changing
// the other.
// currently implemented only for SNPs
//
// theta: "mutation" parameter
// delta: genotyping error prob
// FF[r][ptr]: the prob of copying ind ptr at position r
// n: the individual to be updated
// c: the chromosome to be updated 
// rho: the recom param between each site
//
double ClassPop::FastForwardsAlg ( int n, int c, double theta, double delta, vector <double> * FF, vector<double> & rho)
{
    // int n0,c0;
    // these represent the haplotype which (n,c) copied
    // int m0,d0;
    // these represent the same - used as looping variables
  int nchr = 2*Nind - 1;
  int ptr = 0;
  double sum; // sum of FF[r]

    // set up FF[0]
  for (int n0 = 0; n0 < Nind; ++n0) {
    for (int c0 = 0; c0 < 2; ++c0) {  
      for(int t0 = 0; t0 < SS; ++t0) {
	
	if (n0 == n && c0 == c) {
	  FF[0][ptr] = 0.0; // can't copy itself
	} else {
	  FF[0][ptr] = WEIGHTS[t0];
	}
	
	int r = 0;
      
	int observedallele0 = pop[n].get_orig_allele(0,r);
	int observedallele1 = pop[n].get_orig_allele(1,r);
      
	// the allele currently imputed on the other strand
	int imputedallele = pop[n].get_haplotype(1-c,r);
	// the allele begin copied
	int fromallele = pop[n0].get_haplotype(c0,r);
	
	FF[0][ptr] *= PrObserve(observedallele0, observedallele1, imputedallele, fromallele, nchr, theta, delta, t0, Qptr, r);
	sum += FF[0][ptr];
	ptr++; 
      }
    }
  }
	
  //normalize FF to avoid numerical problems 
  ptr=0;
  for (int n0 = 0; n0 < Nind;  ++n0) {
    for (int c0 = 0; c0 < 2;  ++c0) {
      for (int t0 = 0; t0 < SS; ++t0) {
	  FF[0][ptr++]/=sum;	
      }
    }
  }

  double logsum = log(sum);// log of sum of final FF[r] (if hadn't been normalised)

  for(int r = 1; r < loci_type.size(); ++r) {
    ptr=0;
    int observedallele0 = pop[n].get_orig_allele(0,r);
    int observedallele1 = pop[n].get_orig_allele(1,r);
    int imputedallele = pop[n].get_haplotype(1-c,r);
    sum = 0; // sum of FF[r]

    double exp1 = exp(-rho[r-1]/nchr);

    for (int n0 = 0; n0 < Nind;  ++n0) {
      for (int c0 = 0; c0 < 2;  ++c0) {
	int fromallele = pop[n0].get_haplotype(c0,r);
	for(int t0 = 0; t0 < SS; ++t0) {
	  if (  n0 != n  || c0!=c ) {
	    FF[r][ptr]=((1-exp1) * (WEIGHTS[t0]/nchr) + 
			FF[r-1][ptr] * exp1) * PrObserve(observedallele0, observedallele1, imputedallele, fromallele, nchr, theta, delta, t0, Qptr, r ) ;
	    sum+=FF[r][ptr];
	  }
	  else {
	    FF[r][ptr]= 0;
	  }
	  ptr++;
	  // now calculate FF[m] as a sum over FF[m-1]
	
//  	  if (  n0 != n  || c0!=c ) {
//  	    int fromallele = pop[n0].get_haplotype(c0,r);
//  	    int ptr2 = 0;
//  	    for (int m0 = 0; m0 < Nind;  ++m0) {
//  	      for (int d0 = 0; d0 < 2;  ++d0) {
//  		for (int s0 = 0; s0 < SS; ++s0) {
//  		  if( m0 != n || d0!=c ){
//  		    FF[r][ptr] += FF[r-1][ptr2]
//  		      * TransitionProb(nchr,m0,d0,s0,n0,c0,t0,rho[r-1])
//  		      * PrObserve(observedallele0, observedallele1, imputedallele, fromallele, nchr, theta, delta, t0, Qptr, r );
//  		  }		
//  	      ptr2++;
//  		}
//  	      }
//  	    }	    
//  	  }
	
	} 
      }
    }
      
    //normalize FF to avoid numerical problems 
    ptr=0;
    for (int n0 = 0; n0 < Nind;  ++n0) {
      for (int c0 = 0; c0 < 2;  ++c0) {
	for (int t0 = 0; t0 < SS; ++t0) {
	  FF[r][ptr++]/=sum;	
	}
      }
    }
    logsum+=log(sum);
  }
  
  return logsum;
}

void ClassPop::FastBackwardsAlg ( int n, int c, double theta, double delta, vector<double > * FF, vector<double> & rho, double temperature, vector <int> & CopiedInd, vector<int> & CopiedChr, vector<int> & CopiedTime)
  // backwards part
{
  int nchr = 2*Nind - 1;

  for(int r = loci_type.size()-1; r >= 0; r--) {
    // simulate n0 and c0 from FF[r]


 //     if(temperature !=1){
//        double ftemperature;
//        if(temperature<0.01)
//  	ftemperature = 0.01;
//        else
//  	ftemperature = temperature;
//        for(int u=0;u<FF[r].size();u++){
//  	FF[r][u] = exp(log(FF[r][u])/ftemperature);
//        }
//      }

    int temp = rint2(FF[r]);
    //cout << temp << endl;

    int ptr =0;
    
    int ncopy, ccopy, tcopy;

    for(ncopy = 0; ncopy < Nind; ncopy++) {
      for(ccopy = 0; ccopy < 2; ccopy++) {
	for(tcopy = 0; tcopy < SS; tcopy++) {
	  if(ptr == temp)	goto BACKENDLOOP1;
	  ptr++;
	}
      }
    }   
  BACKENDLOOP1:
    
    CopiedInd[r] = ncopy;
    CopiedChr[r] = ccopy;
    CopiedTime[r] = tcopy;

    // cout << ncopy << "," << ccopy << endl;
    int observedallele0 = pop[n].get_orig_allele(0,r);
    int observedallele1 = pop[n].get_orig_allele(1,r);    
    int imputedallele = pop[n].get_haplotype(1-c,r);

    int fromallele = pop[ncopy].get_haplotype(ccopy,r);

    vector<double> tempprob = vector<double>(2,0.0);

    tempprob[0] = PrHitTarg(r,nchr,tcopy,fromallele,fromallele,Qptr) * PrObserveGivenTruth(observedallele0, observedallele1, imputedallele, fromallele, delta);

    tempprob[1] = PrHitTarg(r,nchr,tcopy,fromallele,1-fromallele,Qptr) * PrObserveGivenTruth(observedallele0, observedallele1, imputedallele, 1-fromallele, delta);

    //    cout << "tempprob:" << tempprob[0] << "," << tempprob[1] << ":";  
    
    double sum = tempprob[0]+tempprob[1];
    tempprob[0] /=sum;
    tempprob[1] /=sum;

    //  if(temperature!=1){
//        tempprob[0] = exp(log(tempprob[0])/temperature);
//        tempprob[1] = exp(log(tempprob[1])/temperature);
//      } 
    
    int choice = rint2(tempprob);
    
    //cout << choice << endl;

    if(choice==0){
      pop[n].update_haplotype(c,r,fromallele);
      //cout << observedallele0 << observedallele1 << imputedallele << fromallele << endl;
    }
    else{
      pop[n].update_haplotype(c,r,1-fromallele);
      //cout << observedallele0 << observedallele1 << imputedallele << 1-fromallele << endl;
    }

    // modify FF[r-1] to reflect the realised value of ncopy,ccopy,tcopy
    if(r>0){
      double expsave = exp(-rho[r-1]/nchr); // compute this outside loop for efficiency
      ptr=0;
      double sum=0;
      for (int m0 = 0; m0 < Nind;  ++m0) {
	for (int d0 = 0; d0 < 2;  ++d0) {
	  for(int s0 = 0; s0 <SS; ++s0) {
	    if( m0 != n || d0!=c ){
	      FF[r-1][ptr] *= TransitionProb(nchr,ncopy,ccopy,tcopy,m0,d0,s0,rho[r-1], expsave);
	    }  
	    sum+=FF[r-1][ptr];
	    ptr ++;
	  }
	}
      }
      ptr=0;
      for (int m0 = 0; m0 < Nind;  ++m0) {
	for (int d0 = 0; d0 < 2;  ++d0) {
	  for(int s0 = 0; s0 < SS; ++s0) {
	    FF[r-1][ptr]/=sum;  
	    ptr ++;
	  }
	}
      }
    }
  }

}




// Prob of observing pair (observed0,observed1) given
// truth is pair (true0,true1)
double PrObserveGivenTruth(int observed0, int observed1, int true0, int true1, double delta)
{
  double prob;

  if(observed0==observed1){
    if(true0 == true1){
      if(observed0==true0)
	prob = (1-delta) * (1-delta);
      else
	prob = delta * delta;
    } else {
      prob = delta * (1-delta);
    }
  } else {
    if(true0 == true1)
      prob = 2*delta*(1-delta);
    else
      prob = delta*delta + (1-delta)*(1-delta);
  }

  return prob;
}

// Prob of observing the pair (observed0,observed1)
// when the opposite strand is imputedallele, and the copied
// allele is fromallele at time time
// currently implemented for SNPs only
double PrObserve(int observed0,int observed1,int imputedallele,int fromallele, int nchr, double theta, double delta, int time, vector<ArrayQ *> & Q, int locus)
{
  // case 1:other true allele is fromallele (prob = n-1/(n-1+theta))
  double prob1=PrObserveGivenTruth(observed0,observed1,imputedallele,fromallele,delta);
  // case 2:other true allele is not fromallele (prob = theta/(n-1+theta))
  double prob2=PrObserveGivenTruth(observed0,observed1,imputedallele,1-fromallele,delta);
  
  return PrHitTarg(locus,nchr,time,fromallele,fromallele,Q)*prob1+PrHitTarg(locus,nchr,time,fromallele,1-fromallele,Q)*prob2;

  //  return ((nchr-1.0)/(nchr-1.0+theta)+0.5*(theta/(nchr-1.0+theta))) * prob1 
  //    + 0.5 * (theta/(nchr-1.0+theta)) * prob2;
}

// Diploid Forwards Algorithm
double ClassPop::DiploidForwardsAlg (vector<vector<double> > & FF,
                                    int n,
                                    vector<double> & totalrho )
{
    // int n0,c0,t0,n1,c1,t1;
    // these represent the individual, chromosome, and time at which the
    // haplotype n is most closely related to the others
    
  // NOTE: NOT YET DONE TO TAKE KNOWN PHASES OR MISSING DATA INTO ACCOUNT

  //  
  // Sum0[whichsum] stores sum over n0,c0,t0 of old FF, 
  // Sum0[1-whichsum] stores same for new FF,
  // Sum1[whichsum], Sum1[1-whichsum] does same for n1,c1,t1
  // Sum[whichsum], Sum[1-whichsum] stores total sums for FF
  vector<vector<double> > Sum0(2, vector<double> (2*Nind*SS, 0.0));
  vector<vector<double> > Sum1(2, vector<double> (2*Nind*SS, 0.0));
  vector<double> Sum(2,0);
  int whichsum=0; 

  double logprob = 0;

  // set up FF[0]
  int ptr0, ptr1;
  ptr0 = 0;

  vector<double>::iterator ptr = FF[0].begin();

  for (int n0 = 0; n0 < Nind; ++n0) {
    for (int c0 = 0; c0 < 2; ++c0) {
      for (int t0 = 0; t0 < SS; ++t0) {
	ptr1 = 0;
	for (int n1 = 0; n1 < Nind; ++n1) {
	  for (int c1 = 0; c1 < 2; ++c1) {
	    for (int t1 = 0; t1 < SS; ++t1) {
	      if (n0 == n || (n1 == n && c1 == 1)) {
		*ptr = 0.0;
	      } else {
		*ptr = WEIGHTS[t0]*WEIGHTS[t1];
	      }
	      int r = 0;
	      int from0 = pop[n0].get_haplotype(c0, r);
	      int from1 = pop[n1].get_haplotype(c1, r);
	      int targ0 = pop[n].get_allele(0, r);
	      int targ1 = pop[n].get_allele(1, r);
	      if ( n1 == n ) {
		*ptr *=
		  PrHitTarg(r, 2*Nind-2, t0, from0, targ0, Qptr) *
		  PrHitTarg(r, 2*Nind-1, t1, targ0, targ1, Qptr) +
		  PrHitTarg(r, 2*Nind-2, t0, from0, targ1, Qptr) *
		  PrHitTarg(r, 2*Nind-1, t1, targ1, targ0, Qptr);
	      } else {
		*ptr *=
		  PrHitTarg(r, 2*Nind-2, t0, from0, targ0, Qptr) *
		  PrHitTarg(r, 2*Nind-1, t1, from1, targ1, Qptr) +
		  PrHitTarg(r, 2*Nind-2, t0, from0, targ1, Qptr) *
		  PrHitTarg(r, 2*Nind-1, t1, from1, targ0, Qptr);
	      }
	      Sum0[whichsum][ptr1++] += *ptr;
	      Sum1[whichsum][ptr0] += *ptr;
	      Sum[whichsum] += *ptr++;
	    }
	  }
	}
	ptr0++;
      }
    }
  }

  logprob += log(Sum[whichsum]);
  
  for(int r = 1; r < loci_type.size(); ++r) {
   
    // compute the sums for F(r-1) over (n0,c0,t0), (n1,c1,t1), and both
	
    //  ptr = ptr0 = 0; 
//      for (int n0 = 0; n0 < Nind;  ++n0) {
//        for (int c0 = 0; c0 < 2;  ++c0) {
//          for (int t0 = 0; t0 < SS; ++t0) {
//  	  ptr1 = 0;
//  	  for (int n1 = 0; n1 < Nind;  ++n1) {
//  	    for (int c1 = 0; c1 < 2;  ++c1) {
//  	      for (int t1 = 0; t1 < SS; ++t1) {
//  		Sum0[ptr1++] += FF[r-1][ptr];
//  		Sum1[ptr0] += FF[r-1][ptr];
//  		Sum += FF[r-1][ptr++];
//  	      }
//  	    }
//  	  }
//  	  ptr0++;
//  	}
//        }
//      }
    
    vector<double>:: iterator ptrnew = FF[r].begin();
    vector<double>:: iterator ptrold = FF[r-1].begin();

    whichsum = 1-whichsum;

    Sum[whichsum]=0;
    fill(Sum0[whichsum].begin(), Sum0[whichsum].end(), 0);
    fill(Sum1[whichsum].begin(), Sum1[whichsum].end(), 0);
      
    double nchrminus2 = Nind + Nind - 2;
    double nchrminus1 = nchrminus2 + 1;    
    int targ0 = pop[n].get_allele(0, r);
    int targ1 = pop[n].get_allele(1, r);
    ptr0 = 0;

    double exp1 = exp(-totalrho[r-1]*(1.0/nchrminus2));
    double exp2 = exp(-totalrho[r-1]*(1.0/nchrminus1));

    for (int n0 = 0; n0 < Nind;  ++n0) {
      for (int c0 = 0; c0 < 2;  ++c0) {
	int from0 = pop[n0].get_haplotype(c0, r);
        for (int t0 = 0; t0 < SS; ++t0) {
	  ptr1 = 0;
	  for (int n1 = 0; n1 < Nind;  ++n1) {
	    for (int c1 = 0; c1 < 2;  ++c1) {
	      int from1 = pop[n1].get_haplotype(c1, r);
	      for (int t1 = 0; t1 < SS; ++t1) {
		if (  n0 != n  && ( n1 !=n || c1==0 )) {		  
		  // factor 1/Sum in below normalizes the FF to avoid over/underflow
		  double hitprob;
		  if ( n1 == n ) {
		    hitprob = PrHitTarg(r, 2*Nind-2, t0, from0, targ0, Qptr) *
		      PrHitTarg(r, 2*Nind-1, t1, targ0, targ1, Qptr) +
		      PrHitTarg(r, 2*Nind-2, t0, from0, targ1, Qptr) *
		      PrHitTarg(r, 2*Nind-1, t1, targ1, targ0, Qptr);
		  } else {
		    hitprob = 
		  PrHitTarg(r, 2*Nind-2, t0, from0, targ0, Qptr) *
		  PrHitTarg(r, 2*Nind-1, t1, from1, targ1, Qptr) +
		  PrHitTarg(r, 2*Nind-2, t0, from0, targ1, Qptr) *
		  PrHitTarg(r, 2*Nind-1, t1, from1, targ0, Qptr);
		  }
		  
		  *ptrnew = (1/Sum[1-whichsum]) *
		    (*ptrold * exp1 * exp2 
		     + Sum0[1-whichsum][ptr1] * exp2 * (1- exp1) * WEIGHTS[t0]/nchrminus2
		     + Sum1[1-whichsum][ptr0] * exp1 * (1- exp2) * WEIGHTS[t1]/nchrminus1
		     + Sum[1-whichsum] * (1 - exp2) * (1 - exp1) * WEIGHTS[t0] * WEIGHTS[t1]/(nchrminus1 * nchrminus2)) * hitprob ;
		  Sum0[whichsum][ptr1] += *ptrnew; 
		  Sum1[whichsum][ptr0] += *ptrnew;
		  Sum[whichsum] += *ptrnew;	 
		} else {
		  *ptrnew = 0;
		}		
		ptr1++;
		ptrnew++;
		ptrold++;
	      }
	    }
	  }
	  ptr0++;
	}
      }
    }
    logprob += log(Sum[whichsum]);
  }

  //cout << logprob << endl;

 
  //Output FF (for debugging)
  //  ptr = 0;
//    for (int n0 = 0; n0 < Nind;  ++n0) {
//        for (int c0 = 0; c0 < 2;  ++c0) {
//          for (int t0 = 0; t0 < SS; ++t0) {
//  	  for (int n1 = 0; n1 < Nind;  ++n1) {
//  	    for (int c1 = 0; c1 < 2;  ++c1) {
//  	      for (int t1 = 0; t1 < SS; ++t1) {
//  		cout << n0 << "," << c0 << "," << t0 << " ; " << n1 << "," << c1 << "," << t1 << " : " << FF[loci_type.size()-1][ptr++] << endl;
//  	      }
//  	    }
//  	  }
//  	}
//        }
//    }

  return logprob;
}
	  

// simulate the copied individual, chr, and time, using FF
//
void ClassPop::DiploidBackwardsAlg ( int n, double theta, double delta, vector<vector<double> > & FF, vector<double> & rho, double temperature, vector <vector<int> > & CopiedInd, vector<vector<int> > & CopiedChr, vector<vector<int> > & CopiedTime)
  // backwards part
{
  int nchr = Nind + Nind - 2;

  for(int r = loci_type.size()-1; r >= 0; r--) {
    // POTENTIAL CHANGING OF FF by a "temperature", for SA
    //     if(temperature !=1){
    //        double ftemperature;
    //        if(temperature<0.01)
    //  	ftemperature = 0.01;
    //        else
    //  	ftemperature = temperature;
    //        for(int u=0;u<FF[r].size();u++){
    //  	FF[r][u] = exp(log(FF[r][u])/ftemperature);
    //        }
    //      }

    // simulate n0,c0,t0,n1,c1,t1 from FF[r]
    int choice = rint2(FF[r]);
    
    int tcopy1 = choice % SS;
    choice /= SS;
    int ccopy1 = choice % 2;
    choice /= 2;
    int ncopy1 = choice % Nind;
    choice /= Nind;
    int tcopy0 = choice % SS;
    choice /= SS;
    int ccopy0 = choice % 2;
    choice /= 2;
    int ncopy0 = choice % Nind;

    int ptr=0;
    //  int ncopy0,ccopy0,tcopy0,ncopy1,ccopy1,tcopy1;
//      for(ncopy0 = 0; ncopy0 < Nind; ncopy0++) {
//        for(ccopy0 = 0; ccopy0 < 2; ccopy0++) {
//  	for(tcopy0 = 0; tcopy0 < SS; tcopy0++) {
//  	  for(ncopy1 = 0; ncopy1 < Nind; ncopy1++) {
//  	    for(ccopy1 = 0; ccopy1 < 2; ccopy1++) {
//  	      for(tcopy1 = 0; tcopy1 < SS; tcopy1++) {
//  		if(ptr == choice)	goto BACKENDLOOP1;
//  		ptr++;
//  	      }
//  	    }
//  	  }   
//  	}
//        }
//      }
//    BACKENDLOOP1:
    
    CopiedInd[0][r] = ncopy0;
    CopiedChr[0][r] = ccopy0;
    CopiedTime[0][r] = tcopy0;
    CopiedInd[1][r] = ncopy1;
    CopiedChr[1][r] = ccopy1;
    CopiedTime[1][r] = tcopy1;

    int observedallele0 = pop[n].get_orig_allele(0,r);
    int observedallele1 = pop[n].get_orig_allele(1,r);    
    int imputedallele0 = pop[n].get_haplotype(0,r);
    int imputedallele1 = pop[n].get_haplotype(1,r);

    int fromallele0 = pop[ncopy0].get_haplotype(ccopy0,r);
    int fromallele1 = pop[ncopy1].get_haplotype(ccopy1,r);


    // decide whether observedallele0 copied fromallele0, or vice versa
    vector<double> tempprob = vector<double>(2,0.0);

    tempprob[0] = PrHitTarg(r,nchr,tcopy0,fromallele0,imputedallele0,Qptr) * PrHitTarg(r,nchr+1,tcopy1,fromallele1,imputedallele1,Qptr);

    tempprob[1] = PrHitTarg(r,nchr,tcopy0,fromallele0,imputedallele1,Qptr) * PrHitTarg(r,nchr+1,tcopy1,fromallele1,imputedallele0,Qptr);


    //  if(temperature!=1){    
    //double sum = tempprob[0]+tempprob[1];
    //tempprob[0] /=sum;
    //tempprob[1] /=sum;
    //        tempprob[0] = exp(log(tempprob[0])/temperature);
    //        tempprob[1] = exp(log(tempprob[1])/temperature);
//      } 
    
    choice = rint2(tempprob);    
    //cout << choice << endl;
    if(choice==1)
      pop[n].flip_phase(r);

    ptr = 0;
    // modify FF[r-1] to reflect the realised value of ncopy0/1,ccopy0/1,tcopy0/1
    if(r>0){
      double expsave1 = exp ( -rho[r-1]/nchr ); // compute exponential outside loop
    // for efficiency
      double expsave2 = exp (-rho[r-1]/(nchr+1));
      ptr=0;
      double sum=0;
      for (int n0 = 0; n0 < Nind;  ++n0) {
	for (int c0 = 0; c0 < 2;  ++c0) {
	  for(int t0 = 0; t0 <SS; ++t0) {
	    for (int n1 = 0; n1 < Nind;  ++n1) {
	      for (int c1 = 0; c1 < 2;  ++c1) {
		for(int t1 = 0; t1 <SS; ++t1) {
		  if (  n0 != n  && ( n1 !=n || c1==0 )) {
		    FF[r-1][ptr] *= TransitionProb(nchr,n0,c0,t0,ncopy0,ccopy0,tcopy0,rho[r-1],expsave1) * TransitionProb(nchr+1,n1,c1,t1,ncopy1,ccopy1,tcopy1,rho[r-1],expsave2) ;
		  }  
		  sum+=FF[r-1][ptr];
		  ptr ++;
		}
	      }
	    }
	  }
	}
      }

      ptr=0;
      for (int n0 = 0; n0 < Nind;  ++n0) {
	for (int c0 = 0; c0 < 2;  ++c0) {
	  for(int t0 = 0; t0 <SS; ++t0) {
	    for (int n1 = 0; n1 < Nind;  ++n1) {
	      for (int c1 = 0; c1 < 2;  ++c1) {
		for(int t1 = 0; t1 <SS; ++t1) {
		  FF[r-1][ptr]/=sum;  
		  ptr ++;
		}
	      }
	    }
	  }
	}
      }

    }
  }

}
