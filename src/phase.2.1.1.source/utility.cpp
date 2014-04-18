
/*
 * $Id: utility.cpp,v 1.13 2003/06/14 00:24:05 stephens Exp $ 
 */

#include "utility.hpp"
#include "constants.hpp"

#include <numeric>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;


int proc_args ( int argc, char ** argv,
                map<string, string> & filenames,
                map<string, int>    & cmdoptions,
		map<string, double>    & d_cmdoptions,
                int    & Niter,
                int    & Nthin,
                int    & Nburn){
  // Filenames    
  //  filenames["clark"]   = "temp.clark";
//    filenames["check"]   = "temp.check";
//    filenames["truth"]   = "temp.truth";
  //filenames["start"]   = "temp.start";
  filenames["endfile"]     = "temp.end";
  filenames["monitorprobs"] = "temp.monitor";
  filenames["deltaFile"] = "eg.delta";
  filenames["pedigreefile"] = "eg.pedigree";
  
  //    
//    cmdoptions["check"]          = 0;
//    cmdoptions["outputforclark"] = 0;
//    cmdoptions["outputtruth"]    = 0;
  cmdoptions["format"]         = 0;
  cmdoptions["usefuzzy"] = 0;
  cmdoptions["usesimplehotspot"] = false;
  cmdoptions["NindToUseForRho"] = 100;
  cmdoptions["useped"] = false;
  cmdoptions["RecomModel"] = 0;  // 0 for general; 1 for simple hotspot; 2 for simple hotspot, fixed pos, 3 for constant (est), 4 for constant (fixed)
  cmdoptions["nhotspot"] = 0;

  //cmdoptions["truthcomp"]      = 0;
  // Which algorithm to use, 
  cmdoptions["EM"] = 0; // set to 1, using -E, to use EM-type method
  cmdoptions["listresolve"] = 0; // set to 1, using -L, to use list method
  
  cmdoptions["method"] = 'R'; // set using -M: 'E' for EM, 'S', 'R', 'C'
  cmdoptions["mcmclistresolve"] = 1;
  cmdoptions["inferrho"] = 0; // this option (-h for "hotspotting") is for
  // use when haplotypes are assumed as input  
  cmdoptions["outputsample"] = 0; // this option (-s) outputs a sample from the posterior
  cmdoptions["hierarchical"] = 0;
  cmdoptions["blocks"] = 0;
  cmdoptions["randomise"] = 0; // set to 1 using -LEr -LSr etc
  cmdoptions["casecontrol"] = 0; // set to 1 using -c
  cmdoptions["NaiveGibbs"]          = 0;    // set to 1 to use Naive Gibbs alg
  //cmdoptions["SD"]             = 1;
  //cmdoptions["usehudson"]      = 0;
  cmdoptions["knowninfo"]    = 0; // default is no known info
  cmdoptions["initmethod"]        =0; // default is to initialise at random
  cmdoptions["Nskip"]           = 0;
  //cmdoptions["useMS"]          = 0;
  cmdoptions["idpresent"]      = 1;
  //cmdoptions["Nset"]           = -1;
  cmdoptions["NDatasets"]      = 1;
  cmdoptions["verbose"]      = 0;
  cmdoptions["inputdelta"]      = 0;
  cmdoptions["testing"]      = 0;
  cmdoptions["AncUpdate"]      = 0;
  cmdoptions["GibbsIter"] = 1000;
  cmdoptions["nperm"] = 0;
  cmdoptions["fRecom"] = 0;
  cmdoptions["NSegments"] = 1;
  cmdoptions["Nrep"] = 1; // number of independent runs to do (set with -x)
  cmdoptions["finalrepmult"] = 1;
  cmdoptions["maxnloci"] = 8;

  d_cmdoptions["phase_threshold"] = 0.9;
  d_cmdoptions["allele_threshold"] = 0.9;
  d_cmdoptions["rhostart"] = 0.0004; //default value, if used
  d_cmdoptions["theta"] = 0.0; // default value: means Wattersons used
  d_cmdoptions["minfreq"] = 0.01;
  d_cmdoptions["output_minfreq"] = 0.01;
  d_cmdoptions["betastart"] = 1.0;
  d_cmdoptions["betaend"] = 0.0;
  d_cmdoptions["left"] = 0; // left and right end of hotspot, if used
  d_cmdoptions["right"] = 0;

  d_cmdoptions["left2"] = 0;
  d_cmdoptions["right2"] = 0;
  d_cmdoptions["maxlambda"] = log(1000.0);
  d_cmdoptions["Hotspotsd"] = 2000; // standard deviation of hotspot width
  d_cmdoptions["MinHotspotSize"] = 200; // minimum size of hotspot
  d_cmdoptions["Hotspotrate"] = 1.0/50000; // one hotspot per 50kb.
  d_cmdoptions["meanRhoMean"] = log(0.0004); // average value of log(rho)
  d_cmdoptions["sdRhoMean"] = 1e100; // expect rho a priori to be within factor 100 of mean
  
  int nperm;
  ifstream prior;

  // Seed for random number generator
  long rngseed = 23423;
  while( ( argc > 1 ) && ( argv[1][0] == '-' ) ) {
    switch(argv[1][1]) {

    case 'a':
	cmdoptions["blocks"] = 1;
	break;

case 'A': 
      cmdoptions["AncUpdate"] = 1;
      break;

    case 'b':
      d_cmdoptions["betastart"] = atof( &argv[1][2] );
      break;

    case 'B':
      d_cmdoptions["betaend"] = atof( &argv[1][2] );
      break;

    case 'c':       // flag to say input file contains case-control info 
      // (before individual id). followed by number of permutations to be done
      // if 0 permutations, use -1.
      cmdoptions["casecontrol"] = 1;
      nperm = atoi ( &argv[1][2] );
      if(nperm > 0 ){
	cmdoptions["nperm"] = nperm;
      }
      if(nperm == 0)
	cmdoptions["nperm"] = 100; // default number of permutations

      if(nperm < 0){
	cmdoptions["nperm"] = 0;
      }
      break;
      
//      case 'C':          /* -Cfilename to take input from hudson's
//  			  program (need -h option set) and output
//  			  it to filename for input to clark's
//  			  program */
      
//        cmdoptions["outputforclark"] = 1;
//        filenames["clark"] = argv[1] + 2;
//        break;
      
    case 'd':
      filenames["deltaFile"] = argv[1] + 2;
      cmdoptions["inputdelta"] = 1; 
      break;

    case 'D':
      cmdoptions["NDatasets"]= atoi(&argv[1][2]);
      if ( cmdoptions["NDatasets"] < 1) {
	cerr << "Error: Number of Datasets should be >0";
      }
      break;

    case 'E':
      cmdoptions["EM"] = 1;
      break;
   
    case 'F': // sets the minimum freq for hap to be included in list under PL
      d_cmdoptions["minfreq"] = atof(&argv[1][2]);
      break;
    
    case 'f':          //f=0 is default format
      // f=1 gives genotypes on single line
      // f=2 gives genotypes on single line, with 'H' to represent hets,
      // and single character for each homozygote
      cmdoptions["format"] = atoi(&argv[1][2]);
      break;
      
    case 'g':          // set number of iterations for naive Gibbs
      // (when naive gibbs used to find starting point)
      cmdoptions["GibbsIter"] = atoi( &argv[1][2] );
      break;
      
    case 'G': // use *only* naive Gibbs - not just to find starting point
      cmdoptions["NaiveGibbs"] = 1;
      break;
     

    case 'h': // infer rho
      cmdoptions["inferrho"] = 1;
      break;

    case 'H': // use hierarchical ligation
      cmdoptions["hierarchical"] = 1;
      break;

   //   case 'H':          // use simulated data from Hudson's program
//        cmdoptions["usehudson"] = 1;
//        cmdoptions["NDatasets"] = atoi( &argv[1][2] );
//        if ( cmdoptions["NDatasets"] < 1 ) { 
//  	cerr << "Error: specify number of data sets to use in"
//  	     << " -H option\n";
//        }
//        break;
      
    //  case 'i':          /* -i0 to initialise to random, -i1 to EM,
//                            -i2 to Gibbs2 */
//        cmdoptions["startmethod"] = atoi( &argv[1][2] );
//        break;
    
    case 'i':   // -i0 to initialise all phases randomly (default), 
      // -i1initfile to initialise all phases to phases in initfile
      // -i2 to initialise all phases to 0
      cmdoptions["initmethod"] = atoi( &argv[1][2] );
      if ( cmdoptions["initmethod"] == 1 )  {
	filenames["initfile"] = argv[1] + 3;
	cerr << "Using initial phase information in " 

	     << filenames["initfile"] << endl;
      }  
      break;
      

    case 'j':          /* Jump first few datasets */
      cmdoptions["Nskip"] = atoi( &argv[1][2] );
      break;

    case 'k':          /* Specify known phases */
      cmdoptions["knowninfo"] = 1;      
      filenames["knownfile"] = argv[1] + 2;
      cerr << "Using known phase information in " 
	   << filenames["knownfile"] << endl;
      if(filenames["knownfile"] == "999"){
	cmdoptions["knowninfo"] = 999;
	cmdoptions["maxnloci"] = 0;
      } 
      break;      

    case 'l': // specify maximum number of loci in each block
      // (specify to be 0 if no maximum)
      cmdoptions["maxnloci"] = atoi( &argv[1][2] );
      break;

    case 'L':
      cmdoptions["listresolve"] = 1;
      cmdoptions["method"] = (int) (argv[1][2]);
      cmdoptions["randomise"] = (int) (argv[1][3]);
      if(cmdoptions["randomise"]=='r')
	cmdoptions["randomise"] = 1;
      else
	cmdoptions["randomise"] = 0;
      break;
      
    case 'M': // used to define method
      // -MC = "Classic" method
      // -MR = modified method with recombination
      // -MR1 n = n simple hotspots (variable pos hotspot)

      // -MR2 n left right left2 right2= simple hotspot model,fixed hotspot pos
      // -MR3 = constant recom rate 
      // -MR4 = constant recom rate at user-specified value

      // -MS = modified method without recombination
      // -ME = modified method with Dirichlet prior
      // -MQ = "Quick" method, which uses S except on final run, when uses R

      cmdoptions["method"] = (int) (argv[1][2]);
      if(cmdoptions["method"]!= 'R' && cmdoptions["method"]!= 'Q' )
	d_cmdoptions["rhostart"] = 0;
      if(cmdoptions["method"] == 'C')
	cmdoptions["mcmclistresolve"] = 0;
      if(cmdoptions["method"]=='R' || cmdoptions["method"] == 'Q'){
	cmdoptions["RecomModel"] = atoi (&argv[1][3]);
	if(cmdoptions["RecomModel"] == 1 || cmdoptions["RecomModel"] == 2){
	  argv++; argc--;
	  cmdoptions["nhotspot"] = atoi (&argv[1][0] );
	  if((cmdoptions["RecomModel"] == 2) && (cmdoptions["nhotspot"] != 1) && (cmdoptions["nhotspot"] != 2) ){
	    cout << "ERROR: number of hotspots for -MR2 must be 1 or 2" << endl;
	    exit(1);
	  }
	  if((cmdoptions["RecomModel"] == 1) && (cmdoptions["nhotspot"] != 1) ){
	    cout << "ERROR: number of hotspots for -MR1 must be 1" << endl;
	    exit(1);
	  }
	} 
	if(cmdoptions["RecomModel"] == 2){
	  argv++; argc--;
	  d_cmdoptions["left"] = atof ( &argv[1][0] );
	  cout << "Left = " << d_cmdoptions["left"] << endl;
	  argv++; argc--;
	  d_cmdoptions["right"] = atof (&argv[1][0] ); 
	  cout << "Right = " << d_cmdoptions["right"] << endl;
	  if(cmdoptions["nhotspot"]==2){
	    argv++; argc--;
	    d_cmdoptions["left2"] = atof ( &argv[1][0] );
	    cout << "Left2 = " << d_cmdoptions["left2"] << endl;
	    argv++; argc--;
	    d_cmdoptions["right2"] = atof (&argv[1][0] ); 
	    cout << "Right2 = " << d_cmdoptions["right2"] << endl;
	  }
	}

      }
	
    //   cmdoptions["randomise"] = (int) (argv[1][3]);
//       if(cmdoptions["randomise"]=='r')
// 	cmdoptions["randomise"] = 1;
//       else
// 	cmdoptions["randomise"] = 0;
      break;
      

   
    case 'm':          // specify file to save monitor probs to
      filenames["monitorprobs"] = argv[1]+2;
      break;
      
    case 'n':          // no ID numbers 
      cmdoptions["idpresent"] = 0;
      break;
      
    case 'N':          // max number to use for estimating rho (default 100)
      cmdoptions["NindToUseForRho"] = atoi( &argv[1][2] );
      break;
    
    case 'O':
      d_cmdoptions["output_minfreq"] = atof(& argv[1][2] );
      break;
   
    case 'p':        
      d_cmdoptions["phase_threshold"] = atof( &argv[1][2] );
      break;
      
    //  case 'r':          /* Randomise input haplotypes (used only
//  			  for testing) */
//        cmdoptions["randomise"] = 1;
//        break;
    
    case 'P':
      cmdoptions["useped"] = atoi(&argv[1][2]);
      filenames["pedigreefile"] =  argv[1]+3;
      break;

      
    case 'q':        
      d_cmdoptions["allele_threshold"] = atof( &argv[1][2] );
      break;
     
    case 'Q':
      cmdoptions["usefuzzy"] = 1;
      break;

    case 'r':
      filenames["priorfile"] = argv[1]+2;
      cerr << "Using priors specified in "  << filenames["priorfile"] << endl;
      prior.open(filenames["priorfile"].c_str());
      prior >> d_cmdoptions["meanRhoMean"];
      prior >> d_cmdoptions["sdRhoMean"];
      prior >> d_cmdoptions["maxlambda"];
      prior >> d_cmdoptions["Hotspotsd"];
      prior >> d_cmdoptions["MinHotspotSize"]; 
      prior >> d_cmdoptions["Hotspotrate"];
      cerr << "Prior mean for rho = " << d_cmdoptions["meanRhoMean"] << endl;
      cerr << "Expect rho to be within factor " << d_cmdoptions["sdRhoMean"] << endl;
      cerr << "Maximum value of lambda = " << d_cmdoptions["maxlambda"] << endl;
      cerr << "Standard deviation of hotspot width prior = " << d_cmdoptions["Hotspotsd"] << endl;
      cerr << "Minimum hotspot width = " << d_cmdoptions["MinHotspotSize"] << endl;
      cerr << "Average number of basepairs per hotspot = " << d_cmdoptions["Hotspotrate"] << endl;

      d_cmdoptions["rhostart"] = d_cmdoptions["meanRhoMean"];

      d_cmdoptions["maxlambda"] = log(d_cmdoptions["maxlambda"]);
      d_cmdoptions["Hotspotrate"] = 1.0/d_cmdoptions["Hotspotrate"];
      d_cmdoptions["meanRhoMean"] = log(d_cmdoptions["meanRhoMean"]);
      d_cmdoptions["sdRhoMean"] = 0.5 * log(d_cmdoptions["sdRhoMean"]);
      
      break;

    case 'R':
      cmdoptions["fRecom"] = 1;
      d_cmdoptions["rhostart"] = atof ( &argv[1][2] );
      break;

    case 's':          // 
      cerr << "Warning: the -s option can produce very large files" << endl;
      cmdoptions["outputsample"] = 1;
      break;
      
    case 'S':          // set SEED
      rngseed = atoi( argv[1] + 2 );          
      break;
      

    case 't': // set theta
      d_cmdoptions["theta"] = atof(&argv[1][2]);
      break;
      
    case 'T':           /* to set to "testing" output */
      cmdoptions["testing"] = 1;      
      break;
      
    case 'v': // set to be verbose  - updates you on allocating memory etc
      cmdoptions["verbose"] = 1;
      break;

    case 'y': // simple hotspot model
      cmdoptions["usesimplehotspot"] = true;
      break;

    case 'x': // set to run several indep runs
      cmdoptions["Nrep"] = atoi(argv[1]+2);
      break;

    case 'X': // set final iterations to be X times as many
      cmdoptions["finalrepmult"] = atoi(argv[1]+2);
      break;

    case 'z':           // specify file name for end phase
      cerr << "Error: the -z option is not implemented in PHASE v 2.0" << endl;
      exit(1);
      //filenames["endfile"] = argv[1] + 2;
      break;
      
    default: 
      cerr << "Error: option " << argv[1] << "unrecognized" << endl;
      exit(1);
    }
    ++argv;
    --argc;
  }
  if ( cmdoptions["usehudson"] == 1 && cmdoptions["useMS"] == 1 )  {
    cerr << "Error: Can't use -M and -H options at the same time" 
	 << endl;
    exit (1);
  }
  
  if(cmdoptions["EM"] == 1 && cmdoptions["NaiveGibbs"] == 1){
    cerr << "Error: Can't use -E and -G options at the same time" << endl;
    exit (1);
  }

  if ( argc < 3) {
    cerr << "usage is PHASE <filename.inp> <filename.out>"
	 << "<number of iterations> <thinning interval> <burn-in>"
	 << endl;
    exit (1);
  }
  filenames["input"]  = argv[1];
  filenames["output"] = argv[2];
  
  if ( argc > 3 ) 
    Niter = atoi( argv[3] );
  if ( argc > 4 ) 
    Nthin = atoi( argv[4] );
  if ( argc > 5 ) 
    Nburn = atoi( argv[5] );
  else
    Nburn = Niter;
  
  // Set random number generator SEED
  //setall ( rngseed, 43243 );
  
  //srandom ( rngseed );
  init_genrand(rngseed);

  return 0;
}



double ranf(){
  double u = genrand_real2();
  //  cout << "Random: " << u << endl;
  return u; //genrand_real2();

  //  return random()/(RAND_MAX+1.0);
}

// generate a random integer according to a user-defined density
int rint2 ( const vector<double> & prob, double psum )
{
    double csum = prob[0];
    double u = ranf();
    
    if(psum == 0.0){ // return a uniform random number if all zeros
      return (int) floor(prob.size() * u);
    }
    else
      {
	if ( psum > 0.0 ) {
	  u *= psum;
	  for (int i = 0; i < prob.size() - 1; ++i) {
            if ( u < csum ) return i;
            csum += prob[i+1];
	  }
	} else {
	  vector<double> cumprob ( prob.size(), 0.0 );
	  // Calculate cdf
	  std::partial_sum ( prob.begin(), prob.end(), cumprob.begin());
	  u *= cumprob[prob.size()-1];
	  for (int i = 0; i < cumprob.size() - 1; ++i) {
            if ( u < cumprob[i] ) return i;
	  }
	}
	return prob.size() - 1;
      }
   
}

//
// normal random generator
// mean mu and sd sigma
// (Ripley (1987), alg 3.17, P82)
//
double rnorm(double mu,double sigma)
{
	double u=0,v=0,x=0,z=0;

	loopstart:
	u=ranf();
	v=0.8578*(2*ranf()-1);
	x=v/u;
	z=0.25*x*x;
	if(z<(1-u)) goto loopend;
	if(z>(0.259/u+0.35)) goto loopstart;
	if(z>(-log(u))) goto loopstart;
	loopend:
	return mu+x*sigma;
}

double logdnorm(double x, double mu, double sigma)
{
  static double logsqrt2pi = 0.9189385;
  return -(logsqrt2pi +log(sigma) + 0.5*(x-mu)*(x-mu)/(sigma*sigma));
}


//
// gamma random generator
// from Ripley, 1987, P230
//
double rgamma(double n,double lambda)
{

	double x=0.0;
	if(n<1)
	{
		const double E=2.71828182;
		const double b=(n+E)/E;
		double p=0.0;
		one: 
		p=b*ranf();
		if(p>1) goto two;
		x=exp(log(p)/n);
		if(x>-log(ranf())) goto one;
		goto three;
		two: 
		x=-log((b-p)/n);
		if (((n-1)*log(x))<log(ranf())) goto one;
		three:;	
	}
	else if(n==1.0)
//
// exponential random variable, from Ripley, 1987, P230
//	
	{
		double a=0.0;
		double u,u0,ustar;
	ten:
		u=ranf();
		u0=u;
	twenty:
		ustar=ranf();
		if(u<ustar) goto thirty;
		u=ranf();
		if(u<ustar) goto twenty;
		a++;
		goto ten;
	thirty:
		return (a+u0)/lambda;
	}
	else
	{
		double static nprev=0.0;
		double static c1=0.0;
		double static c2=0.0;
		double static c3=0.0;
		double static c4=0.0;
		double static c5=0.0;
		double u1;
		double u2;
		if(n!=nprev)
		{
			c1=n-1.0;
			double aa=1.0/c1;
			c2=aa*(n-1/(6*n));
			c3=2*aa;
			c4=c3+2;
			if(n>2.5) c5=1/sqrt(n);
		}
		four:
		u1=ranf();
		u2=ranf();
		if(n<=2.5) goto five;
		u1=u2+c5*(1-1.86*u1);
		if ((u1<=0) || (u1>=1)) goto four;
		five:
		double w=c2*u2/u1;
		if(c3*u1+w+1.0/w < c4) goto six;
		if(c3*log(u1)-log(w)+w >=1) goto four;
		six:
		x=c1*w;		
		nprev=n;
	}	

	return x/lambda;
}

//
// dirichlet random generator
// set b to be ~Dirichlet(a)
//
void rdirichlet(const double * a, const int k, double * b)
{
  int i;
	double sum=0.0;
	for(i=0;i<k;i++)
	{
		b[i]=rgamma(a[i],1);
		sum += b[i];
	}
	for(i=0;i<k;i++)
	{
		b[i] /= sum;
	}
}

//
// Random permutation of 0 to n-1, in perm
//
void rperm(vector<int> & perm,int n)
{
  int i,s,temp,t;
  for(i=0;i<n;i++)
    perm[i]=i;
  for(i=0;i<n;i++)
    {
      t=n-i-1; // t runs from n-1 down to 0
      s=(int) floor((t+1)*ranf()); // s unif on 0 to t
      temp=perm[s]; // swap s and t
      perm[s]=perm[t];
      perm[t]=temp;
    }
}

// {{{ Log
/* 
   $Log: utility.cpp,v $
   Revision 1.13  2003/06/14 00:24:05  stephens
   Adding files, and committing a lot of changes
   that have resulted in version 2.0 of PHASE

   Revision 1.12  2002/02/27 18:56:45  stephens
   Commiting the current source, which is essentially that released as
   PHASE 1.0. The number of quadrature points is reduced to 2. Some code has
   been added to cope with recombination, but this is not included in the release
   and is still under development.

   Revision 1.11  2001/10/21 18:33:29  stephens
   fixed bug in GibbsUpdate

   Revision 1.9  2001/10/16 20:44:17  stephens
   Added way to input position of each marker, by preceding line
   of SSMMMSS with a line beginning with P followed by position of
   each marker.

   Revision 1.8  2001/10/16 18:12:06  stephens
   More minor bugs fixed. (including -m bug)

   Revision 1.7  2001/10/12 23:49:27  stephens
   Various updates, particularly cosmetic, but some bug fixes.
   this version tested on MS and SNP data without missing alleles
   gives very similar answers to the original phase.

   Revision 1.6  2001/10/09 23:03:25  stephens
   Modified to deal with very large files - put in checks for underflow
   and overflow, and changed computation of Q, DiffProb etc to only compute those elements that are necessary.

   Revision 1.5  2001/06/19 17:02:26  stephens
   Changes to computation of arrayFF to make more efficient
   Added facility to store "original phenotype" in indnode,
   in preparation for allowing genotyping error.

   Revision 1.4  2001/05/21 20:17:29  nali
   No real changes

   Revision 1.3  2001/04/24 19:36:37  nali
   rint2 is enhanced so that it is not necessary to provide the sum of the
   (not normalized) probabilities.

   Revision 1.2  2001/04/19 19:50:08  nali
   Minor changes to the other files.

   Revision 1.1  2001/04/17 22:08:44  nali

   Major revisement in overall structure. Created new class ClassPop and
   almost all global functions now became member functions of ClassPop, most
   of them private.

   "mult" removed in update_phase_NR. No other changes in terms of algorithm.
   Haven't check the results yet.

   proc_args() is moved to utility.cpp, which also defines a couple of other
   global functions.

   Revision 1.8  2001/04/16 18:47:06  stephens
   changed output to two files (haplotypes and phases)

   Revision 1.7  2001/04/12 01:53:26  stephens
   Added -D option for reading in several datasets
   Debugged -H option for inputting data from Hudson's simulations
   Reduced Nthin to 1 to reduce time for examples, and
   added a test example for the -H option and -D options

   Revision 1.6  2001/04/07 07:18:10  nali
   Revised to concentrate on implementing SD method.
   Discarded support for EM algorithm, etc..

   Revision 1.5  2001/02/28 04:54:32  nali
   Make EM working and haplotype list right

   Revision 1.4  2001/02/27 08:17:58  nali
   No real changes.

   Revision 1.3  2001/02/16 17:13:17  nali
   New functions: InputRandom, InputHusdonData, etc..
   New member function: make_haplist.

   Revision 1.2  2001/02/12 05:28:42  nali
   Now it can read in and output phenotype data.
   No further processing yet.

   Revision 1.1  2001/02/10 01:13:52  nali
   New files added.

*/
// }}}
