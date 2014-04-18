/*
 * $Id: phase.cpp,v 1.34 2003/06/14 00:24:05 stephens Exp $
 */

#include "classpop.hpp"
#include "constants.hpp"
#include "errcheck.hpp"
#include "utility.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>

using namespace std;
int NHAP = 0; // global variables used for debugging
int NIND = 0;
int NHL = 0;
bool TOOLONG; // used in testing for methods that don't want too many loci

double MCMCResolve(ClassPop & allpop, int Niter, int Nthin, int Nburn, vector<double> & vecDelta,map<string, int> & cmdoptions, map<string, double> & d_cmdoptions, string filename, int & segment, bool collectdata, vector<int> blocks = vector<int>(0) )
{
  double loglik; 
  
  int nloci = allpop.get_nloci();
  int nchr = 2*allpop.get_nind();
  
  vector<double> tempprob(3,0);
  tempprob[0] = 0.4;
  tempprob[1] = 0.3;
  tempprob[2] = 0.3;
  int updateloci = cmdoptions["maxnloci"] ;

  if((nloci <= updateloci) || (cmdoptions["maxnloci"]<=0))
    updateloci = nloci;
  else if(nloci <= 2*cmdoptions["maxnloci"])
    updateloci = nloci / 2;    
  else{
    updateloci = cmdoptions["maxnloci"] - rint2(tempprob,1.0);
    if(updateloci<2)
      updateloci = 2;
  }

  if((nloci>updateloci) && (updateloci>0)){    
    ClassPop LHSpop(allpop,0,updateloci);
    ClassPop RHSpop(allpop,updateloci,nloci);
 
    if(cmdoptions["hierarchical"] == 1){
      // following is for hierarchical ligation, instead of progressive
      LHSpop = ClassPop(allpop,0,nloci/2);
      RHSpop = ClassPop(allpop,nloci/2,nloci);
    }
    allpop = ClassPop(); // remove allpop to save memory

    MCMCResolve(LHSpop, Niter, Nthin, Nburn, vecDelta,cmdoptions,d_cmdoptions,filename,segment, false);
    cerr << segment << " segment operations done" << endl;
    segment ++;
    MCMCResolve(RHSpop, Niter, Nthin, Nburn, vecDelta,cmdoptions,d_cmdoptions,filename,segment, false);    
    cerr << segment << " segment operations done" << endl;
    segment ++;
    
    allpop = ClassPop(LHSpop,RHSpop,(1.0/nchr)*d_cmdoptions["minfreq"]);

    allpop.ResetCounts(); // clear counts taken so far of all the individual phases, for final MCMC run
    
    loglik = allpop.MCMCListResolvePhase(cmdoptions, Niter,Nthin, Nburn, vecDelta,d_cmdoptions, filename, false, collectdata);
    //allpop = ClassPop(LHSpop,RHSpop,cmdoptions["method"], Niter, vecDelta,d_cmdoptions["rhostart"],d_cmdoptions["betastart"], d_cmdoptions["betaend"],(1.0/nchr)*d_cmdoptions["minfreq"]);
    
    //} 
  } else {
    allpop.ResetCounts(); // clear counts taken so far of all the individual phases, for final MCMC run
    bool initialise = true;
    if(cmdoptions["blocks"])
      initialise = false;
    loglik = allpop.MCMCListResolvePhase(cmdoptions, Niter, Nthin, Nburn, vecDelta,d_cmdoptions, filename, initialise, collectdata); 	
  }
  
  return loglik;
}
	  	
int main ( int argc, char** argv)
{

    // Processing command line options
    // Default options     
    map<string, string> filenames;
    map<string, int>    cmdoptions;
    map<string, double> d_cmdoptions; // double command options

    int Niter = 100;           // Number of MCMC iterations 
    int Nthin = 1;             // Thinning interval
    int Nburn = 100;           // Burn-in iterations
    int status = proc_args ( argc, argv, filenames, cmdoptions, d_cmdoptions,
             Niter, Nthin, Nburn);
    
    // Read in data file
    ifstream input (filenames["input"].c_str());
    assure ( input, filenames["input"] );
    // Open outputfile
    ofstream output (filenames["output"].c_str());
    assure ( output, filenames["output"] );
    
    
    string freqfilename = filenames["output"]+"_freqs";
    ofstream freqfile (freqfilename.c_str());

    string monitorfilename = filenames["output"]+"_monitor";
    ofstream monitorfile (monitorfilename.c_str());
    assure ( monitorfile, monitorfilename );

    string hbgfilename = filenames["output"]+"_hbg";
    ofstream hbgfile(hbgfilename.c_str());
    hbgfile.close();
    string probfilename = filenames["output"]+"_probs";
    ofstream probfile(probfilename.c_str());
    probfile.close();

    ifstream ifstr_known;
    ifstream ifstr_init;
    ifstream deltaFile;
    
    for (int dataset = 0; dataset < cmdoptions["Nskip"]; 
         dataset++) {
        // Set loci types
        
        // initialization and read in data
        ClassPop allpop;        

        cerr << "Reading in data" << endl;
        allpop.read_data ( input, cmdoptions, d_cmdoptions, filenames );        
        cerr << endl << "Finished reading" << endl;
	       	
    }
    for (int dataset = 0; dataset < cmdoptions["NDatasets"]; 
         dataset++) {
        // Set loci types
        
        // initialization and read in data
        ClassPop allpop;        

        cerr << "Reading in data" << endl;
        allpop.read_data ( input, cmdoptions, d_cmdoptions, filenames );        
        cerr << endl << "Finished reading" << endl;
       	
        vector<double> vecDelta ( allpop.get_nloci(), 0.0 );

	// input delta from file if necessary
	if(cmdoptions["inputdelta"]==1){
	  if(filenames["deltaFile"] == "1"){
	    for(int r=0; r<vecDelta.size(); r++)
	      vecDelta[r] = 1;
	  } else {
	    deltaFile.open (filenames["deltaFile"].c_str());
	    assure ( deltaFile, filenames["deltaFile"] );
	    for(int r=0; r<vecDelta.size(); r++)
	      deltaFile >> vecDelta[r];
	  }
	}

       	
	if(cmdoptions["knowninfo"] == 1){
	  ifstr_known.open ( filenames["knownfile"].c_str()); 
	  assure (ifstr_known, filenames["knownfile"]);
	}

       	
	if(cmdoptions["initmethod"] == 1){
	  ifstr_init.open ( filenames["initfile"].c_str());
	  assure (ifstr_init, filenames["initfile"]);
	}

	allpop.initialize ( ifstr_known, ifstr_init, cmdoptions["knowninfo"], 
			    cmdoptions["initmethod"], cmdoptions["theta"], vecDelta, cmdoptions["format"]);
	
	if(cmdoptions["knowninfo"] == 1)
	  ifstr_known.close();
	if(cmdoptions["initmethod"] == 1 )
	  ifstr_init.close();

      	
//  if ( cmdoptions["startinfo"] == 2 ||
//               cmdoptions["startinfo"] == 3 ) {
//  	  ifstr_start.open ( filenames["startfile"].c_str());
//  	  assure (ifstr_start, filenames["startfile"]);
//  	  allpop.initialize ( ifstr_start,
//  			      cmdoptions["startinfo"], cmdoptions["format"] );
//  	  ifstr_start.close ();
//          } else {
//  	  // No start file or known file necessary
//  	  allpop.initialize ( ifstr_start,
//  			      cmdoptions["startinfo"], cmdoptions["format"]  );
//          }
	
        allpop.output(cout,true,true);
	int nloci = allpop.get_nloci();
	if(nloci < 25)
	  TOOLONG = false;
	else
	  TOOLONG = true;

	//  ClassPop partpop(allpop,3,10);
//  	partpop.output(cout,true,true);
//  	allpop = partpop;
	
	double bestscore;
	double score;
	ClassPop bestpop;
	for(int count = 0; count < cmdoptions["Nrep"]; count++){
	    
	  // decide which method to use to do the resolving; return
	  // a measure of the goodness of fit in "score"

	  if(cmdoptions["inferrho"]==1){ // just do MCMC on input
	    //haplotypes to infer rho
	    double MeanRhoStart = d_cmdoptions["rhostart"];
	    //allpop.InitialiseRho(MeanRhoStart);
	    double sigmamean = 1;
	    double sigmamult = 0.1;
	    cout << "simple hotspot = " << cmdoptions["usesimplehotspot"] << endl;
	    allpop.InferRho(Niter, sigmamean, sigmamult, cmdoptions["verbose"], d_cmdoptions);
	  } 
	  else if(cmdoptions["useped"]==2){

	  	    allpop.FastHapMapResolve(Niter,Nburn);

	  }
	  // List method
	  else if(cmdoptions["listresolve"]==1){
	    double likelihood;  
	    if(!TOOLONG){
	      ClassPop partpop1(allpop,0,nloci/2);
	      ClassPop partpop2(allpop,nloci/2,nloci);
	      likelihood = partpop1.ListResolvePhase(cmdoptions["method"], Niter, vecDelta,d_cmdoptions["rhostart"],cmdoptions["randomise"],true);
	      //output << "Likelihood = " << likelihood << endl;
	      partpop1.output(cout,true,true);
	      
	      likelihood = partpop2.ListResolvePhase(cmdoptions["method"], Niter, vecDelta,d_cmdoptions["rhostart"],cmdoptions["randomise"],true);
	      //output << "Likelihood = " << likelihood << endl;
	      partpop2.output(cout,true,true);
	      
	      allpop = ClassPop(partpop1,partpop2,0.000001);
	      
	      //likelihood = allpop.ListResolvePhase(cmdoptions["method"], Niter, vecDelta,d_cmdoptions["rhostart"],cmdoptions["randomise"],true);
	      likelihood = allpop.ListResolvePhase(cmdoptions["method"], Niter, vecDelta,d_cmdoptions["rhostart"],cmdoptions["randomise"],false);
	      cout << "Likelihood = " << likelihood << endl;
	      allpop.output(cout, true,true);
	    } 
	  } 
	  else if(cmdoptions["mcmclistresolve"]==1){
	    //allpop.InitialiseRho(d_cmdoptions["rhostart"]);
	    int segment = 0;
	    cout << "Resolving with method " << (char) cmdoptions["method"] << endl;
	    score = MCMCResolve(allpop,Niter,Nthin,Nburn,vecDelta,cmdoptions,d_cmdoptions,filenames["output"],segment,true);
	  }
	  else{
	    // no recombination method
	    if(cmdoptions["fRecom"]==0){
	      
	      cerr << "Finding good starting point..." << endl;
	      allpop.GibbsResolvePhase(cmdoptions["GibbsIter"],DIRPRIOR);
	      
	      cerr << "Starting to resolve phase..." << endl;
	      
	      score = allpop.resolve_phase_NR ( Nburn, Nthin, Niter,
						output, monitorfile, vecDelta, 
						cmdoptions["verbose"], cmdoptions["AncUpdate"], cmdoptions["NaiveGibbs"]);	     
	    }
	    else {
	      allpop.resolve_phase_R ( d_cmdoptions["theta"], 
				       d_cmdoptions["rhostart"],
				       Nburn, Nthin, Niter, 
				       output, monitorfile, vecDelta, cmdoptions["verbose"]);
	    }
	  }
	  
	  //cerr << "SCORE FOR THIS TRIAL =" << score << endl;
	  
	  if ((count ==0) || (score > bestscore)){
	    //cerr << "This beats current score - switching bestpop" << endl;
	    bestpop = allpop;
	    bestscore = score;
	  }
	}
	
	allpop = bestpop;
	
	if(cmdoptions["testing"] ==1){

	  hbgfile.open(hbgfilename.c_str(),ios::app);
	  probfile.open(probfilename.c_str(),ios::app);
	  
	  //allpop.RestoreSavedState();
	  //allpop.output_all_haps (output  ,true , false, false, d_cmdoptions["phase_threshold"]);
	  
	  allpop.output_all_haps (output  ,true , false, true, true, d_cmdoptions["phase_threshold"]);
	  
	  allpop.output_all_correct_probs(probfile);

 	  // print out best guesses according to PL phase 
	  // ie current values of haps.
	  allpop.output_all_haps( hbgfile, true , false, false, true, d_cmdoptions["phase_threshold"]);
	  
	  hbgfile.close();
	  probfile.close();


	}
	else
	  {
	    bool printnames = true;
	    bool printknownphase = true;
	    bool printbestguess = true;
	    bool printmissing = true;

	    allpop.OutputHapList(freqfile, d_cmdoptions["minfreq"]);

	    output << "*************************************************************" << endl;
	    output << "****                                                     ****" << endl;
	    output << "****            Output from PHASE v2.1.1                 ****" << endl;
	    output << "****  Code by M Stephens, with contributions from N Li   ****" << endl;
            output << "****                                                     ****" << endl;
	    output << "*************************************************************" << endl;
	    output << endl;

	    output << "BEGIN COMMAND_LINE " << endl;
	    for(int i = 0; i<argc; i++)
	      output << argv[i] << " ";
	    output << endl;
	    output << "END COMMAND_LINE " << endl;
	    output << endl;

	    output << "BEGIN OUTFILE_LIST" << endl;
	    output << freqfilename << " : haplotype frequency estimates" << endl;
	    output << filenames["output"] + "_pairs : most likely haplotype pairs for each individual" << endl;
	    if(cmdoptions["method"] == 'R' || cmdoptions["method"] == 'Q')
	      output << filenames["output"] + "_recom : estimates of recombination parameters" << endl;
	    if(cmdoptions["casecontrol"])
	      output << filenames["output"] + "_signif : p-value for testing cases vs controls" << endl;
	    if(cmdoptions["outputsample"])
	      output << filenames["output"] + "_sample : sample from posterior distribution of haplotype reconstruction" << endl;
	    
	    output << monitorfilename << " : file for monitoring convergence" << endl;
	    output << "END OUTFILE_LIST" << endl;
	    
	    output << endl;
	     
	    output << "BEGIN INPUT_SUMMARY" << endl;
	    output << "Number of Individuals: " << allpop.get_nind() << endl;
	    output << "Number of Loci: " << allpop.get_nloci() << endl;
	    output << "Positions of loci: ";
	    for(int r=0; r< allpop.get_nloci(); r++)
	      output << allpop.get_position(r) << " ";
	    output << endl;
	    output << "END INPUT_SUMMARY" << endl;
	    
	    output << endl;

	   
	   
	    allpop.MakeHapList(true);
	    	    
	    // output << "List of haplotypes found in best reconstruction," << endl;
// 	    output << "with counts" << endl;
// 	    output << endl;
// 	    allpop.OutputHapList(output,0,false);
// 	    output << endl << endl;
	    
	    output << "List of haplotypes found in best reconstruction, with counts." << endl;
	    output << "(See file " << freqfilename << " for haplotype population frequency estimates)" << endl;
	    output << endl;
	    output << "BEGIN LIST_SUMMARY" << endl;
	    allpop.OutputHapList(output,0,false);
	    output << "END LIST_SUMMARY" << endl;
	    output << endl << endl;


	    output << "Summary of best reconstruction" << endl;
	    output << "(numbers refer to the list of haplotypes given above)" << endl;
	    output << endl;
	    output << "BEGIN BESTPAIRS_SUMMARY" << endl;
	    allpop.OutputHaplistSummary(output);
	    output << "END BESTPAIRS_SUMMARY" << endl;
	    output << endl << endl;
	    
	    output << "Haplotype estimates for each individual, with uncertain phases enclosed in ()" << endl;
	    output << "and uncertain genotypes enclosed in []:" << endl;
	    output << endl;
	    
	    output << "BEGIN BESTPAIRS1" << endl;
	    allpop.output_all_haps (output  ,printknownphase , printnames, printbestguess, printmissing , d_cmdoptions["phase_threshold"], d_cmdoptions["allele_threshold"]);
	    output << "END BESTPAIRS1" << endl;	    
	    
	    output << endl << endl;
	    output << "Haplotype estimates for each individual, with uncertain phases enclosed in ()" << endl;
	    output << "and uncertain genotypes enclosed in []" << endl;
	    output << "with phase known positions indicated by " << SPACEHOLDER << endl;
	    
	    output << endl;
	    
	     output << "BEGIN BESTPAIRS2" << endl;
	    allpop.output_all_haps (output  ,false , printnames, printbestguess, false, d_cmdoptions["phase_threshold"], d_cmdoptions["allele_threshold"]);
	    output << "END BESTPAIRS2" << endl;	 
	    
	    output << endl << endl;
	    
	    output << "Phase probabilities at each site" << endl;
	    output << "with phase known positions indicated by " << SPACEHOLDER << endl;
	    output << "and missing data positions indicated by " << '?' << endl;
	    
	    output << endl;
	    output << "BEGIN PHASEPROBS" << endl;
	    allpop.OutputPhaseProbs( output, false);
	    output << "END PHASEPROBS" << endl;

	    // these used to output final phases, not valid with PHASE v2.0, so removed
	    //ofstream endfile (filenames["endfile"].c_str());
	    //assure( endfile, filenames["endfile"] );
	    //allpop.output_all_phases(endfile, true);
	  }
    }
    

    
    // Close file streams
    input.close();
    output.close();
    monitorfile.close();
    freqfile.close();
       
    return 0;
}

// {{{ Log

/* 
   $Log: phase.cpp,v $
   Revision 1.34  2003/06/14 00:24:05  stephens
   Adding files, and committing a lot of changes
   that have resulted in version 2.0 of PHASE

   Revision 1.33  2002/02/27 18:56:45  stephens
   Commiting the current source, which is essentially that released as
   PHASE 1.0. The number of quadrature points is reduced to 2. Some code has
   been added to cope with recombination, but this is not included in the release
   and is still under development.

   Revision 1.32  2001/10/21 18:33:29  stephens
   fixed bug in GibbsUpdate

   Revision 1.30  2001/10/16 20:44:17  stephens
   Added way to input position of each marker, by preceding line
   of SSMMMSS with a line beginning with P followed by position of
   each marker.

   Revision 1.27  2001/10/09 23:03:25  stephens
   Modified to deal with very large files - put in checks for underflow
   and overflow, and changed computation of Q, DiffProb etc to only compute those elements that are necessary.

   Revision 1.26  2001/10/02 18:03:00  stephens
   just some small additions. haven't committed these for a long time
   because cvs was not working propersly. (now know why: they remounted
   you from /user3 to /home)

   Revision 1.25  2001/06/19 17:02:25  stephens
   Changes to computation of arrayFF to make more efficient
   Added facility to store "original phenotype" in indnode,
   in preparation for allowing genotyping error.

   Revision 1.24  2001/05/30 06:02:18  stephens
   Updated to be considerably more efficient, via introduction of
   various new methods, including introduction of ArrayDiffProb
   and ArrayDiffCount to improve computation al efficiency for SNP
   data. Speedtest.inp now runs in about 7secs.
   Also corrected several bugs. Output now looks more promising
   and major bugs appear to have been eliminated. Convergence of chain
   can now be monitored more easily by the output in temp.monitor, which
   gives the pseudo-likelihood every Nthin repetitions.

   Revision 1.23  2001/05/21 20:17:29  nali
   No real changes

   Revision 1.22  2001/04/24 19:31:31  nali
   Move data input out of the constructor. Member functions read_data and initialize
   have to be called explicitly. Put everything related to hudson data set into a single
   function so that it might go away one day.

   Revision 1.21  2001/04/20 00:32:45  nali
   Put reading loci types into the contructor of ClassPop

   Revision 1.20  2001/04/19 19:50:08  nali
   Minor changes to the other files.

   Revision 1.19  2001/04/17 22:08:44  nali

   Major revisement in overall structure. Created new class ClassPop and
   almost all global functions now became member functions of ClassPop, most
   of them private.

   "mult" removed in update_phase_NR. No other changes in terms of algorithm.
   Haven't check the results yet.

   proc_args() is moved to utility.cpp, which also defines a couple of other
   global functions.

   Revision 1.18  2001/04/16 18:47:06  stephens
   changed output to two files (haplotypes and phases)

   Revision 1.17  2001/04/12 01:53:26  stephens
   Added -D option for reading in several datasets
   Debugged -H option for inputting data from Hudson's simulations
   Reduced Nthin to 1 to reduce time for examples, and
   added a test example for the -H option and -D options

   Revision 1.16  2001/04/11 19:35:46  stephens
   changed Nthin to 10, so that some updates are actually done!
   Results from test sets seem to suggest a bug in the updating.

   Revision 1.15  2001/04/11 19:20:54  stephens
   Added OutputPhase to output the phase of all individuals every iteration.
   Appears to be a bug, since all phases are output as 0 at each step
   in our example.

   Revision 1.14  2001/04/09 16:28:35  nali
   Most part of original ResolvePhase (without recombination) implemented.
   Now compiles and runs.

   Revision 1.13  2001/04/07 07:18:10  nali
   Revised to concentrate on implementing SD method.
   Discarded support for EM algorithm, etc..

   Revision 1.12  2001/04/04 06:26:48  nali
   OutFreq compiled.

   Revision 1.11  2001/02/28 04:54:32  nali
   Make EM working and haplotype list right

   Revision 1.10  2001/02/27 08:16:49  nali
   Make use of the new haplotype list class.

   Revision 1.9  2001/02/21 18:36:23  nali
   indnode.cpp

   Revision 1.8  2001/02/20 16:15:29  nali
   Use assert ( )  for debugging.

   Revision 1.7  2001/02/16 17:13:17  nali
   New functions: InputRandom, InputHusdonData, etc..
   New member function: make_haplist.

   Revision 1.6  2001/02/16 02:41:48  nali

   Remove the character representation of phenotypes.

   Revision 1.5  2001/02/14 23:39:03  nali
   New InputData() function.

   Revision 1.4  2001/02/14 01:55:23  nali

   Move printing phenotypes into a member function call.

   Revision 1.3  2001/02/13 20:40:16  nali

   Two bugs fixed.
   1. Needs a copy constructor if a vector of objects is to be created.
      Since the vector is initialized by first creating a temporary object
      by the default constructor and copying it to every object in the
      vector. Especially when new operator is used in the constructor.
   2. For string stream, char str[10], where the length is necessary.
      Then ostrstream ostr(str, 10), where 10 is the length of char array.

   Revision 1.2  2001/02/12 05:28:42  nali
   Now it can read in and output phenotype data.
   No further processing yet.

   Revision 1.1  2001/02/10 01:13:52  nali
   New files added.

*/

// }}}
