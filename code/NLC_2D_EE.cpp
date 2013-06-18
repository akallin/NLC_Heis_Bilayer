/*******************************************************************
NLCE for the Heisenberg Bilayer!
------------------------------------------------------------------
Linked Cluster Expansion program for the 1D TFIM model
Roger Melko, Ann Kallin, Katie Hyatt June 2012
baesd on a Lanczos code from November 2007, modified 2013
********************************************************************/

#include <utility>
#include <fstream> 
#include <vector> 
#include <math.h>
using namespace std;

#include <iostream>
#include <limits.h>
#include <stdio.h>
#include <time.h>
#include <iomanip>
#include <sstream>

#include "Lanczos_07.h"
#include "GenHam.h"
#include "lapack.h"
#include "simparam.h"
#include "graphs.h"
#include "entropy.h"
  

int main(int argc, char** argv){

  int CurrentArg = 1;
  string InputFile;
  string OutputFile = "output_2d.dat";
  double alpha = 0;

  // flags to set the input file (need to do that), output file (not used), and the Renyi S to be measured
  while (CurrentArg < argc)
    {
      if (argv[ CurrentArg ] == string("-i") || argv[ CurrentArg ] == string("--input"))
        {
	  InputFile = string(argv[ CurrentArg + 1 ]);
        }
      if (argv[ CurrentArg ] == string("-o") || argv[ CurrentArg ] == string("--output"))
        {
	  OutputFile = string(argv[ CurrentArg + 1 ]);
        }
      if (argv [ CurrentArg ] == string("-s") || argv[ CurrentArg ] == string("-S"))
	{
	  alpha = atof( argv [ CurrentArg + 1]);
	}
      CurrentArg++;
    }

    double energy;

    PARAMS prm;  //Read parameters from param.dat  : see simparam.h
    double J;
    double Jperp;

    J=prm.JJ_;
    Jperp=prm.hh_;

    //eigenvector
    vector<l_double> eVec;

    //vector of entropies that gets passed to the entropy function
    // contains both corner and line terms
    vector< pair<double, double> > entVec;

    vector< graph > fileGraphs; //graph objects
    
    // List of weights for the different graphs
    vector<double> WeightEnergy; // WeightMagnetization;
    vector< vector<double> > WeightLineEntropy, WeightCornerEntropy;

    // Running sum of "the property"
    double RunningSumEnergy(0);// RunningSumMagnetization(0);
    vector<double> RunningSumLineEntropy, RunningSumCornerEntropy;

    // Self explanatory (good naming convention, Katie)
    ReadGraphsFromFile(fileGraphs, InputFile);

    //ofstream fout(OutputFile.c_str());
    //fout.precision(10);
    cout.precision(10);
    
    J=1;     

    const int numJperpVals = 1;
    //28 values
    double Jperpvals[numJperpVals] = {1};
    //   double Jperpvals[numhVals] = {3.044};
    
    // The Renyi entropies to measure (if it's not set in commandline)
    vector <double> alphas;
    if(alpha==0){
      for(double a1=1; a1<3.03; a1+=1){
	alphas.push_back(a1);
      }
    }
    else{alphas.push_back(alpha);}
    
    //now that we know the # of renyis, resize entropy vectors
    int numRenyis = alphas.size();
    RunningSumLineEntropy.resize(numRenyis);
    RunningSumCornerEntropy.resize(numRenyis);
    WeightLineEntropy.resize(numRenyis);
    WeightCornerEntropy.resize(numRenyis);
    entVec.resize(numRenyis);
    
    // Loop Over All Values of Jperp -----------------
    for(int Jp=0; Jp<numJperpVals; Jp++){
      entVec.clear();
      entVec.resize(numRenyis);
      
      Jperp = Jperpvals[Jp];
      
      //------------ One Site Graph ----------------
      //  Doesn't exist for Heisenberg!
      // WeightEnergy.push_back(-h); //Energy weight for zero graph (one site)
      // WeightEnergy.push_back(0);
      // WeightMagnetization.push_back(0);  

      //------------ Loop through entropies -------- 
      //for(int a=0; a<numRenyis; a++){
      //	WeightLineEntropy[a].push_back(0);
      //	WeightCornerEntropy[a].push_back(0);
      //}//TEST THIS TO MAKE SURE WE'RE NOT GETTING 2 ZEROS OR SOMETHING !!!___!_!_!_!_!_!_!_!

      //RunningSumEnergy = WeightEnergy.back();      
      //RunningSumMagnetization = WeightMagnetization.back();
      //      for(int a=0; a<numRenyis; a++){
      //	RunningSumLineEntropy[a] = WeightLineEntropy[a].back();
      //	RunningSumCornerEntropy[a] = WeightCornerEntropy[a].back();
      //      }
      
      //------------ Two Site Graph ----------------
      //WeightEnergy.push_back(-sqrt(1+4*h*h)+2*h);
      //WeightEntropy.push_back(TwoSiteEntropy(h,alpha));
      //End of 2 site system stuff!!

      //RunningSumEnergy+=WeightEnergy.back();
      //RunningSumEntropy+=WeightEntropy.back();

      //------------ All the *real* graphs-----------
      for (int i=0; i<fileGraphs.size(); i++){
  	
	if(fileGraphs.at(i).NumberSites<=3){
	  
	  if(fileGraphs.at(i).NumberSites==1){
	    energy=0;
	    for(int a=0; a<numRenyis; a++){
	      entVec[a].first=0;
	      entVec[a].second=0;
	    }
	  }
	  else if(fileGraphs.at(i).NumberSites==2){
	    energy = -0.75*J;
	    for(int a=0; a<numRenyis; a++){
	      entVec[a].first=log(2.0);  //check these values
	      entVec[a].second=0;        //check 'em
	    }
	  }
	  else{
	    energy = -1.0*J;
	    double eig1(5./6.);
	    double eig2(1./6.);
	    for(int a=0; a<numRenyis; a++){
	      if(fabs(1.0-(alphas[a]))<0.000001){entVec[a].first=2.*(-eig1*log(eig1)-eig2*log(eig2));}
	      else{entVec[a].first=-(2./(1.-alphas[a]))*log(pow(eig1,alphas[a])+pow(eig2,alphas[a]));}//check!!!
	      entVec[a].second=0;
	    }
	  }
	}

	else{
	  //---Generate the Hamiltonian---
	  GENHAM HV(fileGraphs.at(i).NumberSites,J,Jperp,fileGraphs.at(i).AdjacencyList); 
	
	  LANCZOS lancz(HV.Vdim);  //dimension of Hilbert space
	  HV.SparseHamJQ();  //generates sparse matrix Hamiltonian for Lanczos

	  //---------- Diagonalize and get Eigenvector -------
	  energy = lancz.Diag(HV, 1, prm.valvec_, eVec); // Hamiltonian, # of eigenvalues to converge, 1 for -values only, 2 for vals AND vectors
	  //cout << "\n ENERGY = " << energy << endl;
	  //cout << "Graph " << i <<" energy: " << energy << endl;
	  Entropy2D(alphas, eVec, entVec, fileGraphs.at(i).RealSpaceCoordinates, HV.Basis);

	}

	//---------- Energy/Entropy NLC Calculation --------
	WeightEnergy.push_back(energy);
	//Entropy1D(alpha, eVec, entVec, mag);
	
	//Loop Here!!!  ALSO MAKE NOTE THAT LINE IS FIRST AND CORNER IS SECOND !_!_!_!_!_!_!_!_!_!_!_!_!_!
	for(int a=0; a<numRenyis; a++){
	  WeightLineEntropy[a].push_back(entVec[a].first);
	  WeightCornerEntropy[a].push_back(entVec[a].second);
	}
	
	//cout<<"Energy " <<i<<" = "<<WeightEnergy.back()<<endl;
	//cout<<"Entropy "<<i<<" = "<<WeightEntropy.back()<<endl;

	for (int j = 0; j<fileGraphs.at(i).SubgraphList.size(); j++){
	  WeightEnergy.back() -= fileGraphs.at(i).SubgraphList[j].second * WeightEnergy[fileGraphs.at(i).SubgraphList[j].first];
	  //  WeightMagnetization.back() -= fileGraphs.at(i).SubgraphList[j].second * WeightMagnetization[fileGraphs.at(i).SubgraphList[j].first];
	  for(int a=0; a<numRenyis; a++){
	    WeightLineEntropy[a].back() -= fileGraphs.at(i).SubgraphList[j].second * WeightLineEntropy[a][fileGraphs.at(i).SubgraphList[j].first];
	    WeightCornerEntropy[a].back() -= fileGraphs.at(i).SubgraphList[j].second * WeightCornerEntropy[a][fileGraphs.at(i).SubgraphList[j].first];
	  }	  
	}
	
	RunningSumEnergy += WeightEnergy.back()*fileGraphs.at(i).LatticeConstant;
	//	RunningSumMagnetization += WeightMagnetization.back()*fileGraphs.at(i).LatticeConstant;
	
	for(int a=0; a<numRenyis; a++){
	  RunningSumLineEntropy[a] += WeightLineEntropy[a].back()*fileGraphs.at(i).LatticeConstant;
	  RunningSumCornerEntropy[a] += WeightCornerEntropy[a].back()*fileGraphs.at(i).LatticeConstant;
	}
	
	if(fileGraphs.size()-1 > i){ if(fileGraphs.at(i).NumberSites != fileGraphs.at(i+1).NumberSites){ 
	    cout <<"Order " <<setw(3)<< fileGraphs.at(i).NumberSites << "    RunningSumEnergy="
		 <<setw(15)<< RunningSumEnergy << "    LineEnt_1= " << setw(15) << RunningSumLineEntropy[0] <<  "    LineEnt_2= " << setw(15) << RunningSumLineEntropy[1] <<  "    LineEnt_3= " << setw(15) << RunningSumLineEntropy[2] << endl;}
	}
	//  << " Magnetization= " << RunningSumMagnetization << endl;
      }
      
      //FIND A GOOD WAY TO OUTPUT THE DATA!_!_!_!_!_!_!_!_!_!_!_!_!

      //   cout<<"S_"<<setw(4)<< alpha<<" h= " <<setw(6)<<h<<" Energy= "<<setw(15)<<RunningSumEnergy<<" LineEnt= "<<setw(15)<<RunningSumLineEntropy
      //	  <<" CornerEnt= "<<setw(15)<<RunningSumCornerEntropy<<" Magnetization= "<<setw(15)<<RunningSumMagnetization<<endl;
     
      for(int a=0; a<alphas.size(); a++){ 
	//	 cout << "h= " << setw(6) << h << " Ener= "<<setw(15)<<RunningSumEnergy<< " Mag= " << setw(15) << RunningSumMagnetization  
	   cout << "Jp= " << setw(6) << Jperp << " Ener= "<<setw(15)<<RunningSumEnergy<< "   " 
	      << "  S_ " << setw (5) << alphas[a] << "  Line= "<< setw(16) << RunningSumLineEntropy[a] 
	      << " Corn=" << setw(17) << RunningSumCornerEntropy[a] << endl;
      }
      cout << endl;

      WeightEnergy.clear();
      // WeightMagnetization.clear();
      RunningSumEnergy=0;
      // RunningSumMagnetization=0;

      for(int a=0; a<numRenyis; a++){
	WeightLineEntropy[a].clear();
	WeightCornerEntropy[a].clear();
	RunningSumLineEntropy[a]=0;
	RunningSumCornerEntropy[a]=0;
      }
    }

    // fout.close();
    return 0;
}
