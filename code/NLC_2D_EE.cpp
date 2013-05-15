/*******************************************************************
Linked Cluster Expansion program for the 1D TFIM model
Roger Melko, Ann Kallin, Katie Hyatt June 2012
baesd on a Lanczos code from November 2007
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
#include <blitz/array.h>
#include <sstream>

BZ_USING_NAMESPACE(blitz)

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
  bool LF = false;
  double alpha = 0;
  // flags to set the input file (need to do that), output file (not used), and low field
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
      if (argv[ CurrentArg ] == string("-LF") || argv[ CurrentArg ] == string("--lowfield"))
	{
	  LF = true;
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
    double h;

    J=prm.JJ_;
    h=prm.hh_;

    //eigenvector
    Array<l_double,1> eVec;

    //vector of entropies that gets passed to the entropy function
    // contains both corner and line terms
    vector< pair<double, double> > entVec;

    vector< graph > fileGraphs; //graph objects
    
    // List of weights for the different graphs
    vector<double> WeightEnergy, WeightMagnetization;
    vector< vector<double> > WeightLineEntropy, WeightCornerEntropy;

    // Running sum of "the property"
    double RunningSumEnergy(0), RunningSumMagnetization(0);
    vector<double> RunningSumLineEntropy, RunningSumCornerEntropy;

    // the magnetization from the 1D calculation (gets measured in entropy function)
    //double mag;

    // Self explanatory (good naming convention, Katie)
    ReadGraphsFromFile(fileGraphs, InputFile);

    //ofstream fout(OutputFile.c_str());
    //fout.precision(10);
    cout.precision(10);
    
    J=1;     

    const int numhVals = 62;
    //28 values
     double hvals[numhVals] = {0.2,0.3,0.4,0.5,0.6,0.7,0.800001,0.9,1.0,
				1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,
				2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.95,3.0,
				3.044,3.05,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,
				4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0,
				5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10};
     //   double hvals[numhVals] = {3.044};

    // The Renyi entropies to measure (if it's not set in commandline)
    vector <double> alphas;
    if(alpha==0){
      for(double a1=0.5; a1<3.03; a1+=0.25){
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

    // the magnetization file name for each h value & the value in it
    string magFile;
    double magOne(0);

    for(int hh=0; hh<numhVals; hh++){
      entVec.clear();
      entVec.resize(numRenyis);

      h = hvals[hh];

      //Make the magnetization file name for a given h value
      if(LF) magOne=1;
      ostringstream s;
      s<<"./MagFRRRRRiles/mag"<<h<<".input";
      magFile = s.str();
      s.clear();
      
      //If the file exists read it in (only gets used in GenHam in low field case)
      ifstream magIn(magFile.c_str());
      if(magIn){
	magIn >> magOne;
      }
      magOne = abs(magOne);
      // Otherwise magOne stays at 1 or the value used for the last h!
      magIn.close();

      //One Site Graph
      WeightEnergy.push_back(-h); //Energy weight for zero graph (one site)
      WeightMagnetization.push_back(0);  
      //Loop through entropies
      for(int a=0; a<numRenyis; a++){
	WeightLineEntropy[a].push_back(0);
	WeightCornerEntropy[a].push_back(0);
      }//TEST THIS TO MAKE SURE WE'RE NOT GETTING 2 ZEROS OR SOMETHING !!!___!_!_!_!_!_!_!_!

      RunningSumEnergy = WeightEnergy.back();      
      RunningSumMagnetization = WeightMagnetization.back();
      for(int a=0; a<numRenyis; a++){
	RunningSumLineEntropy[a] = WeightLineEntropy[a].back();
	RunningSumCornerEntropy[a] = WeightCornerEntropy[a].back();
      }
      
      //Two Site Graph
      //WeightEnergy.push_back(-sqrt(1+4*h*h)+2*h);
      //WeightEntropy.push_back(TwoSiteEntropy(h,alpha));
      //End of 2 site system stuff!!

      //RunningSumEnergy+=WeightEnergy.back();
      //RunningSumEntropy+=WeightEntropy.back();

      //All the *real* graphs
      for (int i=1; i<fileGraphs.size(); i++){ //skip the zeroth graph
  	
	//---Generate the Hamiltonian--- 
	GENHAM HV(fileGraphs.at(i).NumberSites,J,h,fileGraphs.at(i).AdjacencyList,LF,magOne); 
	

	LANCZOS lancz(HV.Vdim);  //dimension of Hilbert space
	HV.SparseHamJQ();  //generates sparse matrix Hamiltonian for Lanczos

	
	//---Diagonalize and get Eigenvector---
	energy = lancz.Diag(HV, 1, prm.valvec_, eVec); // Hamiltonian, # of eigenvalues to converge, 1 for -values only, 2 for vals AND vectors
	//cout << "Graph " << i <<" energy: " << energy << endl;
	
	//---Energy/Entropy NLC Calculation---
	WeightEnergy.push_back(energy);
	//Entropy1D(alpha, eVec, entVec, mag);
	Entropy2D(alphas, eVec, entVec, fileGraphs.at(i).RealSpaceCoordinates);
	WeightMagnetization.push_back(Magnetization(eVec));

	//Loop Here!!!  ALSO MAKE NOTE THAT LINE IS FIRST AND CORNER IS SECOND !_!_!_!_!_!_!_!_!_!_!_!_!_!
	for(int a=0; a<numRenyis; a++){
	  WeightLineEntropy[a].push_back(entVec[a].first);
	  WeightCornerEntropy[a].push_back(entVec[a].second);
	}
	//cout<<"Entropy "<<i<<" = "<<WeightEntropy.back()<<endl;

	for (int j = 0; j<fileGraphs.at(i).SubgraphList.size(); j++){
	  WeightEnergy.back() -= fileGraphs.at(i).SubgraphList[j].second * WeightEnergy[fileGraphs.at(i).SubgraphList[j].first];
	  WeightMagnetization.back() -= fileGraphs.at(i).SubgraphList[j].second * WeightMagnetization[fileGraphs.at(i).SubgraphList[j].first];
	  for(int a=0; a<numRenyis; a++){
	    WeightLineEntropy[a].back() -= fileGraphs.at(i).SubgraphList[j].second * WeightLineEntropy[a][fileGraphs.at(i).SubgraphList[j].first];
	    WeightCornerEntropy[a].back() -= fileGraphs.at(i).SubgraphList[j].second * WeightCornerEntropy[a][fileGraphs.at(i).SubgraphList[j].first];
	  }	  
	}

	// cout<<"h="<<h<<" J="<<J<<" graph #"<<i<<" energy "<<setprecision(12)<<energy<<endl;
	// cout<<"WeightHigh["<<i<<"] = "<<WeightHigh.back()<<endl;
	RunningSumEnergy += WeightEnergy.back()*fileGraphs.at(i).LatticeConstant;
	RunningSumMagnetization += WeightMagnetization.back()*fileGraphs.at(i).LatticeConstant;
	for(int a=0; a<numRenyis; a++){
	  RunningSumLineEntropy[a] += WeightLineEntropy[a].back()*fileGraphs.at(i).LatticeConstant;
	  RunningSumCornerEntropy[a] += WeightCornerEntropy[a].back()*fileGraphs.at(i).LatticeConstant;
	}
	//	cout <<"S_ " <<alpha <<" h= "<< h<< " RunningSumEnergy " << i <<" "<< RunningSumEnergy << " Entropy= " << RunningSumEntropy 
	//  << " Magnetization= " << RunningSumMagnetization << endl;
      }
      
      //FIND A GOOD WAY TO OUTPUT THE DATA!_!_!_!_!_!_!_!_!_!_!_!_!

      //   cout<<"S_"<<setw(4)<< alpha<<" h= " <<setw(6)<<h<<" Energy= "<<setw(15)<<RunningSumEnergy<<" LineEnt= "<<setw(15)<<RunningSumLineEntropy
      //	  <<" CornerEnt= "<<setw(15)<<RunningSumCornerEntropy<<" Magnetization= "<<setw(15)<<RunningSumMagnetization<<endl;
     
      for(int a=0; a<alphas.size(); a++){ 
	 cout << "h= " << setw(6) << h << " Ener= "<<setw(15)<<RunningSumEnergy<< " Mag= " << setw(15) << RunningSumMagnetization 
	      << "  S_ " << setw (5) << alphas[a] << "  Line= "<< setw(16) << RunningSumLineEntropy[a] 
	      << " Corn=" << setw(17) << RunningSumCornerEntropy[a] << endl;
      }
      cout << endl;

      if(LF){ 
      ofstream magOut(magFile.c_str());
      magOut << abs(RunningSumMagnetization);
      magOut.close();
      }

      WeightEnergy.clear();
      WeightMagnetization.clear();
      RunningSumEnergy=0;
      RunningSumMagnetization=0;

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
