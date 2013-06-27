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
  
void ReadMeasurementFile( vector< double > & Measurements, const string & file);

int main(int argc, char** argv){

  int CurrentArg = 1;
  string InputFile;
  string OutputFile = "output_2d.dat";

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
      CurrentArg++;
    }


  //-------------------------------
  // Read Data in from File
  //-------------------------------

  vector<int> Identifier;
  vector<int> Nsites;
  vector<double> Jperp;
  vector<double> Energy;
  vector< vector< pair<double, double> > > RenyiLine, RenyiCorner;
  int Alphas(0);

  ifstream input(InputFile.c_str());
  vector< string > rawLines;
  int currentLine;
  int numLines;
  
  const int memberCount = 6;

  // get the number of lines in the file.
  while ( !input.eof() )
    {
        rawLines.resize(rawLines.size() + 1);
        getline(input, rawLines.back()) ; 
    }

    input.close();
    numLines = rawLines.size()-1;
    
    cout<<"number of lines is " <<numLines<<endl;
    
    Identifier.resize(numLines);
    Nsites.resize(numLines);
    Jperp.resize(numLines);
    Energy.resize(numLines);
    RenyiLine.resize(numLines);
    RenyiCorner.resize(numLines);

    stringstream ss (stringstream::in | stringstream::out);
    

    // Check the number of alphas
    ss << rawLines.at(0);
    string tempstring;
    while(!ss.eof()){
      ss >> tempstring;
      Alphas++;
    }
    cout << "number of graphs is " << Alphas <<endl;
    Alphas = (Alphas - 10)/4;
    cout << "number of Renyis is " << Alphas << endl;

    ss.str("");
    ss.clear();
    tempstring="";


    // Read in everything!

    for  (unsigned int currentLine = 0; currentLine < numLines; currentLine++)
      {
	
	ss << rawLines.at(currentLine);
	
	//------ Read in the first line of the graph ------//
	ss >> tempstring;
        ss >> Identifier[currentLine];
	ss >> tempstring;
        ss >> Nsites[currentLine];
	ss >> tempstring;
	ss >> Jperp[currentLine];
	ss >> tempstring;
        ss >> Energy[currentLine];

	RenyiLine[currentLine].resize(Alphas);
	RenyiCorner[currentLine].resize(Alphas);

	ss >> tempstring;
	for(int i=0; i<Alphas; i++){
	  ss >> RenyiLine[currentLine][i].first;
	  ss >> RenyiLine[currentLine][i].second;	  
	}
	ss >> tempstring;
	for(int i=0; i<Alphas; i++){
	  ss >> RenyiCorner[currentLine][i].first;
	  ss >> RenyiCorner[currentLine][i].second;	  
	}
	

        ss.str("");
        ss.clear();
         
    }   

    // print out results
    /*
      for(int i=0; i<numLines; i++){
      cout << "Identifier " << Identifier[i] << ", # sites " << Nsites[i] << ", Jperp " << Jperp[i] 
      << ", Energy " << Energy[i] << " Line ents ";
      for(int j=0; j<Alphas; j++){ cout << RenyiLine[i][j].first << " " << RenyiLine[i][j].second <<", ";}
      cout << "Corner ents ";
      for(int j=0; j<Alphas; j++){ cout << RenyiCorner[i][j].first << " " << RenyiCorner[i][j].second <<", ";}
      cout << endl;
    }
    */


  //----------------------


    // List of weights for the different graphs
    vector<double> WeightEnergy;
    vector< vector<double> > WeightLineEntropy, WeightCornerEntropy;

    // Running sum of "the property"
    double RunningSumEnergy(0);
    vector<double> RunningSumLineEntropy, RunningSumCornerEntropy;

  

    //ofstream fout(OutputFile.c_str());
    //fout.precision(10);
    cout.precision(10);
    
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
      

      //------------ All the *real* graphs-----------
      for (int i=0; i<fileGraphs.size(); i++){
  	

	//---------- Energy/Entropy NLC Calculation --------
	WeightEnergy.push_back(energy);
	
	//Loop Here!!!  ALSO MAKE NOTE THAT LINE IS FIRST AND CORNER IS SECOND !_!_!_!_!_!_!_!_!_!_!_!_!_!
	for(int a=0; a<numRenyis; a++){
	  WeightLineEntropy[a].push_back(entVec[a].first);
	  WeightCornerEntropy[a].push_back(entVec[a].second);
	}
	


	for (int j = 0; j<fileGraphs.at(i).SubgraphList.size(); j++){
	  WeightEnergy.back() -= fileGraphs.at(i).SubgraphList[j].second * WeightEnergy[fileGraphs.at(i).SubgraphList[j].first];

	  for(int a=0; a<numRenyis; a++){
	    WeightLineEntropy[a].back() -= fileGraphs.at(i).SubgraphList[j].second * WeightLineEntropy[a][fileGraphs.at(i).SubgraphList[j].first];
	    WeightCornerEntropy[a].back() -= fileGraphs.at(i).SubgraphList[j].second * WeightCornerEntropy[a][fileGraphs.at(i).SubgraphList[j].first];
	  }	  
	}
	
	RunningSumEnergy += WeightEnergy.back()*fileGraphs.at(i).LatticeConstant;
	
	for(int a=0; a<numRenyis; a++){
	  RunningSumLineEntropy[a] += WeightLineEntropy[a].back()*fileGraphs.at(i).LatticeConstant;
	  RunningSumCornerEntropy[a] += WeightCornerEntropy[a].back()*fileGraphs.at(i).LatticeConstant;
	}
	
	if(fileGraphs.size()-1 > i){ if(fileGraphs.at(i).NumberSites != fileGraphs.at(i+1).NumberSites){ 
	    cout <<"Order " <<setw(3)<< fileGraphs.at(i).NumberSites << "    RunningSumEnergy="
		 <<setw(15)<< RunningSumEnergy << "    LineEnt_1= " << setw(15) << RunningSumLineEntropy[0] 
		 <<  "    LineEnt_2= " << setw(15) << RunningSumLineEntropy[1] <<  "    LineEnt_3= " << setw(15) 
		 << RunningSumLineEntropy[2] << endl;
	  }
	}
      }
      

      //FIND A GOOD WAY TO OUTPUT THE DATA!_!_!_!_!_!_!_!_!_!_!_!_!     
      for(int a=0; a<alphas.size(); a++){ 
	   cout << "Jp= " << setw(6) << Jperp << " Ener= "<<setw(15)<<RunningSumEnergy<< "   " 
	      << "  S_ " << setw (5) << alphas[a] << "  Line= "<< setw(16) << RunningSumLineEntropy[a] 
	      << " Corn=" << setw(17) << RunningSumCornerEntropy[a] << endl;
      }
      cout << endl;

      WeightEnergy.clear();
      RunningSumEnergy=0;

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
