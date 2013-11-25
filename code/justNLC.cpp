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

//#include "Lanczos_07.h"
//#include "GenHam.h"
//#include "lapack.h"
//#include "simparam.h"
#include "graphs.h"
//#include "entropy.h"
  
int main(int argc, char** argv){

  int CurrentArg = 1;
  string InputFileGraphs;
  string InputFileMeasurements;
  string OutputFile = "output_2d.dat";

  // flags to set the input file (need to do that), output file (not used), and the Renyi S to be measured
  while (CurrentArg < argc)
    {
      if (argv[ CurrentArg ] == string("-g") || argv[ CurrentArg ] == string("--graphs"))
        {
	  InputFileGraphs = string(argv[ CurrentArg + 1 ]);
        }
      if (argv[ CurrentArg ] == string("-m") || argv[ CurrentArg ] == string("--measurements"))
        {
	  InputFileMeasurements = string(argv[ CurrentArg + 1 ]);
        }
      if (argv[ CurrentArg ] == string("-o") || argv[ CurrentArg ] == string("--output"))
        {
	  OutputFile = string(argv[ CurrentArg + 1 ]);
        }
      CurrentArg++;
    }

  //-------------------------------
  // Read Graphs in from File
  //-------------------------------

  vector < graph > fileGraphs;
  ReadGraphsFromFile(fileGraphs, InputFileGraphs);

  //-------------------------------
  // Read Measurement Data in from File
  //-------------------------------

  vector<int> GraphMap; // if they're not in the right order this tells you the order.
  vector<int> Identifier;
  vector<int> Nsites;
  vector<double> Jperp;
  vector<double> Energy;
  vector< vector< pair<double, double> > > RenyiLine, RenyiCorner;
  int Alphas(0);

  ifstream input(InputFileMeasurements.c_str());
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
    cout << "number of entries is " << Alphas <<endl;
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
        // Make GraphMap
        if(Identifier[currentLine]>=GraphMap.size()){ GraphMap.resize(Identifier[currentLine]+1);}
        GraphMap[Identifier[currentLine]] = currentLine;

        //Continue Reading everything in
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
  // Calculate the weights
  //----------------------
    cout << "Weights\n";

    vector<double> WEnergy;
    vector< vector< pair<double, double> > > WLine, WCorner;
    
    WEnergy = Energy;
    WLine = RenyiLine;
    WCorner = RenyiCorner;
    int subGraph;

    //------------ Loop through the graphs and subgraphs
    for (int graph=0; graph<numLines; graph++){
        
        for (int sub = 0; sub < fileGraphs.at(graph).SubgraphList.size(); sub++){
            WEnergy[graph] -= fileGraphs.at(graph).SubgraphList[sub].second * WEnergy[GraphMap[fileGraphs.at(graph).SubgraphList[sub].first]];

            for(int a=0; a<Alphas; a++){
                WLine[graph][a].second -= fileGraphs.at(graph).SubgraphList[sub].second * WLine[GraphMap[fileGraphs.at(graph).SubgraphList[sub].first]][a].second;
                WCorner[graph][a].second -= fileGraphs.at(graph).SubgraphList[sub].second * 
                    WCorner[GraphMap[fileGraphs.at(graph).SubgraphList[sub].first]][a].second;
            }	  
        }      
    }
    
  //----------------------------------------------------------
  // Add all the weights (to get the results for every order)
  //----------------------------------------------------------
    cout << "Sum\n";
    
    //find the maximum order
    int maxSize(0);
    int maxOrder(0);
    for(int i=0; i<numLines; i++){ if(Nsites[i]>maxOrder){maxOrder=Nsites[i];} }
    maxSize = maxOrder-1;

    vector<double> NLCEnergy;
    vector< vector< pair<double, double> > > NLCLine, NLCCorner;

    NLCEnergy.resize(maxSize);
    NLCLine.resize(maxSize);
    NLCCorner.resize(maxSize);

    int entry;
    for (int order=2; order<=maxOrder; order++){
        entry = order-2;
        NLCLine[entry].resize(Alphas);
        NLCCorner[entry].resize(Alphas);

        // change this part to not start at graph 0 //
        for(int graph=0; graph<Energy.size(); graph++){
            if(Nsites[graph]<=order){ 
                NLCEnergy[entry] += WEnergy[graph]*fileGraphs.at(graph).LatticeConstant;
                for(int a=0; a<Alphas; a++){
                    NLCLine[entry][a].second += WLine[graph][a].second * fileGraphs.at(graph).LatticeConstant;
                    NLCCorner[entry][a].second += WCorner[graph][a].second * fileGraphs.at(graph).LatticeConstant;
                }
            }
        }
    }
      
  //----------------------------------------------------------
  // Output results
  //----------------------------------------------------------
    cout << "Output\n";

    //ofstream fout(OutputFile.c_str());
    //fout.precision(10);
    cout.precision(10);

    cout << endl;

    for(int o=2; o<=maxOrder; o++){
      cout <<"Order " <<setw(3)<< o << " Jp " << setw(5) << Jperp[1] << " Energy " << setw(17) << NLCEnergy[o-2]
	   <<" LineEntropies ";

      for(int a=0; a<Alphas; a++){
	cout << setw (5) << RenyiLine[o-2][a].first << setw(17) <<  NLCLine[o-2][a].second ;
      }
      
      cout << " CornerEntropies ";
	
      for(int a=0; a<Alphas; a++){
	cout << setw (5) << RenyiCorner[o-2][a].first << setw(17) <<  NLCCorner[o-2][a].second ;
      }
      cout << endl;
    }
    cout << endl;
    
    
    // fout.close();
    return 0;
}
