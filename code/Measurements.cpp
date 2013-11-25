/*******************************************************************
  NLCE for the Heisenberg Bilayer! 2013
  -----------------------------------------------------------------
  Linked Cluster Expansion program for the 1D TFIM model
  Roger Melko, Ann Kallin, Katie Hyatt June 2012
  based on a Lanczos code from November 2007, modified 2013
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

unsigned long int swap_layers(unsigned long int , int );

int main(int argc, char** argv){

    cout << currentDateTime() << " Beginning Program " << endl;
    int CurrentArg = 1;
    string InputFile = "3x3square.dat"; 
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

    PARAMS prm;  //Read parameters from simparam.h
    double J;
    double Jperp;

    J=prm.JJ_;
    Jperp=prm.Jperp_;

    //eigenvector
    vector<l_double> eVec;

    //vector of entropies that gets passed to the entropy function
    // contains both corner and line terms
    vector< pair<double, double> > entVec;

    vector< graph > fileGraphs; //graph objects

    // List of weights for the different graphs
    double WeightEnergy; // WeightMagnetization;
    vector<double> WeightLineEntropy, WeightCornerEntropy;

    // Self explanatory (good naming convention, Katie)
    ReadGraphsFromFile(fileGraphs, InputFile);
    cout << "graph file is: " << InputFile << endl;

    //ofstream fout(OutputFile.c_str());
    //fout.precision(10);
    cout.precision(10);

    // The Jperp values to use
    vector <double> Jperpvals;
    Jperpvals.push_back(2.522);

    // The Renyi entropies to measure (if it's not set in commandline)
    vector <double> alphas;
    if(alpha==0){
        for(double a1=0.05; a1<5.05; a1+=0.025){
            alphas.push_back(a1);
        }
    }
    else{alphas.push_back(alpha);}

    //now that we know the # of renyis, resize entropy vectors
    int numRenyis = alphas.size();
    WeightLineEntropy.resize(numRenyis);
    WeightCornerEntropy.resize(numRenyis);
    entVec.resize(numRenyis);

    // Loop Over All Values of Jperp -----------------
    for(int Jp=0; Jp<Jperpvals.size(); Jp++){
        entVec.clear();
        entVec.resize(numRenyis);

        Jperp = Jperpvals[Jp];


        // Loop over all the fileGraphs ------------ 
        for (int i=0; i<fileGraphs.size(); i++){

            // Special Case for 1-site, 2-site, 3-site graphs
            // Change this for the bilayer case

                //Now this is a 2 site graph... will this work?  What if the layers are uncoupled?
                if(fileGraphs.at(i).NumberSites==1){
                    //energy depends on Jperp
                    energy= -0.75*Jperp;
                    for(int a=0; a<numRenyis; a++){
                        entVec[a].first=0;
                        entVec[a].second=0;
                    }
                }
            
            // All other graphs
                else{
                //---Generate the Hamiltonian---
                GENHAM HV(fileGraphs.at(i).NumberSites,fileGraphs.at(i).AdjacencyList); 

                // Doesn't do anything really
                LANCZOS lancz(HV.Vdim,J,Jperp,fileGraphs.at(i).NumberSites);  //dimension of Hilbert space 
                //HV.SparseHamJQ();  //generates sparse matrix Hamiltonian for Lanczos

                cout << currentDateTime() << " Starting Lanczos " << endl;

                //---------- Diagonalize and get Eigenvector -------
                energy = lancz.Diag(HV, 1, prm.valvec_, eVec); // Hamiltonian, # of eigenvalues to converge, 1 for -values only, 2 for vals AND vectors

                /* //examine the symmetry of the eigenvector    
                   long unsigned int vi; 
                   for(int v = 0; v<eVec.size(); v++){
                   cout << v << " " << HV.BasPos[swap_layers(HV.Basis[v],fileGraphs.at(i).NumberSites)] <<  " " << eVec[v] << "    " << 
                   eVec[v] + eVec[HV.BasPos[swap_layers(HV.Basis[v],fileGraphs.at(i).NumberSites)]]<< " "  
                   << eVec[v] - eVec[HV.BasPos[swap_layers(HV.Basis[v],fileGraphs.at(i).NumberSites)]] << endl;
                   }
                 */

                cout << currentDateTime() << " energy = " << energy << "  eVec size = " << eVec.size() << endl;
                HV.BasPos.resize(0);
                Entropy2D(alphas, eVec, entVec, fileGraphs.at(i).RealSpaceCoordinates, HV.Basis);
                }
                //---------- Push Back the results --------
                WeightEnergy=energy;
                //Loop Here!!!  ALSO MAKE NOTE THAT LINE IS FIRST AND CORNER IS SECOND !_!_!_!_!_!_!_!_!_!_!_!_!_!
                for(int a=0; a<numRenyis; a++){
                    WeightLineEntropy[a] = entVec[a].first;
                    WeightCornerEntropy[a] = entVec[a].second;
                }



                //Output the Data!!   

            cout <<"Graph " << setw(3) << fileGraphs.at(i).Identifier <<  " Sites " <<setw(2)<< fileGraphs.at(i).NumberSites 
                << "  Jp= " << setw(6) << Jperp << "   E=" <<setw(20)<< WeightEnergy << " Line";

            for(int a=0; a<alphas.size(); a++){ 
                cout  << setw(6) <<  alphas[a] << setw(20) << WeightLineEntropy[a] ;
            }
            cout << " Corner";
            for(int a=0; a<alphas.size(); a++){ 
                cout  <<setw (6) << alphas[a] << setw(20) << WeightCornerEntropy[a];
            }
            cout << endl;

            WeightEnergy=0;

            for(int a=0; a<numRenyis; a++){
                WeightLineEntropy[a]=0;
                WeightCornerEntropy[a]=0;
            }
        }
    }

    // fout.close();
    return 0;
}

unsigned long int swap_layers(unsigned long int init, int N){

    unsigned long int final;
    unsigned long int temp;
    
    temp = (1<<N)-1;

    final = ((init & temp)<<N) +  ((init & (temp<<N))>>N); 

    return final;
}
