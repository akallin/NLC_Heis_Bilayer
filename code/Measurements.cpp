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
    double WeightEnergy; // WeightMagnetization;
    vector<double> WeightLineEntropy, WeightCornerEntropy;

    // Self explanatory (good naming convention, Katie)
    ReadGraphsFromFile(fileGraphs, InputFile);

    //ofstream fout(OutputFile.c_str());
    //fout.precision(10);
    cout.precision(10);

    J=1;     

    const int numJperpVals = 1;
    double Jperpvals[numJperpVals] = {1};

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
                        if(fileGraphs.at(i).RealSpaceCoordinates[0].size()==2){
                            entVec[a].first=log(2.0); 
                        }

                        else{ entVec[a].first=0; }
                        entVec[a].second=0;       
                    }
                }
                else{
                    energy = -1.0*J;
                    double eig1(5./6.);
                    double eig2(1./6.);
                    for(int a=0; a<numRenyis; a++){
                        if(fileGraphs.at(i).RealSpaceCoordinates[0].size()==3){
                            if(fabs(1.0-(alphas[a]))<0.000001){entVec[a].first=2.*(-eig1*log(eig1)-eig2*log(eig2));}
                            else{entVec[a].first=(2./(1.-alphas[a]))*log(pow(eig1,alphas[a])+pow(eig2,alphas[a]));}
                        }
                        else{
                            entVec[a].first=0;
                        }
                        entVec[a].second=0;
                    }
                }
            }

            else{
                //---Generate the Hamiltonian---
                GENHAM HV(fileGraphs.at(i).NumberSites,J,Jperp,fileGraphs.at(i).AdjacencyList); 

                LANCZOS lancz(HV.Vdim);  //dimension of Hilbert space
                //HV.SparseHamJQ();  //generates sparse matrix Hamiltonian for Lanczos

                //---------- Diagonalize and get Eigenvector -------
                energy = lancz.Diag(HV, 1, prm.valvec_, eVec); // Hamiltonian, # of eigenvalues to converge, 1 for -values only, 2 for vals AND vectors
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
                << "  Jp= " << setw(6) << Jperp << "   E=" <<setw(16)<< WeightEnergy << " Line";

            for(int a=0; a<alphas.size(); a++){ 
                cout  << setw(4) <<  alphas[a] << setw(16) << WeightLineEntropy[a] ;
            }
            cout << " Corner";
            for(int a=0; a<alphas.size(); a++){ 
                cout  <<setw (4) << alphas[a] << setw(16) << WeightCornerEntropy[a];
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
