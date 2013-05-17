//c++ class for creating a general Hamiltonian in c++ vectors
//Roger Melko, November 2007
#ifndef GenHam_H
#define GenHam_H

#include <iostream>
#include <vector>
using namespace std;

class GENHAM{

public:
    int Vdim; //dimenson of reduced Hilbert space
  
    vector<vector<long> > PosHam;
    vector<vector< long double > > ValHam;
    //vector<double> DiagHam;
  
    vector<long> Basis;
    vector<long> BasPos;
    vector<int> ConnectCount;
  
    vector< vector< double > > Ham;
  
  GENHAM(const int N_ ,const long double J_, const long double Jh_, vector < pair<int,int> > BBond_, bool Low_, double mag_); 

  //  GENHAM(const int N_ ,const long double J_, const long double Jh_, vector < pair<int,int> > BBond_); 
    void printg();

    vector< double > apply( const vector< double > & );
  
    void SparseHamJQ();

private:
    int Nsite; //number sites
    bool LowField; //High or Low Field expansion
  
    vector< pair < int,int> > Bond;

    long double JJ; //heisenberg exchange value
    long double hh; //next-nearest neighbor exchange value
  long double JJ2; // J2 value (for heis bilayer)
  double Mag;
  
    double HdiagPart(const long, int);
    double HdiagPart(const long);
    double HOFFdBondX(const int, const long);
    double HOFFdBondY(const int, const long);

};

#endif
