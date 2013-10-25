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

        //The sparse matrix represenation of the Hamiltonian
        vector<vector<long> > PosHam;
        vector<vector< long double > > ValHam;

        //The Hilbert Space
        vector<long> Basis;
        vector<long> BasPos;

        //The constructor
        GENHAM(const int N_ ,const long double J_, const long double Jh_, 
                vector < pair<int,int> > BBond_); 

        void printg();

        //The function that makes the Hamiltonian
        void SparseHamJQ();
        
        //The bond list (needs to be public to calc diag terms on the fly)
        vector< pair < int,int> > Bond;

    private:
        int Nsite; //number sites
        bool LowField; //High or Low Field expansion

        
        long double JJ; //heisenberg exchange value
        long double JJ2; // J2 value (for heis bilayer)

        double HdiagPart(const long, int);
        double HdiagPart(const long);
        double HOFFdBondX(const int, const long);

};
#endif
