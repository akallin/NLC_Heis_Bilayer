//Lanczos_07.h
//c++ class for performing a Lanczos diagonalization
//Roger Melko, November 2007 - Blitz++ removed by K. Hyatt, 2012.  
//Modified for NLCE May 2013 - R. Melko

#ifndef LANCZOS_07
#define LANCZOS_07


#include <iostream>
#include <limits.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <iomanip>
#include <vector>
using namespace std;

#include "GenHam.h"

//typedef long double l_double;  //precision for lanczos
typedef double l_double;  //precision for lanczos

class LANCZOS{

  public:
    int Dim; //dimension 

   //Methods
   LANCZOS(const int);
   double Diag(const GENHAM &, const int, const int, vector< l_double > &);
   void tred3(vector< vector<double> >& , vector<double>& , vector<double>& e, const int );

  private:
   int STARTIT;  //how many iterations to "always" perform
   l_double CONV_PREC; //convergence precision

   vector<l_double> V0;  
   vector<l_double> V1;    //Ground state vector
   vector<l_double> V2;

   void apply( vector<l_double> &, const GENHAM &, const vector<l_double>&);  //apply H to |V>
   void Normalize(vector<l_double>& );
   int tqli2(vector<l_double>& , vector<l_double>& , int , vector< vector<l_double > > & , const int );

}; //LANCZOS class


#endif
