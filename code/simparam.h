#ifndef PARAMSIM_H
#define PARAMSIM_H


//Class to read in the simulation parameters from a file
// Adapted to the JQ ED code

class PARAMS
{
  public:
    double JJ_; //the heisenberg exchange
    double Jperp_; //the bilayer heisenberg exchange

    int valvec_; //  1 for -values only, 2 for vals AND vectors
    // FULL_DIAG?

    PARAMS(){
      //initializes commonly used parameters from a file
    
      JJ_ = 1.0;
      Jperp_ = 0.5;
      valvec_ = 2;
    
    }//constructor

}; //PARAMS

#endif
