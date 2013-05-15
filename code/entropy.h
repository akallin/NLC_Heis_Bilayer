//function to calculate the RDM and entropies for a 1d 
#ifndef entropy_H
#define entropy_H

void getEE( vector <double>& alpha1, vector <double>& CornLineEnts, vector< vector<long double> >& SuperMat );


inline double TwoSiteEntropy(double h, double alpha)
{
  double CommonEnt;
  double DiffEnt;
  double unLog;
  CommonEnt = 0.5 + (1. + sqrt(1. + 4.*h*h))/(8.*h*h);
  DiffEnt = h*sqrt(1.+2.*h*h+sqrt(1.+4.*h*h))/2./sqrt(2.0);
  
  unLog = pow(CommonEnt + DiffEnt,alpha) + pow(abs(CommonEnt - DiffEnt),alpha);
  
  // cout << CommonEnt << "  " << DiffEnt << endl;
  if(alpha==1.0){
    // cout << CommonEnt - DiffEnt << endl;
    return -(CommonEnt + DiffEnt)*log(CommonEnt + DiffEnt) 
      - abs(CommonEnt - DiffEnt)*log(abs(CommonEnt - DiffEnt));
  }
  else{
    return (1./(1.-alpha))*log(unLog);
  }
}


inline double Magnetization( Array<l_double,1>& eigs ){
  // The dimension is number of eigenvalues
  long int Dim = eigs.size();

  // Get number of sites from the dimension
  int Nsite = log2(Dim); 

  // -8-8-8-8- Measure the Magnetization!!! -8-8-8-8-
  long double magnetization(0);
  double norm(0);
  int itemp(0);
  for(int i=0; i<Dim; i++){ 
    for (int sp=0; sp<Nsite; sp++){
      itemp += (i>>sp)&1; 
    }
    
    //cout << "   i=" <<i << ",  " << abs(itemp*2-Nsite) << "  eig(i)=" <<  eigs(i);
    //magnetization += sqrt((itemp*2-Nsite)*(itemp*2-Nsite))*eigs(i)*eigs(i);
    magnetization += (itemp*2.0 - Nsite*1.0)*eigs(i)*eigs(i);
    itemp=0;
    norm += eigs(i)*eigs(i);
  }
  //cout << "norm="<< norm << endl; norm=0;
  //cout << "mag= " <<magnetization<<endl;
  return magnetization;
  // -8-8-8-8- End of Magnetization -8-8-8-8-
}
		      
inline void Entropy2D(vector <double>& alpha1, Array<l_double,1>& eigs, vector< pair<double,double> >& ents, vector< vector< int > >& RScoords)
 {
  // Get the graph dimensions from the realspace coordinates
  int xMax = RScoords.size();
  int yMax = RScoords[0].size();

  // The dimension is number of eigenvalues
  long int Dim = eigs.size();

  // Get number of sites from the dimension
  int Nsite = log2(Dim); 

  // The dimensions of region A
  int xSize(0), ySize(0), Adim(0), Bdim(0);

  // A rectangular matrix containing the eigenvalues, used to get the RDM
  vector< vector< long double > > SuperMat;

  // Some temp variables;
  int tempState(-1);         // The current full basis state we're looking at
  int tempSpin(-1);          // The number of the spin that's currently being extracted
  int spinState(-1);         // The state of that spin
  int aState(0), bState(0);  // The basis states for reg A and B extracted from the full basis
  vector <double> tempEnt;

  // make the entropy vector a nonzero size
  tempEnt.resize(alpha1.size());
  for(int a=0; a<ents.size(); a++){
    ents[a].first = 0;
    ents[a].second = 0;
    tempEnt[a] = 0;
  }

  // ------ Line Terms!! ------
  
  // -*-*-*- Horizontal -*-*-*-
  xSize = xMax;
  // Iterate over the horizontal cuts
  for(int ySize=1; ySize<=yMax/2; ySize++){
    // Get the dimensions of region A and B;
    Adim = 1<<xSize*ySize; 
    Bdim = Dim/Adim; 

    // Initialize the matrix of eigenvalues
    SuperMat.resize(Adim);
    for(int q=0; q<Adim; q++){ SuperMat[q].resize(Bdim); }

    // Loop over all the basis states
    for(int i=0; i<Dim; i++){      
      // extractifying the region A and region B states
      tempState = i;
      
      // Loop over region A
      aState=0; // Initialize the state in region A
      for(int y=0; y<ySize; y++){
	for(int x=0; x<xMax; x++){

	  // Figure out the spin number given the x,y coords
	  tempSpin = RScoords[x][y];

	  // Extract the state of tempSpin
	  spinState = ((tempState&(1<<tempSpin))>>tempSpin);

	  // Add the spin state to region A
	  aState += spinState;

	  // Shift the bits by 1 (for the next site)
	  aState = aState<<1;
	}
      }	
      // Unshift aState by 1 (because there was one extra)
      aState = aState>>1;

      // Loop over region B (note y starts at ySize)
      bState=0; // Initialize the state in region B
      for(int y=ySize; y<yMax; y++){
	for(int x=0; x<xMax; x++){

	  // Figure out the spin number given the x,y coords
	  tempSpin = RScoords[x][y];

	  // Extract the state of tempSpin
	  spinState = ((tempState&(1<<tempSpin))>>tempSpin);

	  // Add the spin state to region B
	  bState += spinState;

	  // Shift the bits by 1 (for the next site)
	  bState = bState<<1;
	}
      }	
      // Unshift bState by 1 (because there was one extra)
      bState = bState>>1;

      SuperMat[aState][bState] = eigs(i);
    }
    
    // ------ GET ENTROPY!!! ------
    getEE(alpha1, tempEnt, SuperMat);
    for(int a=0; a<alpha1.size(); a++){
      ents[a].first += tempEnt[a];
      ents[a].second += -(xMax-1)*tempEnt[a];
      if(ySize<(yMax+1)/2){ents[a].first += tempEnt[a];  ents[a].second += -(xMax-1)*tempEnt[a];}
    }
    


    //cout << "Adim " << Adim << "  Bdim " << Bdim << "  Hent=" << getEE(alpha,SuperMat) <<endl;
      
  
  // In the future we can just multiply all renyis by 2 except the middle one for an even system.
  }

  // -*-*-*- Vertical -*-*-*-
  ySize = yMax;
  // Iterate over the vectical cuts
  for(int xSize=1; xSize<=xMax/2; xSize++){
    // Get the dimensions of region A and B;
    Adim = 1<<xSize*ySize; 
    Bdim = Dim/Adim; 

    // Initialize the matrix of eigenvalues
    SuperMat.resize(Adim);
    for(int q=0; q<Adim; q++){ SuperMat[q].resize(Bdim); }
  
    // Loop over all the basis states
    for(int i=0; i<Dim; i++){      
      // extractifying the region A and region B states
      tempState = i;
      
      // Loop over region A
      aState=0; // Initialize the state in region A
      for(int y=0; y<yMax; y++){
	for(int x=0; x<xSize; x++){

	  // Figure out the spin number given the x,y coords
	  tempSpin = RScoords[x][y];

	  // Extract the state of tempSpin
	  spinState = ((tempState&(1<<tempSpin))>>tempSpin);

	  // Add the spin state to region A
	  aState += spinState;

	  // Shift the bits by 1 (for the next site)
	  aState = aState<<1;
	}
      }	
      // Unshift aState by 1 (because there was one extra)
      aState = aState>>1;

      // Loop over region B (note y starts at ySize)
      bState=0; // Initialize the state in region B
      for(int y=0; y<yMax; y++){
	for(int x=xSize; x<xMax; x++){

	  // Figure out the spin number given the x,y coords
	  tempSpin = RScoords[x][y];

	  // Extract the state of tempSpin
	  spinState = ((tempState&(1<<tempSpin))>>tempSpin);

	  // Add the spin state to region B
	  bState += spinState;

	  // Shift the bits by 1 (for the next site)
	  bState = bState<<1;
	}
      }	
      // Unshift bState by 1 (because there was one extra)
      bState = bState>>1;

      SuperMat[aState][bState] = eigs(i);
    }
     
    // ------ GET ENTROPY!!! ------
    getEE(alpha1, tempEnt, SuperMat);
    for(int a=0; a<alpha1.size(); a++){
      ents[a].second += -(yMax-1)*tempEnt[a];
      if(xSize<(xMax+1)/2){ ents[a].second += -(yMax-1)*tempEnt[a]; }
    }
  }

  // -*-*-*-*-*- Corner Terms!! -*-*-*-*-*-
  // Iterate over the corner cuts
  for(int ySize=1; ySize<yMax; ySize++){
    for(int xSize=1; xSize<xMax; xSize++){
      // Get the dimensions of region A and B;
      Adim = 1<<xSize*ySize; 
      Bdim = Dim/Adim; 

      // Initialize the matrix of eigenvalues
      SuperMat.resize(Adim);
      for(int q=0; q<Adim; q++){ SuperMat[q].resize(Bdim); }
      
      // Loop over all the basis states
      for(int i=0; i<Dim; i++){      
	// extractifying the region A and region B states
	tempState = i;
	
	// Loop over region A
	aState=0; // Initialize the state in region A
	for(int y=0; y<ySize; y++){
	  for(int x=0; x<xSize; x++){

	    // Figure out the spin number given the x,y coords
	    tempSpin = RScoords[x][y];

	    // Extract the state of tempSpin
	    spinState = ((tempState&(1<<tempSpin))>>tempSpin);
	    
	    // Add the spin state to region A
	    aState += spinState;
	    
	    // Shift the bits by 1 (for the next site)
	    aState = aState<<1;
	  }
	}	
	// Unshift aState by 1 (because there was one extra)
	aState = aState>>1;

	// Loop over region B (loop over whole state, but do nothing when in region A)
	bState=0; // Initialize the state in region B
	for(int y=0; y<yMax; y++){
	  for(int x=0; x<xMax; x++){
	    if(y<ySize && x<xSize){ continue; }
	    
	    // Figure out the spin number given the x,y coords
	    tempSpin = RScoords[x][y];

	    // Extract the state of tempSpin
	    spinState = ((tempState&(1<<tempSpin))>>tempSpin);

	    // Add the spin state to region B
	    bState += spinState;

	    // Shift the bits by 1 (for the next site)
	    bState = bState<<1;
	  }
	}	
	// Unshift bState by 1 (because there was one extra)
	bState = bState>>1;
	
	SuperMat[aState][bState] = eigs(i);
      }
      
      // ------ GET ENTROPY!!! ------
      getEE(alpha1, tempEnt, SuperMat);
      for(int a=0; a<alpha1.size(); a++){
	ents[a].second += 2.*tempEnt[a];
      }
    }
  }
}

void getEE(vector <double> & alpha2, vector<double > & CornLineEnts, vector< vector<long double> >& SuperMat ){
  
  // The Density Matrix
  Array <double,2> DM;
  long double temp(0);
  int Dim(0);

  // Using SuperMat to get the density matrix
  // If Adim > Bdim TRANSPOSE!!
  if(SuperMat.size()>SuperMat[0].size()){
    Dim = SuperMat[0].size();
    DM.resize(Dim,Dim);
    for(int i=0; i<Dim; i++){
      for(int j=0; j<Dim; j++){
	temp=0;
	for(int k=0; k<SuperMat.size(); k++){
	  temp += SuperMat[k][i]*SuperMat[k][j];
	}
	DM(i,j) = temp;
      }
    }
  }
  // Otherwise, use Adim
  else{
    Dim = SuperMat.size();
    DM.resize(Dim,Dim);
    for(int i=0; i<Dim; i++){
      for(int j=0; j<Dim; j++){
	temp=0;
	for(int k=0; k<SuperMat[0].size(); k++){
	  temp += SuperMat[i][k]*SuperMat[j][k];
	}
	DM(i,j) = temp; 
      }
    }
  }

  // Eigenvalues of the RDM get put in dd
  vector<double> dd;
  
  //Diagonalizing the RDM
  while(dd.size()>0){dd.erase(dd.begin());}
  diagWithLapack_R(DM,dd);
 
  double EE(0);
  double vN(0), renyi(0); 
  temp=0;
  

  for(int a=0; a<alpha2.size(); a++){
    EE=0;
    vN=0;
    renyi=0;
    temp=0;

    // Loop over the eigenvalues
    for(int s=0; s<dd.size(); s++){
      
      if(dd[s]<0){dd[s]=0;} // All eigs should be positive.  If not it's rounding error.
      // If dd[s] is a verrrrry small number, it's probably zero.
      temp=log(dd[s]); 
      if(!(temp>-1000000000)){temp=0;}
      
      vN += -dd[s]*temp;
 
      // Same problem. If they're too small they get set to 0.
      if(abs(dd[s])<1e-15){dd[s]=0;}
      
      renyi+=pow(dd[s],alpha2[a]);
    }    
    
    if(abs(alpha2[a]-1.0)<0.000001){EE = vN;}
    else{EE = 1./(1.-alpha2[a])*log(renyi);}
    
    CornLineEnts[a] = EE;
  }
}

#endif
