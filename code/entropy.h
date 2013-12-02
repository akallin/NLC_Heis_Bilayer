//function to calculate the RDM and entropies for a 1d 
#ifndef entropy_H
#define entropy_H
#include<time.h>
#include<Eigen/Dense>



using Eigen::MatrixXd;
using Eigen::VectorXd;

//void getEE( vector <double>& alpha1, vector <double>& CornLineEnts, vector< vector<double> >& SuperMat );
void getEE( vector <double>& alpha1, vector <double>& CornLineEnts, Eigen::MatrixXd&  SuperMat );
long unsigned int regionDim_NA_N( unsigned int na, unsigned int n, vector<unsigned long int> &Abasis, vector<unsigned long int> &AbasPos );

const std::string currentDateTime();

inline void Entropy2D(vector <double>& alpha1, vector<l_double>& eigs, vector< pair<double,double> >& ents, 
        vector< vector< int > >& RScoords, vector <long unsigned int> Basis){



    // Get the graph dimensions from the realspace coordinates
    int xMax = RScoords.size();
    int yMax = RScoords[0].size();
    int Nsite = xMax*yMax;

    // The dimension is number of eigenvalues
    long unsigned int Dim = eigs.size(); // long needed if we go to... 34 sites (34_C_17 > 2^31)

    // The dimensions of region A/B. 
    unsigned int xSize(0), ySize(0); 
    long unsigned int Adim, Bdim; // could require long
    vector <long unsigned int> Abasis, Bbasis, AbasPos, BbasPos;

    // A rectangular matrix containing the eigenvalues, used to get the RDM
    //vector< vector<double > > SuperMat;
    MatrixXd DuperMat;

    // Some temp variables;
    long unsigned int tempState(-1);         // The current full basis state we're looking at
    int tempSpin(-1);          // The number of the spin that's currently being extracted
    int spinState(-1);         // The state of that spin
    long unsigned int aState(0), bState(0);  // The basis states for reg A and B extracted from the full basis
    vector <double> tempEnt;

    // make the entropy vector a nonzero size
    tempEnt.resize(alpha1.size());
    for(int a=0; a<ents.size(); a++){
        //line ents
        ents[a].first = 0;

        //corner ents
        ents[a].second = 0;

        tempEnt[a] = 0;
    }

    // --------------- Line Terms!! ---------------

    // -*-*-*-*-*-*-*- Horizontal -*-*-*-*-*-*-*-
    xSize = xMax;
    // Iterate over the horizontal cuts
    for(int ySize=1; ySize<=yMax/2; ySize++){
        cout << currentDateTime() << " Line term H" << ySize;

        // Get the dimensions of region A and B;
        // states don't necessarily have Sz=0 in their regions 
        // if NA > N/2 or NB > N/2
        Adim = regionDim_NA_N(xSize*ySize, Nsite, Abasis, AbasPos);
        cout << "    Adim = " << Adim;
        Bdim = regionDim_NA_N(Nsite - xSize*ySize,Nsite, Bbasis, BbasPos);
        cout << "  Bdim = " << Bdim << endl;
        //cout << "Adim = " << Adim << "  Bdim = " << Bdim << endl;

        // Initialize the matrix of eigenvalues
        cout << currentDateTime() << " Initialize Supermat" << endl;
        //SuperMat.clear();
        //SuperMat.resize(Adim);

        DuperMat.resize(Adim,Bdim);
        for(int ad=0;ad<Adim;ad++){
            for(int bd=0;bd<Bdim;bd++){
                DuperMat(ad,bd)=0;
            }
        }
        //for(long unsigned int q=0; q<Adim; q++){ SuperMat[q].resize(Bdim); }
        cout << currentDateTime() << " .... Supermat created" << endl;

        // Loop over all the basis states
        cout << currentDateTime() << " Looping over basis states" << endl;

        for(long unsigned int i=0; i<Dim; i++){      
            // extractifying the region A and region B states
            tempState = Basis[i];

            // Loop over region A
            aState=0; // Initialize the state in region A
            for(int y=0; y<ySize; y++){
                for(int x=0; x<xMax; x++){

                    // Layer 1 --------------------------------------
                    // Figure out the spin number given the x,y coords
                    tempSpin = RScoords[x][y];

                    // Extract the state of tempSpin
                    spinState = ((tempState&(1<<tempSpin))>>tempSpin);

                    // Add the spin state to region A
                    aState += spinState;

                    // Shift the bits by 1 (for the next site)
                    aState = aState<<1;

                    // Layer 2 --------------------------------------
                    // Figure out the spin number given the x,y coords
                    tempSpin += Nsite;

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

                    // Layer 1 --------------------------------------
                    // Figure out the spin number given the x,y coords
                    tempSpin = RScoords[x][y];

                    // Extract the state of tempSpin
                    spinState = ((tempState&(1<<tempSpin))>>tempSpin);

                    // Add the spin state to region B
                    bState += spinState;

                    // Shift the bits by 1 (for the next site)
                    bState = bState<<1;

                    // Layer 2 --------------------------------------
                    // Figure out the spin number given the x,y coords
                    tempSpin += Nsite;

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

            //SuperMat[AbasPos[aState]][BbasPos[bState]] = eigs[i];
            DuperMat(AbasPos[aState],BbasPos[bState]) = eigs[i];
        }
        cout << currentDateTime() <<"   ... Supermat filled" << endl;

        // ------ GET ENTROPY!!! ------
        getEE(alpha1, tempEnt, DuperMat);
        for(int a=0; a<alpha1.size(); a++){
            ents[a].first += tempEnt[a];
            ents[a].second += -(xMax-1)*tempEnt[a];
            if(ySize<(yMax+1)/2){ents[a].first += tempEnt[a];  ents[a].second += -(xMax-1)*tempEnt[a];}
        }
    }



    // -*-*-*-*-*-*-*- Vertical -*-*-*-*-*-*-*-*-
    ySize = yMax;
    
    //if it's a square no need to do the same measurements twice
    if(xMax == yMax){
        for(int a=0; a<alpha1.size(); a++){
            ents[a].first *=2;
            ents[a].second *=2; 
        }
    }
    else{
        // Iterate over the vectical cuts
        for(int xSize=1; xSize<=xMax/2; xSize++){

            cout << currentDateTime() << " Line term V" << xSize ;

            // Get the dimensions of region A and B;
            Adim = regionDim_NA_N(xSize*ySize, Nsite, Abasis, AbasPos);
            Bdim = regionDim_NA_N(Nsite - xSize*ySize, Nsite, Bbasis, BbasPos);
            cout << "   Adim = " << Adim << "  Bdim = " << Bdim << endl;

            // Initialize the matrix of eigenvalues
            cout << currentDateTime() << " Initialize Supermat" << endl;
            //SuperMat.clear();
            //SuperMat.resize(Adim);
            //for(long unsigned int q=0; q<Adim; q++){ SuperMat[q].resize(Bdim); }
            DuperMat.resize(Adim,Bdim);
            for(int ad=0;ad<Adim;ad++){
                for(int bd=0;bd<Bdim;bd++){
                    DuperMat(ad,bd)=0;
                }
            }
            cout << currentDateTime() << " .... Supermat created" << endl;

            // Loop over all the basis states
            for(long unsigned int i=0; i<Dim; i++){      
                // extractifying the region A and region B states
                tempState = Basis[i];

                // Loop over region A
                aState=0; // Initialize the state in region A
                for(int y=0; y<yMax; y++){
                    for(int x=0; x<xSize; x++){

                        // Layer 1 --------------------------------------
                        // Figure out the spin number given the x,y coords
                        tempSpin = RScoords[x][y];

                        // Extract the state of tempSpin
                        spinState = ((tempState&(1<<tempSpin))>>tempSpin);

                        // Add the spin state to region A
                        aState += spinState;

                        // Shift the bits by 1 (for the next site)
                        aState = aState<<1;

                        // Layer 2 --------------------------------------
                        // Figure out the spin number given the x,y coords
                        tempSpin += Nsite;

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

                // Loop over region B (note x starts at xSize)
                bState=0; // Initialize the state in region B
                for(int y=0; y<yMax; y++){
                    for(int x=xSize; x<xMax; x++){

                        // Layer 1 --------------------------------------
                        // Figure out the spin number given the x,y coords
                        tempSpin = RScoords[x][y];

                        // Extract the state of tempSpin
                        spinState = ((tempState&(1<<tempSpin))>>tempSpin);

                        // Add the spin state to region B
                        bState += spinState;

                        // Shift the bits by 1 (for the next site)
                        bState = bState<<1;

                        // Layer 2 --------------------------------------
                        // Figure out the spin number given the x,y coords
                        tempSpin += Nsite;

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

                if(AbasPos[aState]<0 || BbasPos[bState]<0){ cout << "SUPER ERROR!" << endl; exit(1);}
                //SuperMat[AbasPos[aState]][BbasPos[bState]] = eigs[i];
                DuperMat(AbasPos[aState],BbasPos[bState]) = eigs[i];

            }
            cout << currentDateTime() <<"   ... Supermat filled" << endl;

            // ------ GET ENTROPY!!! ------
            getEE(alpha1, tempEnt, DuperMat);
            for(int a=0; a<alpha1.size(); a++){

                ents[a].first += tempEnt[a];
                ents[a].second += -(yMax-1)*tempEnt[a];
                if(xSize<(xMax+1)/2){ents[a].first += tempEnt[a]; ents[a].second += -(yMax-1)*tempEnt[a]; }
            }
        }
    }

    // -*-*-*-*-*- Corner Terms!! -*-*-*-*-*-
    // Iterate over the corner cuts
    for(int ySize=1; ySize<yMax; ySize++){
        for(int xSize=1; xSize<xMax; xSize++){
           
            if(xMax==yMax){if(xSize<ySize){continue;}} 
            cout << currentDateTime() << " Corner term " << ySize*xSize ;

            // Get the dimensions of region A and B;
            Adim = regionDim_NA_N(xSize*ySize, Nsite, Abasis, AbasPos);
            Bdim = regionDim_NA_N(Nsite - xSize*ySize,Nsite, Bbasis, BbasPos);
            cout << "   Adim = " << Adim << "  Bdim = " << Bdim << endl;

            // Initialize the matrix of eigenvalues
            cout << currentDateTime() << " Initialize Supermat" << endl;
            //SuperMat.clear();
            //SuperMat.resize(Adim);
            //for(long unsigned int q=0; q<Adim; q++){ SuperMat[q].resize(Bdim); }
            DuperMat.resize(Adim,Bdim);
            for(int ad=0;ad<Adim;ad++){
                for(int bd=0;bd<Bdim;bd++){
                    DuperMat(ad,bd)=0;
                }
            }
            cout << currentDateTime() << " .... Supermat created" << endl;

            // Loop over all the basis states
            for(long unsigned int i=0; i<Dim; i++){      
                // extractifying the region A and region B states
                tempState = Basis[i];

                // Loop over region A
                aState=0; // Initialize the state in region A
                for(int y=0; y<ySize; y++){
                    for(int x=0; x<xSize; x++){
                        // Layer 1 --------------------------------------
                        // Figure out the spin number given the x,y coords
                        tempSpin = RScoords[x][y];

                        // Extract the state of tempSpin
                        spinState = ((tempState&(1<<tempSpin))>>tempSpin);

                        // Add the spin state to region A
                        aState += spinState;

                        // Shift the bits by 1 (for the next site)
                        aState = aState<<1;

                        // Layer 2 --------------------------------------
                        // Figure out the spin number given the x,y coords
                        tempSpin += Nsite;

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

                        // Layer 1 --------------------------------------
                        // Figure out the spin number given the x,y coords
                        tempSpin = RScoords[x][y];

                        // Extract the state of tempSpin
                        spinState = ((tempState&(1<<tempSpin))>>tempSpin);

                        // Add the spin state to region B
                        bState += spinState;

                        // Shift the bits by 1 (for the next site)
                        bState = bState<<1;

                        // Layer 2 --------------------------------------
                        // Figure out the spin number given the x,y coords
                        tempSpin += Nsite;

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

                if(AbasPos[aState]<0 || BbasPos[bState]<0){ cout << "SUPER ERROR! " << endl; exit(1);}
                //SuperMat[AbasPos[aState]][BbasPos[bState]] = eigs[i];
                DuperMat(AbasPos[aState],BbasPos[bState]) = eigs[i];

            }
            
            cout << currentDateTime() <<"   ... Supermat filled" << endl;
            // ------ GET ENTROPY!!! ------
            getEE(alpha1, tempEnt, DuperMat);
            for(int a=0; a<alpha1.size(); a++){
                ents[a].second += 2.*tempEnt[a];
                if(xMax==yMax && xSize!=ySize){ents[a].second += 2.*tempEnt[a];}
            }//loop over alphas
        }//loop over xSize
    }//loop over ySize
    
}
//End of Entropy2D


/* -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
   |  regionDim_NA_N - Gives the dimension of the Hilbert space for some region A given:             |
   |  *  na - the number of sites in A                                                               |
   |  *  n  - the total number of sites                                                              |
   -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- */
long unsigned int regionDim_NA_N( unsigned na1, unsigned n1, vector<long unsigned int>& basism, vector<long unsigned int>& basPosm )
{
    unsigned na(2*na1), n(2*n1);
    long unsigned int full_dim = 1<<na; 
    basPosm.clear();
    basPosm.resize(full_dim,-99);
    long unsigned int dimm;
    basism.clear();
    basism.resize(0);


    /* If the region is less than or equal to half  |
       |  of the total number of sites (we have a full |
       |  unrestricted basis) Dim = 2^na              */
    if (na<=n/2){
        dimm = full_dim;	
        for (unsigned long i1=0; i1<full_dim; i1++) 
        {
            //Trivial Basis
            basism.push_back(i1);
            basPosm.at(i1)=basism.size()-1;
        }
    }

    // Otherwise Dim is more complicated!
    else{
        dimm = 0;
        long unsigned int temp(0);
        for (unsigned long int i1=0; i1<full_dim; i1++) 
        {
            temp = 0;
            for (int sp =0; sp<na; sp++)
                temp += (i1>>sp)&1;  //unpack bra & count the up spins

            //can't have more than N/2 up spins or down spins
            if (temp<=((n+1)/2) && temp>=(na-(n+1)/2) ){ 
                basism.push_back(i1);
                basPosm.at(i1)=basism.size()-1;
                dimm++;
            }
        }
    }
    return dimm;
}//End of regionDim_NA_N


/* -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
   getEE
   -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- */
void getEE(vector <double> & alpha2, vector<double > & CornLineEnts, Eigen::MatrixXd& SuperMat ){
    
      cout << currentDateTime() << " Multiplying to get DM" << endl;
      int DMdim;
      MatrixXd DM;
    if(SuperMat.cols() >= SuperMat.rows()){
        DMdim = SuperMat.rows();
        DM.resize(DMdim,DMdim);
        DM = SuperMat * SuperMat.transpose();
    }
    else{
        DMdim = SuperMat.cols();
        DM.resize(DMdim,DMdim);
        DM = SuperMat.transpose() * SuperMat;
    }

      cout << currentDateTime() << "Beginning Diagonalization" << endl;
      VectorXd dd = DM.selfadjointView<1>().eigenvalues();
/*
    cout << currentDateTime() << " Getting EE" << endl;
    // Using SuperMat to get the density matrix
    // If Adim > Bdim TRANSPOSE!!
    if(SuperMat.size()>=SuperMat[0].size()){
        Dim = SuperMat[0].size();
        cout << currentDateTime() << " Dim = " << Dim << "  Dim*Dim = " << Dim*Dim << endl;
        //DM.resize(Dim,Dim);
        cout << currentDateTime() << " creating DM ... " << endl;
        Dim_sq = Dim*Dim;
        DM= new double[Dim_sq];  //This is a c-style array
        cout << currentDateTime() << " DM created " << endl;

        for(long int i=0; i<Dim; i++){
            for(long int j=i; j<Dim; j++){
                temp=0;
                for(int k=0; k<SuperMat.size(); k++){
                    temp += SuperMat[k][i]*SuperMat[k][j];
                }
                DM[j*Dim + i] = temp; 
                DM[i*Dim + j] = temp; //matrix is symmetric
            }
            if(i%265==0){cout<<currentDateTime() <<" i = " << i << endl;}
        }
        cout << currentDateTime() << " DM filled" << endl;
    }
    // Otherwise, use Adim
    else{
        Dim = SuperMat.size();
        //DM.resize(Dim,Dim);
        cout << currentDateTime() << " creating DM ... " << endl;
        Dim_sq = Dim*Dim;
        DM= new double[Dim_sq];  //This is a c-style array
        cout << currentDateTime() << " DM created" << endl;

        for(long int i=0; i<Dim; i++){
            for(long int j=i; j<Dim; j++){
                temp=0;
                for(long int k=0; k<SuperMat[0].size(); k++){
                    temp += SuperMat[i][k]*SuperMat[j][k];
                }
                DM[j*Dim + i] = temp;            
                DM[i*Dim + j] = temp; //matrix is symmetric
            }
            if(i%265==0){cout<<currentDateTime()<< " i = " << i << endl;}
        }
        cout << currentDateTime() << " DM filled" << endl;
    }

    // Eigenvalues of the RDM get put in dd
    vector<double> dd;

    //Diagonalizing the RDM
    while(dd.size()>0){dd.erase(dd.begin());}
    cout << currentDateTime() << " Beginning diagonalization" << endl;
    int drim = Dim;
    diagWithLapack_R(DM,dd,drim,drim);
    cout << currentDateTime() << " Diagonalization complete " << endl;

    //clean up DM, unless you need it anywhere below
    delete [] DM;
*/
    double EE(0);
    double renyi(0); 
    VectorXd logdd = dd.array().abs().log().matrix(); 
    VectorXd renyipow;
    double vN;     
    int temp=0;

    vN = - (dd.dot(logdd));
    
    cout << currentDateTime() << " calculating renyis " << endl;
    for(int a=0; a<alpha2.size(); a++){
        EE=0;
        renyi=0;
        temp=0;


        renyipow = dd.array().pow(alpha2[a]).matrix();
        renyi = renyipow.sum();
        // Loop over the eigenvalues
        /*
        for(int s=0; s<dd.size(); s++){

            if(dd[s]<0){dd[s]=0;} // All eigs should be positive.  If not it's rounding error.
            // If dd[s] is a verrrrry small number, it's probably zero.
            temp=log(dd[s]); 
            if(!(temp>-1000000000)){temp=0;}

         
            // Same problem. If they're too small they get set to 0.
            if(fabs(dd[s])<1e-15){dd[s]=0;}

            renyi+=pow(dd[s],alpha2[a]);
        }    */

        if(fabs(alpha2[a]-1.0)<0.000001){EE = vN;}
        else{EE = 1./(1.-alpha2[a])*log(renyi);}

        CornLineEnts[a] = EE;
    }
    cout << currentDateTime() << " EE complete" << endl;
}//End of getEE


const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d[%X]", &tstruct);

    return buf;
}

#endif
