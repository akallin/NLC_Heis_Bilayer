#include<iostream>
#include<vector>
#include<algorithm>
using namespace std;

int main(){

    unsigned maxOrder(15);
    unsigned ID(0);
    bool corner(false);
    //bool corner(true);

    vector< pair<int, int> > idList;
    unsigned LattConst(0), subLattConst(0);
    signed temp;

    //Loop over the orders
    for(unsigned Nsite = 1; Nsite <= maxOrder; Nsite++){

        for(unsigned y = 1 + corner; y < maxOrder; y++){
            for(unsigned x = y; x <= maxOrder/y; x++){
                if(x*y != Nsite){continue;}

                LattConst = (x!=y)+1;
                cout << ID << " " << Nsite << " " << LattConst << " 0\n";

                idList.push_back(make_pair(x,y));       

                // ---------- Site coordinates (x,y) pairs 
                for(unsigned y1=0; y1<y; y1++){
                    for(unsigned x1=0; x1<x; x1++){
                        cout << x1 << " " << y1 << " ";
                    }
                }
                cout << endl; // site coordinates

                // ---------- Bonds
                // x-bonds
                for(unsigned y1=0; y1<y; y1++){
                    for(unsigned x1=0; x1<x-1; x1++){
                        cout << x1+y1*x << " " << x1+y1*x + 1 << " ";
                    }
                }//x-bonds
                // y-bonds
                for(unsigned y1=0; y1<y-1; y1++){
                    for(unsigned x1=0; x1<x; x1++){
                        cout << x1+y1*x << " " << x1+(y1+1)*x << " ";
                    }
                }//y-bonds
                cout << endl;

                for(unsigned sub=0; sub < idList.size()-1; sub++){ 
                
                    temp = (x+1 - idList[sub].first)*(y+1 - idList[sub].second);
                    subLattConst = max(0,temp);
                    if(idList[sub].first!=idList[sub].second){
                        temp =(y+1 - idList[sub].first)*(x+1 - idList[sub].second);
                        subLattConst += max(0,temp);
                    }
                    if(subLattConst > 0){cout<<sub<<" "<<subLattConst<<" ";}
                    //for subgraphs remember to read the 1site case separately.
                }
                cout << endl;
                ID++;

            }
        }
    }
}
