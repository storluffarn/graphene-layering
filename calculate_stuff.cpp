
#include <iostream>
#include "modtomlin.h"

using namespace std;

int main()
{

    uint choice;

    cout << "Hello Astrid! I'm David's printing program. At the moment my current implemented functions are: " << endl << endl << "1. Calculate the potential energy in a point " << endl << "2. Calculate the forces in a point " << endl << "3. Printing all parameter values used" << endl << endl << "Please chose one of the listed functionalities (enter number): " << endl;
    
    cin >> choice;
    
    tomlin afm(spring,supvel,latcon,align,barr1,barr2,kappa1,kappa2, nu2,nu4,temp,0,0,xmass,qmass,xdamp,qdamp, "","","","","","");
    
    switch (choice)
    {
        case 1 :
        {
            double xcoord, qcoord, spos;
            
            cout << "give x coordinate (mn): " << endl;
            cin >> xcoord;
            cout << "give q coordinate (nm): " << endl;
            cin >> qcoord;
            cout << "give support position (nm): " << endl;
            cin >> spos;

            afm.setposx(xcoord*1e-9);
            afm.setposq(qcoord*1e-9);
            afm.setsuppos(spos*1e-9);

            afm.calcpot();
            double pot = afm.getpot();

            cout << "the potential energy at (" << xcoord << "," << qcoord << ") with a support positon " << spos << " is: " << endl << pot << endl;

            break;

        }
        case 2 :
        {
            double xcoord, qcoord, spos;
            
            cout << "give x coordinate (nm): " << endl;
            cin >> xcoord;
            cout << "give q coordinate (nm): " << endl;
            cin >> qcoord;
            cout << "give support position (nm): " << endl;
            cin >> spos;

            afm.setposx(xcoord*1e-9);
            afm.setposq(qcoord*1e-9);
            afm.setsuppos(spos*1e-9);

            afm.xacc();
            afm.qacc();
            
            double xacc = afm.getaccx();
            double qacc = afm.getaccq();

            cout << "the forces (x then q) energy at (" << xcoord << "," << qcoord << ") with a support positon " << spos << " is: " << endl << xacc << endl << qacc << endl;

            break;
        }
        case 3:
        {
            afm.printins();

            break;
        }
        default:
        {
            cout << "Invalid mode chosen, exiting" << endl;

        break;
        }
    }

}
