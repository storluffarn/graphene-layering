
//
// single run modified tomlinson model
//

// includes

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "modtomlin.h"

// Lazy stuff

using namespace std;

typedef unsigned int uint;

// objects

int main()
{
	//input values
	//double barref = 5.0e-20; // 4.0e-20
	//double kapparef = 0.0612245;
	//
	//double spring = 1e-0* 2.0;				// 2.0 //1e-4 for harmonic
	//double supvel = 1.0;				// 1.0
	//double latcon = 2.5e-10;			// 2.5e-10
	//double barr1 = barref;				// 
	//double barr2 = 0.5 * barref;
	//double kappa1 = kapparef;
	//double kappa2 = 0.5 * kapparef;
	//double align = 1.0+0.00;
	//double nu2 = 0.382653;
	//double nu4 = 5.809e17;
	//
	//double xmass = 1e-23;
	//double xdamp = 1.875e13;

	//double qmass = 3.67143e-24;
	//double qdamp = 4.28571e13;
	//
	//double temp = 250;  //5e48
	
	double modif = 0.5;	// 0.5 gives reliable time step dep. 1.0 should be ok
	double tstep = modif * 3e-14;	
	uint tsteps = 1.0/modif * 2.0e5;	// has to be even beucasue lazyness

	uint ttoa = ceil(latcon/(tstep));	// timesteps to minima

	string tfile = "time.csv";
	string xfile = "xout.csv";
	string qfile = "qout.csv";
	string tomfile = "tomout.csv";

	tomlin afm(spring,supvel,latcon,align,barr1,barr2,kappa1,kappa2,
			   nu2,nu4,temp,tstep,tsteps,xmass,qmass,xdamp,qdamp,
			   tfile,xfile,qfile,tomfile,"","");
	
    
    //afm.setposx(1e-10);
    
    // incom
	//afm.setposx(1.23912011781616e-09);
	//afm.setvelx(0.9995276370265342);
	//afm.setaccx(18619338513200.46);
	//afm.setposq(6.874289088010818e-10);
	//afm.setvelq(1.061307759460043);
	//afm.setaccq(45478933617766.52);
    //afm.setsuppos(2.443094999995714e-09);
    
    // com
	//afm.setposx(1.030613785802745e-09);
	//afm.setvelx(0.9995431014034264);
	//afm.setaccx(18645712950047.95);
	//afm.setposq(7.422256383340498e-10);
	//afm.setvelq(0.5066323077421474);
	//afm.setaccq(21750153167849.19);
    //afm.setsuppos(2.448359999995651e-09);
    

	for ( uint k = 0; k < tsteps; k++ )
	{
		//if (k == 20.0*ttoa)
		//{	
		//	afm.tpause();
		//	cout << "resuming at tstep / t = " << k << " / " << k*tstep << endl;
        //}
		//	
		//if (k == 0.25*tsteps)
        //{
		//	afm.treverse();
		//	cout << "reversing at tstep / t = " << k << " / " << k*tstep << endl;
        //}
		//if (k == 0.7*tsteps) //25 and 70 gives more or less match up
        //{
		//	afm.treverse();
		//	cout << "reversing at tstep / t = " << k << " / " << k*tstep << endl;
        //}
		//if (k == 20*ttoa)
		//{
		//	cout << "pausing at tstep / t = " << k << " / " << k*tstep << endl;
		//	afm.tpause();
		//}
		//else if (k == 0.9*tsteps)
		//{
		//	afm.tpause();
		//	cout << "resuming at tstep / t = " << k << " / " << k*tstep << endl;
		//}

		afm.rk4();
		afm.calckin();
		afm.calcpot();
		afm.calcfric();

		afm.pushvals();
		afm.inctime();
	}
	
	//afm.printins();
	//afm.printouts();
	//afm.printavgs();
	afm.writedata();
}

	
