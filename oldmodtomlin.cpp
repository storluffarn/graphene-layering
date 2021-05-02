
//
// single run modified tomlinson model
//


// includes

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "oldmodtomlin.h"

// Lazy stuff

using namespace std;

typedef unsigned int uint;

// objects

int main()
{
	//input values
	
	double spring = 1.0;
	double supvel = 1.0;
	double latconst = 2.5e-10;

	double adhesiontip = 1.3e-21;		// 1.3e-20
	double stiffness = 0.02;				// 0.2
	double coupling = 0.1;
	double adhesionsub = 8.0e17;
	
	double xmass = 1e-23;
	double xveldamp = 3e12;

	double qmass = 5e-22;
	double qposdamp = 0e21;
	double qveldamp = 3e12;

	double timestep = 1.0e-13;
	uint timesteps = 5.0e4;

	string tfile = "time.csv";
	string xfile = "xout.csv";
	string qfile = "qout.csv";
	string tomfile = "tomout.csv";
	
	tomlin afm(spring,supvel,latconst,adhesiontip,stiffness,
			   coupling,adhesionsub,timestep,timesteps,xmass,
			   qmass,qposdamp,xveldamp,qveldamp,xfile,qfile,tomfile,tfile);

	for ( uint k = 0; k < timesteps; k++ )
	{
		afm.rk4();
		afm.kinetic();
		afm.potential();
		afm.friction();

		afm.storevals();
		afm.inctime();
	}

	cout << "final x: " << afm.getposx() << endl << "final q: " << afm.getposq() << endl 
		 << "final fric: " << afm.getfrc() << endl;

	afm.printvals();
}


