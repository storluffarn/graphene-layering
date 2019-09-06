
//
// extracting minima and maxima from tomlinson
//


// includes

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include "modtomlin.h"

// lazy stuff

using namespace std;

typedef unsigned int uint;

// some very manual debugging...

int pingcount = 0; 
void ping(const int line) {cout << "ping at line " << line << endl; pingcount++;}
void getpingcount() {cout << pingcount << endl;}

// threshold for float zero
static const double zero = 1e-11;

int main()
{
	// standard values
	
	double barref = 4.0e-20; // 4.0e-20
	double kapparef = 0.0612245;
	
	double spring = 1e-0* 2.0;				// 2.0 //1e-4 for harmonic
	double supvel = 1.0;				// 1.0
	double latcon = 2.5e-10;			// 2.5e-10
	double barr1 = barref;				// 
	double barr2 = 0.5 * barref;
	double kappa1 = kapparef;
	double kappa2 = 0.5 * kapparef;
	double align = 1.0+0.00;
	double nu2 = 0.382653;
	double nu4 = 5.809e17;
	
	double xmass = 1e-23;
	double xdamp = 1.875e13;

	double qmass = 3.67143e-24;
	double qdamp = 4.28571e13;
	
	double temp = 200;  // 5e48 lest we'd forget
	
	double tmp = 0.5;	// 0.5 gives reliable timestep dep. 1.0 should be ok
	double tstep = tmp * 3e-14;	
	uint tsteps = 1.0/tmp * 2.0e5;	// has to be even beucasue lazyness

	//double ttoa = latcon/tstep;	// timesteps to minima

	string tfile = "time.csv";
	string xfile = "xout.csv";
	string qfile = "qout.csv";
	string ffile = "tomout.csv";
	
	struct point
	{
		double x;
		double y;
	};	

	uint periods = static_cast <uint> (tstep*tsteps*supvel/latcon);
	
	tomlin afm(spring,supvel,latcon,align,barr1,barr2,kappa1,kappa2,
			   nu2,nu4,temp,tstep,tsteps,xmass,qmass,xdamp,qdamp,
			   tfile,xfile,qfile,ffile);

	ofstream fspos;
	fspos.open("slippos.csv");

	uint runs = 1;

	for ( uint l = 0; l < runs; l++)	
	{	
		vector <point> toppts;
		vector <point> botpts;
		
		toppts.reserve(2*periods); 
		botpts.reserve(2*periods);

		bool climb = true;

		double lastpos = 0;
		double lastfric = 0;
		double lasttime = 0;
		
		for ( uint k = 0; k < tsteps; k++ )
		{
			//cout << "now begnning loop: " << k << endl;

			afm.rk4();
			afm.calcfric();

			double pos = afm.getposx();
			double fric = afm.getfric();
			double time = afm.gettime();
					
			double diff = fric - lastfric;

			if (!climb && diff > zero)				// did we stick?
			{
				ping(__LINE__);
				point pt {lastpos,lastfric};
				botpts.push_back(pt);

				climb = true;
			}
			
			if (climb && diff < zero)				// did we slip?
			{
				ping(__LINE__);
				point pt {lastpos,lastfric};
				toppts.push_back(pt);
				
				climb = false;
			}

			afm.inctime();
			lastpos = pos;
			lastfric = fric;
			lasttime = time;

			afm.pushvals();			// remove production runs, just debugging!
		}
	
		for (auto &el : toppts)
		{
			fspos << el.x << ",";
		}
		fspos << endl;
		
		if ((l+1) % 1 == 0)
		{
			cout << "finished " << l+1 << " iterations" << endl;
		}
	}	

	fspos.close();

	afm.writedata();
}




