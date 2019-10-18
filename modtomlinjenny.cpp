
//
// extracting minima and maxima from tomlinson
//


// includes

#include "slicingalgs.h"		// I fucked up my dependancies, modtomlin.h is needed but is included in slicingalgs.h...
#include <chrono>

// lazy stuff

using namespace std;

int main()
{
	// standard values
	
	double barref = 4.0e-20; // 4.0e-20
	double kapparef = 0.0612245;
	
	double spring = 1e-0* 2.0;			// 2.0 //1e-4 for harmonic
	double supvel = 1.0;				// 1.0
	double latcon = 2.5e-10;			// 2.5e-10
	double barr1 = barref;				// 
	double barr2 = 0.5 * barref;
	double kappa1 = kapparef;
	double kappa2 = 0.5 * kapparef;
	double align = 1.0;					//(1+sqrt(5.0)/2.0); 
	double nu2 = 0.382653;
	double nu4 = 5.809e17;
	
	double xmass = 1e-23;
	double xdamp = 1.875e13;

	double qmass = 3.67143e-24;
	double qdamp = 4.28571e13;
	
	double temp = 1.0*200;  // 5e48 lest we'd forget
	
	double tmp = 0.5;	// 0.5 gives reliable timestep dep. 1.0 should be ok
	double tstep = tmp * 3e-14;	
	uint tsteps = 1.0/tmp * 2.0e5;	// has to be even beucasue lazyness

	time_t t = time(0);
	struct tm * now  = localtime(&t);
	char buffer [80];
	strftime (buffer,80,"%Y%m%d%H%M%S",now);

	string tfile = "time.csv";
	string xfile = "xout.csv";
	string qfile = "qout.csv";
	string ffile = "tomout.csv";
	string avgfile = "avgs.csv";

	//string tfile = "time" + buffer + ".csv";
	//string xfile = "xout" + buffer + ".csv";
	//string qfile = "qout" + buffer + ".csv";
	//string ffile = "tomout" + buffer + ".csv";

	//uint periods = static_cast <uint> (tstep*tsteps*supvel/latcon);
	

	ofstream fspos;
	fspos.open("slips.csv");

	// for averaging
	uint skip = 50;			// probably should be a fraction of mean	
	uint halfmeansize = 1000;			// half interval, such that | mean -- mid -- mean |
	
	// for simple slip
	//uint adj = 10;
	//double ttoa = latcon/(tstep*supvel);	// timesteps to minima
	//uint stride = static_cast <uint> (ttoa / 4.0);
	//uint end = tsteps;
	uint adj = 1;
	uint stride = 50;
	uint end = round(tsteps/skip)-skip;		// this is a bit risky, size should be constant over many runs though...
	
	uint runs = 250;

	for ( uint l = 0; l < runs; l++)
	{	
		tomlin afm(spring,supvel,latcon,align,barr1,barr2,kappa1,kappa2,
			   nu2,nu4,temp,tstep,tsteps,xmass,qmass,xdamp,qdamp,
			   tfile,xfile,qfile,ffile,avgfile);
		vector <uint> slips;

		//chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

		for ( uint k = 0; k < tsteps; k++ )
		{
			//cout << "now begnning loop: " << k << endl;

			afm.rk4();
			afm.calcfric();
					
			// multiple algorithms available, see separate .h file
			
			afm.pushvals();			// remove production runs, just debugging!
			afm.inctime();
		}
		
		//chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

		//double looptime = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	
		afm.noisered(halfmeansize,skip);	
		//afm.writedata();		// ONLY FOR DIAGNOSTICS REMOVE LATER
		
		//t1 = std::chrono::high_resolution_clock::now();
	
		// this rn business is fucking ugly...	
		//halfintervals(adj, end, stride, afm.gettimes(),afm.getfrics(), &slips);
		halfintervals(adj, end, stride, afm.getrntimes(),afm.getrnfrics(), &slips);

		for (auto &el : slips)
		{
			//fspos << afm.getrntime(el) << "," << afm.getrnfric(el) << endl;
			fspos << afm.getrntime(el) << "," << afm.getrnfric(el) << endl;
		}
		
		//t2 = std::chrono::high_resolution_clock::now();
		
		//double slicetime = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

		//cout << "loop time: " << looptime << endl << "slicetime " << slicetime << endl;
		
		if ((l+1) % 1 == 0)
		{
			cout << "finished " << l+1 << " iterations" << endl;
		}
	}	

	fspos.close();

}




