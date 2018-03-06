
//
// numerical search for strengthening in the modified tomlinson model
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

// Lazy stuff

using namespace std;

typedef unsigned int uint;

struct variable 
{
	string name;
	double base;
	double inc;
	double val;
};

// least square method
double linreg(vector <double>* pts, vector <double>* crd)
{   
	uint n = pts->size();

	if (n == 1)
	{   
		double slope = pts->at(0);

		return slope;
	}   

	double yavg = accumulate(pts->begin(),pts->end(),0.0)/pts->size();
	double xavg = accumulate(crd->begin(),crd->end(),0.0)/crd->size();

	double a = 0;
	double b = 0;
	for (uint k = 0; k < n; k++)
	{   
		a += (pts->at(k) - yavg) * (crd->at(k) - xavg);
		b += pow((crd->at(k) - xavg),2);
	}   

	double slope = a / b;

	return slope;
}


int main()
{
	// standard values
	
	double spring = 1.0;
	double supvel = 1.0;
	double latconst = 2.5e-10;

	const double adtip0 = 1.3e-20;		// made dynamical below
	const double stiff0 = 0.2;
	const double coupl0 = 1.0;
	const double adsub0 = 8.0e18;
	
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

	// simulation parameters

	bool termin = false;
	bool finish = false;
	bool climb = true;
	double lastfric = 0;
	double lasttime = 0;
	double loops = 0;
	uint maxsteps = floor(timesteps / 2);
	double maxfric = 1.0;
	double maxpos = 1.0;
	double offset = 5.0;
	double step = 0.5;
	uint range = ceil(2*offset / step);
	double strengtol = 1.0/25.0;
	double flattol = strengtol / 10;
	uint pts = 5;
	uint nflat = 2;
	
	variable adhesiontip = {"adhesiontip",adtip0,step*adtip0,offset*adtip0};
	variable stiffness = {"stiffness",stiff0,step*stiff0,offset*stiff0};
	//variable coupling = {"coupling",coupl0,step*coupl0,offset*coupl0};
	//variable adhesionsub = {"adhesionsub",adsub0,step*adsub0,offset*adsub0};	

	map <string,variable> varlist;
	varlist.insert(make_pair("adtip",adhesiontip));
	varlist.insert(make_pair("stiff",stiffness));
	//varlist.insert(make_pair("coupl",coupling));
	//varlist.insert(make_pair("adsub",adhesionsub));

	// Remember to update below when changing the number of variables. Ugly, but best I can do atm
	
	ofstream fsslope;
	fsslope.open("slopes.dat");

	fsslope << "// To help parse. Parameters used (in order), with name, standard value, increment, and initial value:";
	
	for (auto& el : varlist)
		fsslope << " " << el.second.name << ": " << el.second.base << "," << el.second.inc << "," << el.second.val << ".";
	
	fsslope << " The step size is " << step << " times the standard value, the offset is " << offset << " times the standard value, giving a total of " << range << " steps. The number given is the slope of the strengthening, 0 means the termination before strengthening was confirmed, see source code.";

	fsslope << endl;
	
	vector <double> slopes;
	slopes.reserve(pts);

	for (uint l1 = 0; l1 <= range; l1++)
	{
	for (uint l2 = 0; l2 <= range; l2++)
	{
	//for (uint l3 = 0; l3 <= range; l3++)
	//{
	//for (uint l4 = 0; l4 <= range; l4++)
	//{
	
	double adtip = varlist["adtip"].val;
	double stiff = varlist["stiff"].val;
	//double coupl = varlist["coupl"].val;
	//double adsub = varlist["adsub"].val;
	double coupl = coupl0;
	double adsub = adsub0;

	tomlin afm(spring,supvel,latconst,adtip,stiff, coupl,adsub,
			   timestep,timesteps,xmass,qmass,qposdamp,xveldamp,
			   qveldamp,xfile,qfile,tomfile,tfile);
	
	vector <double> fricslope;
	vector <double> timeslope;
	vector <double> fricflat;
	vector <double> timeflat;

	// Find maxima
	while (!termin || !finish)
	{
		afm.rk4();
		afm.friction();

		afm.storevals();
		afm.inctime();
		double fric = afm.getfrc();
		double time = afm.gettime();
		
		uint slips = 0;
		double diff = fric - lastfric;

		if (!climb && diff > 0)				// are we in stick or slip state?
			climb = true;
		
		if (diff < 0 && climb)				// if slipped
		{
			fricslope.push_back(lastfric);	// push frica nd time
			timeslope.push_back(lasttime);
			
			climb = false;
			slips++;

			if (fricslope.size() > nflat)	// check if flat
			{
				double lastfric = fricslope[slips];
				double seclastfric = fricslope[slips-1];
				double flatness = lastfric - seclastfric;
		
				if (flatness < flattol*seclastfric)
					finish = true;
			}
		}
		
		if (loops == maxsteps || fric > maxfric || afm.getposx() > maxpos || fric < 0)		// stop loop if these happen
			termin = true;

		lastfric = fric;
		lasttime = time;
	}

	if (termin)			// if termination qualified discard and move on
	{
		slopes.push_back(0);
		continue;
	}

	else if (finish)	// if finish without errors
	{
		double diff = fricslope[1] - fricslope[0];
		
		if (diff < strengtol*fricslope[0])		// not enough initial strengthening
		{
			slopes.push_back(0);
			continue;
		}

		ostringstream namestream;
		namestream << "run" << l1 << l2 << ".csv";
		//namestream << "run" << l1 << l2 << l3 << l4 << ".csv";
			
		ofstream fsfric;
		fsfric.open(namestream.str());

		fsfric << "// Parameters used:";
		for (auto& el : varlist)
			fsfric << " " << el.second.name << ": " << el.second.val;
		fsfric << endl;
			
		for (uint i = 0; i < fricslope.size(); i++)
			fsfric << timeslope[i] << "," << fricslope[i]  << endl;

		fsfric.close();

		vector <double> fitptsx;
		vector <double> fitptsy;
		fitptsx.push_back(fricslope[0]);
		fitptsy.push_back(timeslope[0]);
		
		uint m = 1;
		
		while (diff > fricslope[0]*strengtol)
		{
			fitptsx.push_back(fricslope[m]);
			fitptsy.push_back(timeslope[m]);

			m++;
			diff = fricslope[m+1] - fricslope[m];
		}
		
		double slope = linreg(&fitptsx,&fitptsy);
		
		slopes.push_back(slope);
	}
	
	else
	{
		cout << "something went terribly wrong, you shouldn't be here, run!" << endl;
		return 0;
	}

	//adhesionsub.val += adhesionsub.inc;	
	//}
	//coupling.val += coupling.inc;	
	//}
	stiffness.val += stiffness.inc;	
	}
	adhesiontip.val += adhesiontip.inc;	
	}

	for (auto& el : slopes)
		fsslope << el << endl;

	fsslope.close();

}






















