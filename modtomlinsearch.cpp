
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

int pingcount = 0; 
void ping(const int line) {cout << "ping at line " << line << endl; pingcount++;}
void getpingcount() {cout << pingcount << endl;}

double zero = 1e-12;

// least square method
double linreg(vector <double>* x, vector <double>* y)
{   
	if (x->size() != y->size())
	{
		cout << "linear regression failure, incompatible vector sizes" << endl;

		return 0;
	}
	else if (x->size() == 0 || y->size() == 0)
	{
		cout << "linear regression failure, vector size is 0" << endl;
	}

	uint n = x->size();

	if (n == 1)
	{   
		double slope = y->at(0);

		return slope;
	}   

	double xavg = accumulate(x->begin(),x->end(),0.0)/x->size();
	double yavg = accumulate(y->begin(),y->end(),0.0)/y->size();

	double a = 0;
	double b = 0;
	for (uint k = 0; k < n; k++)
	{   
		a += (y->at(k) - yavg) * (x->at(k) - xavg);
		b += pow((x->at(k) - xavg),2);
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
	double qposdamp = 2e21;
	double qveldamp = 3e12;

	double timestep = 1.0e-13;
	uint timesteps = 5.0e4;

	string tfile = "";
	string xfile = "";
	string qfile = "";
	string tomfile = "";

	// simulation parameters

	double maxfric = 1.0e-6;
	double maxpos = 1.0e-6;
	double flattol = 0.1;
	double slopetol = 0.75;
	uint periods = static_cast <uint> (timestep*timesteps*supvel/latconst);
	
	// making grids
	
	double gridsize = 10;

	double adtipbot = 0.1;
	double adtiptop = 10;
	double stiffbot = 0.1;
	double stifftop = 10;
	double couplbot = 0.1;
	double coupltop = 10;
	double adsubbot = 0.1;
	double adsubtop = 10;

	vector <double> vars = {adtip0,stiff0,coupl0,adsub0};
	vector <double> varsbot = {adtipbot,stiffbot,couplbot,adsubbot};
	vector <double> varstop = {adtiptop,stifftop,coupltop,adsubtop};

	vector <vector <double>> grids;

	for (uint k = 0; k < vars.size(); k++)
	{
		vector <double> grid;

		double logmin = log(vars[k]*varsbot[k]);
		double logmax = log(vars[k]*varstop[k]);
		double shiftmax = logmax - logmin;
		double step = shiftmax / (double) gridsize;

		for (uint l = 0; l <= gridsize; l++)
		{
			double knot = exp(step*l + logmin);
			grid.push_back(knot);
		}

		grids.push_back(grid);
	}

	struct variable 
	{
		string name;
		double base;
		double bot;
		double top;
		vector <double> grid;
	};

	variable adhesiontip {"adhesiontip",vars[0],varsbot[0],varstop[0],grids[0]};
	variable stiffness {"stiffness",vars[1],varsbot[1],varstop[1],grids[1]};
	variable coupling {"coupling",vars[2],varsbot[2],varstop[2],grids[2]};
	variable adhesionsub {"adhesionsub",vars[3],varsbot[3],varstop[3],grids[3]};

	vector <variable*> varlist = {&adhesiontip,&stiffness,&coupling,&adhesionsub};

	//variable adhesiontip = {"adhesiontip",adtip0,step*adtip0,adtip0/offset};
	//variable stiffness = {"stiffness",stiff0,step*stiff0,stiff0/offset};
	//variable coupling = {"coupling",coupl0,step*coupl0,coupl0/offset};
	//variable adhesionsub = {"adhesionsub",adsub0,step*adsub0,adsub0/offset};	
	
	// this could also have been done with only structs, accessed by varlist::stiffness::val
	//map <string,variable> varlist;		
	//varlist.insert(make_pair("adtip",adhesiontip));
	//varlist.insert(make_pair("stiff",stiffness));
	//varlist.insert(make_pair("coupl",coupling));
	//varlist.insert(make_pair("adsub",adhesionsub));

	// Remember to update below when changing the number of variables. Ugly, but best I can do atm	

	ofstream fsslope;
	fsslope.open("slopes.csv");

	fsslope << "// Some comments to help parsing this file. Parameters used, with name, standard value, lower bound, upper bound:";
	
	for (auto& el : varlist)
		fsslope << " " << el->name << " " << el->base << "," << el->bot << "," << el->top << ".";
	
	fsslope << " The number of grid points is " << gridsize << ", they are logarithmically distributed from the lower to upper search bound, (see grid files). The first column of numbers given in this file contains the slopes of the strengthening, the remaining columnas are diagnistic and contain the maximum friction, as well as the amplitude of the final stick slip period. In the event that the algorithm terminated before fricition was found, error codes are supplied, these are: 101 t > maxtime, 102 fric > maxfric, 103 x > maxpos x, 104 fric < 0, 105 f''(x) = 0";

	fsslope << endl;
	
	vector <vector <double>> slopes;
	slopes.reserve(2*periods);
	
	struct point
	{
		double x;
		double y;
	};

	for (auto& adtip : adhesiontip.grid)
	{
	uint l1 = &adtip - &adhesiontip.grid[0];

	for (auto& stiff : stiffness.grid)
	{
	uint l2 = &stiff - &stiffness.grid[0];
	
	for (auto& coupl : coupling.grid)
	{
	uint l3 = &coupl - &coupling.grid[0];
	
	for (auto& adsub : adhesionsub.grid)
	{
	uint l4 = &adsub - &adhesionsub.grid[0];

	cout << "this is loop: " << l1 << "," << l2 << "," << l3 << "," << l4 << endl;

	bool termin = false;
	bool finish = false;
	bool climb = true;
	bool fail1 = false;
	bool fail2 = false;
	bool fail3 = false;
	bool fail4 = false;
	bool fail5 = false;

	double lastfric = 0;
	double lasttime = 0;
	double lastslope = 0;
	double loops = 0;
	double slope = 0;
	double maxslope = 0;
	
	vector <point> toppts;
	vector <point> botpts;
	vector <double> intslopes;

	toppts.reserve(2*periods); 
	botpts.reserve(2*periods);
	intslopes.reserve(2*periods);
	
	tomlin afm(spring,supvel,latconst,adtip,stiff,coupl,adsub,
			   timestep,timesteps,xmass,qmass,qposdamp,xveldamp,
			   qveldamp,xfile,qfile,tomfile,tfile);

	
	uint slips = 0;
	// Find maxima
	while (!finish && !termin)
	{
		afm.rk4();
		afm.friction();

		double fric = afm.getfrc();
		double time = afm.gettime();
		
		double diff = fric - lastfric;

		if (!climb && diff > 0)				// are we in stick or slip state?
		{
			point pt {lasttime,lastfric};
			botpts.push_back(pt);

			climb = true;
		}
		
		if (diff < 0 && climb)				// if slipped
		{
			point pt {lasttime,lastfric};
			toppts.push_back(pt);
			
			climb = false;
			slips++;
			
			if (toppts.size() > 1)	
			{
				double p1y = toppts.at(slips-2).y;
				double p2y = toppts.at(slips-1).y;
				double p1x = toppts.at(slips-2).x;
				double p2x = toppts.at(slips-1).x;
				
				slope = (p2y-p1y) / (p2x-p1x);
				intslopes.push_back(slope);
	
				// These are iffy (haha), but should work. A non-strengthening case should either be caught by the first conditional, and not the other, since maxslope will be about constant
				if (abs(slope - lastslope) < zero)
				{
					fail5 = true;
					termin = true;
				}
				if (slope > maxslope)
					maxslope = slope;
				else if (slope < maxslope * flattol && p2y-p1y > 0)
					finish = true;
			}
		}
		
		//if (loops == timesteps || fric > maxfric || afm.getposx() > maxpos || fric < 0)		// stop loop if these happen
		//{
		if (loops == timesteps)
		{
			fail1 = true;
			termin = true;
		}
			//cout << "maxloops: " << loops << endl;
		else if (fric > maxfric)
		{
			fail2 = true;
			termin = true;
		}
			//cout << "maxfric: " << fric << endl;
		else if (afm.getposx() > maxpos)
		{
			fail3 = true;
			termin = true;
		}
			//cout << "maxpos: " << afm.getposx() << endl;
		else if (fric < 0)
		{
			fail4 = true;
			termin = true;
		}
			//cout << "negfric: " << fric << endl;

		//	termin = true;
		//}

		afm.inctime();
		lastfric = fric;
		lasttime = time;
		lastslope = slope;
		loops++;
	}
		
	//ostringstream namestream;
	//namestream << "run" << l1 << "," << l2 << ".csv";
	////namestream << "run" << l1 << l2 << l3 << l4 << ".csv";
	//	
	//ofstream fsfric;
	//fsfric.open(namestream.str());

	//fsfric << "Parameters used, ";
	////for (auto& el : varlist)
	////	fsfric << " " << el.first.name  << ": " << el.second.val;
	////fsfric << endl;
	//fsfric << " " << adhesiontip.name << ": " << adtip << " and" ;
	//fsfric << " " << stiffness.name << ": " << stiff << ". Extrema follows below" << endl;

	//uint pts = max(botpts.size(),toppts.size());
	//uint ittop = 0;
	//uint itbot = 0;

	//for (uint i = 0; i < pts; i++)
	//{
	//	if (ittop != toppts.size())
	//	{
	//		fsfric << toppts.at(ittop).x << "," << toppts.at(ittop).y  << endl;
	//		ittop++;
	//	}
	//	if (itbot != botpts.size())
	//	{
	//		fsfric << botpts.at(itbot).x << "," << botpts.at(itbot).y  << endl;
	//		itbot++;
	//	}
	//}

	//fsfric.close();

	if (termin)			// if termination qualified discard and move on
	{
		cout << "loop terminated" << endl;
		
		//cout << "loop terminated" << endl;
		if (fail1)
			slopes.push_back({101,0,0});
		else if (fail2)
			slopes.push_back({102,0,0});
		else if (fail3)
			slopes.push_back({103,0,0});
		else if (fail4)
			slopes.push_back({104,0,0});
		else if (fail5)
			slopes.push_back({105,0,0});
		
		continue;
	}

	else if (finish)	// if finish without errors
	{
		cout << "loop finished" << endl;

		ostringstream namestream;
		namestream << "run(" << l1 << "," << l2 << "," << l3 << "," << l4 << ").csv";
		
		vector <double> vals = {adtip,stiff,coupl,adsub};

		ofstream fsfric;
		fsfric.open(namestream.str());

		fsfric << "// Parameters used, ";
		
		for (auto& el : varlist)
		{
			uint k = &el - &varlist[0];
			fsfric << " " << el->name  << ": " << vals.at(k);
		}

		fsfric << ". Extrema follows below." << endl;

		uint pts = max(botpts.size(),toppts.size());
		uint ittop = 0;
		uint itbot = 0;

		for (uint i = 0; i < pts; i++)
		{
			if (ittop != toppts.size())
			{
				fsfric << toppts.at(ittop).x << "," << toppts.at(ittop).y  << endl;
				ittop++;
			}
			if (itbot != botpts.size())
			{
				fsfric << botpts.at(itbot).x << "," << botpts.at(itbot).y  << endl;
				itbot++;
			}
		}

		fsfric.close();

		uint minpos = 0;
		uint maxpos = 0;
		
		cout << "max slope: " << maxslope << endl;
		cout << "slope tol: " << slopetol*maxslope << endl;
		cout << "in internal slopes " << endl;

		for(auto& el : intslopes)
		{
			cout << el << " ";
			uint k = &el - &intslopes[0];
			
			if (el > slopetol*maxslope)
			{
				if (minpos == 0)
					minpos = k;
				else
					maxpos = k;
			}
		}
		cout << endl;

		vector <double> fitptsx;
		vector <double> fitptsy;
		
		for(uint k = minpos; k < maxpos+2; k++)		// why +2? probably due to discarding of points in the costruction loop
		{
			fitptsx.push_back(toppts.at(k).x);
			fitptsy.push_back(toppts.at(k).y);
		}

		double slope = linreg(&fitptsx,&fitptsy);
		
		cout << "used in fitting:" << endl;
		for (auto& el : fitptsx)
			cout << el << " ";
		cout << endl;

		for (auto& el : fitptsy)
			cout << el << " ";
		cout << endl;

		cout << "slope:" << endl << slope << endl;

		slopes.push_back( {slope,toppts.back().y,toppts.back().y-botpts.back().y} );
	}
	else
	{
		cout << "something went terribly wrong, you shouldn't be here, run!" << endl;
		return 0;
	}

	}
	}
	}
	}

	for (auto& el : slopes)
		fsslope << el[0] << "," << el[1] << "," << el[2] << endl;

	fsslope.close();

	ofstream gstream;
	gstream.open("grids.csv");

	gstream << "Order: ";
	for (auto& el : varlist)
		gstream << " " << el->name << "," ;
	gstream << "." << endl;

	for (auto& grid : grids)
	{
		for (auto& el : grid)
		{
			gstream << el << "," ;
		}
		gstream << endl;;
	}
}




