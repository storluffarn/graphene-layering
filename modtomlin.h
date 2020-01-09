
//
// header file for modified tomlinson model
//


// includes

#include <fstream>
#include <random>
#include <numeric> 
#include <cmath>
#include <iostream>
#include <algorithm> 
#include "common_stuff.h"
//#include <algorithm>

// Lazy stuff

using namespace std;

typedef unsigned int uint;

// constants

//static const double pi = atan(1.0)*4;
//static const double kB = 1.38064852e-23;

// objects

class tomlin
{
	// sub classes
	class xobj
	{
		// dynamic variables
		double pos = 0;
		double vel = 0;
		double acc = 0;
		double kin = 0;

		// parameters
		const double mass;
		const double damp;

		// derived constants
		const double rmass = 1.0/mass;

		// containers
		vector <double> poss;
		vector <double> vels;
		vector <double> accs;

		// functions
		void calckin(){kin = 0.5*mass*pow(vel,2);}

		friend tomlin;

		xobj (double mass, double damp, uint tsteps)
			: mass(mass), damp(damp)
		{
			poss.reserve(tsteps);
			vels.reserve(tsteps);
			accs.reserve(tsteps);
		}
	};

	class qobj
	{
		// dynamic variables
		double pos = 0;
		double vel = 0;
		double acc = 0;
		double kin = 0;

		// parameters
		const double mass;
		const double damp;

		// derived constants
		const double rmass = 1.0/mass;

		// containers
		vector <double> poss;
		vector <double> vels;
		vector <double> accs;
		vector <double> slips;

		// functions
		void calckin(){kin = 0.5*mass*pow(vel,2);}

		friend tomlin;

		qobj (double mass, double damp, uint tsteps)
			: mass(mass), damp(damp)
		{
			poss.reserve(tsteps);
			vels.reserve(tsteps);
			accs.reserve(tsteps);
		}
	};

	// model parameters
	const double spring;
	const double supvel;
	const double latcon;
	const double align;
	const double barr1;	
	const double barr2;	
	const double kappa1;
	const double kappa2;
	const double nu2;
	const double nu4;

	// dyanamic varaibles
	double kin = 0;
	double pot = 0;
	double fric = 0;
	double time = 0;
	double suppos = 0;
	double temp;

	// non physical parameters
	const double tstep;
	const uint tsteps;
	bool reverse = false;
	bool pause = false;

	// derived constants
	const double latcona = latcon;
	const double latconb = latcona/align;
	const double rlatcona = 1.0/latcon;
	const double rlatcona2pi = 2*pi*rlatcona;
	const double rlatconb2pi = 2*pi*latcon/align;

	// objects
	xobj x;
	qobj q;

	// containers
	vector <double> pots;
	vector <double> kins;
	vector <double> frics;
	vector <double> times;
	vector <double> poss;	
	vector <double> lessnoiset;
	vector <double> lessnoisef;
	vector <double> lessnoisesuppos;

	// functions
	void xacc();
	void qacc();
	void testxacc();
	void testqacc();
	pair<double,double> langevin();
	pair<double,double> testlangevin();

	// files
	const string tfile;
	const string xfile;
	const string qfile;
	const string tomfile;
	const string avgfile;
	const string pfile;

	// need for speed
	const double oneosix = 1.0/6.0;

	public:

	// functions
	void rk4();
	void eulmar();
	void calcpot();
	void calcfric();
	void inctime();

	void noisered(uint,uint);

	// in line functions
	void calckin(){x.calckin(); q.calckin(); kin = x.kin + q.kin;}
	void treverse(){reverse = !reverse;}
	void tpause(){pause = !pause;}

	// accessors
	void setposx(double c){x.pos = c;}
	void setvelx(double c){x.vel = c;}
	void setposq(double c){q.pos = c;}
	void setvelq(double c){q.vel = c;}
	void settime(double c){time = c;}
	void setsuppos(double c){suppos = c;}

	double getposx(uint t){return x.poss[t];}
	double getvelx(){return x.vel;}
	double getaccx(){return x.acc;}
	double getposq(){return q.pos;}
	double getvelq(){return q.vel;}
	double getaccq(){return q.acc;}
	uint getposxsize(){return x.poss.size();}

	double gettime(double t){return times[t];}	
	double getkin(){return kin;}
	double getpot(){return pot;}
	double getfric(double t){return frics[t];}
	double getsuppos(){return suppos;} 
	double getrnfric(double t) {return lessnoisef[t];}
	double getrntime(double t) {return lessnoiset[t];}
	double getrnsuppos(double t) {return lessnoisesuppos[t];}

	vector <double>* getposs() {return &x.poss;}
	vector <double>* getfrics() {return &frics;}
	vector <double>* gettimes() {return &times;}
	vector <double>* getrnfrics() {return &lessnoisef;}
	vector <double>* getrntimes() {return &lessnoiset;}
	vector <double>* getrnsupposs() {return &lessnoisesuppos;}

	// io stuff
	void pushvals();
	void writedata();

	// debugging
	void printins();
	void printouts();
	void printavgs();
	void justkicks();

	// constructor
	tomlin(double spring, double supvel, double latcon, double align, double barr1, 
		   double barr2, double kappa1, double kappa2, double nu2, double nu4, 
		   double temp, double tstep, uint tsteps, double xmass, double qmass, 
		   double xdamp, double qdamp, string tfile, string xfile, string qfile, 
		   string tomfile, string avgfile, string pfile)
		 : spring(spring), supvel(supvel), latcon(latcon), align(align), barr1(barr1),
	       barr2(barr2), kappa1(kappa1), kappa2(kappa2), nu2(nu2), nu4(nu4), temp(temp),
		   tstep(tstep), tsteps(tsteps), x(xmass,xdamp,tsteps), 
		   q(qmass,qdamp,tsteps), tfile(tfile), xfile(xfile), 
		   qfile(qfile), tomfile(tomfile), avgfile(avgfile), pfile(pfile)
		{
		   kins.reserve(tsteps);
		   pots.reserve(tsteps);
		   frics.reserve(tsteps);
		   times.reserve(tsteps);
		   poss.reserve(tsteps);
		}
};

// out of line function members

void tomlin::calcpot()
{
	pot = 0.5*spring*pow(x.pos-suppos*time,2) + 
		  (barr1 + kappa1*pow(q.pos,2)) * (1-cos(rlatcona2pi*(x.pos-q.pos))) + 
		  (barr2 + kappa2*pow(q.pos,2)) * (1-cos(rlatcona2pi*align*x.pos)) + 
		  nu2*pow(q.pos,2) + nu4*pow(q.pos,4);
}
	
void tomlin::xacc()
{
	double force = 
		spring*(x.pos-suppos) + 
		rlatcona2pi*(barr1+kappa1*pow(q.pos,2)) * sin(rlatcona2pi*(x.pos-q.pos)) +
		rlatcona2pi*(barr2+kappa2*pow(q.pos,2)) * sin(rlatcona2pi*align*x.pos);
	
    x.acc = - (x.rmass*force + x.damp*x.vel);
	//x.acc = - (x.rmass*force); // use for T = 0 case, remembeer to change for q too! 

	//x.acc = - x.rmass*spring*(x.pos - suppos);
	//x.acc = 0;
}

//void tomlin::testxacc()
//{
//	double force = 
//		spring*(x.pos-suppos) + 
//		2*kappa1*(x.pos-q.pos) * ( 1-cos(rlatcona2pi * q.pos) ) + 
//		2*kappa2*(x.pos-q.pos) * ( 1-cos(rlatconb2pi * x.pos) ) + 
//		rlatconb2pi*sin(rlatconb2pi*x.pos) * (barr2 + kappa2*pow(x.pos-q.pos,2)) + 
//		2*nu2*(x.pos-q.pos) + 4*nu4*pow(x.pos-q.pos,3);
//}

void tomlin::qacc()
{
	double force = 
		-rlatcona2pi*(barr1 + kappa1*pow(q.pos,2)) * sin(rlatcona2pi*(x.pos-q.pos)) +
		2.0*kappa1*q.pos * (1.0-cos(rlatcona2pi*(x.pos-q.pos))) + 
		2.0*kappa2*q.pos * (1.0-cos(rlatcona2pi*align*x.pos)) + 
		2.0*nu2*q.pos + 4.0*nu4*pow(q.pos,3);
	
    q.acc = - (q.rmass*force + q.damp*q.vel);
    //q.acc = - (q.rmass*force);

	//q.acc = 0;
}

//void tomlin::testqacc()
//{
//	double force = 
//		2*kappa1*(x.pos-q.pos) * ( 1-cos(rlatcona2pi * q.pos) ) - 
//		2*kappa2*(x.pos-q.pos) * ( 1-cos(rlatconb2pi * x.pos) ) +
//		rlatcona2pi*sin(1-rlatcona2pi*x.pos) * (barr1 + kappa1*(x.pos-q.pos)) -
//		2*nu2*(x.pos-q.pos) - 4*nu4*pow(x.pos-q.pos,4);
//}

void tomlin::calcfric()
{
	fric = spring*(suppos-x.pos);
}

void tomlin::pushvals()
{
	x.poss.push_back(x.pos);
	x.vels.push_back(x.vel);
	x.accs.push_back(x.acc);
	q.poss.push_back(q.pos);
	q.vels.push_back(q.vel);
	q.accs.push_back(q.acc);

	kins.push_back(kin);
	pots.push_back(pot);
	frics.push_back(fric);
	times.push_back(time);
	poss.push_back(suppos);
}

void tomlin::inctime()
{
	time += tstep;
	
	if (reverse)
	{
		suppos -= supvel*tstep;
	}
	else if (pause)
	{
		suppos = 1.0*suppos;
	}
	else
	{
		suppos += supvel*tstep;
	}
}

void tomlin::writedata()
{
	ofstream xstream, qstream, tomstream, tstream, avgstream, pstream;

	xstream.open(xfile);
		xstream << "position, velocity, acceleration" << endl;
		for (uint k = 0; k < tsteps; k++)
			xstream << x.poss[k] << "," << x.vels[k] << "," << x.accs[k] << endl;
	xstream.close();
	
	qstream.open(qfile);
		qstream << "position, velocity, acceleration" << endl;
		for (uint k = 0; k < tsteps; k++)
			qstream << q.poss[k] << "," << q.vels[k] << "," << q.accs[k] << endl;
	qstream.close();
	
	tomstream.open(tomfile);
		tomstream << "kinetic, potential, friction" << endl;
		for (uint k = 0; k < tsteps; k++)
			tomstream << kins[k] << "," << pots[k] << "," << frics[k] << endl;
	tomstream.close();

	tstream.open(tfile);
		tstream << "time, displacement" << endl;
		for (uint k = 0; k < tsteps; k++)
			tstream << times[k] << "," << poss[k] << endl;
	tstream.close();

	avgstream.open(avgfile);
		avgstream << "time, avged friction " << endl;
		uint size = lessnoisef.size();
		for (uint k = 0; k < size; k++)
			avgstream << lessnoiset[k] << "," << lessnoisef[k] << "," << lessnoisesuppos[k] << endl;
	avgstream.close();
	
	pstream.open(pfile);
		pstream << "spring " << spring << endl << "supvel " << suppos << endl 
				<< "align " << align << endl   
				<< "latcona " << latcona << endl << "latconb "  << latconb  << endl 
				<< "barr1 "  << barr1  << endl << "barr2 "  << barr2  << endl 
				<< "kappa1 " << kappa1 << endl << "kappa2 " << kappa2 << endl 
				<< "nu2 "	  << nu2	<< endl << "nu4 "	 << nu4	   << endl 
				<< "tstep "  << tstep  << endl << "tsteps " << tsteps << endl 
				<< "xmass "  << x.mass  << endl << "qmass "  << q.mass  << endl 
				<< "xdamp "  << x.damp  << endl << "qdamp "  << q.damp  << endl;
	pstream.close();
} 

void tomlin::printins()
{
	cout << "the following parameters was used: " << endl
		 << "spring " << spring << endl << "supvel " << suppos << endl 
		 << "latcon " << latcon << endl << "align "  << align  << endl 
		 << "barr1 "  << barr1  << endl << "barr2 "  << barr2  << endl 
		 << "kappa1 " << kappa1 << endl << "kappa2 " << kappa2 << endl 
		 << "nu2 "	  << nu2	<< endl << "nu4 "	 << nu4	   << endl 
		 << "tstep "  << tstep  << endl << "tsteps " << tsteps << endl 
		 << "xmass "  << x.mass  << endl << "qmass "  << q.mass  << endl 
		 << "xdamp "  << x.damp  << endl << "qdamp "  << q.damp  << endl
		 << "tfile "  << tfile  << endl << "xfile "  << xfile  << endl
		 << "qfile "  << qfile  << endl << "tomfile " << tomfile << endl;
}

void tomlin::printouts()
{
	cout << "and the exit values were: " << endl
		 << "time " << time << endl << "suppos " << suppos << endl << "kin " 
		 << kin << endl		 << "pot " << pot << endl << "tot en " <<  kin + pot 
		 << endl << "xpos " << x.pos << endl << "xvel " << x.vel << endl 
		 << "xacc " << x.acc << endl << "qpos " << q.pos << endl << "qvel " 
		 << q.vel << endl << "qacc " << q.acc << endl;

}

void tomlin::printavgs()
{
	double avgkin = accumulate(kins.begin(),kins.end(),0.0) / kins.size();

	cout << "avg kin " << avgkin << endl;
	
	double avgfric = accumulate(frics.begin(),frics.end(),0.0) / frics.size();

	cout << "avg fric " << avgfric << endl;
}

void tomlin::rk4()	// two degree of freedom RK4 algorithm
{
	// preliminaries

	// container for old positions, velocities, and accelerations
	
	struct oldvals {double pos; double vel; double acc;};
	struct k {double vel; double acc;};
	oldvals oldx;
	oldvals oldq;
	double oldt;

	// calculate accelerations in current configuration
	xacc();
	qacc();
	
	auto kicks = langevin();	// random kicks from langevin dynamics, see below
	double xkick = kicks.first;
	double qkick = kicks.second;
	//double xkick = 0;			// used for debugging
	//double qkick = 0;

	// k1

	// save old values, at this point this is the ones from the beginning of 
	// the timestep

	oldt = time;
	oldx.pos = x.pos;
   	oldx.vel = x.vel;
	oldx.acc = x.acc;
	oldq.pos = q.pos;
	oldq.vel = q.vel;
	oldq.acc = q.acc;
	
	// calculate the RK4 k1 constants, as: 
	// * velocity is the start of timestep velocity since we're in k1, so no time 
	//   has elapsed
	// * acceleration is the acceleration as calculated from the start of timestep 
	//   configuration langevin thermostated since no time has elapsed

	k k1x = {x.vel,  x.acc - x.vel * x.damp + xkick};
	k k1q = {q.vel,  q.acc - q.vel * q.damp + qkick};
	
	// k2
	
	// the RK4 k2 procedure is very similar to that of k1. but for k2 time is advanced
	// by half a timestep. this means we euler forward the dynamic variables half a 
	// timestep (advancing from the start of timestep configuration and using k1 as 
	// the present configuration) before we compute the new acceleration from the 
	// obtained half timestep configuration, from with we then calculate the k2 
	// constants identically to how we did for k1

	time = oldt + 0.5*tstep;
	x.pos = oldx.pos+0.5*tstep*k1x.vel;
	x.vel = oldx.vel+0.5*tstep*k1x.acc;
	q.pos = oldq.pos+0.5*tstep*k1q.vel;
	q.vel = oldq.vel+0.5*tstep*k1q.acc;
	
	xacc();
	qacc();

	k k2x = {x.vel,x.acc - x.vel * x.damp + xkick};
	k k2q = {q.vel,q.acc - q.vel * q.damp + qkick};
	
	// k3
	
	// the RK4 k3 prodecure is identical to that of k2, just the k2 configuration is 
	// used as the present configuration in the euler step

	time = oldt + 0.5*tstep;
	x.pos = oldx.pos+0.5*tstep*k2x.vel;
	x.vel = oldx.vel+0.5*tstep*k2x.acc;
	q.pos = oldq.pos+0.5*tstep*k2q.vel;
	q.vel = oldq.vel+0.5*tstep*k2q.acc;
	
	xacc();
	qacc();
	
	k k3x = {x.vel,x.acc - x.vel * x.damp + xkick};
	k k3q = {q.vel,q.acc - q.vel * q.damp + qkick};
	
	// k4
	
	// the RK4 k4 prodecure is identical to that of k3, but the present configuration 
	// used here is one full timestep from the start of timestep configuration

	time = oldt + tstep;
	x.pos = oldx.pos+tstep*k3x.vel;
	x.vel = oldx.vel+tstep*k3x.acc;
	q.pos = oldq.pos+tstep*k3q.vel;
	q.vel = oldq.vel+tstep*k3q.acc;
	
	xacc();
	qacc();
	
	k k4x = {x.vel,x.acc - x.vel * x.damp + xkick};
	k k4q = {q.vel,q.acc - q.vel * q.damp + qkick};

	// obtained approximation -- just sum everything up according to RK4 formula

	x.pos = oldx.pos + tstep*oneosix*(k1x.vel + 2*k2x.vel + 2*k3x.vel + k4x.vel);
	x.vel = oldx.vel + tstep*oneosix*(k1x.acc + 2*k2x.acc + 2*k3x.acc + k4x.acc);
	q.pos = oldq.pos + tstep*oneosix*(k1q.vel + 2*k2q.vel + 2*k3q.vel + k4q.vel);
	q.vel = oldq.vel + tstep*oneosix*(k1q.acc + 2*k2q.acc + 2*k3q.acc + k4q.acc);

    //x.acc = xkick;
	// reset the timestep (time is ticked at the end of the main loop)
	
	time = oldt;
}

pair <double, double> tomlin::langevin()
{
	random_device rd;
	mt19937 gen(rd());

	double mean = 0;
	double xstd = sqrt(2*x.mass*kB*temp*x.damp/tstep);
	double qstd = sqrt(2*q.mass*kB*temp*q.damp/tstep);

	normal_distribution<double> xkick(mean,xstd);
	normal_distribution<double> qkick(mean,qstd);

	pair <double, double> out = {xkick(gen)*x.rmass,qkick(gen)*q.rmass};

	return out;
}

//pair <double, double> tomlin::testlangevin()
//{
//	random_device rd;
//	mt19937 gen(rd());
//
//	double mean = 0;
//	double xstd = sqrt(2*x.mass*kB*temp*x.damp/tstep);
//	double qstd = sqrt(2*0.5*(x.mass+q.mass)*kB*temp*0.5*(x.damp+q.damp)/tstep);
//
//	normal_distribution<double> xkick(mean,xstd);
//	normal_distribution<double> qkick(mean,qstd);
//
//	pair <double, double> out = {xkick(gen)*x.rmass,qkick(gen)*q.rmass};
//
//	return out;
//}

void tomlin::noisered(uint halfmeansize, uint skip)
{
	// building mean for first point

	double sum = 0;
	double avg = 0;

	uint bgnel = 0;
	uint midel = 0;
	uint endel = 0;

	const uint its = floor((double) halfmeansize/skip);
	const double rmeansize = 1.0/(floor(((double) 2*halfmeansize)/skip)+1);
	
	lessnoiset.reserve(ceil(tsteps/skip));
	lessnoisef.reserve(ceil(tsteps/skip));
	lessnoisesuppos.reserve(ceil(tsteps/skip));

	for (uint l = 0; l < its; l++)		
	{
		sum += frics[endel];
		endel += skip;
	}
	
	midel = endel;
	sum += frics[endel];
	endel += skip;
	
	for (uint l = 0; l < its; l++)
	{
		sum += frics[endel];
		endel += skip;
	}
		
	avg = rmeansize * sum;

	lessnoiset.push_back(times[midel]);
	lessnoisef.push_back(avg);
	lessnoisesuppos.push_back(poss[midel]);

	// calculating rest of means

	while (endel < tsteps)
	{
		avg -= rmeansize*frics[bgnel];
		bgnel += skip;
		
		midel += skip;
		
		avg += rmeansize*frics[endel];
		endel += skip;
		
		lessnoiset.push_back(times[midel]);
		lessnoisef.push_back(avg);
		lessnoisesuppos.push_back(poss[midel]);
	}
}


// ---------- DEBUGGING STARTS HERE --------------

void tomlin::eulmar()		// simple euler for comparison (does not work)
{
	struct oldvals {double pos; double vel; double acc;};
	
	oldvals oldx;
	oldx.pos = x.pos;
	oldx.vel = x.vel;
	oldx.acc = x.acc;

	auto kicks = langevin();
	double xkick = kicks.first;
	//double xkick = 0;
	
	xacc();

	x.vel = oldx.vel + x.acc * tstep - x.damp * oldx.vel + xkick;
	x.pos = oldx.pos + x.vel * tstep;
}

void tomlin::justkicks()
{
	auto kicks = langevin();
	double xkick = kicks.first;

	x.acc = xkick;	
}



