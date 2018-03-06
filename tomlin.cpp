//
//
// Toy program for trying out the Tomlinson model

// includes
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

// global variables
static double pi = atan(1)*4;
static double oos = 1.0/6.0;

// classes

class tomlin
{
	public:
	double pos = 0;
	double vel = 0;
	double acc = 0;
	double kin = 0;
	double pot = 0;
	double mass;
	double resmass;
	double timestep;
	double spring;
	double latconst;
	double reslatconst;
	double latpot;
	double etafric;
	double supvel;
	double time = 0;
	double layermass;
	double interlatpot;

	double calcacc(double,double);
	double calcacc2(double,double);
	void updateTime(){time += timestep;}
	
	// function members
	public:
	
	void rk4();
	void kinetic();
	void potential(double);
	void potential2(double);
	double latfric();
	double lambda(){return latpot*4*pi*pi/(spring*latconst*latconst);}
	double eta(){return etafric*sqrt(mass/spring);}
	
	void setpos(double a){pos = a;}
	void setvel(double a){vel = a;}
	double getpos(){return pos;}
	double getvel(){return vel;}
	double getacc(){return acc;}
	double getkin(){return kin;}
	double getpot(){return pot;}

	// constructor
	tomlin(double x1, double x2, double x3, double x4, double x5, double x6, double x7, double x8, double x9, double x10)
		: pos(x1), mass(x2), timestep(x3), spring(x4), latconst(x5), latpot(x6), etafric(x7), supvel(x8), layermass(x9), interlatpot(x10)
		{
			reslatconst = 1.0/latconst;
			resmass = 1.0/mass;
		}
};
	
double tomlin::calcacc(double pos, double vel)
{
	return -spring*(pos-supvel*time)*resmass - 2*pi*latpot*reslatconst*sin(2*pi*pos*reslatconst)*resmass - etafric*vel;
}

void tomlin::potential(double pos)
{
	pot = 0.5*spring*pow(pos-supvel*time,2) + latpot*(1-cos(2*pi*reslatconst*pos));
}

double tomlin::calcacc2(double pos, double vel)
{
	return -spring*(pos-supvel*time)*resmass - 2*pi*latpot*reslatconst*sin(2*pi*pos*reslatconst)*resmass - etafric*vel - pow(pos,0.5); // * interlatpot*2*pi*reslatconst*sin(2*pi*reslatconst*pos);
}

void tomlin::potential2(double pos)
{
	pot = 0.5*spring*pow(pos-supvel*time,2) + latpot*(1-cos(2*pi*reslatconst*pos)) + pow(pos,0.5); // interlatpot*cos(2*pi*reslatconst*pos);
}

void tomlin::kinetic()
{
	kin = 0.5*mass*vel*vel;
}

double tomlin::latfric()
{
	return spring*(supvel*time-pos);
}

void tomlin::rk4()
{
	struct vparticle {double pos; double vel; double acc;};
	struct k {double vel; double acc;};
	vparticle vp;
	
	acc = calcacc2(pos,vel);
	
	vp.pos = pos;
   	vp.vel = vel;
	vp.acc = acc;
	
	k k1 = {vp.vel,vp.acc};
	
	vp.pos = pos+0.5*timestep*k1.vel;
	vp.vel = vel+0.5*timestep*k1.acc;
	vp.acc = calcacc2(vp.pos,vp.vel);

	k k2 = {vp.vel,vp.acc};
	
	vp.pos = pos+0.5*timestep*k2.vel;
	vp.vel = vel+0.5*timestep*k2.acc;
	vp.acc = calcacc2(vp.pos,vp.vel);

	k k3 = {vp.vel,vp.acc};

	vp.pos = pos+timestep*k3.vel;
	vp.vel = vel+timestep*k3.acc;
	vp.acc = calcacc2(vp.pos,vp.vel);

	k k4 = {vp.vel,vp.acc};

	pos += timestep*oos * (k1.vel + 2*k2.vel + 2*k3.vel + k4.vel);
	vel += timestep*oos * (k1.acc + 2*k2.acc + 2*k3.acc + k4.acc);
	
	updateTime();
}

int main()
{
	// particle parameters	
	
	double timestep = 1e-13;
	double mass = 10e-25;
	double initpos = 0e-9;
	double latpot = 0.5e-19;
	double etafric = 3.5e12;
	double latconst = 1e-9;
	double spring = 1;
	double supvel = 1;
	double layermass = 1e-23;
	double interlatpot = 1e-19;

	tomlin tompart(initpos, mass, timestep, spring, latconst, latpot, etafric, supvel, layermass, interlatpot);

	//tompart.setpos(1.5e-9);
	//tompart.setvel(1e-8);
	//tompart.potential(tompart.getpos());
	//cout << tompart.getpot() << endl;
	//cout << tompart.calcacc(tompart.getpos(),tompart.getvel()) << endl;

	// simulation parameters
	int duration = 20000;

	ofstream fstomlin; fstomlin.open("tomlin.txt");

	for ( int k = 0; k < duration; k++)
	{
		tompart.rk4();
		
		tompart.kinetic();
		tompart.potential2(tompart.getpos());
		
		//cout << "getpos " << tompart.getpos() << endl;
		//cout << "getvel " << tompart.getvel() << endl;
		//cout << "getacc " << tompart.getacc() << endl;
		//cout << "getkin " << tompart.getkin() << endl;
		//cout << "getpot " << tompart.getpot() << endl;
		//cout << "lafric " << tompart.latfric() << endl;
		//cout << endl;

		fstomlin << tompart.getpos() << endl;
		fstomlin << tompart.getvel() << endl;
		fstomlin << tompart.getacc() << endl;
		fstomlin << tompart.getkin() << endl;
		fstomlin << tompart.getpot() << endl;
		fstomlin << tompart.latfric() << endl;
		fstomlin << tompart.time << endl;
	}
	fstomlin.close();

	cout << "eta used " << tompart.eta() << endl;
	cout << "lambda used " << tompart.lambda() << endl;
}

	


