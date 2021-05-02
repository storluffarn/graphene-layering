//
//
// Toy program for trying out the Tomlinson model
//
// nd for adimensional
//

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
	double pos;
	double vel;
	double acc;
	double kin;
	double pot;
	double mass;
	double resmass;
	double timestep;
	double spring;
	double latconst;
	double reslatconst;
	double latpot;
	double etafric;
	double supvel;
	double time;
	double tau;
	double restau;
	double eta;
	double lambda;
	
	void updateTime(){time += timestep;}

	// out-of-line functions
	void rk4();
	void kinetic();
	void potential(double);
	double latfric();
	double calcacc(double,double);

	// inline functions
	
	// accessors
	double getpos(){return pos;}
	double getvel(){return vel;}
	double getacc(){return acc;}
	double getkin(){return kin;}
	double getpot(){return pot;}
	double getlambda(){return lambda;}
	double geteta(){return eta;}
	void setpos(double a){pos = a;}
	void setvel(double a){vel = a;}
	void setlambda(double a){lambda = a;}
	void seteta(double a){eta = a;}
	
	// constructor
	tomlin(double x1, double x2, double x3, double x4, double x5, double x6, double x7, double x8)
		: mass(x2), spring(x4), latconst(x5), latpot(x6), etafric(x7), supvel(x8)
		{
			reslatconst = 1.0/latconst;
			resmass = 1.0/mass;
			pos = 2*pi*x1/latconst;
			tau = sqrt(mass/spring);
			restau = 1.0/tau;
			timestep = x3*restau;
			eta = etafric*tau;
			lambda = latpot*4*pi*pi/(spring*latconst*latconst);
		}
};
	
double tomlin::calcacc(double pos, double vel)
{
	return -(pos-supvel*time) - lambda*sin(pos) - eta*vel;
}

void tomlin::potential(double pos)
{
	pot = 0.5*pow(pos-supvel*time,2) + lambda*(1-cos(pos));
}

void tomlin::kinetic()
{
	kin = 0.5*mass*vel*vel;
}

double tomlin::latfric()
{
	return supvel*time-pos;
}

void tomlin::rk4()
{
	struct vparticle {double pos; double vel; double acc;};
	struct k {double vel; double acc;};
	vparticle vp;
	
	acc = calcacc(pos,vel);
	
	vp.pos = pos;
   	vp.vel = vel;
	vp.acc = acc;
	
	k k1 = {vp.vel,vp.acc};
	
	vp.pos = pos+0.5*timestep*k1.vel;
	vp.vel = vel+0.5*timestep*k1.acc;
	vp.acc = calcacc(vp.pos,vp.vel);

	k k2 = {vp.vel,vp.acc};
	
	vp.pos = pos+0.5*timestep*k2.vel;
	vp.vel = vel+0.5*timestep*k2.acc;
	vp.acc = calcacc(vp.pos,vp.vel);

	k k3 = {vp.vel,vp.acc};

	vp.pos = pos+timestep*k3.vel;
	vp.vel = vel+timestep*k3.acc;
	vp.acc = calcacc(vp.pos,vp.vel);

	k k4 = {vp.vel,vp.acc};

	pos += timestep*oos * (k1.vel + 2*k2.vel + 2*k3.vel + k4.vel);
	vel += timestep*oos * (k1.acc + 2*k2.acc + 2*k3.acc + k4.acc);
	
	updateTime();
}

int main()
{

	// particle parameters	
	//double timestep = 1e-15;
	//double mass = 1.e-26;
	//double initpos = 0e-9;
	//double latpot = 1e-18;
	//double friction = 1e-17;
	//double latconst = 1e-9;
	//double spring = 1e0;
	//double supvel = 1e-7;
	
	double timestep = 1e-15;
	double mass = 1e-25;
	double initpos = 0e-9;
	double latpot = 1e-19;
	double friction = 6e12;
	double latconst = 1e-9;
	double spring = 1;
	double supvel = 200*latconst/sqrt(mass/spring);

	tomlin tompart(initpos, mass, timestep, spring, latconst, latpot, friction, supvel);
	
	tompart.setlambda(3);
	tompart.seteta(2);
	//tompart.setpos(1.5e-9);
	//tompart.setvel(1e-8);
	//tompart.potential(tompart.getpos());
	//cout << tompart.getpot() << endl;
	//cout << tompart.calcacc(tompart.getpos(),tompart.getvel()) << endl;
	//cout << tompart.getlambda() << endl;
	//cout << tompart.geteta() << endl;

	// simulation parameters
	int duration = 10000;

	ofstream fstomlin; fstomlin.open("ndtomlin.txt");

	for ( int k = 1; k < duration; k++)
	{
		tompart.rk4();
		
		tompart.kinetic();
		tompart.potential(tompart.getpos());
		
		//cout << "getpos " << tompart.getpos() << endl;
		//cout << "getvel " << tompart.getvel() << endl;
		//cout << "getacc " << tompart.getacc() << endl;
		//cout << "getpot " << tompart.getpot() << endl;
		//cout << "getkin " << tompart.getkin() << endl;
		//cout << "lafric " << tompart.latfric() << endl;
		//cout << tompart.time << endl;
		//cout << endl;

		fstomlin << tompart.getpos() << endl;
		fstomlin << tompart.getvel() << endl;
		fstomlin << tompart.getacc() << endl;
		fstomlin << tompart.getkin() << endl;
		fstomlin << tompart.getpot() << endl;
		fstomlin << tompart.latfric() << endl;

		
	}
	
	fstomlin.close();
	cout << "eta used " << tompart.geteta() << endl;
	cout << "lambda used " << tompart.getlambda() << endl;
}

	


