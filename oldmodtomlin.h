
//
// header file for modified tomlinson model
//


// includes

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

// Lazy stuff

using namespace std;

typedef unsigned int uint;

// global stuff

static const double pi = atan(1)*4;
uint pingcount;

void ping(const int line) {cout << "ping at line " << line << endl; pingcount++;}

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
		const double veldamp;
		
		// derived constants
		const double massrec = 1.0/mass;

		// containers
		vector <double> poss;
		vector <double> vels;
		vector <double> accs;

		// functions
		void kinetic(){kin = 0.5*mass*pow(vel,2);}

		friend tomlin;

		xobj (double mass, double veldamp, uint timesteps)
			: mass(mass), veldamp(veldamp)
		{
			poss.reserve(timesteps);
			vels.reserve(timesteps);
			accs.reserve(timesteps);
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
		const double posdamp;
		const double veldamp;
		
		// derived constants
		const double massrec = 1.0/mass;
		
		// containers
		vector <double> poss;
		vector <double> vels;
		vector <double> accs;
		
		// functions
		void kinetic(){kin = 0.5*mass*pow(vel,2);}

		friend tomlin;
		
		qobj (double mass, double posdamp, double veldamp, uint timesteps)
			: mass(mass), posdamp(posdamp), veldamp(veldamp)
		{
			poss.reserve(timesteps);
			vels.reserve(timesteps);
			accs.reserve(timesteps);
		}
	};
	
	// dyanamic varaibles
	double kin = 0;
	double pot = 0;
	double frc = 0;
	double time = 0;

	// parameters
	const double spring;
	const double supvel;
	const double latconst;
	const double adhesiontip;						// c1
	const double stiffness;							// c2
	const double coupling;							// c3
	const double adhesionsub;						// c4
	
	// non physical parameters
	const double timestep;
	const uint timesteps;
	
	// derived constants
	const double latconrec = 1.0/latconst;
	const double latconrec2pi = 2*pi*latconrec;

	// objects
	xobj x;
	qobj q;

	// containers
	vector <double> pots;
	vector <double> kins;
	vector <double> frcs;
	vector <double> tims;

	// functions
	void xacc();
	void qacc();

	// files
	const string xfile;
	const string qfile;
	const string tomfile;
	const string tfile;
		
	// repo of speed
	const double oneoversix = 1.0/6.0;

	public:

	// functions
	void rk4();
	void potential();
	void friction();
	
	// in line functions
	void kinetic(){kin = x.kin + q.kin;}
	void inctime(){time += timestep;}

	// accessors
	void setposx(double c){x.pos = c;}
	void setvelx(double c){x.vel = c;}
	void setposq(double c){q.pos = c;}
	void setvelq(double c){q.vel = c;}

	void settime(double c){time = c;}

	double getposx(){return x.pos;}
	double getvelx(){return x.vel;}
	double getaccx(){return x.acc;}
	double getposq(){return q.pos;}
	double getvelq(){return q.vel;}
	double getaccq(){return q.acc;}

	double gettime(){return time;}	
	double getkin(){return kin;}
	double getpot(){return pot;}
	double getfrc(){return frc;}

	// io stuff
	void storevals();
	void printvals();

	// constructor
	tomlin(double spring, double supvel, double latconst, double adhesiontip, double stiffness,
		   double coupling, double adhesionsub, double timestep, uint timesteps, double xmass,
		   double qmass, double qposdamp, double xveldamp, double qveldamp, string xfile,
		   string qfile, string tomfile, string tfile)
		: spring(spring), supvel(supvel), latconst(latconst), adhesiontip(adhesiontip), 
		  stiffness(stiffness), coupling(coupling), adhesionsub(adhesionsub), timestep(timestep), 
		  timesteps(timesteps), x(xmass,xveldamp,timesteps), q(qmass,qposdamp,qveldamp,timesteps),
		  xfile(xfile), qfile(qfile), tomfile(tomfile), tfile(tfile)
		{
			kins.reserve(timesteps);
			pots.reserve(timesteps);
			frcs.reserve(timesteps);
			tims.reserve(timesteps);
		}
};

// out of line function members

void tomlin::potential()
{
	pot = 0.5*spring*pow(x.pos-supvel*time,2) + (adhesiontip + stiffness*pow(q.pos,2))*
		  (1-cos(latconrec2pi*(x.pos-coupling*q.pos))) + adhesionsub*pow(q.pos,4);
}
	
void tomlin::xacc()
{
	double force = spring*(x.pos-supvel*time) + 
				   latconrec2pi * (adhesiontip+stiffness*pow(q.pos,2)) * sin(latconrec2pi*(x.pos-coupling*q.pos));
	
	x.acc = - (x.massrec*force + x.veldamp*x.vel);
}

void tomlin::qacc()
{
	double force = 2.0*stiffness * (1.0-cos(latconrec2pi*(x.pos-coupling*q.pos)))*q.pos -
				   latconrec2pi*coupling*(adhesiontip+stiffness*pow(q.pos,2))*sin(latconrec2pi*(x.pos-coupling*q.pos)) +
				   4.0*adhesionsub*pow(q.pos,3);

	q.acc = - (q.massrec*force + q.posdamp*q.pos + q.veldamp*q.vel);
}

void tomlin::friction()
{
	frc = spring*(supvel*time-x.pos);
}

void tomlin::storevals()
{
	x.poss.push_back(x.pos);
	x.vels.push_back(x.vel);
	x.accs.push_back(x.acc);
	q.poss.push_back(q.pos);
	q.vels.push_back(q.vel);
	q.accs.push_back(q.acc);

	kins.push_back(kin);
	pots.push_back(pot);
	frcs.push_back(frc);
	tims.push_back(time);
}

void tomlin::printvals()
{
	ofstream xstream, qstream, tomstream, tstream;

	xstream.open(xfile);
		xstream << "position,velocity,acceleration" << endl;
		for (uint k = 0; k < timesteps; k++)
			xstream << x.poss[k] << "," << x.vels[k] << "," << x.accs[k] << endl;
	xstream.close();
	
	qstream.open(qfile);
		qstream << "position,velocity,acceleration" << endl;
		for (uint k = 0; k < timesteps; k++)
			qstream << q.poss[k] << "," << q.vels[k] << "," << q.accs[k] << endl;
	qstream.close();
	
	tomstream.open(tomfile);
		tomstream << "kinetic,potential,friction" << endl;
		for (uint k = 0; k < timesteps; k++)
			tomstream << kins[k] << "," << pots[k] << "," << frcs[k] << endl;
	tomstream.close();

	tstream.open(tfile);
		tstream << "time,displacement" << endl;
		for (uint k = 0; k < timesteps; k++)
			tstream << tims[k] << "," << supvel*tims[k] << endl;
	tstream.close();
}

void tomlin::rk4()
{
	// preliminaries

	struct oldvals {double pos; double vel; double acc;};
	struct k {double vel; double acc;};
	oldvals oldx;
	oldvals oldq;

	xacc();
	qacc();
	
	// k1

	oldx.pos = x.pos;
   	oldx.vel = x.vel;
	oldx.acc = x.acc;
	oldq.pos = q.pos;
	oldq.vel = q.vel;
	oldq.acc = q.acc;
	
	k k1x = {x.vel,x.acc};
	k k1q = {q.vel,q.acc};
	
	// k2

	x.pos = oldx.pos+0.5*timestep*k1x.vel;
	x.vel = oldx.vel+0.5*timestep*k1x.acc;
	q.pos = oldq.pos+0.5*timestep*k1q.vel;
	q.vel = oldq.vel+0.5*timestep*k1q.acc;
	
	xacc();
	qacc();

	k k2x = {x.vel,x.acc};
	k k2q = {q.vel,q.acc};
	
	// k3

	x.pos = oldx.pos+0.5*timestep*k2x.vel;
	x.vel = oldx.vel+0.5*timestep*k2x.acc;
	q.pos = oldq.pos+0.5*timestep*k2q.vel;
	q.vel = oldq.vel+0.5*timestep*k2q.acc;
	
	xacc();
	qacc();
	
	k k3x = {x.vel,x.acc};
	k k3q = {q.vel,q.acc};
	
	// k4

	x.pos = oldx.pos+timestep*k3x.vel;
	x.vel = oldx.vel+timestep*k3x.acc;
	q.pos = oldq.pos+timestep*k3q.vel;
	q.vel = oldq.vel+timestep*k3q.acc;
	
	xacc();
	qacc();
	
	k k4x = {x.vel,x.acc};
	k k4q = {q.vel,q.acc};

	// obtained approximation

	x.pos = oldx.pos + timestep*oneoversix*(k1x.vel + 2*k2x.vel + 2*k3x.vel + k4x.vel);
	x.vel = oldx.vel + timestep*oneoversix*(k1x.acc + 2*k2x.acc + 2*k3x.acc + k4x.acc);
	q.pos = oldq.pos + timestep*oneoversix*(k1q.vel + 2*k2q.vel + 2*k3q.vel + k4q.vel);
	q.vel = oldq.vel + timestep*oneoversix*(k1q.acc + 2*k2q.acc + 2*k3q.acc + k4q.acc);

	//cout << "k1x.vel: " << k1x.vel << " k1x.acc: " << k1x.acc << endl;
	//cout << "k1q.vel: " << k1q.vel << " k1q.acc: " << k1q.acc << endl;
	//
	//cout << "k2x.vel: " << k2x.vel << " k2x.acc: " << k2x.acc << endl;
	//cout << "k2q.vel: " << k2q.vel << " k2q.acc: " << k2q.acc << endl;
	//
	//cout << "k3x.vel: " << k3x.vel << " k3x.acc: " << k3x.acc << endl;
	//cout << "k3q.vel: " << k3q.vel << " k3q.acc: " << k3q.acc << endl;
	//
	//cout << "k4x.vel: " << k4x.vel << " k4x.acc: " << k4x.acc << endl;
	//cout << "k4q.vel: " << k4q.vel << " k4q.acc: " << k4q.acc << endl;
	//
	//cout << "x.pos: " << x.pos << " x.vel: " << x.vel << " x.acc: " << x.acc << endl;
	//cout << "q.pos: " << q.pos << " q.vel: " << q.vel << " q.acc: " << x.acc << endl;
}

