
//
// looks as though it has to do with benchmarking?
//

#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <vector>

using namespace std;

//definitions

double const timestep = 1e-14;
double const restimestep = 1/timestep;
const int duration = 1500; // in 1e-14s, number of timesteps

const int T = 120;
const double sigma = 3.4e-10;
const double kB = 1.38064852e-23;
const double epsilon = T*kB;
const double umass = 1.660539040e-27;
const double argonmass = umass*39.948;
const int boxside = 6;
const double boxside2 = 10.299*sigma;
const double gridsize = boxside2/boxside;
const double npart = 2;
const double rcutoff = 2.25 * sigma;
const double oneoversix = 1.0/6.0;

struct particle // all relevant particle properties goes here
{
    double pox;
    double vex;
	double acx;
	double kinE;
	double potE;

	public:
	particle(double a, double b)
		: pox(a), vex(b)
	{}
};

struct virtualpart
{
	double pox;
	double vex;
	double acx;

	public:
	virtualpart(double a, double b)
		: pox(a), vex(b)
	{}
};

//initial conditions

double x1 = boxside2 - 10e-10;
double x2 = boxside2 - 0e-14;
double v1 = -100;
double v2 = 100;

particle p1(x1,v1);
particle p2(x2,v2);

double distance(particle* p1, particle* p2)
{
	double r = p2->pox - p1->pox;
	r -= static_cast<int> (r >= 0 ? r/boxside2 + 0.5 : r/boxside2 - 0.5)*boxside2;
	//r -= static_cast<int> (r/boxside2 + 0.5) * boxside2;
	
	return r;	
}

// 1-d LJ-force

double acc(particle* p1, particle* p2)
{
	//double r = fabs(k2.pox - k1.pox);
	//r -= static_cast<int> (r/boxside2 + 0.5) * boxside2; // +0.5 for cast to work // static_cast<int>(n >= 0 ? n + 0.5 : n - 0.5) for nearest int
	
	double r = distance(p1,p2);

	double acc;
	//if (r <= rcutoff)
		acc = - 24*epsilon*(-pow(sigma,6)/pow(r,7) + 2*pow(sigma,12)/pow(r,13) )/argonmass;
	//else 
	//	acc = 0;

	return acc;
}

double LJpot(particle* p1, particle* p2) 
{
	double r = distance(p1,p2);
	double V;

	//if (r <= rcutoff)
		V = 4*epsilon*(pow(sigma/r,12)-pow(sigma/r,6));
	//else
	//	V = 0;

	return V;
}

double qstep(double q, double qdot, double time) // might want to not return vector if we want to use to do pos and vel steps separetly?
{
	return q + qdot * time;
}

//void rk4sub(vector<vector<particle>>* kn, double time, double stepsize, int it)
//{
//	for (auto& k : *kn)
//	{
//		if (it == 2) stepsize -= 0.5;
//
//		double help = timestep*stepsize;
//
//		if (it == 0) {help = 1; stepsize = 0;}
//
//		k[it].pox = qstep(k[0].pox,help*k[it-1].vex,time+stepsize*timestep);
//		k[it].vex = qstep(k[0].vex,help*k[it-1].acx,time+stepsize*timestep);
//	}
//		
//	double getacc = acc((*kn)[it][0],(*kn)[it][1]);
//
//	(*kn)[0][it].acx = getacc; 
//	(*kn)[1][it].acx = -getacc;	
//}
//
//void rk4(double timestep, double time) // no particle input because they are global
//{
//	double stepsize = 0;
//	vector<vector<particle>> kn (5,vector<particle>(2));
//	kn[0] = partlist;
//
//	for (unsigned int k = 1; k < kn[0].size(); k++)
//	{
//		rk4sub(&kn, time, stepsize, k);
//		
//		stepsize += 0.5;
//	}
//	
//	for(auto& k : partlist)
//	{	
//		auto l = &k - &partlist[0];
//
//		k.pox += timestep/6.0 * (kn[l][0].vex + kn[l][1].vex + kn[l][2].vex + kn[l][3].vex);
//		k.vex += timestep/6.0 * (kn[l][0].acx + kn[l][1].acx + kn[l][2].acx + kn[l][3].vex);
//	}
//}

void rk4naive(particle* p1, particle* p2, double timestep)
{
	double pos1, vel1, pos2, vel2, kacc;
	particle vp1(0,0); particle vp2(0,0);
	struct k {double vex; double acx;};
	
	//k1 -- initial accelerations calculated in main loop
	
	k k11 = {p1->vex,p1->acx};
	k k12 = {p2->vex,p2->acx};
	
	//k2
	
	pos1 = p1->pox+0.5*timestep*k11.vex;
	vel1 = p1->vex+0.5*timestep*k11.acx;
	vp1.pox = pos1; vp1.vex = vel1;	
	
	pos2 = p2->pox+0.5*timestep*k12.vex;
	vel2 = p2->vex+0.5*timestep*k12.acx;
	vp2.pox = pos2; vp2.vex = vel2;	
	
	kacc = acc(&vp1, &vp2);
	
	k k21 = {vel1,kacc};
	k k22 = {vel2,-kacc};

	//k3
	
	pos1 = p1->pox+0.5*timestep*k21.vex;
	vel1 = p1->vex+0.5*timestep*k21.acx;
	vp1.pox = pos1; vp1.vex = vel1;	
	
	pos2 = p2->pox+0.5*timestep*k22.vex;
	vel2 = p2->vex+0.5*timestep*k22.acx;
	vp2.pox = pos2; vp2.vex = vel2;	
	
	kacc = acc(&vp1, &vp2);
	
	k k31 = {vel1,kacc};
	k k32 = {vel2,-kacc};
	
	//k4
	
	pos1 = p1->pox+timestep*k31.vex;
	vel1 = p1->vex+timestep*k31.acx;
	vp1.pox = pos1; vp1.vex = vel1;	
	
	pos2 = p2->pox+timestep*k32.vex;
	vel2 = p2->vex+timestep*k32.acx;
	vp2.pox = pos2; vp2.vex = vel2;	
	
	kacc = acc(&vp1, &vp2);
	
	k k41 = {vel1,kacc};
	k k42 = {vel2,-kacc};

	//rk-step
	
	p1->pox += timestep*oneoversix * (k11.vex + 2*k21.vex + 2*k31.vex + k41.vex);
	p1->vex += timestep*oneoversix * (k11.acx + 2*k21.acx + 2*k31.acx + k41.acx);
	p2->pox += timestep*oneoversix * (k12.vex + 2*k22.vex + 2*k32.vex + k42.vex);
	p2->vex += timestep*oneoversix * (k12.acx + 2*k22.acx + 2*k32.acx + k42.acx);
}

// make sub rutines for printing, and energy calculations
//
//


int main()
{	
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	
	//double force;
	double acc0;

	double Etot1;
	double Etot2;
	double Etot;	
	
	acc0 = acc(&p1,&p2);
	p1.acx = acc0;
	p2.acx = -acc0;

	std::ofstream verlet1dfile;
	verlet1dfile.open("periodicrk.txt");

	for(int k = 0; k < duration; k++)
	{
		rk4naive(&p1, &p2, timestep);

		acc0 = acc(&p1,&p2);
		p1.acx = acc0;
		p2.acx = -acc0;
		
		p1.kinE = argonmass*p1.vex*p1.vex/2;
		p1.potE = LJpot(&p1,&p2);
		
		p2.kinE = argonmass*p2.vex*p2.vex/2;
		p2.potE = LJpot(&p1,&p2);
		
		Etot1 = p1.kinE + p1.potE;
		Etot2 = p2.kinE + p2.potE;
		Etot = p1.kinE + p1.potE - (p2.kinE + p2.potE);
		
		verlet1dfile << p1.pox << " " << p1.vex << " " << p1.kinE << " " << p1.potE << "\n";
		verlet1dfile << p2.pox << " " << p2.vex << " " << p2.kinE << " " << p2.potE << "\n";
		verlet1dfile << Etot1 << " " << Etot2 << " " << Etot << "\n";
		verlet1dfile << p1.acx << " " << p2.acx << "\n";
		verlet1dfile << distance(&p1,&p2) << endl;
	}
verlet1dfile.close();
	
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

	double benchtime = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	
	std::cout << "execution in microseconds: " << benchtime << "\n";

	return 0;
}













