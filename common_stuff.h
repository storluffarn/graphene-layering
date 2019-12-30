
// includes
#include <vector>
#include <iostream>
#include <cmath>
//#include <fstream>
//#include <sstream>
//#include <map>

using namespace std;
typedef unsigned int uint;

// some very manual debugging...
int pingcount = 0;
void ping(const int line) {cout << "ping at line " << line << endl; pingcount++;}
void getpingcount() {cout << pingcount << endl;}

// constants
static const double pi = atan(1.0)*4;
static const double kB = 1.38064852e-23;
;

// standard values
const double barref = 5.0e-20; // 4.0e-20 for double slips // 5.0 for single slips, basically a swith for thermolubricity
const double kapparef = 0.0612245;

const double spring = 1e-0* 2.0;          // 2.0 //1e-4 for harmonic
const double supvel = 0.0;                // 1.0
const double latcon = 2.5e-10;            // 2.5e-10
const double barr1 = barref;              // 
const double barr2 = 0.5 * barref;
const double kappa1 = kapparef;
const double kappa2 = 0.5 * kapparef;
const double align = 1.0; //(1+sqrt(5.0))/2.0; 
const double nu2 = 0.382653;
const double nu4 = 5.809e17;
     
const double xmass = 1e-23;
const double xdamp = 1.875e13;

const double qmass = 3.67143e-24;
const double qdamp = 4.28571e13;

const double temp = 300;  // 5e48 lest we'd forget

// derived constants
const double latcona = 2.5e-10;            // 2.5e-10
const double latcona2pi = 2.0*pi/latcona;            // 2.5e-10
const double latconb = latcona / align;
const double latconb2pi = 2.0*pi/latconb;

// data structures
struct point {uint x; uint y;};
struct dpoint {double x; double y;};
struct hpoint {uint x; double y;};



