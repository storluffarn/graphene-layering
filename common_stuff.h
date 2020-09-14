
// includes
#include <vector>
#include <iostream>
#include <cmath>
//#include <fstream>
//#include <sstream>
//#include <map>
#include <algorithm>
#include <gsl/gsl_poly.h>

using namespace std;
typedef unsigned int uint;


// constants
static const double pi = atan(1.0)*4;
static const double kB = 1.38064852e-23;
;

// standard values
const double barref = 5.0e-20; // 4.0e-20 for double slips // 5.0 for single slips, basically a swith for thermolubricity
const double kapparef = 0.0612245;

const double spring = 1e-0* 2.0;          // 2.0 //1e-4 for harmonic
const double supvel = 1.0;                // 1.0
const double latcon = 2.5e-10;            // 2.5e-10
const double barr1 = barref;             // 
const double barr2 = 0.5 * barref;
const double kappa1 = kapparef;
const double kappa2 = 0.5 * kapparef;
const double align = 1.0; //(1+sqrt(5.0))/2.0; 
const double nu2 = 0.382653;
const double nu4 = 5.809e17;
     
const double xmass = 1e-23;
const double xdamp = 1.875e13;

//const double qmass = 3.67143e-24;
//const double qdamp = 4.28571e13;
const double qmass = 3.67e-24;
const double qdamp = 4.28e13;

const double temp = 300;  // 5e48 lest we'd forget

// derived constants
const double latcona = 2.5e-10;            // 2.5e-10
const double latcona2pi = 2.0*pi/latcona;            // 2.5e-10
const double latconb = latcona / align;
const double latconb2pi = 2.0*pi/latconb;


double c1 = -0.5*pi*kappa1/(latcona*nu4);
double c2 = 0.5*(kappa1 + nu2)/nu4;
double c3 = -0.5*pi*barr1/(latcona*nu4);
double qmax = 0;
double void1, void2;
int tmp = gsl_poly_solve_cubic(c1,c2,c3, &qmax, &void1, &void2);

// data structures
struct point {uint x; uint y;};
struct dpoint {double x; double y;};
struct hpoint {uint x; double y;};

// some very manual debugging...
int pingcount = 0;
void ping(const int line) {cout << "ping at line " << line << endl; pingcount++;}
void getpingcount() {cout << pingcount << endl;}

void printvectoruint (vector <uint>* v)
{ 
    cout << "[" ;

    for_each(v->begin(),v->end(), [](uint el) {cout << el << ", ";} );
    
    cout << "]" << endl;
}

void printvectordouble (vector <double>* v)
{ 
    cout << "[" ;

    for_each(v->begin(),v->end(), [](double el) {cout << el << ", ";} );
    
    cout << "]" << endl;
}

void scalevectordouble (vector <double>* v, double c)
{  
    transform(v->begin(),v->end(),v->begin(),[&c](double el){return el/c;});
}
