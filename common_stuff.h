
// includes
//#include <vector>
//#include <cmath>
#include <iostream>
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

// data structures
struct point {uint x; uint y;};
struct dpoint {double x; double y;};
struct hpoint {uint x; double y;};



