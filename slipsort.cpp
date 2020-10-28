
/*
 *
 *  This program could be much more memory efficient by more pointer usage
 *  ...but memory is cheap as they say!
 *
*/

#include <iostream>                                             // basic io stuff
#include <fstream>                                              // writing to file
#include <cmath>                                                // for math functions
#include <algorithm>                                            // for count_if
#include <vector>        // for using the vector type members
#include <functional>
//#include "common_stuff.h"

typedef unsigned int uint;

using namespace std;

class labelclass
{
    // subclasses

    class dataclass						// data to be analysed
    {
    	// private section
    
    	// data members
    	
    	string filename;							// name of file to read
        vector <double> ts;
    	vector <double> xs;						    // data read
    	vector <double> qs;						    // data read
    	uint size;									// numbers written/read
    
    	// function members
    	
    	void readdata();			// function for reading data
    	
        // friendship is magic
        
    	// public section
        
    	public:
    	
    	// constructor
    	dataclass(string s)
    		: filename(s)							// setting initial values
    	{
    		readdata();		    			// read data from file
    		size = xs.size();						// save size of data
            ts.reserve(size);
    		xs.reserve(size);						// save size of data
    		qs.reserve(size);						// save size of data
    	}
    	
    	//function members
    	
    	// inline functions
    
    	vector <double>* getts(){return &ts;}
    	vector <double>* getxs(){return &xs;}
    	vector <double>* getqs(){return &qs;}
    	double gett(uint i){return ts[i];}
    	double getx(uint i){return xs[i];}
    	double getq(uint i){return qs[i];}	
        uint getsize() {return size;}
    	
    };
    
    struct slippoint
    {
        double t;
        double x;
        double q;
        char label;
    };
    
    struct target
    {
        double x;
        double q;
    }; 

    // object members
    
    dataclass data;
    vector <slippoint> slipdata;
    vector <target> targets;
    
    // data members
    
    string filename1;
    uint size;
    uint tsize;

    // function members
    
    char findlabel(uint);
    
    public: 

    // function members
    void labeldata();
    
    // accessors
    vector <slippoint>* getslipdata() {return &slipdata;}; 
    uint getsize() {return size;}; 
    
    // prints
    void printslipdata(); 

    labelclass (string filename1, vector <double>* targetsx, vector <double>* targetsq)
        : data(filename1)
    {
        size = data.getsize();
        tsize = targetsx->size();

        targets.reserve(tsize);
        slipdata.reserve(size);

        for (uint k = 0; k < tsize; k++)
        {
            target xq = {targetsx->at(k),targetsq->at(k)};
            targets.push_back(xq);
        }
        
        for (uint k = 0; k < size; k++)
        {
            slippoint slip = {data.gett(k),data.getx(k),data.getq(k),};
            slipdata.push_back(slip);
        }
    }
};

// subclass definitions
    
void labelclass::dataclass::readdata()			// reads data from file
{
    ifstream indata (filename);
    
    if (indata.is_open())
    {
        double t, f, q;
        char c;
        uint i = 0;

        while (indata >> t >> c >> f >> c >> q && c == ',')
        {
            double x = -0.5*f + 2.448e-9;        // force to position, spring was 2 N/m, support stopped at 2.448 nm
            
            if (t > 1e-10)  // it is expected that data before this decorrelation time is faulty
            {
                ts.push_back(t);
                xs.push_back(x);
                qs.push_back(q);
            }
        }
    }
}

void labelclass::labeldata()
{
    vector <char> labels;
    labels.reserve(size);
    
    for (uint k = 0; k < size; k++)
    {
        double x = data.getxs()->at(k);
        double q = data.getqs()->at(k);
    
        vector <double> offsetsx;
        vector <double> offsetsq;

        for (uint l = 0; l < tsize; l++)
        {
            double tx = targets[l].x;
            double offsetx = abs(x - tx);
            offsetsx.push_back(offsetx/tx);
            
            double tq = targets[l].q;
            double offsetq = abs(q - tq);
            offsetsq.push_back(offsetq/tq);
        }

        //for (auto &el : offsetsx)
        //    cout << el << " ";
        //cout << endl;

        vector <double> norms (tsize);
        for (auto &el : norms)
        {
            uint l = &el - &norms[0];
            el = offsetsx[l] + offsetsq[l];
        }

        transform(norms.begin(),norms.end(),norms.begin(),[](double n){return 0.5*n;});

        auto minit = min_element(norms.begin(),norms.end());
        uint minel = distance(norms.begin(),minit);

        char label = findlabel(minel);
        
        //cout << label << endl;

        labels.push_back(label);
    }

    for (uint k = 0; k < size; k++)
    {
        slipdata[k].label = labels[k];
    }
}

char labelclass::findlabel(uint n)
{
    if (n == 0)
        return 'A';   
    else if (n == 1)
        return 'B';
    else if (n == 2)
        return 'C';
    else if (n == 3)
        return 'D';
    else if (n == 4)
        return 'E';
    else if (n == 5)
        return 'F';
    else if (n == 6)
        return 'G';
    else if (n == 7)
        return 'H';
    else if (n == 8)
        return 'I';
    else if (n == 9)
        return 'J';
    else if (n == 10)
        return 'K';
    else if (n == 11)
        return 'L';
    else if (n == 12)
        return 'M';
    else if (n == 13)
        return 'N';
    else if (n == 14)
        return 'O';
    else if (n == 15)
        return 'P';
    else if (n == 16)
        return 'Q';
    else if (n == 17)
        return 'R';
    else 
        return '!';
}

void labelclass::printslipdata()
{
    cout << "slip data: t, x, q, label" << endl;
    for (auto &slip : slipdata)
    {
        cout << slip.t << " " << slip.x << " " << slip.q << " " << slip.label << endl;
    }
}

int main()
{
    //string filename1 = "./data/to_jenny/histogram/test.csv";
    string filename1 = "./data/to_jenny/histogram/slips.csv";
    string filename2 = "./data/to_jenny/histogram/sliptos.csv";

    vector <double> targetsx = {0.037,1.268,1.503,1.736,1.966,2.207,1.544,1.767,1.996,2.223,2.462,1.800,2.016,2.241,2.459,2.034,2.254,2.479};

    for (auto &el : targetsx)
        el *= 1e-9;

    vector <double> targetsq = {7.525,7.355,7.217,7.057,6.859,6.778,5.255,5.002,4.808,4.576,4.463,2.922,2.584,2.346,2.042,0.338,0.043,0.000};
    
    for (auto &el : targetsq)
        el *= 1e-10;

    labelclass labeleddatafrom(filename1,&targetsx,&targetsq);
    labeleddatafrom.labeldata();
    auto fromdata = labeleddatafrom.getslipdata();

    labelclass labeleddatato(filename2,&targetsx,&targetsq);
    labeleddatato.labeldata();
    auto todata = labeleddatato.getslipdata();

    uint size = labeleddatafrom.getsize();

    ofstream fs;
    fs.open("./data/labelleddata.dat");

    fs << "Slip data in SI units. Data ordering: slip from (t,x,q,label), slip to (x,q)" << endl;

    for (uint k = 0; k < size; k++)
    {
        fs << fromdata->at(k).t << " " << fromdata->at(k).x << " " << fromdata->at(k).q << " " << fromdata->at(k).label << " "
           << todata->at(k).x << " " << todata->at(k).q << " " << todata->at(k).label << endl; 
    }

    fs.close();

    //labeleddata.printslipdata();
}

