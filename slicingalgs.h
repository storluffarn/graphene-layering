
#include "common_stuff.h"
#include "modtomlin.h"

using namespace std;

typedef unsigned int uint;

static const double zero = 1e-11;

// finding slipping points by inverval halfing
void halfintervals (uint adj, uint end, uint stride, vector <double>* xvals, vector <double>* yvals,  vector <uint>* slips)
{

	// part I find slipping region

	vector <uint> slipsish;

	uint curr = 0;
	uint next = stride;
	bool climbing;
	
	while (next < end)
	{
		//uint next = start + round(adj+exp(logstep*loops));

		double x1 = xvals->at(curr);
		double x2 = xvals->at(next);

		double y1 = yvals->at(curr);
		double y2 = yvals->at(next);

		double slope = (y2 - y1) / (x2 - x1);

		if (slope > 0.0)
			climbing = true;

		if (slope < 0.0 && climbing)
		{
			climbing = false;
			slipsish.push_back(curr);
		}
		curr = next;
		next = curr + stride;
	}
	
	//cout << "slipsish has size: " << slipsish.size() << endl;

	//for (auto &el : slipsish)
	//	cout << el << endl;
	//cout << endl;

	// part II find slipping point

	slips->reserve(round(end/( (double) stride )));
	
	for (uint k = 0; k < slipsish.size() - 1; k++)
	{
		curr = slipsish[k];
		next = slipsish[k]+stride;
		double x1 = 0;
		x1 = xvals->at(curr);

		bool left = true;						// slip is to the left
		uint distance = abs((int)next - (int)curr);		// this is stridem but whatever...
		uint half = 0;	
		half = round( half + distance/2.0 );
		double xhalf = xvals->at(half);

		//cout << "curr is: " << curr << " next is: " << next << " distance is: " << distance << endl;

		//cout << "distance greater than adj? " << ((distance > adj) ? 1 : 0) << " " << endl;	
		while (distance > adj)
		{
			if (left)
			{
				if (x1 - xhalf > zero)
				{
					left = true;
				}
				else
				{
					left = false;
					swap(curr,half);
				}
			}
			else
			{
				if (x1 - xhalf > zero)
				{
					left = false;
				}
				else
				{
					left = true;
					swap(curr,half);
				}
			}
			
			distance = abs((int)half-(int)curr);
			
			x1 += xvals->at(curr);
			half = round( half + distance/2.0 );
			xhalf += xvals->at(half);
		}

		if(left)
		{
			slips->push_back(curr);
		}
		else
		{
			slips->push_back(next);
		}
	}

	// part III post processing
	
//	for (uint k = 0; k < slips->size()-1; k++)
//	{
//		uint curr = slips->at(k);
//		uint next = slips->at(k+1);
//		
//		double y1 = afm->getfric(curr);
//		double y2 = afm->getfric(curr+1000*adj);
//	
//		// remove duplicates
//		if (next - curr < 250*adj)
//			slips->erase(slips->begin() + k+1);
//		
//		// remove semi-slips
//		else if (y1 - y2 < zero)
//			slips->erase(slips->begin() + k+1);
//	}
}


