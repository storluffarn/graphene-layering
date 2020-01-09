
#include "modtomlin.h"

using namespace std;

typedef unsigned int uint;

static const double zero = 1e-11;

// finding slipping points by inverval halfing
void halfintervals (uint mode, uint adj, uint end, uint stride, uint pauseat, vector <double>* xvals, vector <double>* yvals,  vector <uint>* slips)
{

	// part I find slipping region

	vector <uint> slipsish;

	uint curr = 0;
	uint next = stride;

    // finding maxima in noisy stick-slip signal
    if (mode == 1)
    {
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
    }
    else if (mode == 2)
    {
        vector<double>::iterator maxel;
        maxel = max_element(yvals->begin(), yvals->end());
        
        double max = *maxel;
        double low = 0.975*max;
        double slipbound = max - low;
        double flatbound = 0.1;

        cout << "max element: " << max << " slipbound used: " << slipbound <<  endl;

	    while (next < end)
	    {
	    	//uint next = start + round(adj+exp(logstep*loops));

	    	//double x1 = xvals->at(curr);
	    	//double x2 = xvals->at(next);

	    	double y1 = yvals->at(curr);
	    	double y2 = yvals->at(next);

	    	double diff = y1 - y2;

	    	if (diff > slipbound)
	    	{
	    		slipsish.push_back(curr);
	    	}

	    	curr = next;
	    	next = curr + stride;
	    }
        
        uint ttol = 1000;
        double slopetol = -2.0;
        vector <uint> outs; 
        bool slipped = false;

        for (uint k = 0; k < slipsish.size()-1; k++)
        {
            uint slipat = slipsish[k];
            //cout << "pause: " << pauseat << endl << "k: " << slipat <<  endl;
            
            if (slipat > pauseat)
	        {	
                double x1 = xvals->at(slipsish[k]);
	    	    double x2 = xvals->at(slipsish[k+1]);

                double y1 = yvals->at(slipsish[k]);
	    	    double y2 = yvals->at(slipsish[k+1]);

	    	    double slope = (y2 - y1) / (x2 - x1);
           
                //cout << "slope: " << slope << endl;

                if (slope < slopetol && !slipped)
                {
                    outs.push_back(slipat);
                    slipped = true;
                    //cout << "pushed: " << slipat << endl;
                }
                else if (slope > slopetol && slipped)
                {
                    slipped = false;
                    if (k == slipsish.size()-2)
                    {
                        outs.push_back(slipsish[k+1]);
                        //cout << "pushed: " << slipsish[k+1] << endl;
                    }
                }
            }
        }
        
        (*slips) = outs;
        //(*slips) = slipsish;
    }
    else
    {
        cout << "no slicing operation mode selected, exiting" << endl;
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


