
#include "common_stuff.h"
#include "modtomlin.h"

using namespace std;

typedef unsigned int uint;

static const double zero = 1e-11;

void halfintervals (uint adj, uint end, uint stride, tomlin* afm, vector <uint>* slips)
{
	// part I find slipping region

	vector <uint> slipsish;

	//uint loops = 0;
	//double logstep = 1.0;		// emperical parameter to tweak log stepping

	//uint start = 0;
	uint curr = 0;
	uint next = stride;

	while (next < end)
	{
		//uint next = start + round(adj+exp(logstep*loops));

		double x1 = afm->getposx(curr);
		double x2 = afm->getposx(next);

		double y1 = afm->getfric(curr);
		double y2 = afm->getfric(next);

		double slope = (y2 - y1) / (x2 - x1);

		if (slope < -0.0)
		{
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
		next = slipsish[k+1];
		double x1 = 0;
		x1 = afm->getposx(curr);

		bool left = true;						// slip is to the left
		uint distance = abs((int)next - (int)curr);		// this is stridem but whatever...
		uint half = 0;	
		half = round( half + distance/2.0 );
		double xhalf = afm->getposx(half);

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
			
			x1 += afm->getposx(curr);
			half = round( half + distance/2.0 );
			xhalf += afm->getposx(half);
			
			// below can be used to us the mean if points instead
			//x1 = 0;
			//for (int l = -avgpartsize; l < avgpartsize + 1; l++)
			//{
			//	x1 += afm->getposx(curr+l);
			//}
			//x1 /= double(avgpartsize);
			//
			//half = round( half + distance/2.0 );

			//xhalf = 0;
			//for (int l = -avgpartsize; l < avgpartsize + 1; l++)
			//{
			//	xhalf += afm->getposx(half+l);
			//}
			//xhalf /= double(avgpartsize);
			
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

//void firstel (vector <dpoint>* avgpts, double* avg, double skip, uint halfmeansize, double rmeansize, uint* kk, uint* midpt, tomlin* afm)
//{
//	double sum = 0;
//	uint k = *kk;
//	uint its = floor((double) halfmeansize/skip);
//
//	for (uint l = 0; l < its; l++)		
//	{
//		sum += afm->getfric(k);
//		k += skip;
//	}
//
//	*midpt = k;
//	sum += afm->getfric(k);
//	k += skip;
//	
//	for (uint l = 0; l < its; l++)
//	{
//		sum += afm->getfric(k);
//		k += skip;
//	}
//		
//	*avg = sum * rmeansize;
//
//	dpoint firstpoint = {afm->gettime(*midpt),*avg};
//
//	avgpts->push_back(firstpoint);
//
//	*kk = k;
//}
//
//void runningmean (uint skip, uint halfmeansize, uint end, tomlin* afm, vector <dpoint>* avgpts)
//{
//	avgpts->reserve(round(end/skip));
//
//	double avg = 0;
//	
//	double rmeansize = 1.0/(floor(((double) 2*halfmeansize)/skip)+1);
//
//	cout << rmeansize << " " << 1.0 / rmeansize << endl;
//	
//	uint bgnel = 0;
//	uint midel = 0;
//	uint endel = 0;
//
//	firstel(avgpts,&avg,skip,halfmeansize,rmeansize,&endel,&midel, afm);
//
//	while (endel < end)
//	{
//		avg -= rmeansize*afm->getfric(bgnel);
//		bgnel += skip;
//		
//		midel += skip;
//		
//		avg += rmeansize*afm->getfric(endel);
//		endel += skip;
//		
//		dpoint avgel = {afm->gettime(midel), avg};
//
//		avgpts->push_back(avgel);
//	}
//	
//	// mostly for debugging
//	ofstream fs;
//	fs.open("avgs.csv");
//
//	for (uint k = 0; k < avgpts->size(); k++)
//	{
//		fs << avgpts->at(k).x << "," << avgpts->at(k).y << endl;
//	}
//
//	fs.close();
//}



