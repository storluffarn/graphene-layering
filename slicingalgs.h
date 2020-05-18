
#include "modtomlin.h"

using namespace std;

typedef unsigned int uint;

static const double zero = 1e-11;

// finding slipping points by inverval halfing
        void halfintervals (uint mode, uint adj, uint end, uint stride, uint pauseat, vector <double>* tvals, vector <double>* qvals, vector <double>* fvals,  vector <uint>* fslips, vector <uint>* qslips, vector <uint>* slips)
{

    vector <uint> slipsish;

    // method I: identify slipping regions by slope, then fins slipping points
    // by interval halving
    if (mode == 1)
    {
        uint curr = 0;
        uint next = stride;
        bool climbing;
        while (next < end)
        {
            //uint next = start + round(adj+exp(logstep*loops));

            double x1 = tvals->at(curr);
            double x2 = tvals->at(next);

            double y1 = fvals->at(curr);
            double y2 = fvals->at(next);

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
        
        slips->reserve(round(end/( (double) stride )));
        
        for (uint k = 0; k < slipsish.size() - 1; k++)
        {
            curr = slipsish[k];
            next = slipsish[k]+stride;
            double x1 = 0;
            x1 = tvals->at(curr);

            bool left = true;						// slip is to the left
            uint distance = abs((int)next - (int)curr);		// this is stridem but whatever...
            uint half = 0;	
            half = round( half + distance/2.0 );
            double xhalf = tvals->at(half);

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
                
                x1 += tvals->at(curr);
                half = round( half + distance/2.0 );
                xhalf += tvals->at(half);
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
    // method II find slipping rebion by difference tolerance, then find 
    else if (mode == 2)
    {
        uint curr = 0;
        uint next = stride;
        bool climbing;
        while (next < end)
        {
            //uint next = start + round(adj+exp(logstep*loops));

            double x1 = tvals->at(curr);
            double x2 = tvals->at(next);

            double y1 = fvals->at(curr);
            double y2 = fvals->at(next);

            double slope = (y2 - y1) / (x2 - x1);

            if (slope > 0.0)
                climbing = true;

            if (slope < 0.0 && climbing)
            {
                climbing = false;
                slipsish.push_back(next);
            }
            curr = next;
            next = curr + stride;
        }

        //cout << "slipsish has size: " << slipsish.size() << endl;

        //for (auto &el : slipsish)
        //	cout << el << endl;
        //cout << endl;
        
        slips->reserve(round(end/( (double) stride )));
        
        for (uint k = 0; k < slipsish.size() - 1; k++)
        {

            uint curr = slipsish[k];
            uint prev = curr - stride;
            uint next = curr + stride;
            uint distance = stride;

            //bool leftofmin;

            double ycurr = fvals->at(curr);
            double yprev = fvals->at(prev);
            double ynext = fvals->at(next);

            int p;

            if (ycurr - ynext > 0)
            {
                //leftofmin = true;
                p = 1;
            }
            else 
            {
                //leftofmin = false;
                p = -1;
            }

            
            while (distance > adj)
            {
                distance = round(distance / 2.0);
                
                next = round(curr + p*distance);

                ycurr = fvals->at(curr);
                ynext = fvals->at(next);
                
                if (ycurr - ynext > 0)
                {
                    p = 1;
                }
                else
                {
                    p = -1;
                }

                curr = next;
            }

            slips->push_back(curr);
        }
    }
    // slipping points by slope, use this for relaxation, needs small stride
    else if (mode == 3)
    {
        // for frics
        
        uint curr = 0;
        uint next = stride;
        vector<double>::iterator maxel;
        maxel = max_element(fvals->begin(), fvals->end());
        
        double max = *maxel;
        double low = 0.96*max;
        double slipbound = max - low;
        double slopetol = -1.25;

        //slipsish.push_back(pauseat);    

        while (next < end)
        {
            //uint next = start + round(adj+exp(logstep*loops));

            double y1 = fvals->at(curr);
            double y2 = fvals->at(next);
            double x1;          
            double x2;

            double diff = y1 - y2;
        
            if (diff > slipbound && curr > stride) // second conditions bypasses this check for first iteration to avoid overflow
            {
                uint near = curr - round(0.75*stride);
                
                x1 = tvals->at(near);
                x2 = tvals->at(curr);
                y1 = fvals->at(near);
                y2 = fvals->at(curr);
                
                double slope = (y2 - y1) / (x2 - x1);

                //cout << "xslip " << x1 << " " << y1 << " "  << slope << endl;

                if (slope < slopetol)
                {
                    fslips->push_back(curr);
                }
            }

            curr = next;
            next = curr + stride;
        }
        
        //fslips->push_back(curr); // this looks like an ugly hack, but it's really physically motivated -- sort of not a slip, but reports last minimum

        //// for q pos

        //curr = 0;
        //stride = round(0.5*stride);
        //next = stride;
        //maxel = max_element(qvals->begin(), qvals->end());
        //
        //max = *maxel;
        //low = 0.95*max;
        //slipbound = max - low;
        //slopetol = -1.0;
        //double slopetol2 = 2.0;

        ////slipsish.push_back(pauseat);    
        //
        //while (next < end)
        //{
        //    //uint next = start + round(adj+exp(logstep*loops));
        //    //cout << curr << endl;
        //    double y1 = qvals->at(curr);
        //    double y2 = qvals->at(next);
        //    double x1;          
        //    double x2;

        //    double diff = abs(y1 - y2);

        //    if (diff > slipbound && curr > stride)
        //    {
        //        uint near = curr - round(1.0*stride);
        //        
        //        x1 = tvals->at(near);
        //        x2 = tvals->at(curr);
        //        y1 = qvals->at(near);
        //        y2 = qvals->at(curr);
        //        
        //        double slope = abs((y2 - y1) / (x2 - x1));
        //        double tmp = (qvals->at(next) - y2) / (tvals->at(next) - x2);
        //        //cout << "qslip "  << x1 << " "  << slope << endl;

        //        if (slope > slopetol2)
        //        {
        //            qslips->push_back(curr);
        //        }
        //    }

        //    curr = next;
        //    next = curr + stride;
        //}
        
        //qslips->push_back(curr); 

        vector <uint> outs = *fslips;
        vector <uint> tmp = *qslips;        // this is just to collect the f and q slips, later we will eliminate dublicates
        
        // UNCOMMENT THESE TO ADD Q SLIPS
        //outs.insert(outs.end(),tmp.begin(),tmp.end());
        //sort(outs.begin(),outs.end());
        
        //for (auto& el : outs)
        //    cout << el << ",";
        //    cout << endl;
        
        //outs.erase(remove_if(outs.begin(), outs.end(), [stride](uint k, uint l){return l - k < 5*stride;}), outs.end());
        outs.erase(unique(outs.begin(), outs.end(), [stride](uint k, uint l){return l - k < 5*stride;}), outs.end());
        
        //for (auto& el : outs)
        //    cout << el << ",";
        //    cout << endl;
        
        *fslips = outs;    
    }
    else
    {
        cout << "no slicing operation mode selected, exiting" << endl;
    }
}


