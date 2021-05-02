
#include "tomlinclass.cpp"

using namespace std;

typedef unsigned int uint;

static const double zero = 1e-11;

// finding slipping points by inverval halfing
        void halfintervals (uint mode, uint adj, uint end, uint stride, uint pauseat, vector <double>* tvals, vector <double>* qvals, vector <double>* fvals,  vector <uint>* fslips, vector <uint>* qslips, vector <uint>* slips, vector <uint>* sliptos)
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
        double slopetol = -5.0;

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
                uint near = curr + round(0.75*stride);
                
                x2 = tvals->at(near);
                x1 = tvals->at(curr);
                y2 = fvals->at(near);
                y1 = fvals->at(curr);
                
                double slope = (y2 - y1) / (x2 - x1);

                //cout << "foundx " << x1 << " " << y1 << " "  << slope << endl;

                if (slope < slopetol)
                {
                    fslips->push_back(curr);
                    //uint to = curr + 10*stride;
                    //cout << "pushed " << x1 << " with " << y1 << endl;
                    //if (to < end) // this adds 'slip to' functionality
                    //{
                    //    fslips->push_back(to);
                    //    //cout << "slipped at" << endl <<
                    //    //    tvals->at(curr) << " with " <<
                    //    //    fvals->at(curr) << " to" << endl <<
                    //    //    tvals->at(to) << " with " << 
                    //    //    fvals->at(to) << endl;
                    //}
                }
            }

            curr = next;
            next = curr + stride;
        }
        
        //fslips->push_back(curr); // this looks like an ugly hack, but it's really physically motivated -- sort of not a slip, but reports last minimum

        // for q pos

        curr = 0;
        uint qstride = round(0.5*stride); // needs to be shorter to capture spikes
        next = qstride;
        maxel = max_element(qvals->begin(), qvals->end());
        
        max = *maxel;
        low = 0.925*max;
        slipbound = max - low;
        double slopetol2 = 2.5;

        //slipsish.push_back(pauseat);    
        
        while (next < end)
        {
            //uint next = start + round(adj+exp(logstep*loops));
            //cout << curr << endl;
            double y1 = qvals->at(curr);
            double y2 = qvals->at(next);
            double x1;          
            double x2;

            double diff = abs(y1 - y2);

            if (diff > slipbound && curr > qstride)
            {
                uint near = curr + round(0.75*qstride);
                
                x2 = tvals->at(near);
                x1 = tvals->at(curr);
                y2 = qvals->at(near);
                y1 = qvals->at(curr);
                
                double slope = abs((y2 - y1) / (x2 - x1));
                //double tmp = (qvals->at(next) - y2) / (tvals->at(next) - x2);
                //cout << "foundq "  << x1 << " "  << slope << endl;

                if (slope > slopetol2)
                {
                    qslips->push_back(curr);
                    //uint to = curr + 10*qstride; 
                    //if (to < end)
                    //{
                    //    qslips->push_back(to);
                    //    cout << "slipped at" << endl <<
                    //        tvals->at(curr) << " with " <<
                    //        qvals->at(curr) << " to" << endl <<
                    //        tvals->at(to) << " with " << 
                    //        qvals->at(to) << endl;
                    //        
                    //}
                }
            }
            curr = next;
            next = curr + qstride;
        }
        
        //qslips->push_back(curr); 

        vector <uint> aslips = *fslips;
        vector <uint> tmp = *qslips;        // this is just to collect the f and q slips, later we will eliminate dublicates
        
        // UNCOMMENT THESE TO ADD Q SLIPS
        aslips.insert(aslips.end(),tmp.begin(),tmp.end());
        sort(aslips.begin(),aslips.end());
        
        //for (auto& el : aslips)
        //    cout << el << ",";
        //    cout << endl;
        
        //aslips.erase(remove_if(aslips.begin(), aslips.end(), [stride](uint k, uint l){return l - k < 5*stride;}), aslips.end());
        aslips.erase(unique(aslips.begin(), aslips.end(), [stride](uint k, uint l){return l - k < 2*stride;}), aslips.end()); // this isn't really working with slip-to enabled
        
        //for (auto& el : aslips)
        //    cout << el << ",";
        //    cout << endl;

        vector <uint> tos = aslips;
        for (auto& el : tos)
            el += 10*stride;

        *slips = aslips;
        *sliptos = tos;
    }
    else
    {
        cout << "no slicing operation mode selected, exiting" << endl;
    }
}

