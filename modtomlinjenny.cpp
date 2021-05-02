
//
// extracting minima and maxima from tomlinson
//
// TODO
// * during production runs, all quantities are pushed, even though most aren't of interest...
// * the whole treatments of storing slips to file is an affront to god and all that is holy, but you know, it works ...for now
// * data handling should go in a separate file during development, can be marge for production runs
//
//
// NOTE 
//
// All q stuff removed from slicing algs and noise red at the moment. Sorry, there was really no nice way of doing it, it's a MESS of commented and uncommented stuff probably unintelligable to anyone but me, actually, to me too... Best bet to restore is probably to search for all instances of q and see if they should be there or not. Search has to be done in all files...

// includes

#include <fstream>
#include <vector>
#include "slicingalgs.cpp"		// I fucked up my dependancies, modtomlin.h is needed but is included in slicingalgs.h... trying to save thins with common_stuff.h, still bad
#include <chrono>
#include <string>

// lazy stuff

using namespace std;

void jennydist(double tstep, uint res, vector <pair<double,double>>* outs)
{
	// converting notation
	double latcona = latcon;	// not nice, should have consistant naming...
	double latconb = latcon / align;
	
	// some frequent constants
	double tmp = sqrt( pow(2*latcona2pi*kappa1 + latconb2pi*kappa2,2) + 
					   12*nu4*(spring - 2*nu2 - 2*kappa1 - 2*kappa2) );
	double xt = ( 2*latcona2pi*kappa1 + latconb2pi*kappa2 + tmp ) / (12*nu4) 
				  + latcona / 4.0;

	double um = barr2 + 0.25*kappa2*pow(xt,2);
	double km = spring - 2*nu2 - nu4*pow(xt,2) + latcona2pi*kappa1*xt - 2*kappa1 - 2*kappa2;
	double omegak = 1.0/( latconb2pi * sqrt(um/km) );
	double xc = acos(-pow(omegak,2))/latconb2pi;
	double rc = ( latconb2pi*xc + 1/pow(omegak,2)*sin(latconb2pi*xc) ) / latconb2pi;
	double vstar = 2*spring*supvel*xdamp*xmass*latconb/(km*kB*temp)*pow(omegak,2) / 
				   sqrt(1-pow(omegak,4));

	cout << "jennies constants: " << endl
		 << "tmp " << tmp << endl
	     << "xt " << xt << endl
	     << "um " << um << endl
	     << "km " << km << endl
	     << "ok " << omegak << endl
	     << "xc " << xc << endl
	     << "rc " << rc << endl 
	     << "vs " << vstar << endl;
	
	for (uint k = 0; k < res; k++)
	{
		double t = k*tstep;

		double r = (spring*supvel*t - 2*latcona2pi*barr1)/km;
		double f = 1.0 - r/rc;
		double du = 2.0/3.0 * um * pow(2*latconb2pi*rc,1.5) *
						pow(omegak,3)/pow(1-pow(omegak,4),0.25)*pow(f,1.5);
		double fstar = pow(du/(kB*temp),2.0/3.0);

		double p = 1.5*sqrt(fstar)/vstar*exp(-pow(fstar,1.5)-1.0/vstar*exp(-pow(fstar,1.5)));
	
		if (k == 0)
		{
			cout << "t " << t << endl 
				 << "r " << r << endl 
				 << "f " << f << endl 
				 << "du " << du << endl
				 << "fs " << fstar << endl 
				 << "p " << p << endl;
				 
		}

		pair<double,double> out = {supvel*t,p};
		outs->push_back(out);
	}
}

double findqmax()
{
    double p1 = 144*pi*kappa1*nu2*nu4/latcona;
    double p2 = 864*pi*pow(nu4,2.0)*barr1/latcona;
    double p3 = 144*pi*nu4*pow(kappa1,2)/latcona;
    double p4 = 16*pow(pi*kappa1/latcona,3.0);
    double p5 = 4*pow(pi*kappa1/latcona,2.0);
    double p6 = 24*nu4*(kappa1 + nu2);

    double p0 = p2 + p4 - p1 - p3 + sqrt(pow(p2+p4-p1-p3,2) + 4*pow(p6-p5,2));

    double qmax = latcona2pi*( kappa1 * pow( pow(0.5*p0,1.0/3.0)/(12.0*nu4) - (p6-p5)/(6.0*nu4*pow(4.0*p0,1.0/3.0)) + pi*kappa1/(6.0*latcona*nu4) ,2.0) + barr1 ); 

    return qmax;
}

double sangconstantq()
{
    double qmax = findqmax();

	double p0 = latcona*sqrt(xmass/(barr1 + kappa1*pow(qmax,2)));
	double pk = sqrt(xmass/spring);
	double ok = 1.0/latcona2pi*sqrt(spring/(barr1 + kappa1*pow(qmax,2.0)));
	double ke = spring/(1+pow(ok,2));
	double xc = acos(-pow(ok,2))/latcona2pi;
	double rc = (latcona2pi*xc + 1.0/pow(ok,2)*sin(latcona2pi*xc)) / latcona2pi;
	double vs = 2.0*supvel*xdamp*(barr1 + kappa1*pow(qmax,2.0))*pow(p0*ok,2)/(kB*temp*latcona*sqrt(1.0-pow(ok,4)));
	double df = pi*(barr1+kappa1*pow(qmax,2.0))/latcona * pow(1.5*kB*temp/(barr1+kappa1*pow(qmax,2.0)),2.0/3.0) * pow(1.0-pow(ok,4),1.0/6.0) / (1 + pow(ok,2));
	double fc = ke*(rc-latcona/2.0);
	
	double avgfric = fc - df*pow(abs(log(vs)),2.0/3.0);
    
	return avgfric;
}

double kramer(int n)
{
    double qmax = findqmax();
    
    double phi = latcona*(n-0.25) + latconb;
    double lambda = (2*nu2 + 4*nu4*pow(qmax,2) + 2*kappa1) / (latcona2pi*(barr1 + kappa1*pow(qmax,2)));
    double phila = phi - lambda/pow(latcona2pi,2);
    double qc = phila - sqrt(-lambda*(phila + phi) + 2)/latcona2pi;
    double omega = sqrt( latcona2pi*(barr1+kappa1*pow(qmax,2))/xmass * ( lambda + latcona2pi*sqrt(1 - pow(lambda*qc,2))) );
    double du = 4*(nu2 + kappa1)*(phila)*(qc - phila) + 
                8*nu4*phila*( 2*pow(qc,3) + 3*pow(qc,2)*phila + 4*qc*pow(phila,2) - 2*pow(phila,3) ) +
                4*kappa1*phila*(phila-qc)*sqrt(1-pow(lambda*qc,2));

    double tau = pow(omega,2)/(2*pi*xdamp)*exp(-du/(kB*temp));

    cout << "kB*T " << kB*temp<< endl << 
            "n " << n << endl <<
            "qmax "  << qmax << endl <<
            "phi " << phi << endl << 
            "lambda " << lambda << endl << 
            "phila " << phila << endl << 
            "qc " << qc << endl << 
            //"omega " << omega << endl << 
            "du " << du << endl << 
            "tau " << tau << endl <<
            //"1/tau " << 1.0/tau << endl << 
            "part1 " <<  pow(omega,2)/(2*pi*xdamp) << " part2 "  << exp(-du/(kB*temp)) << endl << endl;
    return tau;
}

double sangjenny()
{
	double J = 3.0;
	double km = spring + 2.0*(nu2 + 6.0*nu4*pow(0.5*J*latconb,2.0) + kappa1 + kappa2);
	double um = barr2 + kappa2*pow(0.5*J*latconb,2.0);
	double om = 1.0/latconb2pi*sqrt(km/um);
	double ke = km/(1.0+pow(om,2.0));
	double xc = acos(-pow(om,2.0))/latconb2pi;
	double rc = (latconb2pi*xc + 1.0/pow(om,2.0)*sin(latconb2pi*xc))/latconb2pi;
	double vs = 2.0*spring*supvel*xdamp*xmass*latconb*pow(om,2.0)/(km*kB*temp*sqrt(1.0-pow(om,4.0)));
	double df = pi*um/latconb * pow(1.5*kB*temp/um,2.0/3.0) * pow(1.0-pow(om,4.0),1.0/6.0) / (1 + pow(om,2));
	double fc = ke*(rc-0.5*latconb);
		
	double avgfric = fc - df*pow(abs(log(vs)),2.0/3.0);

	return avgfric;
}

double sangavgfric ()
{
	double p0 = latcona*sqrt(xmass/barr1);
	double pk = sqrt(xmass/spring);
	double ok = p0/(2*pi*pk);
	double ke = spring/(1+pow(ok,2));
	double xc = acos(-pow(ok,2))/latcona2pi;
	double rc = (latcona2pi*xc + 1.0/pow(ok,2)*sin(latcona2pi*xc)) / latcona2pi;
	double vs = 2.0*supvel*xdamp*barr1*pow(p0*ok,2)/(kB*temp*latcona*sqrt(1.0-pow(ok,4)));
	double df = pi*barr1/latcona * pow(1.5*kB*temp/barr1,2.0/3.0) * pow(1.0-pow(ok,4),1.0/6.0) / (1 + pow(ok,2));
	double fc = ke*(rc-latcona/2.0);
		
	double avgfric = fc - df*pow(abs(log(vs)),2.0/3.0);

	return avgfric;
}

double sangsup(double spos)
{
	double p0 = latcona*sqrt(xmass/barr1);
	double pk = sqrt(xmass/spring);
	double ok = p0/(2*pi*pk);
	double ke = spring/(1+pow(ok,2));
	double xc = acos(-pow(ok,2))/latcona2pi;
	double rc = (latcona2pi*xc + 1.0/pow(ok,2)*sin(latcona2pi*xc)) / latcona2pi;
	double vs = 2.0*supvel*xdamp*barr1*pow(p0*ok,2)/(kB*temp*latcona*sqrt(1.0-pow(ok,4)));
   
    double pf = pow(2.0/3.0 * barr1/(kB*temp) * 
                pow(4*pi*rc/latcon, 1.5) * 
                pow(ok,3)/pow(1-pow(ok, 4), 0.25)
                , 2.0/3.0); // prefactor to f in Sang 
        
    double fs = pow(pf,2.0/3.0) * (1 - spos/rc);

    double prob = 1.5*sqrt(fs)/vs * 
                  exp(-pow(fs,1.5)- 1.0/vs - 
                  exp(-pow(fs,1.5)));

	return prob;
    //return 1.0;
}

uint getlabel(vector <double>* targetsx, vector <double>* targetsq, double x, double q)
{
    uint targsize = targetsx->size();
    vector <double> offsetsx;
    vector <double> offsetsq;

    for (uint l = 0; l < targsize; l++)
    {
        double tx = targetsx->at(l);
        double offsetx = abs(x - tx);
        offsetsx.push_back(offsetx/tx);
        //cout << x << " " << tx << endl;
        
        double tq = targetsq->at(l);
        double offsetq = abs(q - tq);
        offsetsq.push_back(offsetq/tq);
    }

    vector <double> norms;
    for (uint l = 0; l < targsize; l++)
    {
        double w1 = 0.9;
        double w2 = offsetsx[l]*(1.0-w1)/offsetsq[l] + 1.0;
        double norm = (w1*offsetsx[l] + w2*offsetsq[l])/2;
        norms.push_back(norm);
    }

    auto minit = min_element(norms.begin(),norms.end());
    uint minel = distance(norms.begin(),minit);

    return minel;
}

int main()
{ 
	// MODEL PARAMTERES SET IN common_stuff.h
	
	double tmp = 0.5;	// 0.5 gives reliable timestep dep. 1.0 should be ok
	double tstep = tmp * 3e-14;	
	uint tsteps = 1.0/tmp * 4e5;	// has to be even beucasue lazyness

	string tfile = "time.csv";
	string xfile = "xout";
	string qfile = "qout";
	string ffile = "tomout.csv";
	string avgfile = "avgs.csv";
	string pfile = "params.dat";
    string sfile = "sliplabels.csv";

	time_t t = time(0);
	struct tm * now  = localtime(&t);
	char buffer [80];
	strftime (buffer,80,"%Y%m%d%H%M%S",now);

    //xfile.append(buffer);
    //qfile.append(buffer);
    
    xfile.append(".csv");
    qfile.append(".csv");

    //vector <double> targetsx = {1.055321,1.281183,1.514248,1.748637,1.782399,1.982108,2.010424,2.036907,2.210648,2.239958,2.258750,2.272320,1.096358,1.289138,1.373806,1.520053,1.647441,1.753444,1.792153,1.869407,1.986442,2.016174,2.046863,2.109832,2.149017,2.215485,2.244193,2.263649,2.372409,2.395439};   // with saddle points
    vector <double> targetsx = {1.055321,1.281183,1.514248,1.748637,1.782399,1.982108,2.010424,2.036907,2.210648,2.239958,2.258750,2.272320};
    transform(targetsx.begin(),targetsx.end(),targetsx.begin(),[](double n){return 1e-9*n;});

    //vector <double> targetsq = {0.768823,0.747456,0.732154,0.717806,0.514427,0.702379,0.493754,0.279007,0.682051,0.474256,0.251936,0.021785,0.804315,0.697276,0.827876,0.676202,0.848242,0.658147,0.434054,0.595167,0.639930,0.406424,0.177963,0.348435,0.623097,0.617785,0.383319,0.143356,0.118994,0.382620i};      // with saddle points
    vector <double> targetsq = {0.768823,0.747456,0.732154,0.717806,0.514427,0.702379,0.493754,0.279007,0.682051,0.474256,0.251936,0.021785};
        transform(targetsq.begin(),targetsq.end(),targetsq.begin(),[](double n){return 1e-9*n;});

    //for (auto &el : targetsx)
    //  cout << spring*(tstep*163224 - el) << endl;

	//uint periods = static_cast <uint> (tstep*tsteps*supvel/latcon);

	// for averaging
	uint skip = 50;			// probably should be a fraction of mean	
	uint halfmeansize = 1000;			// half interval, such that | mean -- mid -- mean |
	
	// for simple slip
        //uint adj = 10;
	uint ttoa = ceil(latcon/(tstep*supvel));	// timesteps to minima
	//uint stride = ttoa / 4.0);
	//uint end = tsteps;
	uint adj = 1;
	uint stride = 50;      // this was 50 efore i started fiddling...
	uint end = round(tsteps/skip)-skip;		// this is a bit risky, size should be constant over many runs though...
	
    uint pauseat = 20;
	uint runs = 1;
		
    //chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    
    //vector <vector <double>> publicslips;
    vector <vector <double>> publicslips1; // 1 indicates using the third method in slicing algs.h, yes the notation is bad, horrific even, still, such is the nature of reality
    vector<pair<double,uint>> publicslips2;

    #pragma omp parallel
    { 
	    //vector <vector <double>> privateslips;
        vector <vector <double>> privateslips1;
        vector<pair<double,uint>> privateslips2;
        
        #pragma omp for nowait
	    for ( uint l = 0; l < runs; l++)
	    {	 
	        vector <uint> slips;
	        vector <uint> sliptos;
	        vector <uint> fslips;
	        vector <uint> qslips;
	    	
            tomlin afm(spring,supvel,latcon,align,barr1,barr2,kappa1,kappa2, nu2,nu4,temp,tstep,tsteps,xmass,qmass,xdamp,qdamp, tfile,xfile,qfile,ffile,avgfile,pfile,sfile);
            
            // incom
            //afm.setposx(1.23912011781616e-09);
	        //afm.setvelx(0.9995276370265342);
	        //afm.setaccx(18619338513200.46);
	        //afm.setposq(6.874289088010818e-10);
	        //afm.setvelq(1.061307759460043);
	        //afm.setaccq(45478933617766.52);
            //afm.setsuppos(2.443094999995714e-09);
            //uint startt = 157832;
	    
            // com
            afm.setposx(1.030613785802745e-09);
	        afm.setvelx(0.9995431014034264);
	        afm.setaccx(18645712950047.95);
	        afm.setposq(7.422256383340498e-10);
	        afm.setvelq(0.5066323077421474);
	        afm.setaccq(21750153167849.19);
            afm.setsuppos(2.448359999995651e-09);
            uint startt = 163224;

            afm.settargts(&targetsx, &targetsq);
            
	    	for ( uint k = 0; k < tsteps; k++ )
	    	{
	    		//cout << "now begnning loop: " << k << endl;

                //if (k == pauseat*ttoa)
                //{
                //    afm.tpause();
                //}

	    		afm.rk4();
	    		afm.calcfric();
	    				
	    		// multiple algorithms available, see separate .h file
	    		
	    		afm.pushvals();			
	    		afm.inctime();
	    	}

	    	afm.noisered(halfmeansize,skip);	
	    	halfintervals(3, adj, end, stride,round(pauseat*ttoa/stride), afm.getrntimes(), afm.getrnqposs(), afm.getrnfrics(), &fslips, &qslips, &slips, &sliptos);
	    	afm.findsortslips();    // use this for integrated slip detection	
	        afm.writedata();		// ONLY FOR DIAGNOSTICS REMOVE LATER (it won't make sense)
	    	//t1 = std::chrono::high_resolution_clock::now();
	    
	    	// this rn business is fooken ugly...	
	    	//halfintervals(adj, end, stride, afm.gettimes(),afm.getfrics(), &slips);
	    	//halfintervals(2, adj, end, stride,round(pauseat*ttoa/stride), afm.getrntimes(),afm.getrnfrics(), &slips);

	    	//for (auto &el : slips)
	    	//{
	    	//	//if(afm.getrntime(el) > pauseat*ttoa*tstep)
            //    //{
            //        //fspos << afm.getrntime(el) << "," << afm.getrnfric(el) << endl;
	    	//	    fspos << afm.getrntime(el) << "," << afm.getrnfric(el) << "," << afm.getrnsuppos(el) << endl;		// rn is for reduced noise btw, you'll thank me later furure me, thank you passed me // future me
            //    //}
	    	//}
            
            //printvectoruint(&slips);

            //if (slips.size() > 1)
            //{
            //    for (uint k = 0; k < slips.size() - 1; k++)
            //    {
            //        uint el1 = slips[k];
            //        uint el2 = slips[k+1];
            //        double diff = afm.getrntime(el2) - afm.getrntime(el1);  
                    //cout << el1 << " " << el2 << " " << diff << endl;
	    	//        privatedecay.push_back(diff);		// rn is for reduced noise btw, you'll thank me later furure me, thank you passed me // future me
            //    }
            //}
            //else if (slips.size() == 1)
            //    fspos << slips[0] << endl;
            //else 
            //    fspos << 0 << endl;

	    	//t2 = std::chrono::high_resolution_clock::now();
	    	
	    	//double slicetime = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

	    	//cout << "loop time: " << looptime << endl << "slicetime " << slicetime << endl;
	    	
	    	//if ((l+1) % 5 == 0)
	    	//{
	    	//	cout << "finished " << l+1 << " out of " << runs <<  " iterations" << endl;
	    	//}
            
		    for (auto& el : slips)
            {    
                double t = afm.getrntime(el);
                double q = afm.getrnqpos(el);
                double f = afm.getrnfric(el);

                vector <double> tmp = {t,f,q};

                privateslips1.push_back(tmp);
            }
            for (uint k = 0; k < afm.getsliptimes()->size(); k++)
            {
                pair <double, uint> tmp = {afm.getsliptimes()->at(k), afm.getsliplabels()->at(k)};  // could have implemented those accessor in a more appealing way...

                privateslips2.push_back(tmp);
            }
            
            //privateslips.push_back(slipdata);
            //privateslips.push_back(sliptimes);
	    }
        #pragma omp critical
        publicslips1.insert(publicslips1.end(), privateslips1.begin(), privateslips1.end());
        publicslips2.insert(publicslips2.end(), privateslips2.begin(), privateslips2.end());
    } 

	//chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

	//double looptime = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    //cout << looptime << endl;
    
    //fspos << "logging slipping events, each line is a realization, data comes in order: friction, time, x and then q position, velocity and acceleration (x.pos, x.vel, ..., q.acc)" << endl; 

    //for (auto& vec : publicslips)
    //{
    //    for (auto& el : vec)
    //    {
    //        //fspos << publictimes[k] << "," << publicslips[k] << endl;
    //        fspos << el << ",";
    //    }
    //    fspos << endl;
    //}
    
    ofstream fsslips;
    fsslips.open("slips.csv");
    
    for (uint k = 0; k < publicslips1.size(); k++)
    { 
        double x = publicslips1[k][1]/spring + supvel*publicslips1[k][0];
        double q = publicslips1[k][2];

        uint label = getlabel(&targetsx,&targetsq,x,q); // not working at the moment

        fsslips << publicslips1[k][0] << "," 
                << publicslips1[k][1] << "," 
                << publicslips1[k][2] << ","
                << label << endl;
    }
    fsslips.close();
    
    ofstream fsslips2;
    fsslips.open("sliplabels.csv");
    
    for (uint k = 0; k < publicslips2.size(); k++)
    { 
        fsslips << publicslips2[k].first << "," 
                << publicslips2[k].second << endl;
    }
    fsslips2.close();

	double ending = 2e-9;
	uint diststeps = 1000;
	double diststep = ending / diststeps;;
	
}


