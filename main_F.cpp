#include <iostream>
#include <fstream>
#include <stdio.h>
#include <unistd.h>
#include <vector>               // std:: vector
//#include <array>                // std::array
#include <cmath>                // sin, cos etc.
#include <complex>
#include <string>
#include <sstream>              // for input/ouput with strings
#include <algorithm>            // std::min_element
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <fftw3.h>              // FFTW
#include <omp.h>                        // OpenMP

#include "boson_gpe_F.hpp"



using namespace std;



#ifndef M_PI
#define M_PI        3.14159265358979323846264338327950288   // pi
#endif


// To make more clear what components of the spin appear
#ifndef X
#define X 0
#endif

#ifndef Y
#define Y 1
#endif

#ifndef Z
#define Z 2
#endif



// Define global variables of list of variables in boson_gpe.hpp
#define VAR(aa,bb,cc) bb aa;
varlist
#undef VAR
char folder[1024];
//



#ifndef TYPE_DEF
#define TYPE_DEF
typedef double my_double;
typedef complex<my_double> comp;
typedef vector< comp > state_type;  // The type of container used to hold the state, vector of spins.
#endif






// Prints execution time in hours, minutes and seconds
void print_time (int t)
{
    int hours = (int)t/3600;
    int minutes = (int)(t-3600*hours)/60;
    int seconds = (int)(t-3600*hours-60*minutes);
    printf("Execution time: %i hours, %i minutes, %i seconds\n",hours,minutes,seconds);
}//









int main(int argc, char* argv[])
{
    
    long int runtime = (long)time(NULL);
    
// Initialises global variables with values as given in command line
    
    if (argc==varnum+1)
    {
        #define STR(aa,bb) snprintf(aa,1024,argv[bb]);
        #define VAR(aa,bb,cc) aa=(bb)atof(argv[cc]);
        varlist
        strlist
        #undef VAR
        #undef STR
    }
    else { cout << "Error: Missing variable" << endl; return 0;}
    
    char discr[1024];      // Append to file name
    //snprintf(discr,1024,"");       // In case nothing to append
    snprintf(discr,1024,"_N%i_dt%g_h%g",N,dt,h);
    
// OMP and FFTW settings
    
    int numprocs = omp_get_num_procs();
    cout << "Number of threads available: " << numprocs
    << endl;
    
    if (tnum!=0&&tnum<numprocs) {
        numprocs = tnum;
    }
    omp_set_num_threads(numprocs);
    cout << "Using OMP with " << numprocs << " threads\n";
    
    if (!fftw_init_threads()) cout << "Error with fftw_thread\n";
    fftw_plan_with_nthreads(numprocs);
    
    
    
    
// Initialize random number generator
    
    long int random_seed = seed;
    if (seed<0) random_seed = (long)time(NULL)+1234*(long)getpid();
    //long int random_seed2 = (rdtsc() % (long)pow(2,30)); //different type of random seed, use this for different seeds within same program
    
    gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set (r,random_seed);
    
    
    
    
// Containers for outputing CR and F data, and DPF and DPrho data
    
    vector<double> timesCR;
    vector< vector<double> > outCR;
    vector< vector< vector<comp> > > outCRpdep(Nbins);
    vector< vector<double> > physMomenta(Nbins);    // Holds physical momenta and number of momenta per bin
    
    vector<double> timesF;
    vector< vector< vector<comp> > > outFpdep(Nbins);
	
    vector<double> timesDPrho;
    vector< vector< vector<comp> > > outDPrhopdep(Nbins);
	
    vector<double> timesDPF;
    vector< vector< vector<comp> > > outDPFpdep(Nbins);
    
    
    
// Settings
    
    vector<double> times;                       // put time results here

    
    
    // Integrate with split-step method and iterate. Print separate files for different iterations
    
    for (int i=0; i<iter; i++)
    {
        
        int filenumber = nfile*iter+i;
        
        // Initialize state and add noise
        boson_gpe bosonGPE(r);         // initialize gpe class
        bosonGPE.addTWAnoise();             // Add TWA noise to initial conditions
        
        
        // Export physical momenta
        if (cphysMom==1 && i==0)
        {
            char buffer_mom[bsize];
            snprintf(buffer_mom,1024,"%s/physMom_Nbins%i_V%i.txt",folder,Nbins,N);
            ofstream output_mom(buffer_mom);
            output_mom.precision(12);
            
            if (output_mom.is_open())
            {
                output_mom << "# Momentum Bin Number - Physical Momentum - Momenta per bin\n";
                
                int binfullcount=0;
                vector<double> momtemp;
                for (int n=0; n<Nbins; n++)
                {
                    momtemp = bosonGPE.getPhysMom(n);
                    if (momtemp[1]!=0) { output_mom << binfullcount << '\t' << momtemp[0] << '\t' << momtemp[1] << endl;  binfullcount++;}
                }
            }
            
            output_mom.close();
        }
        
        
        // Export initial sp
        if (spstep>=0)
        {
            char buffer_sp0[bsize];
            snprintf(buffer_sp0,1024,"%s/sp_gpe%s_%06.f_%i.txt",folder,discr,0.0,filenumber);
            ofstream output_sp0(buffer_sp0);
            output_sp0.precision(12);
        
            if (output_sp0.is_open())
            {
                bosonGPE.write_spec(output_sp0);      // print initial spectrum
            }
        
            output_sp0.close();
        }
        
        
        // Export intial op
        ofstream output_op;
        
        if (opstep>=0)
        {
            char buffer_op[bsize];
            snprintf(buffer_op,1024,"%s/op_gpe%s_%i.txt",folder,discr,filenumber);
            output_op.open(buffer_op);
            output_op.precision(12);
                
            if (output_op.is_open())
            {
                output_op << "# Time - Magnetization - Total N - Kinetic Energy - Potential Energy\n";
                    
                #define VAR(aa,bb,cc) << "# " << #aa << "=" << aa << "\n"
                output_op varlist;
                #undef VAR          // This adds to the op file a list of the parameters used for the simulation
                    
                if (opstep>=0)
                {
                    bosonGPE.write_op(output_op,0);       // print magnetization, energy, etc at initial time
                }
            }
        }
		
        // Export initial DPspec
        if (DPspstep>=0)
        {
            char buffer_DPsp0[bsize];
            snprintf(buffer_DPsp0,1024,"%s/DPsp_gpe%s_%06.f_%i.txt",folder,discr,0.0,filenumber);
            ofstream output_DPsp0(buffer_DPsp0);
            output_DPsp0.precision(12);
        
            if (output_DPsp0.is_open())
            {
				output_DPsp0 << "# Momentum - Phase-Phase - Dichte-Dichte - Re(Dichte-Phase) - Im(Dichte-phase)\n";
				
                bosonGPE.write_DPspec(output_DPsp0);      // print initial density-phase spectrum
            }
        
            output_DPsp0.close();
        }
        
        
        // Save initial phi(p=0)
        if (phi0step>=0)
        {
            bosonGPE.save_phi0t(0);     // saves phi(t,p=0) for calculation of unequal time correlators
        }
        
        
        // Save physical momenta for CR and F output
        if (i==0 && ( ( CRmodeResp==2 || CRmodeResp==3 ) || ( Fstep>0 && twF>=0 ) || ( DPFstep>0 && twDPF>=0 ) || ( DPrhostep>0 && twDPrho>=0 ) ) )
        {
            for (int n=0; n<Nbins; n++) physMomenta[n] = bosonGPE.getPhysMom(n);
        }
        
        
        
        //////
        //      Integrate GPE
        //////
        double next_sp = spstep;
        double next_op = opstep;
        double next_phi0 = phi0step;
        double next_cr = tw;
        double next_F = twF;
		double next_DPsp = DPspstep;
		double next_DPrho = twDPrho;
		double next_DPF = twDPF;
        int with_h = 0;     // If 0: system not triplicated yet, 1: yes
		int with_hDP = 0;     // If 0: system not triplicated yet, 1: yes (phase-dichte)
        int F_saved = 0;    // If 0: phi(twF) not saved yet, 1: yes
		int DPF_saved = 0;    // If 0: phase(twDPF) and dichte(twDPF) not saved yet, 1: yes
        int CRCount = 0;    // Counts number of times CR has been saved within a given iteration
        int FCount = 0;    // Counts number of times F has been saved within a given iteration
		int DPrhoCount = 0;    // Counts number of times dichte-phase rho has been saved within a given iteration
		int DPFCount = 0;    // Counts number of times dichte-phase F has been saved within a given iteration
        
        // In case tw=0, triple system and output at t=0
        if (0>=tw-0.01*dt && tw>=0 && with_h!=1)
        {
            bosonGPE.triple_system();
            bosonGPE.rotation_hR();
            bosonGPE.rotation_hI();
            with_h = 1;
        }// Triple system
        if (0>=next_cr-0.01*dt && CRstep>=0 && tw>=0)
        {
            next_cr += CRstep;
            
            if (i!=0)
            {
                if (CRmodeResp==0 || CRmodeResp==1)
                {
                    vector<double> temp(bosonGPE.save_CR());
                    for (int k=0; k<6; k++) outCR[CRCount][k] += temp[k];
                }
                else if (CRmodeResp==2 || CRmodeResp==3)
                {
                    vector< vector<comp> > tempp(bosonGPE.save_CR_pdep());
                    for (int n=0; n<Nbins; n++) { outCRpdep[n][CRCount][0] += tempp[n][0]; outCRpdep[n][CRCount][1] += tempp[n][1]; }
                }
            }
            else
            {
                if (CRmodeResp==0 || CRmodeResp==1)
                {
                    outCR.push_back(bosonGPE.save_CR());
                    timesCR.push_back(0);     // Save times for first iteration
                }
                else if (CRmodeResp==2 || CRmodeResp==3)
                {
                    vector< vector<comp> > tempp(bosonGPE.save_CR_pdep());
                    for (int n=0; n<Nbins; n++) outCRpdep[n].push_back(tempp[n]);
                    timesCR.push_back(0);     // Save times for first iteration
                }
            }
            
            CRCount++;
        }// Save CR
        
        
        
        // In case twF=0, save state of system and output at t=0
        if (0>=twF-0.01*dt && twF>=0 && F_saved!=1)
        {
            bosonGPE.save_phi_at_twF();
            F_saved = 1;
        }
        if (0>=next_F-0.01*dt && Fstep>=0 && twF>=0)
        {
            next_F += Fstep;
            
            if (i!=0)
            {
                vector< vector<comp> > tempp(bosonGPE.save_F_pdep());
                for (int n=0; n<Nbins; n++) { outFpdep[n][FCount][0] += tempp[n][0]; outFpdep[n][FCount][1] += tempp[n][1]; }
            }
            else
            {
                vector< vector<comp> > tempp(bosonGPE.save_F_pdep());
                for (int n=0; n<Nbins; n++) outFpdep[n].push_back(tempp[n]);
                timesF.push_back(0);     // Save times for first iteration
            }
            
            FCount++;
        }// Save F
		
		
		
        // In case twDPrho=0, triple system and output at t=0
        if (0>=twDPrho-0.01*dt && twDPrho>=0 && with_hDP!=1)
        {
            bosonGPE.triple_system_DP();
            bosonGPE.rotation_hD();
            bosonGPE.rotation_hP();
            with_hDP = 1;
        }// Triple system
        if (0>=next_DPrho-0.01*dt && DPrhostep>=0 && twDPrho>=0)
        {
            next_DPrho += DPrhostep;
            
            if (i!=0)
            {
                vector< vector<comp> > tempp(bosonGPE.save_DPrho_pdep());
                for (int n=0; n<Nbins; n++) { for (int c=0; c<outDPrhopdep[n][DPrhoCount].size(); c++) outDPrhopdep[n][DPrhoCount][c] += tempp[n][c]; }
            }
            else
            {
                vector< vector<comp> > tempp(bosonGPE.save_DPrho_pdep());
                for (int n=0; n<Nbins; n++) outDPrhopdep[n].push_back(tempp[n]);
                timesDPrho.push_back(0);     // Save times for first iteration
            }
            
            DPrhoCount++;
        }// Save DPrho
		
		
		
        // In case twDPF==0, save dichte-phase of system and output at t=0
        if (0>=twDPF-0.01*dt && twDPF>=0 && DPF_saved!=1)
        {
            bosonGPE.save_DP_at_twDPF();
            DPF_saved = 1;
        }
        if (0>=next_DPF-0.01*dt && DPFstep>=0 && twDPF>=0)
        {
            next_DPF += DPFstep;
            
            if (i!=0)
            {
                vector< vector<comp> > tempp(bosonGPE.save_DPF_pdep());
                for (int n=0; n<Nbins; n++) { outDPFpdep[n][DPFCount][0] += tempp[n][0]; outDPFpdep[n][DPFCount][1] += tempp[n][1]; outDPFpdep[n][DPFCount][2] += tempp[n][2]; }
            }
            else
            {
                vector< vector<comp> > tempp(bosonGPE.save_DPF_pdep());
                for (int n=0; n<Nbins; n++) outDPFpdep[n].push_back(tempp[n]);
                timesDPF.push_back(0);     // Save times for first iteration
            }
            
            DPFCount++;
        }// Save DPF
        

        
        // Time evolution
        for (double t=dt; t<=tmax; t+=dt)
        {
            // Dynamics
            bosonGPE.dynamics();    // calculates phi(t+dt)
            if (with_h==1) { bosonGPE.dynamics_phihR(); bosonGPE.dynamics_phihI(); }      // calculates phi_h(t+dt)
			if (with_hDP==1) { bosonGPE.dynamics_phihD(); bosonGPE.dynamics_phihP(); }      // calculates phi_hDP(t+dt)
            
            
            // Triple system with external field and rotate with h
            
            if (t>=tw-0.01*dt && tw>=0 && with_h!=1)
            {
                bosonGPE.triple_system();
                bosonGPE.rotation_hR();
                bosonGPE.rotation_hI();
                with_h = 1;
            }
            
            if (t>=twF-0.01*dt && twF>=0 && F_saved!=1)
            {
                bosonGPE.save_phi_at_twF();
                F_saved = 1;
            }
			
            if (t>=twDPrho-0.01*dt && twDPrho>=0 && with_hDP!=1)
            {
                bosonGPE.triple_system_DP();
                bosonGPE.rotation_hD();
                bosonGPE.rotation_hP();
                with_hDP = 1;
            }
			
            if (t>=twDPF-0.01*dt && twDPF>=0 && DPF_saved!=1)
            {
                bosonGPE.save_DP_at_twDPF();
                DPF_saved = 1;
            }
            
			
            
            // Output
            
            if (t>=next_phi0-0.01*dt && phi0step>=0)   // substract -0.01 in order to avoid:  1 != 0.5+0.5
            {
                next_phi0 += phi0step;
                bosonGPE.save_phi0t(t);     // saves phi(t,p=0) for calculation of unequal time correlators
            }
            
            if (t>=next_op-0.01*dt && opstep>=0)
            {
                next_op += opstep;
                
                if (output_op.is_open())
                {
                    bosonGPE.write_op(output_op,t);       // print magnetization, energy etc at endtime of integration
                }
            }
        
            if (t>=next_sp-0.01*dt && spstep>=0)
            {
                next_sp += spstep;
                
                char buffer_sp[bsize];
                snprintf(buffer_sp,1024,"%s/sp_gpe%s_%06.f_%i.txt",folder,discr,t,filenumber);
                ofstream output_sp(buffer_sp);
                output_sp.precision(12);
                
                if (output_sp.is_open())
                {
                    bosonGPE.write_spec(output_sp);       // print spectrum at endtime of integration
                }
                
                output_sp.close();
            }
			
            if (t>=next_DPsp-0.01*dt && DPspstep>=0)
            {
                next_DPsp += DPspstep;
                
                char buffer_DPsp[bsize];
                snprintf(buffer_DPsp,1024,"%s/DPsp_gpe%s_%06.f_%i.txt",folder,discr,t,filenumber);
                ofstream output_DPsp(buffer_DPsp);
                output_DPsp.precision(12);
                
                if (output_DPsp.is_open())
                {
                    bosonGPE.write_DPspec(output_DPsp);       // print density-phase spectrum at endtime of integration
                }
                
                output_DPsp.close();
            }
            
            
            // Output rho
            
            if (t>=next_cr-0.01*dt && CRstep>=0 && tw>=0)
            {
                next_cr += CRstep;
                
                if (i!=0)
                {
                    if (CRmodeResp==0 || CRmodeResp==1)
                    {
                        vector<double> temp(bosonGPE.save_CR());
                        for (int k=0; k<6; k++) outCR[CRCount][k] += temp[k];
                    }
                    else if (CRmodeResp==2 || CRmodeResp==3)
                    {
                        vector< vector<comp> > tempp(bosonGPE.save_CR_pdep());
                        for (int n=0; n<Nbins; n++) { outCRpdep[n][CRCount][0] += tempp[n][0]; outCRpdep[n][CRCount][1] += tempp[n][1]; }
                    }
                }
                else
                {
                    if (CRmodeResp==0 || CRmodeResp==1)
                    {
                        outCR.push_back(bosonGPE.save_CR());
                        timesCR.push_back(t);     // Save times for first iteration
                    }
                    else if (CRmodeResp==2 || CRmodeResp==3)
                    {
                        vector< vector<comp> > tempp(bosonGPE.save_CR_pdep());
                        for (int n=0; n<Nbins; n++) outCRpdep[n].push_back(tempp[n]);
                        timesCR.push_back(t);     // Save times for first iteration
                    }
                }
                
                
                
                CRCount++;
            }
            
            
            
            // Output F
            
            if (t>=next_F-0.01*dt && Fstep>=0 && twF>=0)
            {
                next_F += Fstep;
                
                if (i!=0)
                {
                    vector< vector<comp> > tempp(bosonGPE.save_F_pdep());
                    for (int n=0; n<Nbins; n++) { outFpdep[n][FCount][0] += tempp[n][0]; outFpdep[n][FCount][1] += tempp[n][1]; }
                }
                else
                {
                    vector< vector<comp> > tempp(bosonGPE.save_F_pdep());
                    for (int n=0; n<Nbins; n++) outFpdep[n].push_back(tempp[n]);
                    timesF.push_back(t);     // Save times for first iteration
                }
                
                FCount++;
            }
			
			
			
			
            // Output DPrho
			
	        if (t>=next_DPrho-0.01*dt && DPrhostep>=0 && twDPrho>=0)
	        {
	            next_DPrho += DPrhostep;
            
	            if (i!=0)
	            {
	                vector< vector<comp> > tempp(bosonGPE.save_DPrho_pdep());
	                for (int n=0; n<Nbins; n++) { for (int c=0; c<outDPrhopdep[n][DPrhoCount].size(); c++) outDPrhopdep[n][DPrhoCount][c] += tempp[n][c]; }
	            }
	            else
	            {
	                vector< vector<comp> > tempp(bosonGPE.save_DPrho_pdep());
	                for (int n=0; n<Nbins; n++) outDPrhopdep[n].push_back(tempp[n]);
	                timesDPrho.push_back(t);     // Save times for first iteration
	            }
            
	            DPrhoCount++;
	        }
			
			
			
            // Output DPF
            
            if (t>=next_DPF-0.01*dt && DPFstep>=0 && twDPF>=0)
            {
                next_DPF += DPFstep;
                
                if (i!=0)
                {
                    vector< vector<comp> > tempp(bosonGPE.save_DPF_pdep());
                    for (int n=0; n<Nbins; n++) { outDPFpdep[n][DPFCount][0] += tempp[n][0]; outDPFpdep[n][DPFCount][1] += tempp[n][1]; outDPFpdep[n][DPFCount][2] += tempp[n][2]; }
                }
                else
                {
                    vector< vector<comp> > tempp(bosonGPE.save_DPF_pdep());
                    for (int n=0; n<Nbins; n++) outDPFpdep[n].push_back(tempp[n]);
                    timesDPF.push_back(t);     // Save times for first iteration
                }
                
                DPFCount++;
            }
            
            
            
            
            
            
            
        }//End of time evolution
        
        
        
        if (timesCR.size() != outCR.size() &&  ( CRmodeResp==0 || CRmodeResp==1 ) ) cout << "Sizes of timesCR and outCR are different.\n";
        if (timesCR.size() != outCRpdep[0].size() && ( CRmodeResp==2 || CRmodeResp==3 ) ) cout << "Sizes of timesCR and outCRpdep[0] are different.\n";
        if (timesF.size() != outFpdep[0].size() && twF>=0 && Fstep>=0 ) cout << "Sizes of timesF and outFpdep[0] are different.\n";
		if (timesDPrho.size() != outDPrhopdep[0].size() && twDPrho>=0 && DPrhostep>=0 ) cout << "Sizes of timesDPrho and outDPrhopdep[0] are different.\n";
		if (timesDPF.size() != outDPFpdep[0].size() && twDPF>=0 && DPFstep>=0 ) cout << "Sizes of timesDPF and outDPFpdep[0] are different.\n";

        if (output_op.is_open())
        {
            output_op.close();
        }
        
        
        // Print phi(t,p=0)
        if (phi0step>=0)
        {
            char buffer_phi[bsize];
            snprintf(buffer_phi,1024,"%s/phi0%s_%i.txt",folder,discr,filenumber);
            ofstream output_phi(buffer_phi);
            
            if (output_phi.is_open())
            {
                output_phi << "# Time - real(phi(t,p=0)) - imag(phi(t,p=0)) - |phi(t,p=0)|^2\n";
                bosonGPE.write_phi0t(output_phi);
            }
            
            output_phi.close();
        }
        
        
    }//End of iterations
    
    
    // Normalize CR output and output it
    
    if (tw>=0 && CRstep>=0)
    {
        // Normalize
        double Nfac = 1.0/iter;
        if (CRmodeResp==0 || CRmodeResp==1)
        {
            for (int k=0; k<timesCR.size(); k++)
            {
                if (outCR[0].size() != 6) cout << "Problem with size of CR output vector" << endl;
                for (int m=0; m<6; m++)
                {
                    outCR[k][m] *= Nfac;
                }
            }
        }
        else if (CRmodeResp==2 || CRmodeResp==3)
        {
            for (int n=0; n<Nbins; n++)
            {
                for (int k=0; k<timesCR.size(); k++)
                {
                    if (outCRpdep[n][0].size() != 2) cout << "Problem with size of CRpdep output vector" << endl;
                    outCRpdep[n][k][0] *= Nfac;
                    outCRpdep[n][k][1] *= Nfac;
                }
            }
        }
        
        // Naming
        char modeOut[bsize];
        switch (CRmodeResp) {
                case 0:
                snprintf(modeOut,1024,"_0mode");
                break;
                
                case 1:
                snprintf(modeOut,1024,"_Integrated");
                break;
                
                case 2:
                snprintf(modeOut,1024,"_Pdep");
                break;
                
                case 3:
                snprintf(modeOut,1024,"_PdepUni");
                break;
                
            default:
                cout << "Error: Wrong type of mode to output for linear response.\n";
                break;
        }
        
        // Output
        if (CRmodeResp==0 || CRmodeResp==1)
        {
            char buffer_cr[bsize];
            snprintf(buffer_cr,1024,"%s/CR_gpe%s%s_tw%06.f_runs%i_file%i.txt",folder,modeOut,discr,tw,iter,nfile);
            ofstream output_cr(buffer_cr);
            output_cr.precision(12);
            
            if (output_cr.is_open())
            {
                output_cr << "# Waiting time tw - Time t - < phi(k=0) >_hR real - imag - < phi(k=0) >_hI real - imag - < phi(k=0) > real - imag \n";
                
                #define VAR(aa,bb,cc) << "# " << #aa << "=" << aa << "\n"
                output_cr varlist;
                #undef VAR
                
                for (int i=0; i<timesCR.size(); i++)
                {
                    output_cr << tw << '\t' << timesCR[i] << '\t' << outCR[i][0] << '\t' << outCR[i][1] << '\t' << outCR[i][2] << '\t' << outCR[i][3]
                    << '\t' << outCR[i][4] << '\t' << outCR[i][5] << endl;
                }
            }
            output_cr.close();
        }
        else if (CRmodeResp==2 || CRmodeResp==3)
        {
            int binfullcount=0;
            for (int n=0; n<Nbins; n++)
            {
                if (physMomenta[n][1]!=0 && ( physMomenta[n][0]<=Pcutoff || Pcutoff<0 ) ) // if bin not empty and momentum smaller than Pcutoff (careful: bin_momentum != phys_momentum)
                {
                    char buffer_cr[bsize];
                    snprintf(buffer_cr,1024,"%s/CR_gpe%s%s_p%i_tw%06.f_runs%i_file%i.txt",folder,modeOut,discr,binfullcount,tw,iter,nfile);
                    ofstream output_cr(buffer_cr);
                    output_cr.precision(12);
                    
                    if (output_cr.is_open())
                    {
                        //output_cr << "# Waiting time tw - Time t - Re( < i[phi(t,p),phi^dagger(tw,p)] >) - Im(.) - Re( < i[phi(t,p),phi(tw,p)] >) - Im(.) \n";
                        output_cr << "# Time (t-t_w) - Re( < i[phi(t,p),phi^dagger(tw,p)] >) - Im(.) \n";
                        
                        #define VAR(aa,bb,cc) << "# " << #aa << "=" << aa << "\n"
                        output_cr varlist;
                        #undef VAR
                        
                        for (int i=0; i<timesCR.size(); i++)
                        {
                            //output_cr << tw << '\t' << timesCR[i] << '\t' << outCRpdep[n][i][0].real() << '\t' << outCRpdep[n][i][0].imag()
                            //<< '\t' << outCRpdep[n][i][1].real() << '\t' << outCRpdep[n][i][1].imag() << endl;
                            output_cr << timesCR[i]-tw << '\t' << outCRpdep[n][i][0].real() << '\t' << outCRpdep[n][i][0].imag() << endl;
                        }
                    }
                    output_cr.close();
                    
                    binfullcount++;
                }
            }
        }
    }
    
    
    
    
    
    
    // Normalize F output and output it
    
    if (twF>=0 && Fstep>=0)
    {
        // Normalize
        double Nfac = 1.0/iter;
        for (int n=0; n<Nbins; n++)
        {
            for (int k=0; k<timesF.size(); k++)
            {
                if (outFpdep[n][0].size() != 2) cout << "Problem with size of Fpdep output vector" << endl;
                outFpdep[n][k][0] *= Nfac;
                outFpdep[n][k][1] *= Nfac;
            }
        }
        
        // Output
        int binfullcount=0;
        for (int n=0; n<Nbins; n++)
        {
            if (physMomenta[n][1]!=0 && ( physMomenta[n][0]<=Pcutoff || Pcutoff<0 ) ) // if bin not empty and momentum smaller than Pcutoff (careful: bin_momentum != phys_momentum)
            {
                char buffer_F[bsize];
                snprintf(buffer_F,1024,"%s/F_gpe%s_p%i_tw%06.f_runs%i_file%i.txt",folder,discr,binfullcount,twF,iter,nfile);
                ofstream output_F(buffer_F);
                output_F.precision(12);
                
                if (output_F.is_open())
                {
                    output_F << "# Time (t-twF) - Re( < phi(t,p) phi^dagger(twF,p) >) - Im(.) \n";
                    
                    #define VAR(aa,bb,cc) << "# " << #aa << "=" << aa << "\n"
                    output_F varlist;
                    #undef VAR
                    
                    for (int i=0; i<timesF.size(); i++)
                    {
                        output_F << timesF[i]-twF << '\t' << outFpdep[n][i][0].real() << '\t' << outFpdep[n][i][0].imag() << endl;
                    }
                }
                output_F.close();
                
                binfullcount++;
            }
        }
        
        
    }
	
	
	
	
	
	
	
    // Normalize DPrho output and output it
    
    if (twDPrho>=0 && DPrhostep>=0)
    {
        // Normalize
        double Nfac = 1.0/iter;
        for (int n=0; n<Nbins; n++)
        {
            for (int k=0; k<timesDPrho.size(); k++)
            {
				for (int c=0; c<outDPrhopdep[n][k].size(); c++)
				{
					outDPrhopdep[n][k][c] *= Nfac;
				}
            }
        }
        
        // Output
        int binfullcount=0;
        for (int n=0; n<Nbins; n++)
        {
            if (physMomenta[n][1]!=0 && ( physMomenta[n][0]<=Pcutoff || Pcutoff<0 ) ) // if bin not empty and momentum smaller than Pcutoff (careful: bin_momentum != phys_momentum)
            {
                char buffer_DPrho[bsize];
                snprintf(buffer_DPrho,1024,"%s/DPrho_gpe_N%i_dt%g_hDP%g_p%i_tw%06.f_runs%i_file%i.txt",folder,N,dt,hDP,binfullcount,twDPrho,iter,nfile);
                ofstream output_DPrho(buffer_DPrho);
                output_DPrho.precision(12);
                
                if (output_DPrho.is_open())
                {
                    output_DPrho << "# Time (t-t_w) - Re( < i[phase,phase^dagger] >) - Im(.) - Re( < i[dichte,dichte^dagger] >) - Im(.) - Re( < i[phase,dichte^dagger] >) - Im(.) - Re( < i[dichte,phase^dagger] >) - Im(.) \n";
                    
                    #define VAR(aa,bb,cc) << "# " << #aa << "=" << aa << "\n"
                    output_DPrho varlist;
                    #undef VAR
                    
                    for (int i=0; i<timesDPrho.size(); i++)
                    {
                        output_DPrho << timesDPrho[i]-twDPrho;
						for (int c=0; c<outDPrhopdep[n][i].size(); c++)
						{
							output_DPrho << '\t' << outDPrhopdep[n][i][c].real() << '\t' << outDPrhopdep[n][i][c].imag();
						}
						output_DPrho << endl;
                    }
                }
                output_DPrho.close();
                
                binfullcount++;
            }
        }
    }
	
	
	
	
	
	
    // Normalize DPF output and output it
    
    if (twDPF>=0 && DPFstep>=0)
    {
        // Normalize
        double Nfac = 1.0/iter;
        for (int n=0; n<Nbins; n++)
        {
            for (int k=0; k<timesDPF.size(); k++)
            {
                if (outDPFpdep[n][0].size() != 3) cout << "Problem with size of DPFpdep output vector" << endl;
                outDPFpdep[n][k][0] *= Nfac;
                outDPFpdep[n][k][1] *= Nfac;
				outDPFpdep[n][k][2] *= Nfac;
            }
        }
        
        // Output
        int binfullcount=0;
        for (int n=0; n<Nbins; n++)
        {
            if (physMomenta[n][1]!=0 && ( physMomenta[n][0]<=Pcutoff || Pcutoff<0 ) ) // if bin not empty and momentum smaller than Pcutoff (careful: bin_momentum != phys_momentum)
            {
                char buffer_DPF[bsize];
                snprintf(buffer_DPF,1024,"%s/DPF_gpe%s_p%i_tw%06.f_runs%i_file%i.txt",folder,discr,binfullcount,twDPF,iter,nfile);
                ofstream output_DPF(buffer_DPF);
                output_DPF.precision(12);
                
                if (output_DPF.is_open())
                {
                    output_DPF << "# Time (t-twF) - Re( < phase(t,p) phase^dagger(twDPF,p) >) - Im(.) - Re( < dichte(t,p) dichte^dagger(twDPF,p) >) - Im(.) - Re( < phase(t,p) dichte^dagger(twDPF,p) >) - Im(.) \n";
                    
                    #define VAR(aa,bb,cc) << "# " << #aa << "=" << aa << "\n"
                    output_DPF varlist;
                    #undef VAR
                    
                    for (int i=0; i<timesDPF.size(); i++)
                    {
                        output_DPF << timesDPF[i]-twDPF << '\t' << outDPFpdep[n][i][0].real() << '\t' << outDPFpdep[n][i][0].imag()
								 						<< '\t' << outDPFpdep[n][i][1].real() << '\t' << outDPFpdep[n][i][1].imag()
														<< '\t' << outDPFpdep[n][i][2].real() << '\t' << outDPFpdep[n][i][2].imag() << endl;
                    }
                }
                output_DPF.close();
                
                binfullcount++;
            }
        }
        
        
    }
    
    
    
    
    
    
    

    
    runtime = (long)time(NULL) - runtime;       // Execution time in seconds
    print_time ((int)runtime);                   // Print execution time in hours, minutes, seconds
    
    gsl_rng_free(r);
    fftw_cleanup_threads();
    
    return 0;
    
    
    

}// MAIN










/*
 //[ Measures the total pseudo-cycles since the processor was powered on. Used to generate random seeds.
 unsigned long long rdtsc(){
 unsigned int lo,hi;
 __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
 return ((unsigned long long)hi << 32) | lo;
 }//]
 */








