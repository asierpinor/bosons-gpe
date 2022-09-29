#include <iostream>
#include <fstream>
//#include <stdio.h>
//#include <unistd.h>
#include <vector>               // std:: vector
#include <complex>
#include <cmath>                // sin, cos etc.
#include <fftw3.h>              // FFTW
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <string>
#include <sstream>              // for input/ouput with strings

#include "boson_gpe_F.hpp"

using namespace std;

#ifndef X
#define X 0
#endif //X

#ifndef Y
#define Y 1
#endif //Y

#ifndef Z
#define Z 2
#endif //Z

#ifndef M_PI
#define M_PI        3.14159265358979323846264338327950288   // pi
#endif

#ifndef TYPE_DEF
#define TYPE_DEF
typedef double my_double;
typedef complex<my_double> comp;
typedef vector< comp > state_type;  // The type of container used to hold the state, vector of spins.
#endif



boson_gpe::boson_gpe(gsl_rng * rndm) : r(rndm)
{
    dim = 0;
    if (Nx>1) { dim++; Nlist.push_back(Nx); }
    if (Ny>1) { dim++; Nlist.push_back(Ny); }
    if (Nz>1) { dim++; Nlist.push_back(Nz); }
    if (dim==0) { dim = DIM; Nlist.resize(dim,N); }
    
    Ntotal=1;
    for (int i=0; i<dim; i++) { Ntotal*=Nlist[i]; }
    volume = pow(ls,dim)*Ntotal;
    I = comp(0.0,1.0);
	
    NforR2Cfourier = ( Ntotal/Nlist[dim-1] )*( Nlist[dim-1]/2 + 1 );
    NlistforR2C = Nlist;
    NlistforR2C[dim-1] = Nlist[dim-1]/2+1;
    
    counter = 0;
    nTimes = (int)ceil(tmax/dt);
    allTimes.resize(nTimes+2,-1);
    phi0t.resize(nTimes+2,comp(0.0,0.0));   // Initializes vectors holding times and phi(t,p=0). The +2 is just profilaxis and the -1 is for later uses.
    
    bin_log();
    
    phiX = (comp *) fftw_malloc(sizeof(comp) * (Ntotal));
    phiK = (comp *) fftw_malloc(sizeof(comp) * (Ntotal));
    
    int * NlistAsArray = & Nlist[0];
    
    switch (dim) {
        case 1:
            fft_xk = fftw_plan_dft_1d(Nlist[0],(fftw_complex *)phiX, (fftw_complex*)phiK, FFTW_FORWARD, FFTW_MEASURE);
            fft_kx = fftw_plan_dft_1d(Nlist[0],(fftw_complex *)phiK, (fftw_complex*)phiX, FFTW_BACKWARD, FFTW_MEASURE);
            break;
            
        case 2:
            fft_xk = fftw_plan_dft_2d(Nlist[0],Nlist[1],(fftw_complex *)phiX, (fftw_complex*)phiK, FFTW_FORWARD, FFTW_MEASURE);
            fft_kx = fftw_plan_dft_2d(Nlist[0],Nlist[1],(fftw_complex *)phiK, (fftw_complex*)phiX, FFTW_BACKWARD, FFTW_MEASURE);
            break;
            
        case 3:
            fft_xk = fftw_plan_dft_3d(Nlist[0],Nlist[1],Nlist[2],(fftw_complex *)phiX, (fftw_complex*)phiK, FFTW_FORWARD, FFTW_MEASURE);
            fft_kx = fftw_plan_dft_3d(Nlist[0],Nlist[1],Nlist[2],(fftw_complex *)phiK, (fftw_complex*)phiX, FFTW_BACKWARD, FFTW_MEASURE);
            break;
            
        default:
            fft_xk = fftw_plan_dft(dim,NlistAsArray,(fftw_complex *)phiX, (fftw_complex*)phiK, FFTW_FORWARD, FFTW_MEASURE);
            fft_kx = fftw_plan_dft(dim,NlistAsArray,(fftw_complex *)phiK, (fftw_complex*)phiX, FFTW_BACKWARD, FFTW_MEASURE);
            break;
            
    }
    
    //for (int i=0; i<Ntotal; i++) phiX[i] = b[i];  // fill field
    for (int i=0; i<Ntotal; i++) phiX[i] = 0;  // fill field
    
    fill_kTable();
    
    // Bogoliubov arrays
    bogHappened=0;      // Arrays not allocated yet.
    
    // Binary coding to select potential term
    gTerms = 0;
    if (g!=0) gTerms += 1;
    if (gn!=0) gTerms += 2;
    if (glong!=0) gTerms += 4;
    
    
    
    if (glong!=0)
    {
		cout << "WARNING: Using glong!=0 may collide with density-phase output\n";
        
        Vint = (double *) fftw_malloc(sizeof(double) * Ntotal);
        VintK = (comp *) fftw_malloc(sizeof(comp) * NforR2Cfourier);
        density = (double *) fftw_malloc(sizeof(double) * Ntotal);
        densityK = (comp *) fftw_malloc(sizeof(comp) * NforR2Cfourier);
        
        fft_Vint = fftw_plan_dft_r2c(dim,NlistAsArray,Vint,(fftw_complex*)VintK,FFTW_ESTIMATE);
        fft_density = fftw_plan_dft_r2c(dim,NlistAsArray,density,(fftw_complex*)densityK,FFTW_ESTIMATE);
        fft_convolution = fftw_plan_dft_c2r(dim,NlistAsArray,(fftw_complex*)densityK,density,FFTW_ESTIMATE);
        
        fill_Vint();
        
        fftw_execute(fft_Vint);
        
    }
    
	
	// -------------
	//	Commutators
	// -------------
	
    // Fundamental Fields: Fill external field array hext and define Fourier transforms
    if (tw>=0 && tw<tmax)
    {
        if (CRmodeResp==1)
        {
            hextR = new comp [Ntotal];
            hextI = new comp [Ntotal];
            for (int n=0; n<Ntotal; n++) { hextR[n] = h*comp(2*(int)gsl_rng_uniform_int(r,2)-1,0); hextI[n] = h*comp(2*(int)gsl_rng_uniform_int(r,2)-1,0); }
        }
        else if (CRmodeResp==2)
        {
            hextR = new comp [Ntotal];
            hextI = new comp [Ntotal];
            hextR_k = new comp [Ntotal];
            hextI_k = new comp [Ntotal];
            for (int n=0; n<Ntotal; n++) { hextR[n] = h*comp(2*(int)gsl_rng_uniform_int(r,2)-1,0); hextI[n] = h*comp(2*(int)gsl_rng_uniform_int(r,2)-1,0); }
            
            fft_hextR = fftw_plan_dft(dim,NlistAsArray,(fftw_complex *)hextR, (fftw_complex*)hextR_k, FFTW_FORWARD, FFTW_MEASURE);
            fft_hextI = fftw_plan_dft(dim,NlistAsArray,(fftw_complex *)hextI, (fftw_complex*)hextI_k, FFTW_FORWARD, FFTW_MEASURE);
            fftw_execute(fft_hextR); fftw_execute(fft_hextI);
            fftw_destroy_plan(fft_hextR); fftw_destroy_plan(fft_hextI);
            
            commSpectrum = new comp [Nbins*4];
            
            
            // Cut-off momenta higher than Pcutoff to get better convergence
            if (Pcutoff>0)
            {
                double fakt = 4.0/(ls*ls);
                double ksq;
                
                vector<int> indices(dim,0);
                
                for (int n=0; n<Ntotal; n++)
                {
                    get_indices(indices,n);
                    
                    double ksq = 0; // Physical momentum squared
                    for (int i=0; i<dim; i++) ksq += sin(M_PI*indices[i]/Nlist[i])*sin(M_PI*indices[i]/Nlist[i]);
                    ksq *= fakt;
                    
                    if (sqrt(ksq)>Pcutoff) { hextR_k[n] = 0; hextI_k[n] = 0; }
                    
                }
                
                // Fourier transfrom back
                fftw_plan fft_hextR_back = fftw_plan_dft(dim,NlistAsArray,(fftw_complex *)hextR_k, (fftw_complex*)hextR, FFTW_BACKWARD, FFTW_MEASURE);
                fftw_plan fft_hextI_back = fftw_plan_dft(dim,NlistAsArray,(fftw_complex *)hextI_k, (fftw_complex*)hextI, FFTW_BACKWARD, FFTW_MEASURE);
                fftw_execute(fft_hextR_back); fftw_execute(fft_hextI_back);
                fftw_destroy_plan(fft_hextR_back); fftw_destroy_plan(fft_hextI_back);
                
                // Normalize
                double norm = 1.0/((double)Ntotal);
                for (int n=0; n<Ntotal; n++) { hextR[n] *= norm;   hextI[n] *= norm; }
                
            }
            
        }
        else if (CRmodeResp==3)
        {
            commSpectrum = new comp [Nbins*4];
        }
        
        phiX_hR = (comp *) fftw_malloc(sizeof(comp) * (Ntotal));
        phiK_hR = (comp *) fftw_malloc(sizeof(comp) * (Ntotal));
        fft_xk_hR = fftw_plan_dft(dim,NlistAsArray,(fftw_complex *)phiX_hR, (fftw_complex*)phiK_hR, FFTW_FORWARD, FFTW_MEASURE);
        fft_kx_hR = fftw_plan_dft(dim,NlistAsArray,(fftw_complex *)phiK_hR, (fftw_complex*)phiX_hR, FFTW_BACKWARD, FFTW_MEASURE);
        
        phiX_hI = (comp *) fftw_malloc(sizeof(comp) * (Ntotal));
        phiK_hI = (comp *) fftw_malloc(sizeof(comp) * (Ntotal));
        fft_xk_hI = fftw_plan_dft(dim,NlistAsArray,(fftw_complex *)phiX_hI, (fftw_complex*)phiK_hI, FFTW_FORWARD, FFTW_MEASURE);
        fft_kx_hI = fftw_plan_dft(dim,NlistAsArray,(fftw_complex *)phiK_hI, (fftw_complex*)phiX_hI, FFTW_BACKWARD, FFTW_MEASURE);
        
    }
	
	
    // Density-Phase: Fill external field array hextD and hextP and define Fourier transforms for DPrho computation
    if (twDPrho>=0 && twDPrho<tmax)
    {
        hextD = new comp [Ntotal];
        hextP = new comp [Ntotal];
        hextD_k = new comp [Ntotal];
        hextP_k = new comp [Ntotal];
        for (int n=0; n<Ntotal; n++) { hextD[n] = hDP*comp(2*(int)gsl_rng_uniform_int(r,2)-1,0); hextP[n] = hDP*comp(2*(int)gsl_rng_uniform_int(r,2)-1,0); }
        
        fft_hextD = fftw_plan_dft(dim,NlistAsArray,(fftw_complex *)hextD, (fftw_complex*)hextD_k, FFTW_FORWARD, FFTW_MEASURE);
        fft_hextP = fftw_plan_dft(dim,NlistAsArray,(fftw_complex *)hextP, (fftw_complex*)hextP_k, FFTW_FORWARD, FFTW_MEASURE);
        fftw_execute(fft_hextD); fftw_execute(fft_hextP);
        fftw_destroy_plan(fft_hextD); fftw_destroy_plan(fft_hextP);
        
        
        // Cut-off momenta higher than Pcutoff to get better convergence
        if (Pcutoff>0)
        {
            double fakt = 4.0/(ls*ls);
            double ksq;
            
            vector<int> indices(dim,0);
            
            for (int n=0; n<Ntotal; n++)
            {
                get_indices(indices,n);
                
                double ksq = 0; // Physical momentum squared
                for (int i=0; i<dim; i++) ksq += sin(M_PI*indices[i]/Nlist[i])*sin(M_PI*indices[i]/Nlist[i]);
                ksq *= fakt;
                
                if (sqrt(ksq)>Pcutoff) { hextD_k[n] = 0; hextP_k[n] = 0; }
                
            }
            
            // Fourier transfrom back
            fftw_plan fft_hextD_back = fftw_plan_dft(dim,NlistAsArray,(fftw_complex *)hextD_k, (fftw_complex*)hextD, FFTW_BACKWARD, FFTW_MEASURE);
            fftw_plan fft_hextP_back = fftw_plan_dft(dim,NlistAsArray,(fftw_complex *)hextP_k, (fftw_complex*)hextP, FFTW_BACKWARD, FFTW_MEASURE);
            fftw_execute(fft_hextD_back); fftw_execute(fft_hextP_back);
            fftw_destroy_plan(fft_hextD_back); fftw_destroy_plan(fft_hextP_back);
            
            // Normalize
            double norm = 1.0/((double)Ntotal);
            for (int n=0; n<Ntotal; n++) { hextD[n] *= norm;   hextP[n] *= norm; }
            
        }
            
        
        phiX_hD = (comp *) fftw_malloc(sizeof(comp) * (Ntotal));
        phiK_hD = (comp *) fftw_malloc(sizeof(comp) * (Ntotal));
        fft_xk_hD = fftw_plan_dft(dim,NlistAsArray,(fftw_complex *)phiX_hD, (fftw_complex*)phiK_hD, FFTW_FORWARD, FFTW_MEASURE);
        fft_kx_hD = fftw_plan_dft(dim,NlistAsArray,(fftw_complex *)phiK_hD, (fftw_complex*)phiX_hD, FFTW_BACKWARD, FFTW_MEASURE);
        
        phiX_hP = (comp *) fftw_malloc(sizeof(comp) * (Ntotal));
        phiK_hP = (comp *) fftw_malloc(sizeof(comp) * (Ntotal));
        fft_xk_hP = fftw_plan_dft(dim,NlistAsArray,(fftw_complex *)phiX_hP, (fftw_complex*)phiK_hP, FFTW_FORWARD, FFTW_MEASURE);
        fft_kx_hP = fftw_plan_dft(dim,NlistAsArray,(fftw_complex *)phiK_hP, (fftw_complex*)phiX_hP, FFTW_BACKWARD, FFTW_MEASURE);
        
    }
    
	
	
	// -----------------
	//	Anticommutators
	// -----------------
	
    // For F(t,twF) output
    if (twF>=0)
    {
        phi_at_twF = (comp *) fftw_malloc(sizeof(comp) * (Ntotal));
    }
	
    // For density-phase F(t,twDPF) output
    if (twDPF>=0)
    {
        phase_at_twDPF = (comp *) fftw_malloc(sizeof(comp) * (Ntotal));
		dichte_at_twDPF = (comp *) fftw_malloc(sizeof(comp) * (Ntotal));
    }
	
	
	
	// ----------------
	//	Density-Phase
	// ----------------
	
	// Density-phase output: spectrum, F or rho
	if (DPspstep>0 || twDPF>=0 || twDPrho>=0)
	{
        dichte = (double *) fftw_malloc(sizeof(double) * Ntotal);
        dichteK = (comp *) fftw_malloc(sizeof(comp) * NforR2Cfourier);
        phase = (double *) fftw_malloc(sizeof(double) * Ntotal);
        phaseK = (comp *) fftw_malloc(sizeof(comp) * NforR2Cfourier);
        fft_dichte = fftw_plan_dft_r2c(dim,NlistAsArray,dichte,(fftw_complex*)dichteK,FFTW_ESTIMATE);
		fft_phase = fftw_plan_dft_r2c(dim,NlistAsArray,phase,(fftw_complex*)phaseK,FFTW_ESTIMATE);
		
		if (twDPrho>=0)
		{
	        dichte_hD = (double *) fftw_malloc(sizeof(double) * Ntotal);
	        dichteK_hD = (comp *) fftw_malloc(sizeof(comp) * NforR2Cfourier);
	        phase_hD = (double *) fftw_malloc(sizeof(double) * Ntotal);
	        phaseK_hD = (comp *) fftw_malloc(sizeof(comp) * NforR2Cfourier);
	        fft_dichte_hD = fftw_plan_dft_r2c(dim,NlistAsArray,dichte_hD,(fftw_complex*)dichteK_hD,FFTW_ESTIMATE);
			fft_phase_hD = fftw_plan_dft_r2c(dim,NlistAsArray,phase_hD,(fftw_complex*)phaseK_hD,FFTW_ESTIMATE);
		
	        dichte_hP = (double *) fftw_malloc(sizeof(double) * Ntotal);
	        dichteK_hP = (comp *) fftw_malloc(sizeof(comp) * NforR2Cfourier);
	        phase_hP = (double *) fftw_malloc(sizeof(double) * Ntotal);
	        phaseK_hP = (comp *) fftw_malloc(sizeof(comp) * NforR2Cfourier);
	        fft_dichte_hP = fftw_plan_dft_r2c(dim,NlistAsArray,dichte_hP,(fftw_complex*)dichteK_hP,FFTW_ESTIMATE);
			fft_phase_hP = fftw_plan_dft_r2c(dim,NlistAsArray,phase_hP,(fftw_complex*)phaseK_hP,FFTW_ESTIMATE);
		}
	}
    

    
    
}




boson_gpe::~boson_gpe()
{
   	
    fftw_destroy_plan(fft_xk);
    fftw_destroy_plan(fft_kx);
    fftw_free(phiX); fftw_free(phiK);
    
    if (glong!=0)
    {
        fftw_free(Vint); fftw_free(VintK);
        fftw_free(density); fftw_free(densityK);
        fftw_destroy_plan(fft_Vint);
        fftw_destroy_plan(fft_density);
        fftw_destroy_plan(fft_convolution);
    }
    
	
	// Commutators
	
    if (tw>=0 && tw<tmax)
    {
        if (CRmodeResp==1)
        {
            delete [] hextR;
            delete [] hextI;
        }
        else if (CRmodeResp==2)
        {
            delete [] hextR;
            delete [] hextI;
            delete [] hextR_k;
            delete [] hextI_k;
            delete [] commSpectrum;
        }
        else if (CRmodeResp==3)
        {
            delete [] commSpectrum;
        }
        
        fftw_destroy_plan(fft_xk_hR);
        fftw_destroy_plan(fft_kx_hR);
        fftw_free(phiX_hR); fftw_free(phiK_hR);
        
        fftw_destroy_plan(fft_xk_hI);
        fftw_destroy_plan(fft_kx_hI);
        fftw_free(phiX_hI); fftw_free(phiK_hI);
    }
	
    if (twDPrho>=0 && twDPrho<tmax)
    {
        delete [] hextD;
        delete [] hextP;
        delete [] hextD_k;
        delete [] hextP_k;
        //delete [] commSpectrum;
        
        fftw_destroy_plan(fft_xk_hD);
        fftw_destroy_plan(fft_kx_hD);
        fftw_free(phiX_hD); fftw_free(phiK_hD);
        
        fftw_destroy_plan(fft_xk_hP);
        fftw_destroy_plan(fft_kx_hP);
        fftw_free(phiX_hP); fftw_free(phiK_hP);
    }
	
	
	
	// Anticommutators
    
    if (twF>=0)
    {
        fftw_free(phi_at_twF);
    }
	
    if (twDPF>=0)
    {
        fftw_free(phase_at_twDPF);
		fftw_free(dichte_at_twDPF);
    }
	
	
	
	// Density-Phase
	
    if (DPspstep>0 || twDPF>=0 || twDPrho>=0)
    {
        fftw_free(dichte); fftw_free(dichteK);
		fftw_free(phase); fftw_free(phaseK);
        fftw_destroy_plan(fft_dichte);
		fftw_destroy_plan(fft_phase);
		
		if (twDPrho>=0)
		{
	        fftw_free(dichte_hD); fftw_free(dichteK_hD);
			fftw_free(phase_hD); fftw_free(phaseK_hD);
	        fftw_destroy_plan(fft_dichte_hD);
			fftw_destroy_plan(fft_phase_hD);
			
	        fftw_free(dichte_hP); fftw_free(dichteK_hP);
			fftw_free(phase_hP); fftw_free(phaseK_hP);
	        fftw_destroy_plan(fft_dichte_hP);
			fftw_destroy_plan(fft_phase_hP);
		}
    }
    
	
	
    delete [] momPerBin;
    delete [] binAvFactor;
    delete [] whichBin;
    delete [] kBin;
    delete [] spectrum;
    delete [] kTable;
    
    if (bogHappened!=0)
    {
        delete [] bogSpectrum;
        delete [] bogK;
    }
}





// Gets indices of each dimension assuming row-order
void boson_gpe::get_indices (vector<int> & indices, int n)
{
    int temp = 0;
    int rest = n;
    int block = Ntotal;
    
    while (temp<dim)
    {
        block = block/Nlist[temp];
        indices[temp] = rest/block;     //should be able to do this more efficiently
        rest -= indices[temp]*block;
        temp++;
    }
}// get_indices






// Gets indices of each dimension assuming row-order
void boson_gpe::get_indices_r2c (vector<int> & indices, int n)
{
    int temp = 0;
    int rest = n;
    int block = NforR2Cfourier;
    
    while (temp<dim)
    {
        block = block/NlistforR2C[temp];
        indices[temp] = rest/block;     //should be able to do this more efficiently
        rest -= indices[temp]*block;
        temp++;
    }
}// get_indices






// Computes the position in the array given the individual indices of each dimension
int boson_gpe::get_array_position (vector<int> indices)
{
    int pos = indices[0];
    for (int i=1; i<dim; i++) pos = pos*Nlist[i] + indices[i];
    
    return pos;
}// get_array_position






// Fill long-range interaction matrix (distances taken with respect to zero)
void boson_gpe::fill_Vint ()
{
    vector<int> indices(dim,0);
    double dist, tempdist;
    
    Vint[0]=0;  // No self-interaction term
    
    for (int n=1; n<Ntotal; n++)
    {
        get_indices(indices,n);
        dist = 0;
        
        for (int i=0; i<dim; i++) { tempdist = min(indices[i],Nlist[i]-indices[i]);  dist += tempdist*tempdist; }
        dist = sqrt(dist);
        
        Vint[n] = glong/pow(dist,explong);
    }
}// fill_Vint




vector<double> boson_gpe::getPhysMom (int binNumber)
{
    if (binNumber>=Nbins) cout << "Wrong use of getPhysMom: n larger than Nbins.\n";
    
    vector<double> temp;
    
    temp.push_back(kBin[binNumber]);    // Physical momentum of bin n
    temp.push_back(momPerBin[binNumber]);   // Number of momenta in that bin

    return temp;
}







/////////////////////////////////////////////////////////////////////////

///////////////            INITIAL CONDITIONS          //////////////////

/////////////////////////////////////////////////////////////////////////




// Add noise to initial condition
void boson_gpe::addTWAnoise ()
{
    switch (IC) {
            
        case 1:
            addTWAnoise_box();
            break;
            
        case 2:
            addTWAnoise_homoDensity();
            break;
			
	    case 4:
	        addTWAnoise_oneQuadrant();
            break;
			
		case 5:
	        addTWAnoise_randomPhase();
            break;
            
        default:
            addTWAnoise_box();
            break;
    }
}



//[ Add TWA noise to initial conditions
void boson_gpe::addTWAnoise_box ()
{
    
    double fakt = 4.0/(ls*ls);
    double ksq;
    
    vector<int> indices(dim,0);
    
    for (int n=0; n<Ntotal; n++)
    {
        get_indices(indices,n);
        
        double ksq = 0; // Physical momentum squared
        for (int i=0; i<dim; i++) ksq += sin(M_PI*indices[i]/Nlist[i])*sin(M_PI*indices[i]/Nlist[i]);
        ksq *= fakt;
        
        sigma = sqrt(amp*volume/2.0);         // Without quantum one-half. Assumes large amplitude rescaled
        
        if (ksq > Qmax*Qmax || sqrt(ksq) < Qmin)
        {
            sigma = 0.0;    // Neglecting the quantum 1/2
        }
        
        phiK[n] = comp(gsl_ran_gaussian(r,sigma),gsl_ran_gaussian(r,sigma));

    }// Add gaussian random numbers in fourier space.
    
    fftw_execute(fft_kx);
    
    for (int i=0; i<Ntotal; i++)
    {
        phiX[i] = phiX[i] / comp(Ntotal,0.0);   // Volume factor due to Fourier normalization. N^3 or volume?
    }
    
    
    
}//]







//[ Add TWA noise to initial conditions, set density to be homogeneous
void boson_gpe::addTWAnoise_homoDensity ()
{
    
    double fakt = 4.0/(ls*ls);
    double ksq;
    
    vector<int> indices(dim,0);
    
    for (int n=0; n<Ntotal; n++)
    {
        get_indices(indices,n);
        
        double ksq = 0; // Physical momentum squared
        for (int i=0; i<dim; i++) ksq += sin(M_PI*indices[i]/Nlist[i])*sin(M_PI*indices[i]/Nlist[i]);
        ksq *= fakt;
        
        sigma = sqrt(amp*volume/2.0);         // Without quantum one-half. Assumes large amplitude rescaled
        
        if (ksq > Qmax*Qmax || sqrt(ksq) < Qmin)
        {
            sigma = 0.0;    // Neglecting the quantum 1/2
        }
        
        phiK[n] = comp(gsl_ran_gaussian(r,sigma),gsl_ran_gaussian(r,sigma));
        
    }// Add gaussian random numbers in fourier space.
    
    fftw_execute(fft_kx);
    
    for (int i=0; i<Ntotal; i++)
    {
        phiX[i] = phiX[i] / comp(Ntotal,0.0);   // Volume factor due to Fourier normalization. N^3 or volume?
    }
    
    // Count number of particles
    double partNum = 0.0;
    
    #pragma omp parallel for reduction(+:partNum)
    for (int i=0; i<Ntotal; i++)
    {
        partNum += sqrabs(phiX[i]);
    }
	
	double meanDens = sqrt(partNum / Ntotal);
    
    // Normalize to equal density
    #pragma omp parallel for
    for (int i=0; i<Ntotal; i++)
    {
        phiX[i] = phiX[i] / abs(phiX[i]) * meanDens;
    }
    
    
}//]






//[ Add TWA noise to initial conditions, restrict Psi to postive real and imag parts
void boson_gpe::addTWAnoise_oneQuadrant ()
{
    
    double fakt = 4.0/(ls*ls);
    double ksq;
    
    vector<int> indices(dim,0);
    
    for (int n=0; n<Ntotal; n++)
    {
        get_indices(indices,n);
        
        double ksq = 0; // Physical momentum squared
        for (int i=0; i<dim; i++) ksq += sin(M_PI*indices[i]/Nlist[i])*sin(M_PI*indices[i]/Nlist[i]);
        ksq *= fakt;
        
        sigma = sqrt(amp*volume/2.0);         // Without quantum one-half. Assumes large amplitude rescaled
        
        if (ksq > Qmax*Qmax || sqrt(ksq) < Qmin)
        {
            sigma = 0.0;    // Neglecting the quantum 1/2
        }
        
        phiK[n] = comp(gsl_ran_gaussian(r,sigma),gsl_ran_gaussian(r,sigma));
        
    }// Add gaussian random numbers in fourier space.
    
    fftw_execute(fft_kx);
    
    for (int i=0; i<Ntotal; i++)
    {
        phiX[i] = phiX[i] / comp(Ntotal,0.0);   // Volume factor due to Fourier normalization. N^3 or volume?
    }
    
    // Restrict to positive real and imaginary parts
    for (int i=0; i<Ntotal; i++)
    {
        phiX[i] = comp(phiX[i].real(),abs(phiX[i].imag()));
    }
    
    
}//]





//[ Add TWA noise to initial conditions, const density, random phase
void boson_gpe::addTWAnoise_randomPhase ()
{
    double dens = amp;
	double rn;
	
    for (int i=0; i<Ntotal; i++)
    {
		rn = gsl_rng_uniform(r);
        phiX[i] = dens * comp(cos(2*M_PI*rn),sin(2*M_PI*rn));
    }
	
	fftw_execute(fft_xk);
}//]
    
    
    
    
    
    
    
    
    
    /////////////////////////////////////////////////////////////////////////
    
    ///////////////////            DYNAMICS          ////////////////////////
    
    /////////////////////////////////////////////////////////////////////////
    


// Makes two copies of the system, which will be rotated at t=tw with hR and hI, respectively.
void boson_gpe::triple_system ()
{
    for (int n=0; n<Ntotal; n++) { phiX_hR[n] = phiX[n]; phiX_hI[n] = phiX[n]; }
}



// Makes two copies of the system, which will be rotated at t=twDPrho with hD and hP, respectively.
void boson_gpe::triple_system_DP ()
{
    for (int n=0; n<Ntotal; n++) { phiX_hD[n] = phiX[n]; phiX_hP[n] = phiX[n]; }
}





// Save phi(twF,p)
void boson_gpe::save_phi_at_twF ()
{
    for (int n=0; n<Ntotal; n++)
    {
        phi_at_twF[n] = phiK[n];
    }
}



// Save phase(twDPF,p) and dichte(twDPF,p)
void boson_gpe::save_DP_at_twDPF ()
{
	// Fill dichte and phase arrays
	for (int i=0; i<Ntotal; i++)
	{
		dichte[i] = sqrabs(phiX[i]);
		phase[i] = atan2(phiX[i].imag(),phiX[i].real());
	}
	
	// Fourier transform to dichteK and phaseK
	fftw_execute(fft_dichte);
	fftw_execute(fft_phase);
	
	// Save dichteK and phaseK
    for (int n=0; n<Ntotal; n++)
    {
        phase_at_twDPF[n] = phaseK[n];
		dichte_at_twDPF[n] = dichteK[n];
    }
}






// Fill kTable with kinetic exp factor for integration step
void boson_gpe::fill_kTable ()
{
    kTable = new comp [Ntotal];
    
    double fakt = 4.0/(ls*ls);
    double kinfakt = 1.0/(2.0*m);
    
    vector<int> indices(dim,0);
    
    for (int n=0; n<Ntotal; n++)
    {
        get_indices(indices,n);
        
        double ksq = 0; // Physical momentum squared
        for (int i=0; i<dim; i++) ksq += sin(M_PI*indices[i]/Nlist[i])*sin(M_PI*indices[i]/Nlist[i]);
        ksq *= fakt;
        
        kTable[n] = exp(-I*kinfakt*dt*ksq);     // fill with kinetic exp factor for integration step
    }

}



// Calculates phi(t+dt) from GPE with split-step method
void boson_gpe::dynamics ()
{
    
    //Potential part
    switch (gTerms) {
        case 0:
            dynamics_pot_000(); break;
            
        case 1:
            dynamics_pot_001(); break;
            
        case 2:
            dynamics_pot_010(); break;
            
        case 3:
            dynamics_pot_011(); break;
            
        case 4:
            dynamics_pot_100(); break;
            
        case 5:
            dynamics_pot_101(); break;
            
        case 6:
            dynamics_pot_110(); break;
            
        case 7:
            dynamics_pot_111(); break;
            
        default:
            cout << "Problem with gTerms\n";
            break;
    }
    
    //Kinetic part
    dynamics_kin();
    
}




// Multiply time evolution with potential part for (glong,gn,g)=(0,0,0)
void boson_gpe::dynamics_pot_000 () { }



// Multiply time evolution with potential part for (glong,gn,g)=(0,0,*)
void boson_gpe::dynamics_pot_001 ()
{
    double tt;
    double intfakt = g*dt;
    
    #pragma omp parallel for private(tt)
    for(int i=0; i<Ntotal; i++)
    {
        tt = (sqrabs(phiX[i]))*intfakt;
        phiX[i] *= comp(cos(tt),-sin(tt));      // multiply with exp(-I*g*dt*|phi|^2)
    }
}

// Multiply time evolution with potential part for (glong,gn,g)=(0,*,0)
void boson_gpe::dynamics_pot_010 ()
{
    double tt;
    double intfakt6 = gn*dt;
    double exponent = expphi/2-1;
    
    #pragma omp parallel for private(tt)
    for(int i=0; i<Ntotal; i++)
    {
        tt = pow(sqrabs(phiX[i]),exponent)*intfakt6;
        phiX[i] *= comp(cos(tt),-sin(tt));
    }
}


// Multiply time evolution with potential part for (glong,gn,g)=(0,*,*)
void boson_gpe::dynamics_pot_011 ()
{
    double tt;
    double intfakt = g*dt;
    double intfakt6 = gn*dt;
    double exponent = expphi/2-1;
    
    #pragma omp parallel for private(tt)
    for(int i=0; i<Ntotal; i++)
    {
        tt = (sqrabs(phiX[i]))*intfakt + pow(sqrabs(phiX[i]),exponent)*intfakt6;
        phiX[i] *= comp(cos(tt),-sin(tt));
    }
}




// Multiply time evolution with potential part for (glong,gn,g)=(*,0,0)
void boson_gpe::dynamics_pot_100 ()
{
    interaction_longrange(); // Computes convolution and stores it in array density
    
    double tt;
    double intfaktlong = dt;
    
    #pragma omp parallel for private(tt)
    for(int i=0; i<Ntotal; i++)
    {
        tt = density[i]*intfaktlong;
        phiX[i] *= comp(cos(tt),-sin(tt));
    }
}




// Multiply time evolution with potential part for (glong,gn,g)=(*,0,*)
void boson_gpe::dynamics_pot_101 ()
{
    interaction_longrange(); // Computes convolution and stores it in array density
    
    double tt;
    double intfaktlong = dt;
    double intfakt = g*dt;
    
    #pragma omp parallel for private(tt)
    for(int i=0; i<Ntotal; i++)
    {
        tt = density[i]*intfaktlong + sqrabs(phiX[i])*intfakt;
        phiX[i] *= comp(cos(tt),-sin(tt));
    }
}




// Multiply time evolution with potential part for (glong,gn,g)=(*,*,0)
void boson_gpe::dynamics_pot_110 ()
{
    interaction_longrange(); // Computes convolution and stores it in array "density"
    
    double tt;
    double intfaktlong = dt;
    double intfakt6 = gn*dt;
    double exponent = expphi/2-1;
    
    #pragma omp parallel for private(tt)
    for(int i=0; i<Ntotal; i++)
    {
        tt = density[i]*intfaktlong + pow(sqrabs(phiX[i]),exponent)*intfakt6;
        phiX[i] *= comp(cos(tt),-sin(tt));
    }
}




// Multiply time evolution with potential part for (glong,gn,g)=(*,*,*)
void boson_gpe::dynamics_pot_111 ()
{
    interaction_longrange(); // Computes convolution and stores it in array "density"
    
    double tt;
    double intfaktlong = dt;
    double intfakt6 = gn*dt;
    double intfakt = g*dt;
    double exponent = expphi/2-1;
    
    #pragma omp parallel for private(tt)
    for(int i=0; i<Ntotal; i++)
    {
        tt = density[i]*intfaktlong + pow(sqrabs(phiX[i]),exponent)*intfakt6 + sqrabs(phiX[i])*intfakt;
        phiX[i] *= comp(cos(tt),-sin(tt));
    }
}




// Multiply time evolution with kinetic part
void boson_gpe::dynamics_kin ()
{
    double fouriernorm = 1.0/Ntotal;
    fftw_execute(fft_xk);     //forward transformation
    
    #pragma omp parallel for
    for(int i=0; i<Ntotal; i++) phiK[i] *= kTable[i]; // Multiply with exp(-I*dt*k^2/(2*m))
    
    fftw_execute(fft_kx);     //backward transformation
    
    #pragma omp parallel for
    for(int i=0;i<Ntotal;i++) phiX[i] *= fouriernorm; //Normalization
}




// Computes interaction part of GPE for long-range interactions
void boson_gpe::interaction_longrange ()
{
    #pragma omp parallel for
    for (int n=0; n<Ntotal; n++) density[n] = sqrabs(phiX[n]);  // Compute |psi|^2
    
    fftw_execute(fft_density);

    #pragma omp parallel for
    for (int n=0; n<NforR2Cfourier; n++) densityK[n] = densityK[n]*VintK[n];    // Convolute in Fourier space

    fftw_execute(fft_convolution);
    
    #pragma omp parallel for
    for (int n=0; n<Ntotal; n++) density[n] = density[n]/Ntotal; //Normalize
}









// Calculates phi(t+dt) from GPE with split-step method for phi_hR
void boson_gpe::dynamics_phihR ()
{
    if (gTerms!=1) cout << "Problem with gTerms for dynamics with h.\n";
    
    //Potential part
    double tt;
    double intfakt = g*dt;
    
    #pragma omp parallel for private(tt)
    for(int i=0; i<Ntotal; i++)
    {
        tt = (sqrabs(phiX_hR[i]))*intfakt;
        phiX_hR[i] *= comp(cos(tt),-sin(tt));      // multiply with exp(-I*g*dt*|phi|^2)
    }
            
    
    //Kinetic part
    double fouriernorm = 1.0/Ntotal;
    fftw_execute(fft_xk_hR);     //forward transformation
    
    #pragma omp parallel for
    for(int i=0; i<Ntotal; i++) phiK_hR[i] *= kTable[i]; // Multiply with exp(-I*dt*k^2/(2*m))
    
    fftw_execute(fft_kx_hR);     //backward transformation
    
    #pragma omp parallel for
    for(int i=0;i<Ntotal;i++) phiX_hR[i] *= fouriernorm; //Normalization
}




// Calculates phi(t+dt) from GPE with split-step method for phi_hI
void boson_gpe::dynamics_phihI ()
{
    if (gTerms!=1) cout << "Problem with gTerms for dynamics with h.\n";
    
    //Potential part
    double tt;
    double intfakt = g*dt;
    
    #pragma omp parallel for private(tt)
    for(int i=0; i<Ntotal; i++)
    {
        tt = (sqrabs(phiX_hI[i]))*intfakt;
        phiX_hI[i] *= comp(cos(tt),-sin(tt));      // multiply with exp(-I*g*dt*|phi|^2)
    }
    
    
    //Kinetic part
    double fouriernorm = 1.0/Ntotal;
    fftw_execute(fft_xk_hI);     //forward transformation
    
    #pragma omp parallel for
    for(int i=0; i<Ntotal; i++) phiK_hI[i] *= kTable[i]; // Multiply with exp(-I*dt*k^2/(2*m))
    
    fftw_execute(fft_kx_hI);     //backward transformation
    
    #pragma omp parallel for
    for(int i=0;i<Ntotal;i++) phiX_hI[i] *= fouriernorm; //Normalization
}




// Turn on external field hR at t=tw
void boson_gpe::rotation_hR()
{
    switch (CRmodeResp) {
            case 0:
            {
                comp hR (h,0.0);
                comp tt = I*hR/2.0;
                #pragma omp parallel for
                for (int n=0; n<Ntotal; n++) phiX_hR[n] += tt;
                break;
            }
            
            case 1:
            case 2:
            {
                #pragma omp parallel for
                for (int n=0; n<Ntotal; n++) phiX_hR[n] += I*hextR[n]/2.0;
                break;
            }
            
            case 3:
            {
                comp hR (h,0.0);
                phiX_hR[0] += I*hR/2.0;
                break;
            }
            
        default:
            cout << "Error: Wrong type of mode to output for linear response (rotation_hR).\n";
            break;
    }
    
    
    fftw_execute(fft_xk_hR);     //forward transformation
}



// Turn on external field hI at t=tw
void boson_gpe::rotation_hI()
{
    switch (CRmodeResp) {
            case 0:
            {
                comp hI (0.0,h);
                comp tt = I*hI/2.0;
                #pragma omp parallel for
                for (int n=0; n<Ntotal; n++) phiX_hI[n] += tt;
                break;
            }
            
            case 1:
            case 2:
            {
                #pragma omp parallel for
                for (int n=0; n<Ntotal; n++) phiX_hI[n] += -hextI[n]/2.0;
                break;
            }
            
            case 3:
            {
                comp hI (0.0,h);
                phiX_hI[0] += I*hI/2.0;
                break;
            }
            
        default:
            cout << "Error: Wrong type of mode to output for linear response (rotation_hI).\n";
            break;
    }
    
    fftw_execute(fft_xk_hI);     //forward transformation
}






// Calculates phi(t+dt) from GPE with split-step method for phi_hD
void boson_gpe::dynamics_phihD ()
{
    if (gTerms!=1) cout << "Problem with gTerms for dynamics with h.\n";
    
    //Potential part
    double tt;
    double intfakt = g*dt;
    
    #pragma omp parallel for private(tt)
    for(int i=0; i<Ntotal; i++)
    {
        tt = (sqrabs(phiX_hD[i]))*intfakt;
        phiX_hD[i] *= comp(cos(tt),-sin(tt));      // multiply with exp(-I*g*dt*|phi|^2)
    }
            
    
    //Kinetic part
    double fouriernorm = 1.0/Ntotal;
    fftw_execute(fft_xk_hD);     //forward transformation
    
    #pragma omp parallel for
    for(int i=0; i<Ntotal; i++) phiK_hD[i] *= kTable[i]; // Multiply with exp(-I*dt*k^2/(2*m))
    
    fftw_execute(fft_kx_hD);     //backward transformation
    
    #pragma omp parallel for
    for(int i=0;i<Ntotal;i++) phiX_hD[i] *= fouriernorm; //Normalization
}



// Calculates phi(t+dt) from GPE with split-step method for phi_hP
void boson_gpe::dynamics_phihP ()
{
    if (gTerms!=1) cout << "Problem with gTerms for dynamics with h.\n";
    
    //Potential part
    double tt;
    double intfakt = g*dt;
    
    #pragma omp parallel for private(tt)
    for(int i=0; i<Ntotal; i++)
    {
        tt = (sqrabs(phiX_hP[i]))*intfakt;
        phiX_hP[i] *= comp(cos(tt),-sin(tt));      // multiply with exp(-I*g*dt*|phi|^2)
    }
    
    
    //Kinetic part
    double fouriernorm = 1.0/Ntotal;
    fftw_execute(fft_xk_hP);     //forward transformation
    
    #pragma omp parallel for
    for(int i=0; i<Ntotal; i++) phiK_hP[i] *= kTable[i]; // Multiply with exp(-I*dt*k^2/(2*m))
    
    fftw_execute(fft_kx_hP);     //backward transformation
    
    #pragma omp parallel for
    for(int i=0;i<Ntotal;i++) phiX_hP[i] *= fouriernorm; //Normalization
}



// Turn on external field hD at t=twDPrho
void boson_gpe::rotation_hD()
{
    #pragma omp parallel for
    for (int n=0; n<Ntotal; n++) phiX_hD[n] *= comp(cos(hextD[n].real()/2.0),sin(hextD[n].real()/2.0));
    
    fftw_execute(fft_xk_hD);     //forward transformation
}


// Turn on external field hP at t=twDPrho
void boson_gpe::rotation_hP()
{
    //#pragma omp parallel for
    for (int n=0; n<Ntotal; n++)
	{
		//if ( abs(hextP[n].real()) > sqrabs(phiX_hP[n]) ) cout << "Problem with perturbation field: |h|/n = " << abs(hextP[n].real()) / sqrabs(phiX_hP[n]) << "\n";
		if ( abs(hextP[n].real()) < sqrabs(phiX_hP[n]) ) phiX_hP[n] *= (1 - 0.5 * hextP[n].real() / sqrabs(phiX_hP[n]) );
	}
    
    fftw_execute(fft_xk_hP);     //forward transformation
}







/*
// Calculates phi(t+dt) from GPE with split-step method
void boson_gpe::dynamics_with_h ()
{
    
    //Potential part
    switch (gTerms) {
            
        case 1:
            dynamics_pot_001_with_h(); break;
            
        default:
            cout << "Problem with gTerms for dynamics with h.\n";
            break;
    }
    
    //Kinetic part
    dynamics_kin_with_h();
    
    //External field h part
    dynamics_h();
    
}





// Multiply time evolution (with h) with potential part for (glong,gn,g)=(0,0,*)
void boson_gpe::dynamics_pot_001_with_h ()
{
    double tt;
    double intfakt = g*dt;
    
    #pragma omp parallel for private(tt)
    for(int i=0; i<Ntotal; i++)
    {
        tt = (sqrabs(phiX_h[i]))*intfakt;
        phiX_h[i] *= comp(cos(tt),-sin(tt));      // multiply with exp(-I*g*dt*|phi|^2)
    }
}





// Multiply time evolution (with h) with kinetic part
void boson_gpe::dynamics_kin_with_h ()
{
    double fouriernorm = 1.0/Ntotal;
    fftw_execute(fft_xk_h);     //forward transformation
    
    #pragma omp parallel for
    for(int i=0; i<Ntotal; i++) phiK_h[i] *= kTable[i]; // Multiply with exp(-I*dt*k^2/(2*m))
    
    fftw_execute(fft_kx_h);     //backward transformation
    
    #pragma omp parallel for
    for(int i=0;i<Ntotal;i++) phiX_h[i] *= fouriernorm; //Normalization
}





// Add external field part to time evolution
void boson_gpe::dynamics_h()
{
    comp tt = -I*dt/2.0;
    
    #pragma omp parallel for
    for (int n=0; n<Ntotal; n++) phiX_h[n] += tt*hext[n];
}
*/





    
    
    
    
    
    
    
    
    /////////////////////////////////////////////////////////////////////////
    
    ////////////////////            OUTPUT          /////////////////////////
    
    /////////////////////////////////////////////////////////////////////////
    



// Calculates the Bogoliubov modes from the free modes
void boson_gpe::bogtrafo ()
{
    double N0 = magnetization();           // Number of particles in 0-mode
    k0 = sqrt(2*m*g*N0/volume);     // Transition momentum
    
    double fakt = 4.0/(ls*ls);
    int indk, indminusk;
    
    vector<int> indices(dim,0);
    vector<int> indicesMinusK(dim,0);
    
    for (int n=0; n<Ntotal; n++)
    {
        //Indices of k and -k
        get_indices(indices,n);
        for (int i=0; i<dim; i++) indicesMinusK[i] = (Nlist[i]-indices[i])%Nlist[i];
        
        indk = n;
        indminusk = get_array_position(indicesMinusK);
        
        // Physical momentum squared. Should use this or lattice momentum???
        double ksq = 0;
        for (int i=0; i<dim; i++) ksq += sin(M_PI*indices[i]/Nlist[i])*sin(M_PI*indices[i]/Nlist[i]);
        ksq *= fakt;
        
        yk = sqrt(ksq)/k0;
        alphak = 1.0 + yk*yk - yk*sqrt(2.0 + yk*yk);
        uk = 1.0/sqrt(1-alphak*alphak);
        vk = -alphak*uk;
        
        bogK[indk] = uk*phiK[indk] - vk*phiK[indminusk];     // This assumes that the fft_xk has already been executed in write_spec
        
    }// Careful, should omit 0-mode -> change
    
}





double boson_gpe::magnetization ()
{
    comp temp = 0.0;
    
    for(int i=0; i<Ntotal; i++) temp += phiX[i];
    
    return sqrabs(temp)/volume;
}






double boson_gpe::totaln()
{
    double partNum = 0.0;
    
    #pragma omp parallel for reduction(+:partNum)
    for(int i=0; i<Ntotal; i++)
    {
        partNum += sqrabs(phiX[i]);
    }
    
    return partNum;
}






double boson_gpe::kinetic ()
{
    double kin = 0.0;
    
    vector<int> indices(dim,0);
    vector<int> indicesRight(dim,0); // indices of dimension
    
    int pos;
    vector<int> posRight(dim,0); // positions in array
    
    for (int n=0; n<Ntotal; n++)
    {
        get_indices(indices,n);
        pos = n;

        for (int i=0; i<dim; i++)
        {
            indicesRight = indices;
            indicesRight[i] = (indicesRight[i]+1)%Nlist[i];
            
            posRight[i] = get_array_position(indicesRight);
        }
        
        for (int i=0; i<dim; i++) kin += sqrabs( phiX[pos]-phiX[posRight[i]] ); // using grad(psi)*grad(psi^\dagger) with right derivative
    }
    
    return (kin/(2.0*m))/volume;
}





double boson_gpe::potenergy ()
{
    double pot = 0.0;
    
    if (g!=0)
    {
        #pragma omp parallel for reduction(+:pot)
        for(int i=0; i<Ntotal; i++) pot += sqrabs(phiX[i])*sqrabs(phiX[i]);
    }
    
    double exponent = expphi/2;
    double potn = 0.0;
    
    if (gn!=0)
    {
        #pragma omp parallel for reduction(+:potn)
        for(int i=0; i<Ntotal; i++) potn += pow(sqrabs(phiX[i]),exponent);
    }
    
    double potlong = 0.0;
    
    if (glong!=0)
    {
        #pragma omp parallel for
        for (int n=0; n<Ntotal; n++) density[n] = sqrabs(phiX[n]);  // Compute |psi|^2
        
        fftw_execute(fft_density);
        
        #pragma omp parallel reduction(+:potlong)
        {
            vector<int> indices(dim,0);
        
            #pragma omp for
            for (int n=0; n<NforR2Cfourier; n++)
            {
                get_indices_r2c(indices,n);
                
                if ( indices[dim-1]!=NlistforR2C[dim-1]-1 && indices[dim-1]!=0 ) potlong += 2*sqrabs(densityK[n])*VintK[n].real();
                else potlong += sqrabs(densityK[n])*VintK[n].real();
            }
        }
        potlong = potlong/Ntotal;
    }
    
    return ( pot*g*0.5 + potn*gn/exponent + potlong*0.5 )/volume;
}




double boson_gpe::densityVariance()
{
    double meanDens = 0, var = 0;
    
	// Compute mean
    #pragma omp parallel for reduction(+:meanDens)
    for(int i=0; i<Ntotal; i++)
    {
        meanDens += abs(phiX[i]);
    }
	meanDens = (meanDens / Ntotal);
	
	// Compute variance of density
	for (int i=0; i<Ntotal; i++)
	{
		var += sqr(abs(phiX[i])-meanDens) ;
	}
	var *= 1.0/Ntotal;
    
    return sqrt(var)/meanDens;
}

double boson_gpe::angularVariance()
{
    comp meanPsi = 0;
	double var = 0;
    
	// Compute mean Psi
    for(int i=0; i<Ntotal; i++)
    {
        meanPsi += phiX[i];
    }
	meanPsi = meanPsi * (1.0/ Ntotal);
	
	// Compute variance of density
	double temp;
	for (int i=0; i<Ntotal; i++)
	{
		temp = abs(atan2(phiX[i].imag(),phiX[i].real())-atan2(meanPsi.imag(),meanPsi.real()));
		temp = min(temp, abs(atan2(phiX[i].imag(),phiX[i].real())-atan2(meanPsi.imag(),meanPsi.real()) + 2*M_PI) );
		temp = min(temp, abs(atan2(phiX[i].imag(),phiX[i].real())-atan2(meanPsi.imag(),meanPsi.real()) - 2*M_PI) );
		var += sqr(temp) ;
	}
	var *= 1.0/Ntotal;
    
    return sqrt(var);
}




// Correlation function < phi_x phi_x^* >
double boson_gpe::two_pointfct ()
{
    double temp = 0.0;
    
    for(int i=0; i<Ntotal; i++) temp += sqrabs(phiX[i]);
    
    return temp/volume;
}




// Correlation function < phi_x phi_x >
comp boson_gpe::two_pointfct_anomalous ()
{
    comp temp = 0.0;
    
    for(int i=0; i<Ntotal; i++) temp += phiX[i]*phiX[i];
    
    return temp/volume;
}



// Correlation function < (phi_x phi_x^*)^2 >
double boson_gpe::four_pointfct ()
{
    double temp = 0.0;
    
    for(int i=0; i<Ntotal; i++) temp += sqrabs(phiX[i])*sqrabs(phiX[i]);
    
    return temp/volume;
}




// Correlation function < (phi_x phi_x^*)^3 >
double boson_gpe::six_pointfct ()
{
    double temp = 0.0;
    
    for(int i=0; i<Ntotal; i++) temp += sqrabs(phiX[i])*sqrabs(phiX[i])*sqrabs(phiX[i]);
    
    return temp/volume;
}




// Save phi(t,0)
void boson_gpe::save_phi0t (double t)
{
    allTimes[counter] = t;
    phi0t[counter] = phiK[0];   // zero mode
    
    counter++;
}







// Computes and returns commutators commutators [phi,phi^dagger+phi] and [phi,phi^dagger-phi] in linear response at given time
vector<double> boson_gpe::save_CR ()
{
    vector<double> temp(6,0.0);
    
    comp respR(0,0);
    comp respI(0,0);
    
    switch (CRmodeResp) {
            case 0: // 0-mode
            {
                respR = (phiK_hR[0]-phiK[0])/(h*Ntotal);
                respI = -I*(phiK_hI[0]-phiK[0])/(h*Ntotal);
                break;
            }
            
            case 1: // Integrated autoresponse
            {
                //#pragma omp parallel for
                for (int n=0; n<Ntotal; n++)
                {
                    respR += hextR[n] * (phiX_hR[n]-phiX[n]);
                    respI += -I * hextI[n] * (phiX_hI[n]-phiX[n]);
                }
                respR = respR/(h*h*Ntotal);
                respI = respI/(h*h*Ntotal);
                break;
            }
            
        default:
            cout << "Error: Wrong type of mode to output for linear response (rotation_hI).\n";
            break;
    }
    
    //comp respR = (phiK_hR[0]-phiK[0])/(h*Ntotal);
    //comp respI = -I*(phiK_hI[0]-phiK[0])/(h*Ntotal);
    
    temp[0] = respR.real();
    temp[1] = respR.imag();
    temp[2] = respI.real();
    temp[3] = respI.imag();
    temp[4] = phiK[0].real();
    temp[5] = phiK[0].imag();

    return temp;
    
}






// Computes and returns commutators [phi,phi^dagger] and [phi,phi] in linear response at given time for all momenta
vector< vector<comp> > boson_gpe::save_CR_pdep ()
{
    vector<comp> temp0(2,0); // First column will be [phi,phi^dagger] and second [phi,phi]
    vector< vector<comp> > temp(Nbins,temp0);
    
    for (int n=0; n<Nbins; n++) { for (int i=0; i<4; i++) commSpectrum[4*n+i]=0.0; }
    
    if (CRmodeResp!=2 && CRmodeResp!=3) { cout << "Problem with CR output: wrong CRmodeResp for this option.\n"; return temp;}
    
    switch (CRmodeResp) {
            case 2:
            {
                for (int n=0; n<Ntotal; n++)
                {
                    commSpectrum[4*whichBin[n]] += conj(hextR_k[n]) * (phiK_hR[n]-phiK[n]);
                    commSpectrum[4*whichBin[n]+1] += -I * conj(hextI_k[n]) * (phiK_hI[n]-phiK[n]);
                    commSpectrum[4*whichBin[n]+2] += conj(hextR_k[n]) * (phiK_hR[n]-phiK[n]);
                    commSpectrum[4*whichBin[n]+3] += -I * conj(hextI_k[n]) * (phiK_hI[n]-phiK[n]);
                    //commSpectrum[4*whichBin[n]+2] += hextR_k[n] * (phiK_hR[n]-phiK[n]);
                    //commSpectrum[4*whichBin[n]+3] += -I * hextI_k[n] * (phiK_hI[n]-phiK[n]);
                }
                
                double normFac = 1/(h*h*Ntotal);
                for (int n=0; n<Nbins; n++) { for (int i=0; i<4; i++) commSpectrum[4*n+i] *= binAvFactor[n]*normFac; }
                
                break;
            }
            
            case 3:
            {
                for (int n=0; n<Ntotal; n++)
                {
                    commSpectrum[4*whichBin[n]] += phiK_hR[n]-phiK[n];
                    commSpectrum[4*whichBin[n]+1] += -I * (phiK_hI[n]-phiK[n]);
                    commSpectrum[4*whichBin[n]+2] += phiK_hR[n]-phiK[n];
                    commSpectrum[4*whichBin[n]+3] += -I * (phiK_hI[n]-phiK[n]);
                    //commSpectrum[4*whichBin[n]+2] += hextR_k[n] * (phiK_hR[n]-phiK[n]);
                    //commSpectrum[4*whichBin[n]+3] += -I * hextI_k[n] * (phiK_hI[n]-phiK[n]);
                }
                
                double normFac = 1.0/h;
                for (int n=0; n<Nbins; n++) { for (int i=0; i<4; i++) commSpectrum[4*n+i] *= binAvFactor[n]*normFac; }
                
                break;
            }
            
        default:
            cout << "Error: Wrong type of mode to output for p-dep output of linear response.\n";
            break;
    }
    
    
    for (int n=0; n<Nbins; n++)
    {
        temp[n][0] = commSpectrum[4*n] + commSpectrum[4*n+1];
        temp[n][1] = commSpectrum[4*n+2] - commSpectrum[4*n+3];
    }
    
    return temp;
}








// Computes and returns F(t,twF) for all momenta
vector< vector<comp> > boson_gpe::save_F_pdep ()
{
    vector<comp> temp0(2,0); // First column will be phi*phi^dagger and second phi*phi
    vector< vector<comp> > temp(Nbins,temp0);
    

    // Fill arrays
    for (int n=0; n<Ntotal; n++)
    {
        temp[whichBin[n]][0] += phiK[n] * conj(phi_at_twF[n]);
        temp[whichBin[n]][1] += phiK[n] * phi_at_twF[n];   // wrong choice, anomalous should be phi(p)*phi(-p)
    }
    
    // Normalize
    double normFac = 1.0/((double)Ntotal);
    for (int n=0; n<Nbins; n++) { temp[n][0] *= binAvFactor[n]*normFac;     temp[n][1] *= binAvFactor[n]*normFac; }
    
    return temp;
}





// Computes and returns commutators [phase,phase^dagger], [dichte,dichte^dagger] and [phase,dichte^dagger] in linear response at given time for all momenta
vector< vector<comp> > boson_gpe::save_DPrho_pdep ()
{
    vector<comp> temp0(4,0); // First column will be [phase,phase^dagger]  and second [dichte,dichte^dagger]
    vector< vector<comp> > temp(Nbins,temp0);
	
	// Fill dichte and phase arrays
	for (int i=0; i<Ntotal; i++)
	{
		dichte[i] = sqrabs(phiX[i]);
		phase[i] = atan2(phiX[i].imag(),phiX[i].real());
		
		dichte_hD[i] = sqrabs(phiX_hD[i]);
		phase_hD[i] = atan2(phiX_hD[i].imag(),phiX_hD[i].real());
		
		dichte_hP[i] = sqrabs(phiX_hP[i]);
		phase_hP[i] = atan2(phiX_hP[i].imag(),phiX_hP[i].real());
	}
	
	// Fourier transform to dichteK and phaseK
	fftw_execute(fft_dichte); fftw_execute(fft_phase);
	fftw_execute(fft_dichte_hD); fftw_execute(fft_phase_hD);
	fftw_execute(fft_dichte_hP); fftw_execute(fft_phase_hP);
	
    
    // Fill arrays
    for (int n=0; n<NforR2Cfourier; n++)
    {
		vector<int> indices(dim,0);
        get_indices_r2c(indices,n);
		int n_ord = get_array_position(indices);	// compute position of indices in normal order
		
        if ( indices[dim-1]!=NlistforR2C[dim-1]-1 && indices[dim-1]!=0 )
		{
	        temp[whichBin[n_ord]][0] += 2.0 * conj(hextP_k[n]) * (phaseK_hP[n]-phaseK[n]);	// [phase,phase]
	        temp[whichBin[n_ord]][1] += 2.0 * conj(hextD_k[n]) * (dichteK_hD[n]-dichteK[n]);	// [dichte,dichte]
			temp[whichBin[n_ord]][2] += 2.0 * conj(hextD_k[n]) * (phaseK_hD[n]-phaseK[n]);	// [phase,dichte]
			temp[whichBin[n_ord]][3] += 2.0 * conj(hextP_k[n]) * (dichteK_hP[n]-dichteK[n]);	// [dichte,phase] (just for checking)
		}
        else
		{
	        temp[whichBin[n_ord]][0] += conj(hextP_k[n]) * (phaseK_hP[n]-phaseK[n]);	// [phase,phase]
	        temp[whichBin[n_ord]][1] += conj(hextD_k[n]) * (dichteK_hD[n]-dichteK[n]);	// [dichte,dichte]
			temp[whichBin[n_ord]][2] += conj(hextD_k[n]) * (phaseK_hD[n]-phaseK[n]);	// [phase,dichte]
			temp[whichBin[n_ord]][3] += conj(hextP_k[n]) * (dichteK_hP[n]-dichteK[n]);	// [dichte,phase] (just for checking)
		}
	}
	
    // Normalize
    double normFac = 1/(hDP*hDP*Ntotal);
    for (int n=0; n<Nbins; n++)
	{
		temp[n][0] *= binAvFactor[n]*normFac;
		temp[n][1] *= 2.0*binAvFactor[n]*normFac;
		temp[n][2] *= 2.0*binAvFactor[n]*normFac;	// Factor 2 due to perturbing with dichte
		temp[n][3] *= binAvFactor[n]*normFac;
	}
    
    return temp;
}








// Computes and returns phase-dichte F(t,twDPF) for all momenta
vector< vector<comp> > boson_gpe::save_DPF_pdep ()
{
    vector<comp> temp0(3,0); // 1st column: phase*phase^dagger, 2nd column: dichte*dichte^dagger, 3rd column: phase*dichte^dagger
    vector< vector<comp> > temp(Nbins,temp0);
    
	// Fill dichte and phase arrays
	for (int i=0; i<Ntotal; i++)
	{
		dichte[i] = sqrabs(phiX[i]);
		phase[i] = atan2(phiX[i].imag(),phiX[i].real());
	}
	
	// Fourier transform to dichteK and phaseK
	fftw_execute(fft_dichte);
	fftw_execute(fft_phase);

    // Fill arrays
    for (int n=0; n<NforR2Cfourier; n++)
    {
		vector<int> indices(dim,0);
        get_indices_r2c(indices,n);
		int n_ord = get_array_position(indices);	// compute position of indices in normal order
		
        if ( indices[dim-1]!=NlistforR2C[dim-1]-1 && indices[dim-1]!=0 )
		{
	        temp[whichBin[n_ord]][0] += 2.0 * phaseK[n] * conj(phase_at_twDPF[n]);
	        temp[whichBin[n_ord]][1] += 2.0 * dichteK[n] * conj(dichte_at_twDPF[n]);
			temp[whichBin[n_ord]][2] += 2.0 * phaseK[n] * conj(dichte_at_twDPF[n]);
		}
        else
		{
	        temp[whichBin[n_ord]][0] += phaseK[n] * conj(phase_at_twDPF[n]);
	        temp[whichBin[n_ord]][1] += dichteK[n] * conj(dichte_at_twDPF[n]);
			temp[whichBin[n_ord]][2] += phaseK[n] * conj(dichte_at_twDPF[n]);
		}
	}	
    
    // Normalize
    double normFac = 1.0/((double)Ntotal);
    for (int n=0; n<Nbins; n++) { temp[n][0] *= binAvFactor[n]*normFac;     temp[n][1] *= binAvFactor[n]*normFac; temp[n][2] *= binAvFactor[n]*normFac; }
    
    return temp;
}







/*
// Computes autocorelation function C(t,tw)
comp boson_gpe::autocorr_C ()
{
    comp temp = comp(0.0,0.0);
 
    for (int n=0; n<Ntotal; n++)
    {
        temp += phiX[n]*conj(phi_at_tw[n]);
    }
    
    return temp/(volume);
}//autocorr_C




// Computes (integrated) autocorelation response function R(t,tw)
comp boson_gpe::autocorr_R ()
{
    comp temp = comp(0.0,0.0);
    
    for (int n=0; n<Ntotal; n++)
    {
        temp += phiX_h[n]*conj(hext[n]);
    }
    
    return temp/(volume*h*h);
}//autocorr_R
*/





//[ Output of the initial conditions for each run
void boson_gpe::write_IC (int run, ofstream & output)
{
    output << "# Run " << run << ":\n" ;
    
    for (int n=0; n<Ntotal; n++)
    {
        output << phiX[n].real() << '\t' << phiX[n].imag() << '\n';
    }
    
    output << "\n\n";
    
}//]





//[ Output of the distribution function
void boson_gpe::write_spec (ofstream & output)
{
    fftw_execute(fft_xk);
    
    for(int i=0; i<Nbins; i++) spectrum[i]=0.0;
    
    for(int i=0; i<Ntotal; i++) { spectrum[whichBin[i]] += abs(phiK[i])*abs(phiK[i]) ;}    // add psi_k*psi_k^\dagger to corresponding bin
    
    for(int i=0; i<Nbins; i++) spectrum[i]*=binAvFactor[i]/volume;     // average over number of k's in each bin
    
    for (int i=0; i<Nbins; i++)
    {
        if (momPerBin[i]!=0)
        {
            output << kBin[i] << '\t' << spectrum[i] << endl;
        }
    }
    
}//]




//[ Output of the Bogoliubov distribution function
void boson_gpe::write_bogspec (ofstream & output)
{
    if (bogHappened==0)
    {
        bogSpectrum = new my_double [Nbins];
        bogK = new comp [Ntotal];
        bogHappened = 1;
    }
    
    bogtrafo();
    
    for(int i=0; i<Nbins; i++) bogSpectrum[i]=0.0;
    
    for(int i=0; i<Ntotal; i++) { bogSpectrum[whichBin[i]] += abs(bogK[i])*abs(bogK[i]) ;}    // add psi_k*psi_k^\dagger to corresponding bin
    
    for(int i=0; i<Nbins; i++) bogSpectrum[i]*=binAvFactor[i]/volume;     // average over number of k's in each bin
    
    for (int i=0; i<Nbins; i++)
    {
        if (momPerBin[i]!=0)
        {
            output << kBin[i] << '\t' << bogSpectrum[i] << endl;
        }
    }
    
}//]





// Output of the total particle number, energy and magnetization
void boson_gpe::write_op (ofstream & output, double time)
{
    output << time << '\t' << magnetization() << '\t' << totaln() << '\t'
    << kinetic() << '\t' << potenergy() << '\t' << densityVariance() << '\t' << angularVariance() << endl;
}//write_op




// Output of the on-site n-point functions
void boson_gpe::write_np (ofstream & output, double time)
{
    output << time << '\t' << two_pointfct() << '\t' << two_pointfct_anomalous().real() << '\t' << two_pointfct_anomalous().imag() << '\t'
    << four_pointfct() << '\t' << six_pointfct() << endl;
}//write_np




// Output of phi(t,p=0)
void boson_gpe::write_phi0t (ofstream & output)
{
    
    for (int n=0; n<counter; n++)
    {
        output << allTimes[n] << '\t' << phi0t[n].real() << '\t' << phi0t[n].imag() << '\t' << sqrabs(phi0t[n]) << '\n';
    }
    
}//write_phi0t


// Output of linear response data
void boson_gpe::write_CR (ofstream & output, vector<double> & timesCR, vector< vector<double> > & outCR)
{
    for (int i=0; i<timesCR.size(); i++)
    {
        output << tw << '\t' << timesCR[i] << '\t' << outCR[i][0] << '\t' << outCR[i][1] << '\t' << outCR[i][2] << '\t' << outCR[i][3]
        << '\t' << outCR[i][4] << '\t' << outCR[i][5] << endl;
    }
    
}//write_cr



//[ Output of the density-phase distribution functions
void boson_gpe::write_DPspec (ofstream & output)
{
	// Fill dichte and phase arrays
	for (int i=0; i<Ntotal; i++)
	{
		dichte[i] = sqrabs(phiX[i]);
		phase[i] = atan2(phiX[i].imag(),phiX[i].real());
	}
	
	// Fourier transform to dichteK and phaseK
	fftw_execute(fft_dichte);
	fftw_execute(fft_phase);
	
	// Create output arrays and initialize to zero
	double * spectrumPP = new double [Nbins];	// phase-phase correlator
	double * spectrumDD = new double [Nbins];	// dichte-dichte
	comp * spectrumDP = new comp [Nbins];	// dichte-phase
	
	for (int n=0; n<Nbins; n++)
	{
		spectrumPP[n]=0;
		spectrumDD[n]=0;
		spectrumDP[n]=0;
	}
	
    // Average over bins
    //#pragma omp parallel for
    for (int n=0; n<NforR2Cfourier; n++)
    {
		vector<int> indices(dim,0);
        get_indices_r2c(indices,n);
		int n_ord = get_array_position(indices);	// compute position of indices in normal order
        
        if ( indices[dim-1]!=NlistforR2C[dim-1]-1 && indices[dim-1]!=0 )
		{
			spectrumPP[whichBin[n_ord]] += 2.0 * sqrabs(phaseK[n]);
			spectrumDD[whichBin[n_ord]] += 2.0 * sqrabs(dichteK[n]);
			spectrumDP[whichBin[n_ord]] += 2.0 * phaseK[n] * conj(dichteK[n]);
		}
        else
		{
			spectrumPP[whichBin[n_ord]] += sqrabs(phaseK[n]);
			spectrumDD[whichBin[n_ord]] += sqrabs(dichteK[n]);
			spectrumDP[whichBin[n_ord]] += phaseK[n] * conj(dichteK[n]);
		}
    }
    
    for (int i=0; i<Nbins; i++)
	{
		spectrumPP[i] *= binAvFactor[i]/volume;
		spectrumDD[i] *= binAvFactor[i]/volume;
		spectrumDP[i] *= binAvFactor[i]/volume;
	} 
    
	// Output
    for (int i=0; i<Nbins; i++)
    {
        if (momPerBin[i]!=0)
        {
            output << kBin[i] << '\t' << spectrumPP[i] << '\t' << spectrumDD[i] << '\t' << spectrumDP[i].real() << '\t' << spectrumDP[i].imag() << endl;
        }
    }
	
	// Delete arrays
	delete [] spectrumPP;
	delete [] spectrumDP;
	delete [] spectrumDD;
    
}//write_DPspec



/*
// Output of the autocorrelation C and R functions
void boson_gpe::write_CR (ofstream & output, double time)
{
    comp respR = (phiK_hR[0]-phiK[0])/(h*Ntotal);
    comp respI = -I*(phiK_hI[0]-phiK[0])/(h*Ntotal);
    
    output << tw << '\t' << time << '\t' << respR.real() << '\t' << respR.imag() << '\t' << respI.real() << '\t' << respI.imag()
                 << '\t' << phiK[0].real() << '\t' << phiK[0].imag() << endl;
}//write_cr
 */


/*
// Output of the autocorrelation C and R functions
void boson_gpe::write_CR (ofstream & output, double time)
{
    comp autoC = autocorr_C();
    comp autoR = autocorr_R();
    
    output << tw << '\t' << time << '\t' << autoC.real() << '\t' << autoC.imag() << '\t' << autoR.real() << '\t' << autoR.imag() << endl;
}//write_cr
*/










/*

 void boson_gpe::outphis (ofstream & output, double time)
 {
     output << tw << '\t' << time << '\t';
     
     for (int n=0; n<1; n++)
     {
         output << phiX[n].real() << '\t' << phiX_h[n].real() << '\t' << phiX[n].imag() << '\t' << phiX_h[n].imag()
                << '\t' << hext[n].real() << '\t' << hext[n].imag() << '\t';
     }
     
     output << endl;
 }


*/







// Divides the physical momentum space into bins
void boson_gpe::bin_log()
{
    momPerBin = new int [Nbins];
    binAvFactor = new my_double [Nbins];
    whichBin = new int [Ntotal];
    kBin = new my_double [Nbins];
    spectrum = new my_double [Nbins];
    
    for(int i=0; i<Nbins; i++) momPerBin[i]=0.0;
    
    int maxN = Nlist[0];
    for (int i=1; i<dim; i++)
    {
        if (Nlist[i]>maxN) maxN = Nlist[i];
    }// longest side of grid
    
    double smallest = 2*sin(M_PI/maxN)/ls;      // smallest momentum on grid (apart from 0)
    
    double biggest = sqrt(dim*4.0)/ls;          // largest momentum on grid
    
    logstep = (log(biggest)-log(smallest))/(Nbins-1);    // log width of bin
    double smlog = log(smallest);
    double linstep=exp(logstep);
    
    
    kBin[0]=0.0;
    double binleft = smallest; // *exp(0.5*step);
    for(int i=1; i<Nbins; i++) { kBin[i]=binleft; binleft*=linstep; }     // fill kBin with left momentum of each bin
    

    double fakt = 4.0/(ls*ls);
    
    vector<int> indices(dim,0);
    
    for (int n=0; n<Ntotal; n++)
    {
        get_indices(indices,n);
        
        double ksq = 0; // Physical momentum squared
        for (int i=0; i<dim; i++) ksq += sin(M_PI*indices[i]/Nlist[i])*sin(M_PI*indices[i]/Nlist[i]);
        ksq *= fakt;
        
        double logk=0.5*log(ksq);
        int mybin=(int)((logk-smlog)/logstep)+1;
        if (ksq==0.0) mybin=0;
        if (mybin>=Nbins) mybin=Nbins-1;        // calculate corresponding bin for given k
        
        momPerBin[mybin]++;                         // 1 momentum more in bin mybin
        whichBin[n] = mybin;                      // associate bin to given fourier-momentum
        
    }
    
    for(int i=0; i<Nbins; i++)
    {
        if (momPerBin[i]!=0) binAvFactor[i]=1.0/momPerBin[i];
        else binAvFactor[i]=0.0;
    }
     
    
}// bin_log


























/*
// Try convolution in two different ways to check accuracy of FFTW method
void prueba_convol ()
{
    for (int n=0; n<Ntotal; n++) density[n] = sqrabs(phiX[n]);
    
    vector<int> indicesX(dim,0);
    vector<int> indicesY(dim,0);
    vector<int> indicesDiff(dim,0);
    int posInt, temp;
    
    // Direct sum
    vector<double> direct_sum(Ntotal,0);
    
    for (int n=0; n<Ntotal; n++)
    {
        get_indices(indicesX,n);
        
        for (int m=0; m<Ntotal; m++)
        {
            get_indices(indicesY,m);
            
            for (int i=0; i<dim; i++)
            {
                temp = (indicesX[i]-indicesY[i]) % (Nlist[i]);
                if (temp<0) temp += Nlist[i];
                indicesDiff[i] = temp;
            }
            
            posInt = get_array_position(indicesDiff);
            
            direct_sum[n] += Vint[posInt]*density[m];
        }
    }
    
    
    
    // With FFT
    fftw_execute(fft_density);
    
    for (int n=0; n<NforR2Cfourier; n++) densityK[n] = densityK[n]*VintK[n];
    
    fftw_execute(fft_convolution);
    
    
    // Out
    for (int n=0; n<Ntotal; n++) {
        cout << n << '\t' << direct_sum[n] << '\t' << density[n]/Ntotal << endl;
    }
}
*/





















