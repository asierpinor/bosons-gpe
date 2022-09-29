#ifndef BOSON_GPE_H
#define BOSON_GPE_H

using namespace std;



//This makes a list of input variables with their type, as well as a running index {1,2,...,varnum}. See meaning of variables in .sh file
#define varlist VAR(N,int,1) VAR(DIM,int,2) VAR(Nx,int,3) VAR(Ny,int,4) VAR(Nz,int,5) VAR(ls,double,6) \
    VAR(m,double,7) VAR(g,double,8) VAR(gn,double,9) VAR(glong,double,10) VAR(expphi,double,11) VAR(explong,double,12) VAR(h,double,13) VAR(hDP,double,14) \
    VAR(IC,int,15) VAR(amp,double,16) VAR(Qmax,double,17) VAR(Qmin,double,18) VAR(tmax,double,19) VAR(dt,double,20) \
    VAR(iter,int,21) VAR(tw,double,22) VAR(twF,double,23) VAR(twDPrho,double,24) VAR(twDPF,double,25) VAR(bsize,int,26) VAR(tnum,int,27) VAR(seed,long int,28) VAR(nfile,int,29) VAR(Nbins,int,30) \
    VAR(spstep,double,31) VAR(bogstep,double,32) VAR(opstep,double,33) VAR(npstep,double,34) VAR(phi0step,double,35) VAR(Fstep,double,36) \
    VAR(CRstep,double,37) VAR(DPspstep,double,38) VAR(DPFstep,double,39) VAR(DPrhostep,double,40) VAR(CRmodeResp,int,41) VAR(cphysMom,int,42) VAR(Pcutoff,double,43)
#define strlist STR(folder,44)
#define varnum 44

//

// Define extern variables of previous list of variables
#define VAR(aa,bb,cc) extern bb aa;
varlist
#undef VAR
extern char folder[1024];
//



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





class boson_gpe {
    
	// Settings
    vector<int> Nlist, NlistforR2C;
    int dim;            // Lattice dimension: depends on input
    int Ntotal;
    int NforR2Cfourier;
    double volume;
    
    // For initial conditions
    gsl_rng * r;
    double sigma;
    
    // For integration of GPE
    comp I;
    comp * phiX, * phiK;
    fftw_plan fft_xk, fft_kx;
    comp * kTable;
	
	// Interactions
    int gTerms;     // selects potential term based on the values of the g's
    
	// Long-range interactions
    double * Vint, * density;  // Long-range interaction and |psi|^2
    comp * VintK, * densityK;
    fftw_plan fft_Vint, fft_density, fft_convolution;
    
    
	// ------------
	// Commutators
	// ------------
	
    // For dynamics with h
    comp * phiX_hR, * phiK_hR;
    comp * phiX_hI, * phiK_hI;        // Same but with external field hR or hI turned on
    fftw_plan fft_xk_hR, fft_kx_hR;
    fftw_plan fft_xk_hI, fft_kx_hI;
    comp * hextR, * hextI;  // Array of x-dependent external field (both are real in reality, define them like this to avoid r2c output size)
    comp * hextR_k, * hextI_k;
    fftw_plan fft_hextR, fft_hextI;
	
	
    // For dynamics with hDP
    comp * phiX_hD, * phiK_hD;
    comp * phiX_hP, * phiK_hP;
    fftw_plan fft_xk_hD, fft_kx_hD;
    fftw_plan fft_xk_hP, fft_kx_hP;
    comp * hextD, * hextP;  // Array of x-dependent external field (both are real in reality, define them like this to avoid r2c output size)
    comp * hextD_k, * hextP_k;
    fftw_plan fft_hextD, fft_hextP;
    
	
	// ----------------
	// Anticommutators
	// ----------------
    
    // For F output
    int nTimes;                         // Number of times measured
    int counter;                        // Saves running number of times that phi0t has been saved
    vector<double> allTimes;        // save all times in this file
    vector<comp> phi0t;             // save phi(t,p=0) here
    comp * commSpectrum;   // temporary array to hold the bin-averaged p-dependent commutators
    comp * phi_at_twF;      // Save phi(twF,p) at time t=twF for computation of F(t,twF)
	comp * phase_at_twDPF;      // Save phase(twDPF,p) at time t=twF for computation of F_phase(t,twDPF)
	comp * dichte_at_twDPF;      // Save dichte(twDPF,p) at time t=twF for computation of F_phase(t,twDPF)
	
	
	// --------------
	// Density-Phase
	// --------------
	
	// For density-phase output: spectrum, F or rho
    double * dichte, * phase;
    comp * dichteK, * phaseK;
	fftw_plan fft_dichte, fft_phase;
	
	// For density-phase rho output
    double * dichte_hD, * phase_hD;
	double * dichte_hP, * phase_hP;
    comp * dichteK_hD, * phaseK_hD;
	comp * dichteK_hP, * phaseK_hP;
	fftw_plan fft_dichte_hD, fft_phase_hD;
	fftw_plan fft_dichte_hP, fft_phase_hP;
    
	
	
    
    // Bogoliubov modes
    my_double k0;                       // Transition momentum
    comp * bogK;                        // Bogoliubov field modes
    my_double * bogSpectrum;            // Spectrum of Bogoliubov quasiparticles
    my_double alphak, yk, uk, vk;       // Parameters for Bogoliubov transformation
    int bogHappened;                    // To check if arrays have been allocated or not
	
	
	
    // For binning of momenta
    double logstep;
    int * momPerBin;          // number of k's in one bin
    my_double * binAvFactor;        // 1/(number of k's in one bin)
    int * whichBin;                 // associates a bin number to each fourier-momentum (i,j,k)
    my_double * kBin;               // contains left values of k that define bins
    my_double * spectrum;           // will hold final binned spectrum
    
    
    
    
public:
    
    boson_gpe(gsl_rng * rndm);
    ~boson_gpe();
    
    inline double sqr (double a) {return a*a;}
    inline double sqrabs (comp z) {return sqr(z.real())+sqr(z.imag());}    // Absolute value squared of complex number

    void get_indices (vector<int> & indices, int n);    // Gets indices of each dimension assuming row-order
    void get_indices_r2c (vector<int> & indices, int n);    // Gets indices of each dimension assuming row-order for r2c lattice dimension
    int get_array_position (vector<int> indices);    // Computes the position in the array given the individual indices of each dimension
    
    void fill_Vint ();    // Fill long-range interaction matrix (distances taken with respect to zero)
    
    vector<double> getPhysMom (int binNumber);  // Outputs physical momentum and momenta per bin for given bin number n
    
    
///////////////            INITIAL CONDITIONS          //////////////////
    
    
    void addTWAnoise ();    // Choose which noise to add to initial condition
    void addTWAnoise_box ();    // Add TWA noise to initial conditions
    void addTWAnoise_homoDensity ();    // Homogeneous density, only phase noise
	void addTWAnoise_oneQuadrant ();    // Restrict Psi(x) to positive real and imag parts
	void addTWAnoise_randomPhase ();    // Constant density, random phase
    
    
    
///////////////////            DYNAMICS          ////////////////////////
    
    
    void triple_system ();		// Triples system into one evolving as before, one perturbed with hextR and one perturbed with hextI. Happens after waiting time tw.
	void triple_system_DP ();		// Triples system into one evolving as before, one perturbed with hextD and one perturbed with hextP. Happens after waiting time twDPrho.
	
    void save_phi_at_twF ();    // Saves state phi(twF,p) in array phi_at_twF
	void save_DP_at_twDPF ();    // Saves phase(twDPF,p) and dichte(twDPF,p) in arrays phase_at_twDPF and dichte_at_twDPF, respectively.
    
    void fill_kTable ();    // Fill kTable with kinetic exp factor for integration step
    
    void dynamics ();    // Calculates phi(t+dt) from GPE with split-step method
    
    void dynamics_pot_000 ();    // Multiply time evolution with potential part for (glong,gn,g)=(0,0,0)
    void dynamics_pot_001 ();    // Multiply time evolution with potential part for (glong,gn,g)=(0,0,*)
    void dynamics_pot_010 ();    // Multiply time evolution with potential part for (glong,gn,g)=(0,*,0)
    void dynamics_pot_011 ();    // Multiply time evolution with potential part for (glong,gn,g)=(0,*,*)
    void dynamics_pot_100 ();    // Multiply time evolution with potential part for (glong,gn,g)=(*,0,0)
    void dynamics_pot_101 ();    // Multiply time evolution with potential part for (glong,gn,g)=(*,0,*)
    void dynamics_pot_110 ();    // Multiply time evolution with potential part for (glong,gn,g)=(*,*,0)
    void dynamics_pot_111 ();    // Multiply time evolution with potential part for (glong,gn,g)=(*,*,*)
    
    void dynamics_kin ();    // Multiply time evolution with kinetic part
    
    void interaction_longrange ();    // Computes interaction part of GPE for long-range interactions
    
    //void dynamics_with_h ();    // Calculates phi(t+dt) from GPE with split-step method
    //void dynamics_pot_001_with_h ();    // Multiply time evolution (with h) with potential part for (glong,gn,g)=(0,0,*)
    //void dynamics_kin_with_h ();    // Multiply time evolution (with h) with kinetic part
    //void dynamics_h();    // Add external field part to time evolution
    
    void dynamics_phihR ();
    void dynamics_phihI ();
    void rotation_hR ();
    void rotation_hI ();
	
    void dynamics_phihD ();
    void dynamics_phihP ();
    void rotation_hD ();
    void rotation_hP ();
    

    
////////////////////            OUTPUT          /////////////////////////

    
    void bogtrafo ();    // Calculates the Bogoliubov modes from the free modes
    
    double magnetization ();
    double totaln();
    double kinetic ();
    double potenergy ();
	double densityVariance ();
	double angularVariance ();
    
    double two_pointfct ();    // Correlation function < phi_x phi_x^* >
    comp two_pointfct_anomalous ();    // Correlation function < phi_x phi_x >
    double four_pointfct ();    // Correlation function < (phi_x phi_x^*)^2 >
    double six_pointfct ();    // Correlation function < (phi_x phi_x^*)^3 >
    
    void save_phi0t (double t);    // Save phi(t,0)
    vector<double> save_CR ();   // Save linear response info
    vector< vector<comp> > save_CR_pdep ();   // Save momentum-dependent linear response info
    vector< vector<comp> > save_F_pdep ();   // Save momentum-dependent F(t,twF)
	vector< vector<comp> > save_DPrho_pdep ();   // Save momentum-dependent phase-dichte rho(t,twDPrho)
	vector< vector<comp> > save_DPF_pdep ();   // Save momentum-dependent phase-dichte F(t,twDPF)
    
    //comp autocorr_C ();    // Computes autocorelation function C(t,tw)
    //comp autocorr_R ();    // Computes (integrated) autocorelation response function R(t,tw)
    
    
    void write_IC (int run, ofstream & output);     // Output of the initial conditions for each run
    void write_spec (ofstream & output);    // Output of the distribution function
    void write_bogspec (ofstream & output);    // Output of the Bogoliubov distribution function
    void write_op (ofstream & output, double time);    // Output of the total particle number, energy and magnetization
    void write_np (ofstream & output, double time);    // Output of the on-site n-point functions
    void write_phi0t (ofstream & output);     // Output of phi(t,p=0)
    //void write_CR (ofstream & output, double time);    // Output of the autocorrelation C and R functions
    void write_CR (ofstream & output, vector<double> & timesCR, vector< vector<double> > & outCR);
	
	void write_DPspec (ofstream & output);    // Output of the density-phase distribution functions
    
    void outphis (ofstream & output, double time);
    
    void bin_log();     // Divides the physical momentum space into bins
    
    
    
    
};
//]












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


























#endif // BOSON_GPE_H


