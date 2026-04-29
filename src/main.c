#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <fftw3-mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "parameters.h"
#include "lpt.h"
#include "msg.h"
#include "power.h"
#include "comm.h"
#include "pm.h"
#include "cola.h"
#include "fof.h"
#include "write.h"
#include "timer.h"
#include "mem.h"
#include "move.h"
#include "subsample.h"
#include "coarse_grid.h"
#include "solve_growth.h"

//#define LOGTIMESTEP 1

float f_heart(float x, float y, float z) {
    float a;
a = x * x + 9.0f / 4.0f * y * y + z * z - 1;
    return a * a * a - x * x * z * z * z - 9.0f / 80.0f * y * y * z * z * z;
}


float h_heart(float x, float z) {
float y;
    for ( y = 1.0f; y >= 0.0f; y -= 0.001f)
        if (f_heart(x, y, z) <= 0.0f)
            return y;
    return 0.0f;
}


void print_heart() {
float z,x,v,y0,ny,nx,nz,nd,d;
    for ( z = 1.5f; z > -1.5f; z -= 0.05f) {
        for ( x = -1.5f; x < 1.5f; x += 0.025f) {
             v = f_heart(x, 0.0f, z);
            if (v <= 0.0f) {
                 y0 = h_heart(x, z);
                 ny = 0.01f;
                 nx = h_heart(x + ny, z) - y0;
                 nz = h_heart(x, z + ny) - y0;
                 nd = 1.0f / sqrtf(nx * nx + ny * ny + nz * nz);
                 d = (nx + ny - nz) * nd * 0.5f + 0.5f;
                putchar(".:-=+*#%@"[(int)(d * 5.0f)]);
            }
            else
                putchar(' ');
        }
        putchar('\n');
    }
}



int mpi_init(int* p_argc, char*** p_argv);
void fft_init(int threads_ok);
void snapshot_time(const float aout, const int iout, 
		   Particles const * const particles, 
		   Snapshot * const snapshot, 
		   const char fof_filename[], 
		   const char subsample_filename[], 
		   const char cgrid_filename[], const int cgrid_nc,
		   void* const mem1, const size_t size1,
		   const int write_longid,
		   const double fof_linking_factor);

void write_slice(const char filebase[], Particle* p, const int np, const float dz, const float boxsize);
void write_snapshot_slice(const char filebase[], Snapshot const * const snapshot, const float boxsize);

// xiaodong: some functions for the comov_r calculation
// # functions 
//
//
double gb_omegam;
double Simpson(double (* fun)(), double x1, double x2, int nstep){
	int i;
	double delta_x, f, x, rlt;
	nstep *= 2;
	delta_x = (x2 - x1) /(double)nstep;
	rlt = (*fun)(x1) + (*fun)(x2);
	for(i=2; i<=nstep-2; i+=2){
		x = x1 + delta_x * i; f = (* fun)(x); rlt += 2*f;
	}
	for(i=1; i<=nstep-1; i+=2){
		x = x1 + delta_x * i; f = (* fun)(x); rlt += 4*f;
	}
	rlt *= (delta_x /3.);
	return rlt;
}


double Hz(double z){
	double y;
	y = gb_omegam * (1.+z)*(1.+z)*(1.+z) + (1. - gb_omegam) * pow(1.+z, 3.*(1.+gb_w));
	y = sqrt(y);
	return y;
}

double inv_Hz(double z){
	return 299792.458/Hz(z)/100.;
}

double fun_x2(double x){return x*x;}

// structure storing particles who needed to be checked (whether they cross_boundary)
struct chain_1000 {
	long long ids[1000]; int nowi; // index marking the current last particle,
	float xyzvs[1000][3];  struct chain_1000 * next; } ; // chain consists of arrays with length 1000
void free_chain_1000(struct chain_1000 * head){
	struct chain_1000 * nowp=head, * nextp;
	while(nowp!= NULL){ nextp = nowp->next; free(nowp); nowp = nextp; }	
}
void print_chain_1000(struct chain_1000 * head){
	struct chain_1000 * nowp = head;
	while(nowp!=NULL){ 
		int i; printf("(%d-particles): ", nowp->nowi); for(i=0;i<nowp->nowi;i++) printf("%d ", (int)(nowp->ids[i])); nowp=nowp->next;} 
	putchar('\n');
}

int main(int argc, char* argv[])
{
  const int multi_thread= mpi_init(&argc, &argv);
  msg_init();
                                                      timer_set_category(Init);

  //
  // Initialization / Memory allocation
  //						      
  Parameters param;
  read_parameters(argc, argv, &param);



  const int nc_factor= param.pm_nc_factor;
  const double OmegaM= param.omega_m;
  const double OmegaLambda= 1.0 - OmegaM;
  const double Hubble= param.h;
  const double sigma8= param.sigma8;

  // xiaodong: variables for lightcone output begin
    int only_output_1eighth = param.only_output_1eighth;
    double zmax, rmax=0.;
    int nshift=0, shift1, shift2; 
  // xiaodong: variables for lightcone output end

  msg_set_loglevel(param.loglevel);

  confirm_parameters(&param, argv[argc-1]);
  gb_solve_D_E_Qs(); // xiaodong: check that

  fft_init(multi_thread);
  comm_init(param.pm_nc_factor*param.nc, param.nc, param.boxsize);

  const int nsteps= param.ntimestep;
  const double a_final= param.a_final;

#ifndef LOGTIMESTEP
  const double da= a_final/nsteps;
  /* NOTE: To maintain consistency with the original COLA implementation,
   * a_init is forced to equal da (i.e., a_init = a_final/nsteps).
   * The user-provided a_init parameter is currently ignored to avoid
   * time-step inconsistencies. See developing/a_init.md for details.
   */
  const double a_init= da;
#else
  const double a_init= param.a_init;
#endif


  // xiaodong variables begin
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    struct chain_1000 * head_precheck=NULL, * p_precheck=NULL, *head_nowcheck=NULL, * p_nowcheck=NULL;
    int is_first_lightcone_shell = 1, do_check_cross=1;
    float a1,a2,a3,amid, z1,z2,z3,zmid, r1,r2,rmid, r_upbound, r_lowbound, r_nextstep_lowbound;
    float max_possible_v = 1500., max_travel_distance;
    long long n_check_cross_skip_total =0, n_check_cross_output_total =0;
    long long n_write_total = 0;

    char check_cross_filename[6000];
    FILE * fp_check_cross_test=NULL;
    if (param.use_lightcone)
    {
        sprintf(check_cross_filename,"%s_check_cross_test_myrank%d.txt",param.lightcone_basename,myrank);
        if(myrank==0) fp_check_cross_test = fopen(check_cross_filename, "w");
    }
  // xiaodong variables end


  
  //power_init("camb0_matterpower.dat", a_init, sigma8, OmegaM, OmegaLambda);
  power_init(param.power_spectrum_filename, a_init,
	     sigma8, OmegaM, OmegaLambda);


  Memory mem;
  allocate_shared_memory(param.nc, nc_factor, param.np_alloc_factor, &mem);
  lpt_init(param.nc, mem.mem1, mem.size1);
  const int local_nx= lpt_get_local_nx();

  Particles* particles=
     allocate_particles(param.nc, local_nx, param.np_alloc_factor);
  Snapshot* snapshot= allocate_snapshot(param.nc, local_nx,
				  particles->np_allocated, mem.mem2, mem.size2);
  snapshot->boxsize= param.boxsize;
  snapshot->omega_m= OmegaM;
  snapshot->h= Hubble;
  //strncpy(snapshot->filename, param.snapshot_filename, 64);
  snapshot->filename= param.snapshot_filename;
  

  pm_init(nc_factor*param.nc, nc_factor, param.boxsize, param.np_alloc_factor,
	  mem.mem1, mem.size1, mem.mem2, mem.size2);
  fof_init(particles->np_allocated, param.nc, mem.mem1, mem.size1);
  subsample_init(param.subsample_factor, param.random_seed);

  const int nout= param.n_zout;
  float* aout= malloc(sizeof(float)*nout);
  for(int i=0; i<nout; i++) {
    aout[i] = (float)(1.0/(1 + param.zout[i]));
    msg_printf(verbose, "zout[%d]= %lf, aout= %f\n", 
	       i, param.zout[i], aout[i]);
  }


  //
  // Many realizations with different initial conditions
  //
  for(int irealization=0; irealization<param.nrealization; irealization++) {
    MPI_Barrier(MPI_COMM_WORLD);
    msg_printf(verbose, "\n%d/%d realization.\n", 
	       irealization+1, param.nrealization);
    int seed= param.random_seed + irealization;

    int iout= 0;

    // Sets initial grid and 2LPT desplacement

    timer_set_category(LPT);
    lpt_set_displacement(a_init, OmegaM, seed,
			 param.boxsize, particles);
    snapshot->seed= seed;

    if(param.init_filename) {
      // Writes initial condition to file for e.g. Gadget N-body simulation
      char filename[2560]; // TODO: use variable length for filename
      sprintf(filename, "%s_%05d", param.init_filename, seed);
      set_noncola_initial(particles, snapshot);
      write_snapshot(filename, snapshot, param.write_longid);
    }

    timer_set_category(COLA);

#ifndef LOGTIMESTEP
    //particles->a_v= 0.5*da;
    //  I thought 0.5*da is the natural initial time for leap frog integration,
    //  but da gives much better matter power. I don't know why. (Feb 12, 2014)
    particles->a_v= da;
    msg_printf(info, "timestep linear in a\n");
#else
    const double loga_init= log(a_init);
    const double dloga= (log(a_final) - log(a_init))/nsteps;
    particles->a_v= exp(log(a_init) - 0.5*dloga);
    msg_printf(info, "timestep linear in loga\n");
#endif


    // ----------------------------------------------------
    // xiaodong's code for lightcone output begin
    // *** beging of setting
    if (param.use_lightcone)
    {
        // 1: set omegam
        gb_omegam = OmegaM;
        // 2: set maximal radius of the lightcone
        zmax = param.lightcone_zmax;
        // *** end of setting

        rmax = Simpson(&inv_Hz, 0., zmax, 512);
        nshift = ceil(rmax /  param.boxsize) + 1 ;

        shift1 = -nshift; shift2 = nshift-1;
        //shift1 = -nshift-1; shift2 = nshift;
        //shift1 = -nshift-2; shift2 = nshift+1;
        //
        if(only_output_1eighth) shift1 = 0;
        if(myrank==0 || 1){
        printf("  (write_lightcone) myrank = %d  lightcone_basename = %s   da = %.4lf  \n ", myrank, param.lightcone_basename, da );
        printf(" (write_lightcone) myrank = %d  zmax, rmax = %lf %lf \n ", myrank, zmax, rmax);
        printf(" (write_lightcone) nshift = rmax/boxize = %lf / %lf = %d \n ", rmax, param.boxsize, nshift);
        printf(" (write_lightcone) xyzmin/max after shift =  %lf,  %lf \n ", 0+shift1*param.boxsize, param.boxsize * (shift2+1));
        }
        //printf("Integrate of x from 0 to 1: %f\n", Simpson(&fun_x2, 0, 3, 100 ));
        //printf("gb_omegam = %lf \n", gb_omegam);
        //printf("Comov_r at z = 0., 1.0: %lf %lf\n", Simpson(&inv_Hz, 0, 0.5, 100 ),  Simpson(&inv_Hz, 0, 1., 100 ));
        // ----------------------------------------------------
    }

    // xiaodong's code for lightcone output end
    // ----------------------------------------------------


    //
    // Time evolution loop
    //
    //   TODO: allow nstep=1?
    //


    if(nout > 0 && nsteps > 1 && a_final > a_init) {
      msg_printf(normal, "Time integration a= %g -> %g, %d steps\n", 
		 a_init, a_final, nsteps);
    for (int istep=1; istep<=nsteps; istep++) {
	msg_printf(normal, "\n\nTimestep %d/%d\n", istep, nsteps);
      
                                                            timer_start(comm);
        // move particles to other nodes
        move_particles2(particles, param.boxsize, mem.mem1, mem.size1 );

                                                            timer_stop(comm);

        pm_calculate_forces(particles); 

#ifndef LOGTIMESTEP
	double avel0= (istep-0.5)*da;
	double apos0=  istep*da;
	
	double avel1= (istep+0.5)*da;
	double apos1= (istep+1.0)*da;
#else
	float avel0= exp(loga_init + (istep-0.5)*dloga);
	float apos0= exp(loga_init + istep*dloga);
	
	float avel1= exp(loga_init + (istep+0.5)*dloga);
	float apos1= exp(loga_init + (istep+1)*dloga);
#endif
      
	while(iout < nout && avel0 <= aout[iout] && aout[iout] <= apos0) {
	  // Time to write output
	  snapshot_time(aout[iout], iout, particles, snapshot, param.fof_filename, param.subsample_filename, param.cgrid_filename, param.cgrid_nc, mem.mem1, mem.size1, param.write_longid, param.fof_linking_factor);
	  iout++;
	}
	//printf("apos0 apos1 = %f %f \n", avel0, avel1);
	//if(iout >= nout) break;

	// Leap-frog "kick" -- velocities updated
	cola_kick(particles, OmegaM, avel1);

	while(iout < nout && apos0 < aout[iout] && aout[iout] <= avel1) {
	  // Time to write output
	  snapshot_time(aout[iout], iout,  particles, snapshot, param.fof_filename, param.subsample_filename, param.cgrid_filename, param.cgrid_nc, mem.mem1, mem.size1, param.write_longid, param.fof_linking_factor);
	  iout++;
	}
	//if(iout >= nout) break;

	// Leap-frog "drift" -- positions updated
	cola_drift(particles, OmegaM, apos1);

    // ----------------------------------------------------
    // xiaodong's code for lightcone output begin
    /* definition of Particles struct
    typedef float float3[3];
    typedef struct {
    float x[3]; 
    float dx1[3]; // ZA displacement // 
    float dx2[3]; // 2LPT displacement
    float v[3];   // velocity
    long long id;
    } Particle;
    typedef struct {
    Particle* p; float3* force; float a_x, a_v;
    int np_local, np_allocated; long long np_total; float np_average; } Particles;
    */

	// will write the particle if its observed z is between z1 and z2 (i.e. consisitent with current z)

	//a1=  (istep+1. - 0.5 - 0.2 * (float)nsteps / 20.) *da; 	// allowing 20% overlap between two timesteps
	//a2=  (istep+1. + 0.5 + 0.2 * (float)nsteps / 20. ) *da; 	// allowing 20% overlap between two timesteps
    if (param.use_lightcone){
        a1=  (istep+1. - 0.5) * da; // - 0.2 * (float)nsteps / 20.) *da;
        a2=  (istep+1. + 0.5) * da; //+ 0.2 * (float)nsteps / 20. ) *da;
        amid = (a1+a2) / 2.; a3=  a2; //(istep+1. - 2.5 - 0.2 * (float)nsteps / 20.) *da;  // a very small a, large r: use to decide when to start output
        zmid = 1/amid - 1.; z1 = 1./a1 - 1.; z2 = 1./a2 - 1.;  z3 = 1./a3 - 1.;
        rmid = Simpson(&inv_Hz, 0., zmid, 512); 
        r1 = Simpson(&inv_Hz, 0., z1, 512); r2 = Simpson(&inv_Hz, 0., z2, 512); //r3 = Simpson(&inv_Hz, 0., z3, 512);
        // a1 is small; r1 is large; r1 is the upbound, r2 is the low bound

        r_upbound = r1; r_lowbound = r2; // this shell has r_lowbound < r < r_upbound
        // if(istep>5) printf("r_nextstep_lowbound, r_lowbound = %f %f \n", r_nextstep_lowbound, r_lowbound); 
        r_nextstep_lowbound = Simpson(&inv_Hz, 0., 1./( (istep+2.+0.5)*da )-1., 512);
            max_travel_distance = max_possible_v / 3.e5 * (r_lowbound - r_nextstep_lowbound); // 粒子在两个 snapshot 之间极限 travel 距离；
    }

	if(param.use_lightcone && myrank==0){
	  printf(" (write_lightcone) istep = %d || a in (%5.3f %5.3f), z in (%5.3f %5.3f), r in (%10.3f %10.3f)\n", istep, a1,a2, z2,z1, r2,r1);
  	  printf("                   amid, zmid, rmid = %lf, %lf, %lf\n", amid, zmid, rmid);
	}

	// filename
    char filename[3000]; 
    char filename2[3000]; 
    char filename_test[3000];

    if (param.use_lightcone)
    {
        sprintf(filename, "%s_lightcone.istep%d_thread%d",param.lightcone_basename,istep, myrank);
        sprintf(filename2, "%s_lightcone.istep%d_thread%d.npar",param.lightcone_basename,istep, myrank);
        sprintf(filename_test, "%s_lightcone.istep%d_thread%d_rand_test.txt",param.lightcone_basename,istep, myrank);
    }

	// check whether we need to output lightcone
	if( (!param.use_lightcone) || ( z3>zmax || istep == 1)){}
	else{
	  if(myrank==0){
	    if(is_first_lightcone_shell) printf("                   The first lightcone-shell start to output!!! (zmax=%.3f, rmax=%.3f)\n",zmax,rmax);
	    printf("                   call cola_set_snapshot :    ");
	  }
	cola_set_snapshot((double)amid, particles, snapshot);

	//if(head_nowcheck!=NULL) free_chain_1000(head_nowcheck);
	head_nowcheck = (struct chain_1000 *) malloc(sizeof(struct chain_1000)); p_nowcheck = head_nowcheck;
	p_nowcheck->nowi = -1; p_nowcheck->next = NULL;

	// write particles
	printf("                   (myrank=%d) write lightcone particles to:    %s \n", myrank, filename);
	int ipar, ishift1, ishift2, ishift3, nchecked_now=0, nchecked_pre=0, npar_shifted=0;
	int nowi, nowi_pre, prepar_end_flag, n_check_cross_skip=0, n_check_cross_output=0; 
    long long iwrite=0, id_pre, id_now;
    float x1, x2, x3, x1_shift, x2_shift, x3_shift, r, vfac=1./sqrt(avel1), tmpv[3], x1_shift_pre, x2_shift_pre, x3_shift_pre,  r_pre, r_now; 
    FILE *fp, *fptest;
    float rstat_min=1.0e30, rstat_max=0., rstat_mean=0.;

        fp = fopen(filename,"wb");
        fptest = fopen(filename_test,"w");

    fwrite(&rmid, sizeof(float), 1, fp);

    if(myrank==0){ printf("                              will process    %lld   particles \n", snapshot->np_total); }
    p_precheck=head_precheck; nowi_pre=0;  prepar_end_flag=0;
    if(p_precheck!=NULL){
    id_pre=p_precheck->ids[nowi_pre];  
		//print_chain_1000(p_precheck);
	  	//printf(" (test) Set ipar_pre = %d \n",ipar_pre);
	}


	for(ipar=0; ipar<snapshot->np_local; ipar++){
	   int checked = 0;
	   x1 = snapshot->p[ipar].x[0]; x2 = snapshot->p[ipar].x[1]; x3 = snapshot->p[ipar].x[2]; id_now = snapshot->p[ipar].id;
	   tmpv[0] = snapshot->p[ipar].v[0]*vfac; tmpv[1] = snapshot->p[ipar].v[1]*vfac; tmpv[2] = snapshot->p[ipar].v[2]*vfac;
	   for(ishift1=shift1; ishift1<=shift2; ishift1++){
	   for(ishift2=shift1; ishift2<=shift2; ishift2++){
	   for(ishift3=shift1; ishift3<=shift2; ishift3++){
		x1_shift=x1+param.boxsize*ishift1;  x2_shift=x2+param.boxsize*ishift2;  x3_shift=x3+param.boxsize*ishift3;
		r = sqrt(x1_shift*x1_shift + x2_shift*x2_shift + x3_shift*x3_shift); npar_shifted++;
		int check_cross_skip = 0, check_cross_output = 1;
		if(r < rmax && r_lowbound < r && r < r_upbound + max_travel_distance ){
			//fwrite(&x1, sizeof(double), 1, fp);
			// 上一个输出过，这一个也要输出的，就不要再输出了
			// 上一个没输出，这一个也没有输出的，就输出了
			// 如果是第一次输出 lightcone，那没有上面的比较，就不管了
			// 如果上面的例子已经达到 end 了，那就也不管，不比较了
			if(do_check_cross && !is_first_lightcone_shell ){ // check lightcone cross-boundary if not being the first shell to be outputed
				// 找到这个 ipar 在上一步的 check_cross 数组中的位置
				if(id_pre > id_now){} // 不存在，跳过
				else{
					while(id_pre < id_now && !prepar_end_flag){ // 有可能存在，顺序搜索
						nowi_pre ++;
						if(nowi_pre >= p_precheck->nowi){
							if(p_precheck->next == NULL){ prepar_end_flag = 1; break; }
							else{p_precheck=p_precheck->next; nowi_pre = 0;} }	
						id_pre = p_precheck->ids[nowi_pre]; 
						}
					if(id_pre == id_now){ // 找到了
					//if(ipar_pre == snapshot->p[ipar].id){ // 找到了
						nchecked_pre ++;
						x1_shift_pre=p_precheck->xyzvs[nowi_pre][0]+param.boxsize*ishift1;  
						x2_shift_pre=p_precheck->xyzvs[nowi_pre][1]+param.boxsize*ishift2;  
						x3_shift_pre=p_precheck->xyzvs[nowi_pre][2]+param.boxsize*ishift3;  
						r_pre = sqrt(x1_shift_pre*x1_shift_pre + x2_shift_pre*x2_shift_pre + x3_shift_pre*x3_shift_pre);
						r_now = r;
						if(r_pre >r_upbound && r_now < r_upbound && abs(r_pre - r_now)<max_travel_distance*2){
							check_cross_skip = 0; // 上一个输出过了，这一个就不要输出了 
							n_check_cross_skip ++;
							if(myrank==0){
							  fprintf(fp_check_cross_test, "%lld %.3f %.3f %.3f   %.3f %.3f %.3f    %.3f %.3f %.3f 0\n", 
							    snapshot->p[ipar].id, x1_shift_pre,x2_shift_pre,x3_shift_pre, 
							    x1_shift,x2_shift,x3_shift, r_pre,r_now,r_upbound);}
						}
						else if(r_pre<r_upbound && r_now > r_upbound && abs(r_pre - r_now)<max_travel_distance*2){
							check_cross_output = 1; // 上一个没有输出过，这一个也不会输出 -- 要补回来
							n_check_cross_output ++;
							if(myrank==0){
							  fprintf(fp_check_cross_test, "%lld %.3f %.3f %.3f   %.3f %.3f %.3f    %.3f %.3f %.3f 1\n", 
							    snapshot->p[ipar].id, x1_shift_pre,x2_shift_pre,x3_shift_pre, 
							    x1_shift,x2_shift,x3_shift, r_pre,r_now,r_upbound);}
						}
					}
				}
			}
			// 需要输出的情形： 当前粒子在上边界内，且没有被 skip；当前粒子在上边界以外，但被要求 output 补上
			if( (r < r_upbound && !check_cross_skip) || check_cross_output ){ // 
				fwrite(&(snapshot->p[ipar].id), sizeof(long long), 1, fp);
				fwrite(&x1_shift,sizeof(float),1,fp); 
				fwrite(&x2_shift,sizeof(float),1,fp); 
				fwrite(&x3_shift,sizeof(float),1,fp);
				if(rand()%1000<5) { fprintf(fptest, "%10.3f %10.3f %10.3f \n", x1_shift, x2_shift, x3_shift); }
				fwrite(tmpv, sizeof(float), 3, fp);
				rstat_mean=rstat_mean+r; if(r<rstat_min) rstat_min = r; if(r>rstat_max) rstat_max = r;
			}
			iwrite++;
		}
		if(do_check_cross) // 对当前步骤，跟下边界比较，看哪些粒子需要存起来作 boundary-cross 检查
		 if( ((r>r_lowbound && r - max_travel_distance < r_lowbound) || (r<r_lowbound && r+max_travel_distance>r_lowbound)) && !checked){
			p_nowcheck->nowi ++; // 准备存储到 nowi 位置
			if(p_nowcheck->nowi == 1000){ // 如果已经到尾巴，那么开辟新节点
				p_nowcheck->next = (struct chain_1000 *) malloc (sizeof(struct chain_1000));
				p_nowcheck = p_nowcheck->next; p_nowcheck->next=NULL; p_nowcheck->nowi = 0;
			}
			nowi = p_nowcheck->nowi; //  需要存储 位置 以及身份识别的 ipar 
			p_nowcheck->xyzvs[nowi][0] = x1; p_nowcheck->xyzvs[nowi][1] = x2; p_nowcheck->xyzvs[nowi][2] = x3;
			p_nowcheck->ids[nowi] = id_now; 
			checked = 1; nchecked_now ++;   // checked 变量防止 box-shift 之后重复记录
			//printf("p_check->nowi, nchecked_now = %d %d \n", p_nowcheck->nowi, nchecked_now);
			//printf("particle for check cross: ipar, x,y,z, r, r_lowbound, max_travel_distance = %d   %f %f %f   %f %f  %f\n",ipar,x1,x2,x3,r, r_lowbound, max_travel_distance);
        }
	   } // For 1 454
       } // For 2 455
       } // For 3 456
	}
	  rstat_mean /= (float)(iwrite);
	  n_check_cross_skip_total += n_check_cross_skip;
	  n_check_cross_output_total += n_check_cross_output;
	  n_write_total += iwrite;
	  printf("                   (myrank=%d)    %lld   (among %d) particles written with %.3f <= r <= %.3f (avg=%.3f)\n",myrank,iwrite,npar_shifted,rstat_min,rstat_max,rstat_mean);
	  printf("                   (myrank=%d)    #-check_cross_from_pre =   %d (%.3f%% w.r.t. npar_shifted) skip=%d, %.3f%%; output=%d, %.3f%% (w.r.t. iwrite)\n", 
	    myrank, nchecked_pre, nchecked_pre/(double)npar_shifted*100., 
	    n_check_cross_skip, n_check_cross_skip/(double)iwrite*100., 
	    n_check_cross_output, n_check_cross_output/(double)iwrite*100. );
	  printf("                   (myrank=%d) check cross for-next: %d (%.3f%%/%.3f%%  w.r.t. np_local/npar_shifted) (vmax,max_travel=%.1f %.1f)\n", 
	    myrank,nchecked_now,nchecked_now/(double)snapshot->np_local*100.,nchecked_now/(double)npar_shifted*100.,max_possible_v,max_travel_distance);
	  //printf(" (test) now print head_nowcheck: \n"); print_chain_1000(head_nowcheck);
	  if(!head_precheck) {free_chain_1000(head_precheck);} 
	  head_precheck = head_nowcheck; head_nowcheck=NULL;
	  //printf(" (test) now print head_precheck: \n"); print_chain_1000(head_precheck);
	  long long n_neg=-1; fwrite(&n_neg,sizeof(long long),1,fp); fclose(fp); fclose(fptest);
	  fp = fopen(filename2, "wb"); fwrite(&iwrite, sizeof(long long), 1, fp); fclose(fp);
	  is_first_lightcone_shell = 0;
	  if(myrank==0) {
		 //getchar();
	  }
	  MPI_Barrier(MPI_COMM_WORLD); 

	  //if(n_check_cross_skip>1) exit(0);
	} // else 416
    } // istep 323
    } // nout 320

    if (param.use_lightcone)
    {
    printf("                   (myrank=%d)    # written-to-file =   %lld      skip=%lld, %.3f%%; output=%lld, %.3f%%  \n",
	    myrank, n_write_total, n_check_cross_skip_total, n_check_cross_skip_total / (double)n_write_total * 100, 
	    n_check_cross_output_total, n_check_cross_output_total / (double)n_write_total * 100 );
    long long n_write_allnode, n_check_cross_skip_allnode, n_check_cross_output_allnode;
    MPI_Reduce(&n_write_total, &n_write_allnode, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&n_check_cross_skip_total, &n_check_cross_skip_allnode, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&n_check_cross_output_total, &n_check_cross_output_allnode, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if(myrank==0){
      fclose(fp_check_cross_test);
      printf("----------------------------------\n");
      printf("   (SUM-UP ALL NODES)   # written-to-file =   %lld      skip=%lld, %.3f%%; output=%lld, %.3f%%  \n",
	    n_write_allnode, n_check_cross_skip_allnode, n_check_cross_skip_allnode/ (double)n_write_allnode* 100, 
	    n_check_cross_output_allnode, n_check_cross_output_allnode/ (double)n_write_allnode* 100 );
      printf(" (write_lightcone) Done. da = %.3f   ", da );
      printf(" | zmax, rmax = %.3f %.3f  ", zmax, rmax);
      printf(" | nshift = rmax/boxize = %.3f / %.3f = %d \n ", rmax, param.boxsize, nshift);
      printf("----------------------------------\n");

    }
  }
    timer_print();
    // xiaodong's code for lightcone output end
    // ----------------------------------------------------
  }

  //move_particles(particles);
  //write_slice("slice", particles->p, particles->np_local, param.boxsize/32, param.boxsize);
  if(myrank==0){
    print_heart();
//  int i;
//  for(i=0;i<100;i++)
//	printf("  &&&&& 写这个程序我真心醉了  ######...");
//    putchar('\n');
  }

  MPI_Finalize();
  return 0;
}

int mpi_init(int* p_argc, char*** p_argv)
{
  // MPI+OpenMP paralellization: MPI_THREAD_FUNNELED
  // supported by mpich2 1.4.1, but now by openmpi 1.2.8

#ifdef _OPENMP
  int thread_level, hybrid_parallel;
  MPI_Init_thread(p_argc, p_argv, MPI_THREAD_FUNNELED, &thread_level);
  hybrid_parallel = (thread_level >= MPI_THREAD_FUNNELED);

  int myrank; MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if(myrank == 0) {
    if(hybrid_parallel)
      printf("MPI + multi thread supported (MPI_THREAD_FUNNELED).\n");
    else
      printf("Warning: MPI + multi thread not supported. 1 thread per node.\n");
  }
	
  return hybrid_parallel;
#else
  MPI_Init(p_argc, p_argv);
  int myrank; MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if(myrank == 0)
    printf("MPI only without OpenMP\n");
  return 0;
#endif

}

void fft_init(int threads_ok)
{
  // Initialize FFTW3

#ifdef _OPENMP
  if(threads_ok)
    threads_ok= fftwf_init_threads();
  if(!threads_ok)
    msg_printf(warn, "Multi-thread FFTW not supported.\n");
#endif
    
  fftwf_mpi_init();

#ifdef _OPENMP
  if(threads_ok) {
    int nthreads= omp_get_max_threads();
    fftwf_plan_with_nthreads(nthreads);
    msg_printf(info, "Multi-threaded FFTW: %d threads\n", nthreads);
  }
#endif

}


void write_slice(const char filebase[], Particle* p, const int np, const float dz, const float boxsize)
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  char filename[3000];
  sprintf(filename, "%s%d.txt", filebase, myrank);
  FILE* fp= fopen(filename, "w");
  if(fp == 0)
    msg_abort(0020, "Unable to write to %s\n", filename);

  
  for(int i=0; i<np; i++) {
    float x= fmod(p->x[0], boxsize);
    float y= fmod(p->x[1], boxsize);
    float z= fmod(p->x[2], boxsize);
    if(0.0f < z && z < dz) {
      fprintf(fp, "%e %e %e\n",  x, y, z);
    }
    p++;
  }
  fclose(fp);
}

void write_snapshot_slice(const char filebase[], Snapshot const * const snapshot, const float boxsize)
{
  msg_printf(normal, "Writing snapshot slices %s\n", filebase);
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  char filename[3000];
  sprintf(filename, "%s_%d.txt", filebase, myrank);
  FILE* fp= fopen(filename, "w");
  if(fp == 0)
    msg_abort(0020, "Unable to write to %s\n", filename);

  ParticleMinimum const * p= snapshot->p;
  const int np= snapshot->np_local;
  
  for(int i=0; i<np; i++) {
    float x= fmod(p->x[0], boxsize);
    float y= fmod(p->x[1], boxsize);
    float z= fmod(p->x[2], boxsize);
    if(0.0f < z && z < boxsize/32) {
      fprintf(fp, "%e %e %e\n",  x, y, z);
    }
    p++;
  }
  fclose(fp);
}
    
void snapshot_time(const float aout, const int iout, 
		   Particles const * const particles, 
		   Snapshot * const snapshot,
		   const char fof_filename[], 
		   const char subsample_filename[], 
		   const char cgrid_filename[], const int cgrid_nc,
		   void* const mem1, const size_t size1,
		   const int write_longid,
		   const double fof_linking_factor)
{
  // Halo finding and snapshot outputs

  char filebase[6000];
  const int isnp= iout+1;
  char suffix= 'a' + iout;

  // If you prefer other output filename format (not a,b,c, for example
  // Change the sprintf statements below
  // sprintf(filebase, "%s%05d_%03d", snapshot->filename, snapshot->seed, iout);
  // would give snp00001_000
  
                                                       timer_set_category(Snp);
  cola_set_snapshot(aout, particles, snapshot);

  const int nc= snapshot->nc; assert(nc > 0);
  const float boxsize= snapshot->boxsize; assert(boxsize > 0.0f);

  // FoF halo catalogue (text file)
  if(fof_filename) {
    const float ll= (float)(fof_linking_factor*boxsize/nc);// FoF linking length
    fof_find_halos(snapshot, ll);  // move_particles done here
    sprintf(filebase, "%s%05d.%03d.txt", fof_filename, snapshot->seed, iout);
    fof_write_halos(filebase);
  }


                                                       timer_start(write);
  // Gadget snapshot for all particles
  if(snapshot->filename) {
    sprintf(filebase, "%s%05d.%03d", snapshot->filename, snapshot->seed, iout);
    write_snapshot(filebase, snapshot, write_longid);
  }
                                                       timer_stop(write);
						       timer_start(sub);
  // particle subsample (binary)
  if(subsample_filename) {
    sprintf(filebase, "%s%05d.%03d.b", subsample_filename, snapshot->seed, iout);

    // This is an alternative option of subsamplig (regular subsampling)
    // write_subsample(filebase, subsample_factor, snapshot, mem1, size1);

    write_random_sabsample(filebase, snapshot, mem1, size1);
  }

  // coarse mesh (binary)
  if(cgrid_filename) {
    sprintf(filebase, "%s%05d.%03d.b", cgrid_filename, snapshot->seed, iout);
    coarse_grid2(filebase, snapshot, cgrid_nc, mem1, size1);
  }

                                                       timer_stop(sub);
  // text file of a slice
  // sprintf(filebase, "slice%02d", isnp); temp
  // write_snapshot_slice(filebase, snapshot, boxsize);

  const double z_out= 1.0/aout - 1.0;
  msg_printf(normal, "snapshot %d (%03d) written z=%4.2f a=%5.3f\n", 
	     isnp, iout, z_out, aout);
                                                       timer_set_category(COLA);
}
