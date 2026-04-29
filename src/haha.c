main.c:  MPI_Comm_rank(MPI_COMM_WORLD, &myrank); // xiaodong
main.c:  //power_init("camb0_matterpower.dat", a_init, sigma8, OmegaM, OmegaLambda);
main.c:	     sigma8, OmegaM, OmegaLambda);
main.c:  snapshot->omega_m= OmegaM;
main.c:  subsample_init(param.subsample_factor, param.random_seed);
main.c:    MPI_Barrier(MPI_COMM_WORLD);
main.c:    int seed= param.random_seed + irealization;
main.c:    lpt_set_displacement(a_init, OmegaM, seed,
main.c:    // set omegam
main.c:    gb_omegam = OmegaM;
main.c:    //printf("Integrate of x from 0 to 1: %f\n", Simpson(&fun_x2, 0, 3, 100 ));
main.c:    //printf("gb_omegam = %lf \n", gb_omegam);
main.c:    //printf("Comov_r at z = 0., 1.0: %lf %lf\n", Simpson(&inv_Hz, 0, 0.5, 100 ),  Simpson(&inv_Hz, 0, 1., 100 ));
main.c:                                                            timer_start(comm);
main.c:                                                            timer_stop(comm);
main.c:	cola_kick(particles, OmegaM, avel1);
main.c:	cola_drift(particles, OmegaM, apos1);
main.c:	  printf("                   (myrank=%d)    #-check_cross_from_pre =   %d (%.3f%% w.r.t. npar_shifted) skip=%d, %.3f%%; output=%d, %.3f%% (w.r.t. iwrite)\n", 
main.c:	  MPI_Barrier(MPI_COMM_WORLD); 
main.c:    MPI_Reduce(&n_write_total, &n_write_allnode, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
main.c:    MPI_Reduce(&n_check_cross_skip_total, &n_check_cross_skip_allnode, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
main.c:    MPI_Reduce(&n_check_cross_output_total, &n_check_cross_output_allnode, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
main.c:  int myrank; MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
main.c:  int myrank; MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
main.c:    int nthreads= omp_get_max_threads();
main.c:  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
main.c:  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
main.c:    write_random_sabsample(filebase, snapshot, mem1, size1);
mem.c:#include "comm.h"
mem.c:  const int nnode= comm_nnode();
mem.c:    fftwf_mpi_local_size_3d(nc, nc, nc/2+1, MPI_COMM_WORLD,
mem.c:  ptrdiff_t ncomplex_lpt= 12*size_lpt_one;
mem.c:	     (int)(ncomplex_lpt*sizeof(fftwf_complex)/(1024*1024)));
mem.c:    fftwf_mpi_local_size_3d_transposed(Ngrid, Ngrid, Ngrid/2+1, MPI_COMM_WORLD,
mem.c:  ptrdiff_t ncomplex_pm= size_pm_one;
mem.c:	     (int)(ncomplex_pm*sizeof(fftwf_complex)/(1024*1024)));
mem.c:  ptrdiff_t ncomplex1= ncomplex_lpt > ncomplex_pm ? ncomplex_lpt : ncomplex_pm;
mem.c:  size_t size1= sizeof(fftwf_complex)*ncomplex1;
mem.c:    ncomplex1= size_fof/sizeof(fftwf_complex) + 1;
mem.c:  mem->mem1= fftwf_alloc_complex(ncomplex1);
mem.c:  mem->size1= sizeof(fftwf_complex)*ncomplex1;
mem.c:  size_t ncomplex2= (Ngrid/2+1)*Ngrid*local_ny; //ncomplex_pm;
mem.c:  size_t size2= sizeof(fftwf_complex)*(Ngrid/2+1)*Ngrid*local_ny;
mem.c:    ncomplex2= size_snapshot/sizeof(fftwf_complex) + 1;
mem.c:  mem->mem2= fftwf_alloc_complex(ncomplex2);
mem.c:  mem->size2= sizeof(fftwf_complex)*ncomplex2;
move.c:#include "comm.h"
move.c:  const int ThisNode= comm_this_node();
move.c:  const int NNode= comm_nnode();
move.c:  const float x_left= comm_xleft(0);
move.c:  const float x_right= comm_xright(0);
move.c:    int nsend_max= comm_share_int(nsend, MPI_MAX);
move.c:  int np_max= comm_reduce_int(particles->np_local, MPI_MAX);
move.c:  MPI_Reduce(&np, &np_global, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
move.c:  const float x_left= comm_xleft(dix);
move.c:  const float x_right= comm_xright(dix);
move.c:  const int node_to= comm_node(dix);
move.c:  const int node_from= comm_node(-dix);
move.c:  MPI_Barrier(MPI_COMM_WORLD); // debug!!!
move.c:	       &nrecv, 1, MPI_INT, node_from, tag, MPI_COMM_WORLD, &status); 
move.c:               p+np_local, nrecv*sizeof(Particle), MPI_BYTE, node_from, tag,
move.c:	       MPI_COMM_WORLD, &status); tag++;
move.c:  msg_printf(verbose, "Received %d particles from node %d\n", nrecv, node_from);
move_min.c:// This code is automatically generated from move.c
move_min.c:#include "comm.h"
move_min.c:  const int ThisNode= comm_this_node();
move_min.c:  const int NNode= comm_nnode();
move_min.c:  const float x_left= comm_xleft(0);
move_min.c:  const float x_right= comm_xright(0);
move_min.c:    int nsend_max= comm_share_int(nsend, MPI_MAX);
move_min.c:  int np_max= comm_reduce_int(particles->np_local, MPI_MAX);
move_min.c:  MPI_Reduce(&np, &np_global, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
move_min.c:  const float x_left= comm_xleft(dix);
move_min.c:  const float x_right= comm_xright(dix);
move_min.c:  const int node_to= comm_node(dix);
move_min.c:  const int node_from= comm_node(-dix);
move_min.c:  MPI_Barrier(MPI_COMM_WORLD); // debug!!!
move_min.c:	       &nrecv, 1, MPI_INT, node_from, tag, MPI_COMM_WORLD, &status); 
move_min.c:               p+np_local, nrecv*sizeof(ParticleMinimum), MPI_BYTE, node_from, tag,
move_min.c:	       MPI_COMM_WORLD, &status); tag++;
move_min.c:  msg_printf(verbose, "Received %d particles from node %d\n", nrecv, node_from);
msg.c:// Utility functions to write message, record computation time ...
msg.c:  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
msg.c:  MPI_Abort(MPI_COMM_WORLD, errret);
pm.c:// This file contains some standard functions for a PM code. 
pm.c:#include "comm.h"
pm.c:static fftwf_complex* fftdata;
pm.c:static fftwf_complex* density_k;
pm.c:  #pragma omp atomic
pm.c:    fftwf_mpi_local_size_3d_transposed(Ngrid, Ngrid, Ngrid/2+1, MPI_COMM_WORLD,
pm.c:  size_t bytes= sizeof(fftwf_complex)*total_size;
pm.c:    fftdata=  fftwf_alloc_complex(total_size);
pm.c:    assert(size1 >= total_size*sizeof(fftwf_complex));
pm.c:    fftdata= (fftwf_complex*) mem1;
pm.c:    density_k= fftwf_alloc_complex((NgridL/2+1)*NgridL*Local_ny_td);
pm.c:    assert(size2 >= sizeof(fftwf_complex)*(NgridL/2+1)*NgridL*Local_ny_td);
pm.c:    density_k= (fftwf_complex*) mem2;
pm.c:		    MPI_COMM_WORLD, FFTW_MEASURE | FFTW_MPI_TRANSPOSED_OUT);
pm.c:		     MPI_COMM_WORLD, FFTW_MEASURE | FFTW_MPI_TRANSPOSED_IN);
pm.c:  MPI_Comm_rank(MPI_COMM_WORLD, &ThisNode);
pm.c:  MPI_Comm_size(MPI_COMM_WORLD, &NNode);
pm.c:		MPI_COMM_WORLD);
pm.c:		MPI_COMM_WORLD);
pm.c:void PtoMesh(const Particle P[], const int np, float* const density)
pm.c:  msg_printf(verbose, "Calculating PtoMesh\n");
pm.c:  #pragma omp parallel for default(shared)
pm.c:  #pragma omp parallel for default(shared)
pm.c:    // Buffer particles are copied from adjacent nodes, instead
pm.c:void compute_density_k(void)
pm.c:  #pragma omp parallel for default(shared)
pm.c:// Calculate one component of force mesh from precalculated density(k)
pm.c:void compute_force_mesh(const int axes)
pm.c:  fftwf_complex* const FN11= fftdata;
pm.c:  fftwf_complex* const P3D= density_k;
pm.c:#pragma omp parallel for default(shared)
pm.c:  #pragma omp parallel for default(shared)     
pm.c:	       MPI_COMM_WORLD, &status); 
pm.c:		 MPI_COMM_WORLD, &status);
pm.c:  MPI_Reduce(&nsend, &nsend_global, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
pm.c:			MPI_COMM_WORLD, &status);
pm.c:  MPI_Reduce(&sum, &sum_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
pm.c:  PtoMesh(particles->p, np_plus_buffer, (float*) fftdata);
pm.c:  compute_density_k();
pm.c:    compute_force_mesh(axes);
pm.c:                                                            timer_start(comm);
pm.c:                                                            timer_stop(comm);
power.c:static double Omega, OmegaLambda;
power.c:void power_init(const char filename[], const double a_init, const double sigma8, const double omega_m, const double omega_lambda)
power.c:  Omega= omega_m;
power.c:  OmegaLambda= omega_lambda;
power.c:  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
power.c:  MPI_Bcast(&NPowerTable, 1, MPI_INT, 0, MPI_COMM_WORLD);
power.c:	    MPI_COMM_WORLD);
power.c:  MPI_Bcast(&Norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
power.c:    msg_abort(3010, "Input sigma8 %f is far from target sigma8 %f\n",
power.c:  return pow(a / (Omega + (1 - Omega - OmegaLambda) * a + OmegaLambda * a * a * a), 1.5);
power.c:  hubble_a = sqrt(Omega / (a * a * a) + (1 - Omega - OmegaLambda) / (a * a) + OmegaLambda);
read.c:#include "comm.h"
read.c:  const int this_node= comm_this_node();
read.c:  MPI_Bcast(&numfiles, 1, MPI_INT, 0, MPI_COMM_WORLD);
read.c:      msg_printf(verbose, "%d particles read from %s.\n", 
read.c:    MPI_Bcast(&np_snapshot, 1, MPI_INT, 0, MPI_COMM_WORLD);
read.c:    MPI_Bcast(&vfac_back, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
read.c:    MPI_Bcast(x, np_snapshot * 6, MPI_FLOAT, 0, MPI_COMM_WORLD);
read.c:    const float x_left= comm_xleft(0);
read.c:    const float x_right= comm_xright(0);
read.c:  MPI_Bcast(&header, sizeof(GadgetHeader), MPI_BYTE, 0, MPI_COMM_WORLD);
read.c:  snapshot->np_average= (float)(((double) snapshot->np_total)/comm_nnode());
read.c:  snapshot->omega_m= (float) header.omega0;
read.c:  MPI_Reduce(&np_local, &np_total, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
read_param.c:// -- comment
read_param.c:static void remove_comments(char* str)
read_param.c:  // any string after -- is a coment
read_param.c:    msg_abort(1003, "Error: string value not starting from \": %s\n",
read_param.c:  MPI_Comm_rank(MPI_COMM_WORLD, &myrank_);
read_param.c:  MPI_Bcast(param, sizeof(Parameters), MPI_BYTE, 0, MPI_COMM_WORLD);
read_param.c:  param->omega_m = -1.0;
read_param.c:  param->random_seed=  1;
read_param.c:  param->omega_l=-1.;
read_param.c:    remove_comments(buf);
read_param.c:    else if(strcmp(name, "random_seed") == 0)
read_param.c:      param->random_seed= get_int(p);
read_param.c:    else if(strcmp(name, "omega_m") == 0)
read_param.c:      param->omega_m= get_double(p);
read_param.c:    else if(strcmp(name, "omega_l") == 0)
read_param.c:      param->omega_l = get_double(p);
read_param.c:  if(param->omega_l < 0.){
read_param.c:     param->omega_l = 1. - param->omega_m;	  
read_param.c:     printf("Set omega_l as 1-omega_m... %g\n", param->omega_l);
read_param.c:  const int ret1= MPI_Bcast(len, 1, MPI_INT, 0, MPI_COMM_WORLD);
read_param.c:  const int ret2= MPI_Bcast(*pstring, n, MPI_CHAR, 0, MPI_COMM_WORLD);
read_param.c:  const int ret1= MPI_Bcast(len, 1, MPI_INT, 0, MPI_COMM_WORLD);
read_param.c:  const int ret2= MPI_Bcast(*parray, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
read_param_lua.c:  MPI_Comm_rank(MPI_COMM_WORLD, &myrank_);
read_param_lua.c:  MPI_Bcast(param, sizeof(Parameters), MPI_BYTE, 0, MPI_COMM_WORLD);
read_param_lua.c:  param->random_seed= read_int(L, "random_seed");
read_param_lua.c:  param->omega_m= read_double(L, "omega_m");
read_param_lua.c:  const int ret1= MPI_Bcast(len, 1, MPI_INT, 0, MPI_COMM_WORLD);
read_param_lua.c:  const int ret2= MPI_Bcast(*pstring, n, MPI_CHAR, 0, MPI_COMM_WORLD);
read_param_lua.c:  const int ret1= MPI_Bcast(len, 1, MPI_INT, 0, MPI_COMM_WORLD);
read_param_lua.c:  const int ret2= MPI_Bcast(*parray, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
solve_growth.c://double gb_om = 0.3071, gb_oml=0.6929, gb_w=-1.5, gb_H0=70;
solve_growth.c://double gb_om = gb_om, gb_w = gb_w, gb_oml = gb_oml, gb_H0=70.;
solve_growth.c:double gb_om,gb_w, gb_oml, gb_H0;
solve_growth.c:// compute T^2 D1 = Q (d/da) ( Q (d/da) D1 ) = Q dQ/da dD1/da + Q^2  (d/da)^2 D1
solve_growth.c:	return gb_H0 * pow( gb_om*pow(a,-3.0) + gb_oml*pow(a,-3.0*(1.0+gb_w)) + (1.-gb_om-gb_oml)/(a*a), 0.5);
solve_growth.c:	//return gb_H0 * pow( gb_om*pow(a,-3.0) + gb_oml*pow(a,-3.0*(1.0+gb_w) + (1.-gb_om-gb_oml)/(a*a)), 0.5);
solve_growth.c://	y1 = 3. * pow(a,-4.) / y * gb_H0 * (1.+gb_om*(-1.+y-gb_w)+gb_w);
solve_growth.c://	y2 = pow(a,-3.)/y*(1.+(-1.+y)*gb_om);
solve_growth.c://	y1 = -3.*gb_H0*gb_om + 2.*a*gb_H0*(-1.+gb_oml+gb_om);
solve_growth.c://	y2 = sqrt(gb_oml+(gb_om-a*(-1.+gb_oml+gb_om))/y);
solve_growth.c:	y1 = -3.*y * gb_om + 2.*a * y * (gb_oml+gb_om-1.) - 3.*gb_oml*(1.+gb_w);
solve_growth.c:	y2 = gb_oml / y + gb_om - a*(gb_oml+gb_om-1.);
solve_growth.c://	y1 = (3.*gb_H0*gb_H0*gb_om)/(2.*a*a*a) * D - 2*H*a*H*Dprime;
solve_growth.c:	y1 = (3.*gb_H0*gb_H0*gb_om)/(2.*a*a*a) * D - 2*H*a*H*Dprime;
solve_growth.c:	y1 = (3.*gb_H0*gb_H0*gb_om)/(2.*a*a*a) * (E-D*D);
solve_growth.c:	y1 = (3.*gb_H0*gb_H0*gb_om)/(2.*a*a*a) * D - 2*H*a*H*Dprime;
solve_growth.c:	y1 = (3.*gb_H0*gb_H0*gb_om)/(2.*a*a*a) * (E-D*D);
solve_growth.c:	//gb_om = gb_om; gb_w = gb_w; gb_oml = gb_oml;
solve_growth.c:	//printf(" (get_D_E) gb_om, gb_oml, gb_w = %f %f %f \n",gb_om,gb_oml,gb_w);
solve_growth.c:	//gb_om = gb_om; gb_w = gb_w; gb_oml = gb_oml;
solve_growth.c:	printf(" (gb_solve_D_E_Qs) gb_om, gb_oml, gb_w = %f %f %f \n",gb_om,gb_oml,gb_w);
subsample.c:#include "comm.h"
subsample.c:static gsl_rng* Random_generator= 0;
subsample.c:  Random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
subsample.c:  const int this_node= comm_this_node();
subsample.c:  gsl_rng_set(Random_generator, 2*seed + 100*this_node);
subsample.c:   gsl_rng_free(Random_generator);
subsample.c:// newer version below (write_random_subsample()) uses MPI_File_write
subsample.c:void write_random_sabsample1(const char filename[], Snapshot const * const snapshot, void* const mem, size_t mem_size)
subsample.c:    if(gsl_rng_uniform(Random_generator) < SubsampleFactor)
subsample.c:  const int this_node= comm_this_node();
subsample.c:  const int nnode= comm_nnode();
subsample.c:    MPI_Gather(&nsub, 1, MPI_INT, nsub_recv, 1, MPI_INT, 0, MPI_COMM_WORLD);
subsample.c:		   0, MPI_COMM_WORLD);
subsample.c:    subsample.omega_m= snapshot->omega_m;
subsample.c:void write_random_sabsample(const char filename[], Snapshot const * const snapshot, void* const mem, size_t mem_size)
subsample.c:    if(gsl_rng_uniform(Random_generator) < SubsampleFactor) {
subsample.c:  header[1]= snapshot->omega_m*rho_crit*pow(snapshot->boxsize, 3.0)/np;
subsample.c:  header[2]= snapshot->omega_m;
subsample.c:  const int this_node= comm_this_node();
subsample.c:  const int nnode= comm_nnode();
subsample.c:    MPI_Gather(&nsub, 1, MPI_INT, nsub_recv, 1, MPI_INT, 0, MPI_COMM_WORLD);
subsample.c:		   0, MPI_COMM_WORLD);
subsample.c:    subsample.omega_m= snapshot->omega_m;
timer.c:static const char * SubName[]= {"", "fft", "assign", "force_mesh", "pforce", "check", "comm", "evolve", "write", "kd_build", "kd_link", "interpolate", "global", "smalldata"};
write.c:#include "comm.h"
write.c:  sprintf(filename, "%s.%d", filebase, comm_this_node());
write.c:  const double omega_m= snapshot->omega_m;
write.c:  MPI_Reduce(&np_send, &np_total, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
write.c:  MPI_Bcast(&np_total, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
write.c:  const double m= omega_m*rho_crit*pow(boxsize, 3.0)/np_total;
write.c:  header.num_files= comm_nnode();
write.c:  header.omega0= omega_m;
write.c:  header.omega_lambda= 1.0 - omega_m;
write.c:  const double omega_m= snapshot->omega_m;
write.c:  const double m= omega_m*rho_crit*pow(boxsize, 3.0)/np;
write.c:  header.num_files= comm_nnode();
write.c:  header.omega0= omega_m;
write.c:  header.omega_lambda= 1.0 - omega_m;
write.c:  int inode= comm_this_node();
write.c:  const float omega_m= snapshot->omega_m;
write.c:  const float m= omega_m*rho_crit*pow(boxsize, 3.0)/np;  
write.c:  fwrite(&snapshot->omega_m, sizeof(float), 1, fp);
write.c:  // Gather number of particles to compute the offset for writing
write.c:  const int this_node= comm_this_node();
write.c:  const int nnode= comm_nnode();
write.c:    MPI_Gather(&np, 1, MPI_INT, np_local, 1, MPI_INT, 0, MPI_COMM_WORLD);
write.c:		   &np_partial_sum, 1, MPI_INT, 0, MPI_COMM_WORLD);
write.c:  ret= MPI_File_open(MPI_COMM_WORLD, filename,
write.c:    // header: boxsize, m_particle, omega_m, h, a, redshift & np_total
