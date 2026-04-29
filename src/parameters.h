#ifndef PARAMETERS_H
#define PARAMETERS_H

typedef struct {
  int nc;
  int pm_nc_factor;

  double np_alloc_factor;

  int ntimestep;
  int random_seed;
  int nrealization;

  int loglevel;
  
  double a_init;
  double a_final;
  double boxsize;

  double omega_m, sigma8, h;
  double omega_k, omega_l, de_w; // xiaodong add: omegal,  dark energy EoS

  int only_output_1eighth; // xiaodong add; if set as 1, only output 1/8 of the sky

  double* zout; int n_zout;

  /*
  char power_spectrum_file[64];
  char fof_filename[64];
  char snapshot_filename[64];
  char subsample_filename[64];
  char cgrid_filename[64];
  char init_filename[64];
  */
  char* power_spectrum_filename; int strlen_power_spectrum_filename;
  char* fof_filename;        int strlen_fof_filename;
  char* snapshot_filename;   int strlen_snapshot_filename;
  char* subsample_filename;  int strlen_subsample_filename;
  char* cgrid_filename;      int strlen_cgrid_filename;
  char* init_filename;       int strlen_init_filename;
  // 'strlen' includes the null character. It is the number of chars.


  //int subsample_factor;
  double fof_linking_factor;
  double subsample_factor;
  int cgrid_nc;

  int write_longid; 

  // xiaodong: variables for lightcone
  int use_lightcone;  // boolean flag: 1 to enable lightcone, 0 to disable (default)
  double lightcone_zmax;
  char* lightcone_basename;  int strlen_lightcone_basename;
  int use_solve_growth;

} Parameters;

//int read_parameters(const char filename[], Parameters* const param);
int read_parameters(const int argc, char * argv[], 
		    Parameters* const param);

void confirm_parameters(Parameters* const param, const char filename[]);

#endif
