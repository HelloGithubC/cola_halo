

extern double gb_om,gb_w, gb_H0, gb_oml;
extern int gb_use_solve_growth;
//extern int solvegrowth_flag;

extern void gb_solve_D_E_Qs();
extern void get_D_E(double nowa, double *Da, double *Ea, double *QdDda, double * QdEda);
extern void intpl_D_E_Qs(double, double*, double*, double*, double*);
//extern void gb_Tsq_
extern void gb_Tsquare_DE(double, double, double, double, double, float *, float *);
