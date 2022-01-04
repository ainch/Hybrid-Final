#include "xypic.cuh"
extern float dtc;
extern int nfsp;
extern int Gsize, Csize;
extern int ngx,ngy,ncx,ncy;
extern float dx,dy;
extern int MainGas;
extern float Total_Pressure;
extern GCA *vec_C;
extern GGA *vec_G;
extern Fluid *FG;	// fluid species
extern GFC *Fluid_sp;
extern GFG *Fluid_Den, *Fluid_Src;
extern int Conti_xnum, Conti_ynum;
extern Con_RegionX *Conti_x;
extern Con_RegionY *Conti_y;
extern float gputime;
extern cudaEvent_t start, stop;
#ifndef __CUDA_FLUID_CUH__
#define __CUDA_FLUID_CUH__
extern Fluid *dev_FG;
extern GFG *dev_FG_Den, *dev_FG_Src;
void Set_Fluid_cuda();
void Sync_Fluid_GFCtoGFG_forDen(GFC *A, GFG *B);
void Sync_Fluid_GFGtoGFC_forSource(GFG *A, GFC *B);
void Cal_D_Argon(GCA *vecC, GGA *vecG, GFC *data);
float Ar_meta_omega(float value);
void Cal_D_Oxygen(GCA *vecC, GGA *vecG, GFC *data);
void Cal_D_ArO2(GCA *vecC, GGA *vecG, GFC *data);
void calculate_gummel_coef_x();
void calculate_gummel_coef_y();
void Solve_Continuity_eqn();
void Solve_Continuity_eqn_check();
void Solve_Density_x(int isp);
void Solve_Density_y(int isp);
void Calculate_Flux_x(int isp);
void Calculate_Flux_y(int isp);
int tridiag(float *a, float *b, float *c, float *d, float *x, float *gam, int min, int max);
__global__ void calculate_gummel_coef_x(int nfsp, int size, int ngy, int Gsize,int Csize, float dx, Con_RegionX *val, Fluid *info, GFC *Cdata, GFG *Gdata);
__global__ void calculate_gummel_coef_y(int nfsp, int size, int ngy, int Gsize,int Csize, float dy, Con_RegionY *val, Fluid *info, GFC *Cdata, GFG *Gdata);
#endif