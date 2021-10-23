#include "xypic.cuh"
extern int nfsp;
extern int Gsize, Csize;
extern int ngx,ngy,ncx,ncy;
extern float dx,dy;
extern int MainGas;
extern Fluid *FG;	// fluid species
extern GFC *Host_C_F;
extern GGA *dev_GvecSet;
extern GCA *dev_CvecSet;
extern GGA *vec_G;
extern GCA *vec_C;
extern float Total_Pressure;
extern dim3 CONTI_GRID,CONTI_BLOCK;
#ifndef __CUDA_FLUID_CUH__
#define __CUDA_FLUID_CUH__
int Conti_Flag;
Fluid *dev_FG;
GFC *dev_C_F;
GFG *Host_G_F;
GFG *dev_G_F;
int Conti_xnum, Conti_ynum;
dim3 CONTIx_GRID,CONTIx_BLOCK;
dim3 CONTIy_GRID,CONTIy_BLOCK;
Con_RegionX *Conti_x, *dev_Conti_x;
Con_RegionY *Conti_y, *dev_Conti_y;

void Set_Fluid_cuda();
int Cal_XRegion_check();
int Cal_YRegion_check();
void Set_Con_Region(int isp, Con_RegionX *Cx,Con_RegionY *Cy);
void Set_Con_Boundary(int isp, Con_RegionX *Cx,Con_RegionY *Cy);
__global__ void Cal_D_Argon(int nfsp, int ncy, int Csize, float press, GCA *vec_C, GGA *vec_G, GFC *data);
__device__ float Ar_meta_omega(float value);
__global__ void Cal_D_Oxygen(int nfsp, int ncy, int Csize, float press, GCA *vec_C, GGA *vec_G, GFC *data);
__global__ void Cal_D_ArO2(int nfsp, int ncy, int Csize, float press, GCA *vec_C, GGA *vec_G, GFC *data);
__global__ void calculate_gummel_coef_x(int nfsp, int size, int ngy, int Gsize,int Csize, float dx, Con_RegionX *val, Fluid *info, GFC *Cdata, GFG *Gdata);
__global__ void calculate_gummel_coef_y(int nfsp, int size, int ngy, int Gsize,int Csize, float dy, Con_RegionY *val, Fluid *info, GFC *Cdata, GFG *Gdata);
#endif