#include "xypic.h"
extern int DumpFlag;
extern int nfsp;
extern int ncx,ncy,Csize;
extern int ngx,ngy,Gsize;
extern float dx,dy;
extern GCA *vec_C;
extern GGA *vec_G;
extern Fluid *FG;	// fluid species
#ifndef __FLUID_H__
#define __FLUID_H__
extern int Conti_Flag;
extern GFC *Fluid_sp;
extern GFG *Fluid_Den, *Fluid_Src;
extern int Conti_xnum, Conti_ynum;
extern Con_RegionX *Conti_x;
extern Con_RegionY *Conti_y;
void Fluid_Setting();
int Cal_XRegion_check();
int Cal_YRegion_check();
void Set_Con_Region(int isp, Con_RegionX *Cx,Con_RegionY *Cy);
void Set_Con_Boundary(int isp, Con_RegionX *Cx,Con_RegionY *Cy);
#endif