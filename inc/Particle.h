#include "xypic.h"
extern long int seed; // related to Random number
extern float dx,dy;
extern int ngx,ngy,Gsize;
extern int ncx,ncy,Csize;
extern float fncx,fncy,fngx,fngy;
extern float xlength,ylength,zlength;
extern Species *SP;// particle species
extern int NP_LIMIT; //Each of particle limit
extern HCP *PtD;
extern GCA *vec_C;
extern float **MFMalloc(int sizeX,int sizeY);
extern void MFInit(float **M,float C,int sizeX,int sizeY);

#ifndef __PARTICLE_H__
#define __PARTICLE_H__
int h_nvel;
float *vsave;
float F(float v);
float frand();
void maxwellv(float *vx_local, float *vy_local, float *vz_local, float vti);
void SetParticleLoad(int isp, float Ninit, int load_type, float x_left,float x_right, float y_top, float y_bottom, float vti);
#endif