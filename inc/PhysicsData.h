#include "xypic.h"

extern int nsp, nfsp, nBG;
extern Species *SP;// particle species
extern Fluid *FG;	// fluid species
extern BackG *BG;	// background species
extern int TnRct; // Total Number of reaction 
extern int CX_TEC_Flag;
extern CollF *Coll_Flag;

extern float *VFMalloc(int size);
extern void VFInit(float *V,float C,int size);
#ifndef __PHYSICSDATA_H__
#define __PHYSICSDATA_H__
float LOGX_MIN,LOGX_MAX;
float dLOGX,idLOGX;
int N_LOGX;
ArCollD *Ar_Data;
O2CollD *O2_Data;
ArO2CollD *ArO2_Data;
void CrossSection(int* Cxnum, int* CXsize, float *CXx, float *CXy);
float ArSigmaEl(float energy);
float ArSigmaEx(float energy);
float ArSigmaIz(float energy);
float ArSigmaCX(float energy);
float ArSigmaSC(float energy);
float O2sigmaE(float energy);
float O2sigmaR(float energy);
float O2sigmaE1(float energy);
float O2sigmaE2(float energy);
float O2sigmaE3(float energy);
float O2sigmaE4(float energy);
float O2sigmaSD(float energy);
float O2sigmaSS(float energy);
float O2sigmaL1(float energy);
float O2sigmaDA(float energy);
float O2sigmaL2(float energy);
float O2sigmaL3(float energy);
float O2sigmaL4(float energy);
float O2sigmaI(float energy);
float O2sigmaLE(float energy);
float O2sigmaPD(float energy);
float O2sigmaDISSIZ(float energy);
float O2sigmaEO(float energy);
float O2sigmaE1D(float energy);
float O2sigmaE1S(float energy);
float O2sigmaE3P0(float energy);
float O2sigmaE5S0(float energy);
float O2sigmaE3S0(float energy);
float O2sigmaIO(float energy);
float O2sigmaeD(float energy);
float O2sigmaneD(float energy);
float O2sigmanegOD(float energy);
float O2sigmaDR(float energy);
float O2sigmaMN(float energy);
float O2sigmaMNO(float energy);
float O2sigmaCE(float energy);
float O2sigmaCEOO(float energy);
float O2sigmaCEO(float energy);
float O2sigmaCEO2O(float energy);
float O2sigmaFRAG(float energy);
float O2sigmaneE(float energy);
float O2sigmaneuE(float energy);
float O2sigmapoE(float energy);
float O2sigmapooE(float energy);
float O2sigmaO2O2(float energy);
float O2sigmaOO(float energy);
void Argon_CrossSectionSET(CollF *CF);
void Oygen_CrossSectionSET(CollF *CF);
void ArO2_CrossSectionSET(CollF *CF);
#endif