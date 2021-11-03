#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include "def.h"
#include "parson.h"

extern int ***TIMalloc(int sizeX,int sizeY,int sizeZ);
extern float ***TFMalloc(int sizeX,int sizeY,int sizeZ);
extern int **MIMalloc(int sizeX,int sizeY);
extern float **MFMalloc(int sizeX,int sizeY);
extern int *VIMalloc(int size);
extern float *VFMalloc(int size);
extern void TFFree(float ***T,int sizeX,int sizeY);
extern void TIFree(int ***T,int sizeX,int sizeY);
extern void MFFree(float **M,int sizeX);
extern void MIFree(int **M,int sizeX);
extern void TFInit(float ***T,float C,int sizeX,int sizeY,int sizeZ);
extern void TIInit(int ***T,int C,int sizeX,int sizeY,int sizeZ);
extern void MFInit(float **M,float C,int sizeX,int sizeY);
extern void MIInit(int **M,int C,int sizeX,int sizeY);
extern void VFInit(float *V,float C,int size);
extern void VIInit(int *V,int C,int size);
extern void MFCopy(float **M,float **C,int sizeX,int sizeY);
extern void MICopy(int **M,int **C,int sizeX,int sizeY);
extern void VFCopy(float *V,float *C,int size);
extern void VICopy(int *V,int *C,int size);
#ifndef __XYPIC_C_H__
#define __XYPIC_C_H__

#endif
#ifndef __XYPIC_CU_H__
#define __XYPIC_CU_H__

#endif
//https://m.blog.naver.com/PostView.naver?isHttpsRedirect=true&blogId=icysword&logNo=140202424642 
//Stream concept


