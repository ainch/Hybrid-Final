#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include "def.h"

#ifndef __XYPIC_C_H__
#define __XYPIC_C_H__
// C Function declaration
void display_title();
void display_finish();
void InputRead(int argc, char *argv[]);
void start();
void DumpRead(int argc, char *argv[]);


int ***TIMalloc(int sizeX,int sizeY,int sizeZ);
float ***TFMalloc(int sizeX,int sizeY,int sizeZ);
int **MIMalloc(int sizeX,int sizeY);
float **MFMalloc(int sizeX,int sizeY);
int *VIMalloc(int size);
float *VFMalloc(int size);

void TFFree(float ***T,int sizeX,int sizeY);
void TIFree(int ***T,int sizeX,int sizeY);
void MFFree(float **M,int sizeX);
void MIFree(int **M,int sizeX);

void TFInit(float ***T,float C,int sizeX,int sizeY,int sizeZ);
void TIInit(int ***T,int C,int sizeX,int sizeY,int sizeZ);
void MFInit(float **M,float C,int sizeX,int sizeY);
void MIInit(int **M,int C,int sizeX,int sizeY);
void VFInit(float *V,float C,int size);
void VIInit(int *V,int C,int size);

void MFCopy(float **M,float **C,int sizeX,int sizeY);
void MICopy(int **M,int **C,int sizeX,int sizeY);
void VFCopy(float *V,float *C,int size);
void VICopy(int *V,int *C,int size);
#endif
#ifndef __XYPIC_CU_H__
#define __XYPIC_CU_H__
// CU Function declaration
extern void main_cuda();
#endif



