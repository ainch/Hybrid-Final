#include "Tecplot.h"

void Initial_Particle_Save(int size,HCP *PtData){
	FILE *fp;
	char filename[512];
	int i,j,index,isp;
	float **buf_x,**buf_y,**buf_vx,**buf_vy,**buf_vz;
	int *SAVE_NP;

	buf_x = MFMalloc(nsp,size);
	buf_y = MFMalloc(nsp,size);
	buf_vx = MFMalloc(nsp,size);
	buf_vy = MFMalloc(nsp,size);
	buf_vz = MFMalloc(nsp,size);
	SAVE_NP = VIMalloc(nsp);
	MFInit(buf_x,0.0,nsp,size);
	MFInit(buf_x,0.0,nsp,size);
	MFInit(buf_x,0.0,nsp,size);
	MFInit(buf_x,0.0,nsp,size);
	MFInit(buf_x,0.0,nsp,size);
	VIInit(SAVE_NP,0,nsp);
	for(isp=0;isp<nsp;isp++){
		if(SP[isp].np>size){
			SAVE_NP[isp] = size;
			for(i=0;i<SAVE_NP[isp];i++){
				//printf("isp= %d, i = %d\n",isp,i);
				index=frand()*((float)SP[isp].np-1.0);
				buf_x[isp][i]=dx*PtData[isp].x[index];
				buf_y[isp][i]=dy*PtData[isp].y[index];
				buf_vx[isp][i]=PtData[isp].vx[index];
				buf_vy[isp][i]=PtData[isp].vy[index];
				buf_vz[isp][i]=PtData[isp].vz[index];
			}
		}else{
			SAVE_NP[isp]=SP[isp].np;
			for(i=0;i<SAVE_NP[isp];i++){
				//printf("isp= %d, i = %d\n",isp,i);
				buf_x[isp][i]=dx*PtData[isp].x[i];
				buf_y[isp][i]=dy*PtData[isp].y[i];
				buf_vx[isp][i]=PtData[isp].vx[i];
				buf_vy[isp][i]=PtData[isp].vy[i];
				buf_vz[isp][i]=PtData[isp].vz[i];
			}
		}
	}
	sprintf(filename,"Initial_Pt_TEC.dat");
	fp = fopen(filename,"w");
	fprintf(fp, "TITLE = \"2D PIC Movie\"\n");
	fprintf(fp, "VARIABLES = \"X(m)\",\"Y(m)\",\"Vx(m/s)\",\"Vy(m/s)\",\"Vz(m/s)\"\n");
	// GEOMETRY
    fprintf(fp, "GEOMETRY\n");
    fprintf(fp, "F=POINT\n");
    fprintf(fp, "CS=GRID\n");
    fprintf(fp, "X=0.00,Y=0.00,Z=0.00\n");
    fprintf(fp, "C=BLACK\n");
    fprintf(fp, "S=GLOBAL\n");
    fprintf(fp, "L=SOLID\n");
    fprintf(fp, "PL=4\n");
    fprintf(fp, "LT=0.1\n");
    fprintf(fp, "CLIPPING=CLIPTOVIEWPORT\n");
    fprintf(fp, "DRAWORDER=AFTERDATA\n");
    fprintf(fp, "MFC=\"\"\n");
    fprintf(fp, "T=RECTANGLE %g %g\n",xlength,ylength);
    for(i=0;i<CondNUM;i++){
        fprintf(fp, "GEOMETRY\n");
        fprintf(fp, "F=POINT\n");
        fprintf(fp, "CS=GRID\n");
        fprintf(fp, "X=%g,Y=%g,Z=0.00\n",CondX0[i]*dx,CondY0[i]*dy);
        fprintf(fp, "C=BLACK\n");
        fprintf(fp, "S=GLOBAL\n");
        fprintf(fp, "L=SOLID\n");
        fprintf(fp, "PL=4\n");
        fprintf(fp, "LT=0.1\n");
        fprintf(fp, "CLIPPING=CLIPTOVIEWPORT\n");
        fprintf(fp, "DRAWORDER=AFTERDATA\n");
        fprintf(fp, "MFC=\"\"\n");
        fprintf(fp, "T=RECTANGLE %g %g\n",CondX1[i]*dx-CondX0[i]*dx,CondY1[i]*dy-CondY0[i]*dy);
    }
    for(i=0;i<DielNUM;i++){
        fprintf(fp, "GEOMETRY\n");
        fprintf(fp, "F=POINT\n");
        fprintf(fp, "CS=GRID\n");
        fprintf(fp, "X=%g,Y=%g,Z=0.00\n",DielX0[i]*dx,DielY0[i]*dy);
        fprintf(fp, "C=BLACK\n");
        fprintf(fp, "S=GLOBAL\n");
        fprintf(fp, "L=SOLID\n");
        fprintf(fp, "PL=4\n");
        fprintf(fp, "LT=0.1\n");
        fprintf(fp, "CLIPPING=CLIPTOVIEWPORT\n");
        fprintf(fp, "DRAWORDER=AFTERDATA\n");
        fprintf(fp, "MFC=\"\"\n");
        fprintf(fp, "T=RECTANGLE %g %g\n",DielX1[i]*dx-DielX0[i]*dx,DielY1[i]*dy-DielY0[i]*dy);
    }
	for(isp=0;isp<nsp;isp++){
		fprintf(fp, "ZONE T=\"ZONE %d\"\n",isp+1);
		fprintf(fp, " STRANDID=0, SOLUTIONTIME=1\n");
		fprintf(fp, " I=%d, J=1, K=1, ZONETYPE=Ordered\n",SAVE_NP[isp]);
		fprintf(fp, " DATAPACKING=POINT\n");
		fprintf(fp, " DT=(SINGLE SINGLE SINGLE SINGLE SINGLE)\n");
		for (i = 0; i < SAVE_NP[isp]; i++) {
			fprintf(fp,"%0.3e %0.3e %0.3e %0.3e %0.3e\n",buf_x[isp][i],buf_y[isp][i],buf_vx[isp][i],buf_vy[isp][i],buf_vz[isp][i]);
		}
	}
	fclose(fp);	
	MFFree(buf_x,nsp);
	MFFree(buf_y,nsp);
	MFFree(buf_vx,nsp);
	MFFree(buf_vy,nsp);
	MFFree(buf_vz,nsp);
	free(SAVE_NP);
}
void Cross_Section_data_Save(){
	FILE *fp;
	int i,kk;
	int nbar;
	char filename[512];
    
    if(MainGas==ARGON){
		sprintf(filename,"CrossX_Ar_TEC.dat");
		fp = fopen(filename,"w");
		fprintf(fp, "TITLE = \"Argon Cross Section Data\"\n");
		// VARIABLES NAMES
		fprintf(fp, "VARIABLES = \"Energy (eV)\",\n");
		fprintf(fp,"\"R0-e+Ar>e+Ar\",\n");
		fprintf(fp,"\"R1-e+Ar>e+Ar*\",\n");
		fprintf(fp,"\"R2-e+Ar>e+Ar*m\",\n");
		fprintf(fp,"\"R3-e+Ar>2e+Ar^\",\n");
		fprintf(fp,"\"R4-e+Ar*m>2e+Ar^\",\n");
		fprintf(fp,"\"R5-Ar+Ar^>Ar^+Ar\",\n");
		fprintf(fp,"\"R6-Ar+Ar^>Ar+Ar^\",\n");
		// NUMBER 0F ZONE
		fprintf(fp, "ZONE I = %d, F=BLOCK,\n", N_LOGX);
		fprintf(fp, " T=\"1-D DATA\",\n");
		// DATA1
		//
		nbar = 6;
		// E_array
		kk = 0;fprintf(fp, "\t");
		for(i=0;i<N_LOGX;i++){
			fprintf(fp, "%.3e\t",Ar_Data[i].xee);
			kk++;
			if(kk==nbar){
				fprintf(fp, "\n\t");
				kk=0;
			}
		}
		kk = 0;fprintf(fp, "\t");
		for(i=0;i<N_LOGX;i++){
			fprintf(fp, "%.3e\t",Ar_Data[i].cx_0);
			kk++;
			if(kk==nbar){
				fprintf(fp, "\n\t");
				kk=0;
			}
		}
		kk = 0;fprintf(fp, "\t");
		for(i=0;i<N_LOGX;i++){
			fprintf(fp, "%.3e\t",Ar_Data[i].cx_1);
			kk++;
			if(kk==nbar){
				fprintf(fp, "\n\t");
				kk=0;
			}
		}
		kk = 0;fprintf(fp, "\t");
		for(i=0;i<N_LOGX;i++){
			fprintf(fp, "%.3e\t",Ar_Data[i].cx_2);
			kk++;
			if(kk==nbar){
				fprintf(fp, "\n\t");
				kk=0;
			}
		}
		kk = 0;fprintf(fp, "\t");
		for(i=0;i<N_LOGX;i++){
			fprintf(fp, "%.3e\t",Ar_Data[i].cx_3);
			kk++;
			if(kk==nbar){
				fprintf(fp, "\n\t");
				kk=0;
			}
		}
		kk = 0;fprintf(fp, "\t");
		for(i=0;i<N_LOGX;i++){
			fprintf(fp, "%.3e\t",Ar_Data[i].cx_4);
			kk++;
			if(kk==nbar){
				fprintf(fp, "\n\t");
				kk=0;
			}
		}
		kk = 0;fprintf(fp, "\t");
		for(i=0;i<N_LOGX;i++){
			fprintf(fp, "%.3e\t",Ar_Data[i].cx_5);
			kk++;
			if(kk==nbar){
				fprintf(fp, "\n\t");
				kk=0;
			}
		}
		kk = 0;fprintf(fp, "\t");
		for(i=0;i<N_LOGX;i++){
			fprintf(fp, "%.3e\t",Ar_Data[i].cx_6);
			kk++;
			if(kk==nbar){
				fprintf(fp, "\n\t");
				kk=0;
			}
		}
		fclose(fp);
    }else if(MainGas==OXYGEN){
		sprintf(filename,"CrossX_O2_TEC.dat");
		fp = fopen(filename,"w");
		fprintf(fp, "TITLE = \"Oxygen Cross Section Data\"\n");
		// VARIABLES NAMES
		fprintf(fp, "VARIABLES = \"Energy (eV)\",\n");
        fprintf(fp,"\" R0-e+O2>e+O2\",\n");
        fprintf(fp,"\" R1-e+O2>e+O2*\",\n");
        fprintf(fp,"\" R2-e+O2>e+O2*\",\n");
        fprintf(fp,"\" R3-e+O2>e+O2A\",\n");
        fprintf(fp,"\" R4-e+O2>e+O2B\",\n");
        fprintf(fp,"\" R5-e+O2>e+O2*\",\n");
        fprintf(fp,"\" R6-e+O2>OP+O-\",\n");
        fprintf(fp,"\" R7-e+O2>e+2OP\",\n");
        fprintf(fp,"\" R8-e+O2>e+OP+OD\",\n");
        fprintf(fp,"\" R9-e+O2>e+2OD\",\n");
        fprintf(fp,"\"R10-e+O2>2e+O2^\",\n");
        fprintf(fp,"\"R11-e+O2>e+OP+O*\",\n");
        fprintf(fp,"\"R12-e+O2>e+O^+O-\",\n");
        fprintf(fp,"\"R13-e+O2>2e+O^+OP\",\n");
        fprintf(fp,"\"R14-e+O2A>2e+O2+\",\n");
        fprintf(fp,"\"R15-e+O2A>OP+O-\",\n");
        fprintf(fp,"\"R16-e+O2A>e+O2\",\n");
        fprintf(fp,"\"R17-e+O2A>e+O2\",\n");
        fprintf(fp,"\"R18-e+O2A>e+2OP\",\n");
        fprintf(fp,"\"R19-e+O2A>e+OP+OD\",\n");
        fprintf(fp,"\"R20-e+O2A>e+2OD\",\n");
        fprintf(fp,"\"R21-e+O2A>2e+O^+OP\",\n");
        fprintf(fp,"\"R22-e+O2B>2e+O2^\",\n");
        fprintf(fp,"\"R23-e+O2B>OP+O-\",\n");
        fprintf(fp,"\"R24-e+O2B>e+O2\",\n");
        fprintf(fp,"\"R25-e+O2B>e+O2\",\n");
        fprintf(fp,"\"R26-e+O2B>e+2O\",\n");
        fprintf(fp,"\"R27-e+O2B>e+OP+OD\",\n");
        fprintf(fp,"\"R28-e+O2B>e+2OD\",\n");
        fprintf(fp,"\"R29-e+O2B>2e+O^+OP\",\n");
        fprintf(fp,"\"R30-e+O->2e+OP\",\n");
        fprintf(fp,"\"R31-e+O2^>OP+OD\",\n");
        fprintf(fp,"\"R32-e+OP>e+OP\",\n");
        fprintf(fp,"\"R33-e+OP>e+OD\",\n");
        fprintf(fp,"\"R34-e+OP>e+O*\",\n");
        fprintf(fp,"\"R35-e+OP>e+O*\",\n");
        fprintf(fp,"\"R36-e+OP>e+O*\",\n");
        fprintf(fp,"\"R37-e+OP>2e+O^\",\n");
        fprintf(fp,"\"R38-e+OP>e+O*\",\n");
        fprintf(fp,"\"R39-e+OD>2e+O^\",\n");
        fprintf(fp,"\"R40-e+OD>e+OP\",\n");
        fprintf(fp,"\"R41-O-+O2>O-+O2\",\n");
        fprintf(fp,"\"R42-O-+O2>e+OP+O2\",\n");
        fprintf(fp,"\"R43-O-+OP>e+O2\",\n");
        fprintf(fp,"\"R44-O-+O2^>OP+O2\",\n");
        fprintf(fp,"\"R45-O-+O^>2OP\",\n");
        fprintf(fp,"\"R46-O-+O2A>e+OP+O2\",\n");
        fprintf(fp,"\"R47-O2^+OP>O2+O^\",\n");
        fprintf(fp,"\"R48-O2^+O2>O2+O2^\",\n");
        fprintf(fp,"\"R49-O2^+O2>O2^+O2\",\n");
        fprintf(fp,"\"R50-O2^+O2>O^+OP+O2\",\n");
        fprintf(fp,"\"R51-O2^+O2A>O2+O2^\",\n");
        fprintf(fp,"\"R52-O2^+O2B>O2+O2^\",\n");
        fprintf(fp,"\"R53-O^+O2>OP+O2^\",\n");
        fprintf(fp,"\"R54-O^+O2>O^+O2\",\n");
        fprintf(fp,"\"R55-O^+OP>OP+O^\",\n");
        fprintf(fp,"\"R56-O^+O2A>O2^+OP\",\n");
        fprintf(fp,"\"R57-O^+O2B>O2^+OP\",\n");
		// NUMBER 0F ZONE
		fprintf(fp, "ZONE I = %d, F=BLOCK,\n", N_LOGX);
		fprintf(fp, " T=\"1-D DATA\",\n");
		// DATA1
		//
		nbar = 6;
		// E_array
		kk = 0;fprintf(fp, "\t");
		for(i=0;i<N_LOGX;i++){
			fprintf(fp, "%.3e\t",O2_Data[i].xee);
			kk++;
			if(kk==nbar){
				fprintf(fp, "\n\t");
				kk=0;
			}
		}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_0);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_1);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_2);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_3);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_4);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_5);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_6);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_7);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_8);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_9);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_10);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_11);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_12);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_13);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_14);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_15);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_16);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_17);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_18);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_19);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_20);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_21);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_22);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_23);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_24);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_25);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_26);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_27);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_28);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_29);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_30);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_31);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_32);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_33);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_34);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_35);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_36);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_37);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_38);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_39);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_40);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_41);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_42);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_43);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_44);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_45);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_46);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_47);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_48);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_49);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_50);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_51);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_52);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_53);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_54);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_55);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_56);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",O2_Data[i].cx_57);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		fclose(fp);
    }else if(MainGas==ARO2){
		sprintf(filename,"CrossX_ArO2_TEC.dat");
		fp = fopen(filename,"w");
		fprintf(fp, "TITLE = \"Argon+Oxygen Cross Section Data\"\n");
		// VARIABLES NAMES
		fprintf(fp, "VARIABLES = \"Energy (eV)\",\n");
        fprintf(fp,"\" R0-e+Ar>e+Ar\",\n");
        fprintf(fp,"\" R1-e+Ar>e+Ar*\",\n");
        fprintf(fp,"\" R2-e+Ar>e+Ar*m\",\n");
        fprintf(fp,"\" R3-e+Ar>2e+Ar+\",\n");
        fprintf(fp,"\" R4-e+Ar*m>2e+Ar^\",\n");
        fprintf(fp,"\" R5-e+O2>e+O2\",\n");
        fprintf(fp,"\" R6-e+O2>e+O2*\",\n");
        fprintf(fp,"\" R7-e+O2>e+O2*\",\n");
        fprintf(fp,"\" R8-e+O2>e+O2A\",\n");
        fprintf(fp,"\" R9-e+O2>e+O2B\",\n");
        fprintf(fp,"\"R10-e+O2>e+O2*\",\n");
        fprintf(fp,"\"R11-e+O2>OP+O-\",\n");
        fprintf(fp,"\"R12-e+O2>e+2OP\",\n");
        fprintf(fp,"\"R13-e+O2>e+OP+OD\",\n");
        fprintf(fp,"\"R14-e+O2>e+2OD\",\n");
        fprintf(fp,"\"R15-e+O2>2e+O2+\",\n");
        fprintf(fp,"\"R16-e+O2>e+OP+O*\",\n");
        fprintf(fp,"\"R17-e+O2>e+O++O-\",\n");
        fprintf(fp,"\"R18-e+O2>2e+O^+OP\",\n");
        fprintf(fp,"\"R19-e+O2A>2e+O2^\",\n");
        fprintf(fp,"\"R20-e+O2A>OP+O-\",\n");
        fprintf(fp,"\"R21-e+O2A>e+O2\",\n");
        fprintf(fp,"\"R22-e+O2A>e+O2\",\n");
        fprintf(fp,"\"R23-e+O2A>e+2O\",\n");
        fprintf(fp,"\"R24-e+O2A>e+OP+OD\",\n");
        fprintf(fp,"\"R25-e+O2A>e+2OD\",\n");
        fprintf(fp,"\"R26-e+O2A>2e+O^+OP\",\n");
        fprintf(fp,"\"R27-e+O2B>2e+O2^\",\n");
        fprintf(fp,"\"R28-e+O2B>OP+O-\",\n");
        fprintf(fp,"\"R29-e+O2B>e+O2\",\n");
        fprintf(fp,"\"R30-e+O2B>e+O2\",\n");
        fprintf(fp,"\"R31-e+O2B>e+2O\",\n");
        fprintf(fp,"\"R32-e+O2B>e+OP+OD\",\n");
        fprintf(fp,"\"R33-e+O2B>e+2OD\",\n");
        fprintf(fp,"\"R34-e+O2B>2e+O++OP\",\n");
        fprintf(fp,"\"R35-e+O->2e+OP\",\n");
        fprintf(fp,"\"R36-e+O2+>OP+OD \",\n");
        fprintf(fp,"\"R37-e+OP>e+OP\",\n");
        fprintf(fp,"\"R38-e+OP>e+OD\",\n");
        fprintf(fp,"\"R39-e+OP>e+O*\",\n");
        fprintf(fp,"\"R40-e+OP>e+O*\",\n");
        fprintf(fp,"\"R41-e+OP>e+O*\",\n");
        fprintf(fp,"\"R42-e+OP>2e+O^\",\n");
        fprintf(fp,"\"R43-e+OP>e+O*\",\n");
        fprintf(fp,"\"R44-e+OD>2e+O+\",\n");
        fprintf(fp,"\"R45-e+OD>e+O\",\n");
        fprintf(fp,"\"R46-O-+O2>O-+O2\",\n");
        fprintf(fp,"\"R47-O-+O2>e+OP+O2\",\n");
        fprintf(fp,"\"R48-O-+OP>e+O2\",\n");
        fprintf(fp,"\"R49-O-+O2^>OP+O2\",\n");
        fprintf(fp,"\"R50-O-+O^>2OP\",\n");
        fprintf(fp,"\"R51-O-+O2A>e+OP+O2\",\n");
        fprintf(fp,"\"R52-O2^+OP>O2+O^\",\n");
        fprintf(fp,"\"R53-O2^+O2>O2+O2^\",\n");
        fprintf(fp,"\"R54-O2^+O2>O2^+O2\",\n");
        fprintf(fp,"\"R55-O2^+O2>O^+OP+O2\",\n");
        fprintf(fp,"\"R56-O2^+O2A>O2+O2^\",\n");
        fprintf(fp,"\"R57-O2^+O2B>O2+O2^\",\n");
        fprintf(fp,"\"R58-O2^+Ar>O2+Ar^\",\n");
        fprintf(fp,"\"R59-O2^+Ar>O2^+Ar^\",\n");
        fprintf(fp,"\"R60-O^+O2>OP+O2^\",\n");
        fprintf(fp,"\"R61-O^+O2>O^+O2\",\n");
        fprintf(fp,"\"R62-O^+OP>OP+O^\",\n");
        fprintf(fp,"\"R63-O^+O2A>O2^+OP\",\n");
        fprintf(fp,"\"R64-O^+O2B>O2^+OP\",\n");
        fprintf(fp,"\"R65-Ar^+Ar>Ar+Ar^\",\n");
        fprintf(fp,"\"R66-Ar^+Ar>Ar++Ar\",\n");
        fprintf(fp,"\"R67-Ar^+O2>O2+Ar^\",\n");
		// NUMBER 0F ZONE
		fprintf(fp, "ZONE I = %d, F=BLOCK,\n", N_LOGX);
		fprintf(fp, " T=\"1-D DATA\",\n");
		// DATA1
		//
		nbar = 6;
		// E_array
		kk = 0;fprintf(fp, "\t");
		for(i=0;i<N_LOGX;i++){
			fprintf(fp, "%.3e\t",ArO2_Data[i].xee);
			kk++;
			if(kk==nbar){
				fprintf(fp, "\n\t");
				kk=0;
			}
		}
		// 2.value
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_0);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_1);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_2);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_3);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_4);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_5);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_6);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_7);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_8);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_9);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_10);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_11);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_12);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_13);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_14);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_15);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_16);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_17);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_18);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_19);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_20);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_21);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_22);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_23);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_24);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_25);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_26);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_27);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_28);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_29);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_30);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_31);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_32);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_33);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_34);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_35);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_36);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_37);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_38);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_39);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_40);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_41);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_42);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_43);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_44);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_45);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_46);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_47);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_48);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_49);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_50);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_51);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_52);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_53);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_54);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_55);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_56);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_57);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_58);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_59);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_60);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_61);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_62);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_63);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_64);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_65);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_66);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		kk=0;fprintf(fp,"\t");for(i=0;i<N_LOGX;i++){fprintf(fp,"%.3e\t",ArO2_Data[i].cx_67);kk++;if(kk==nbar){fprintf(fp,"\n\t");kk=0;}}
		fclose(fp);
    }
}
void Main_Variable_printorSave(){
    int i,j;
    FILE *fp;
    int nbar,kk,Gbuf;
   //int Boundary;   //Boundary Condition Constant 0~4
	//int CondID;		// Conductor ID , NO Conductor is zero
	//int Face;
	//float Area;
	//float Temp;
	//float BackDen1;
	//float BackVel1;
	//float BackDen2;
	//float BackVel2;
	//float Lap_Pot;
	//float Pois_Pot;
	//float Ex;
	//float Ey;
    fp = fopen("GsizeData_TEC.dat", "w");
    fprintf(fp, "TITLE = \"2D-PIC Gsize data set\"\n");
    fprintf(fp, "VARIABLES = \"X\",\"Y\",\n");
    fprintf(fp, "\"Boundary\",\n");
    fprintf(fp, "\"CondID\",\n");
	fprintf(fp, "\"Face\",\n");
	fprintf(fp, "\"Area\",\n");
    fprintf(fp, "\"Temp\",\n");
    fprintf(fp, "\"BackDen1\",\n");
	fprintf(fp, "\"BackVel1\",\n");
	fprintf(fp, "\"BackDen2\",\n");
	fprintf(fp, "\"BackVel2\",\n");
	fprintf(fp, "\"Lap_Pot\",\n");
	fprintf(fp, "\"Pois_Pot\",\n");
	fprintf(fp, "\"Ex\",\n");
	fprintf(fp, "\"Ey\",\n");
	fprintf(fp, "ZONE I = %d, J = %d, F=BLOCK,\n", ngx, ngy);
    nbar = 6;
    // x_Garray
    kk = 0;for(j=0;j<ngy;j++){fprintf(fp, "\t");for(i=0;i<ngx;i++){
    fprintf(fp, "%.3e\t",
			x_Garray[i]
    );kk++;if(kk==nbar){fprintf(fp, "\n\t");kk=0;}}kk=0;fprintf(fp, "\n");}
    // y_Garray
    kk = 0;for(j=0;j<ngy;j++){fprintf(fp, "\t");for(i=0;i<ngx;i++){
    fprintf(fp, "%.3e\t",
			y_Garray[j]
    );kk++;if(kk==nbar){fprintf(fp, "\n\t");kk=0;}}kk=0;fprintf(fp, "\n");}

    kk = 0;for(j=0;j<ngy;j++){fprintf(fp, "\t");for(i=0;i<ngx;i++){
    fprintf(fp, "%d\t",
			vec_G[i*ngy+j].Boundary
    );kk++;if(kk==nbar){fprintf(fp, "\n\t");kk=0;}}kk=0;fprintf(fp, "\n");}
    kk = 0;for(j=0;j<ngy;j++){fprintf(fp, "\t");for(i=0;i<ngx;i++){
    fprintf(fp, "%d\t",
			vec_G[i*ngy+j].CondID
    );kk++;if(kk==nbar){fprintf(fp, "\n\t");kk=0;}}kk=0;fprintf(fp, "\n");}
    kk = 0;for(j=0;j<ngy;j++){fprintf(fp, "\t");for(i=0;i<ngx;i++){
    fprintf(fp, "%d\t",
			vec_G[i*ngy+j].Face
    );kk++;if(kk==nbar){fprintf(fp, "\n\t");kk=0;}}kk=0;fprintf(fp, "\n");}
    kk = 0;for(j=0;j<ngy;j++){fprintf(fp, "\t");for(i=0;i<ngx;i++){
    fprintf(fp, "%.3e\t",
			vec_G[i*ngy+j].Area
    );kk++;if(kk==nbar){fprintf(fp, "\n\t");kk=0;}}kk=0;fprintf(fp, "\n");}
	kk = 0;for(j=0;j<ngy;j++){fprintf(fp, "\t");for(i=0;i<ngx;i++){
    fprintf(fp, "%.3e\t",
			vec_G[i*ngy+j].Temp
    );kk++;if(kk==nbar){fprintf(fp, "\n\t");kk=0;}}kk=0;fprintf(fp, "\n");}
    kk = 0;for(j=0;j<ngy;j++){fprintf(fp, "\t");for(i=0;i<ngx;i++){
    fprintf(fp, "%.3e\t",
			vec_G[i*ngy+j].BackDen1
    );kk++;if(kk==nbar){fprintf(fp, "\n\t");kk=0;}}kk=0;fprintf(fp, "\n");}
	kk = 0;for(j=0;j<ngy;j++){fprintf(fp, "\t");for(i=0;i<ngx;i++){
    fprintf(fp, "%.3e\t",
			vec_G[i*ngy+j].BackVel1
    );kk++;if(kk==nbar){fprintf(fp, "\n\t");kk=0;}}kk=0;fprintf(fp, "\n");}
	kk = 0;for(j=0;j<ngy;j++){fprintf(fp, "\t");for(i=0;i<ngx;i++){
    fprintf(fp, "%.3e\t",
			vec_G[i*ngy+j].BackDen2
    );kk++;if(kk==nbar){fprintf(fp, "\n\t");kk=0;}}kk=0;fprintf(fp, "\n");}
	kk = 0;for(j=0;j<ngy;j++){fprintf(fp, "\t");for(i=0;i<ngx;i++){
    fprintf(fp, "%.3e\t",
			vec_G[i*ngy+j].BackVel2
    );kk++;if(kk==nbar){fprintf(fp, "\n\t");kk=0;}}kk=0;fprintf(fp, "\n");}
	kk = 0;for(j=0;j<ngy;j++){fprintf(fp, "\t");for(i=0;i<ngx;i++){
    fprintf(fp, "%.3e\t",
			vec_G[i*ngy+j].Lap_Pot
    );kk++;if(kk==nbar){fprintf(fp, "\n\t");kk=0;}}kk=0;fprintf(fp, "\n");}
	kk = 0;for(j=0;j<ngy;j++){fprintf(fp, "\t");for(i=0;i<ngx;i++){
    fprintf(fp, "%.3e\t",
			vec_G[i*ngy+j].Pois_Pot
    );kk++;if(kk==nbar){fprintf(fp, "\n\t");kk=0;}}kk=0;fprintf(fp, "\n");}
	kk = 0;for(j=0;j<ngy;j++){fprintf(fp, "\t");for(i=0;i<ngx;i++){
    fprintf(fp, "%.3e\t",
			vec_G[i*ngy+j].Ex
    );kk++;if(kk==nbar){fprintf(fp, "\n\t");kk=0;}}kk=0;fprintf(fp, "\n");}
	kk = 0;for(j=0;j<ngy;j++){fprintf(fp, "\t");for(i=0;i<ngx;i++){
    fprintf(fp, "%.3e\t",
			vec_G[i*ngy+j].Ey
    );kk++;if(kk==nbar){fprintf(fp, "\n\t");kk=0;}}kk=0;fprintf(fp, "\n");}
	// GEOMETRY
    fprintf(fp, "GEOMETRY\n");
    fprintf(fp, "F=POINT\n");
    fprintf(fp, "CS=GRID\n");
    fprintf(fp, "X=0.00,Y=0.00,Z=0.00\n");
    fprintf(fp, "C=BLACK\n");
    fprintf(fp, "S=GLOBAL\n");
    fprintf(fp, "L=SOLID\n");
    fprintf(fp, "PL=4\n");
    fprintf(fp, "LT=0.1\n");
    fprintf(fp, "CLIPPING=CLIPTOVIEWPORT\n");
    fprintf(fp, "DRAWORDER=AFTERDATA\n");
    fprintf(fp, "MFC=\"\"\n");
    fprintf(fp, "T=RECTANGLE %g %g\n",xlength,ylength);
    for(i=0;i<CondNUM;i++){
        fprintf(fp, "GEOMETRY\n");
        fprintf(fp, "F=POINT\n");
        fprintf(fp, "CS=GRID\n");
        fprintf(fp, "X=%g,Y=%g,Z=0.00\n",CondX0[i]*dx,CondY0[i]*dy);
        fprintf(fp, "C=BLACK\n");
        fprintf(fp, "S=GLOBAL\n");
        fprintf(fp, "L=SOLID\n");
        fprintf(fp, "PL=4\n");
        fprintf(fp, "LT=0.1\n");
        fprintf(fp, "CLIPPING=CLIPTOVIEWPORT\n");
        fprintf(fp, "DRAWORDER=AFTERDATA\n");
        fprintf(fp, "MFC=\"\"\n");
        fprintf(fp, "T=RECTANGLE %g %g\n",CondX1[i]*dx-CondX0[i]*dx,CondY1[i]*dy-CondY0[i]*dy);
    }
    for(i=0;i<DielNUM;i++){
        fprintf(fp, "GEOMETRY\n");
        fprintf(fp, "F=POINT\n");
        fprintf(fp, "CS=GRID\n");
        fprintf(fp, "X=%g,Y=%g,Z=0.00\n",DielX0[i]*dx,DielY0[i]*dy);
        fprintf(fp, "C=BLACK\n");
        fprintf(fp, "S=GLOBAL\n");
        fprintf(fp, "L=SOLID\n");
        fprintf(fp, "PL=4\n");
        fprintf(fp, "LT=0.1\n");
        fprintf(fp, "CLIPPING=CLIPTOVIEWPORT\n");
        fprintf(fp, "DRAWORDER=AFTERDATA\n");
        fprintf(fp, "MFC=\"\"\n");
        fprintf(fp, "T=RECTANGLE %g %g\n",DielX1[i]*dx-DielX0[i]*dx,DielY1[i]*dy-DielY0[i]*dy);
    }
    fclose(fp);
}
void Field_Laplace_Solution_Save(char *Filename,float **Sol){
    int i,j,k;
    FILE *fp;
	char filename[512];
    int nbar,kk,Gbuf;
    sprintf(filename,"%s_Solution_TEC.dat",Filename);
    fp = fopen(filename, "w");
    fprintf(fp, "TITLE = \"2D-PIC CPU PCG SOLUTION\"\n");
    fprintf(fp, "VARIABLES = \"X\",\"Y\",\n");
	for(k=0;k<CondNUMR;k++){
		fprintf(fp, "\"Cond %d_Solution\",\n",k);
	}
	fprintf(fp, "ZONE I = %d, J = %d, F=BLOCK,\n", ngx, ngy);
    nbar = 6;
    // x_Garray
    kk = 0;for(j=0;j<ngy;j++){fprintf(fp, "\t");for(i=0;i<ngx;i++){
    fprintf(fp, "%.3e\t",
			x_Garray[i]
    );kk++;if(kk==nbar){fprintf(fp, "\n\t");kk=0;}}kk=0;fprintf(fp, "\n");}
    // y_Garray
    kk = 0;for(j=0;j<ngy;j++){fprintf(fp, "\t");for(i=0;i<ngx;i++){
    fprintf(fp, "%.3e\t",
			y_Garray[j]
    );kk++;if(kk==nbar){fprintf(fp, "\n\t");kk=0;}}kk=0;fprintf(fp, "\n");}
	for(k=0;k<CondNUMR;k++){
    	kk = 0;for(j=0;j<ngy;j++){fprintf(fp, "\t");for(i=0;i<ngx;i++){
    	fprintf(fp, "%g\t",
				Sol[k][i*ngy+j]
    	);kk++;if(kk==nbar){fprintf(fp, "\n\t");kk=0;}}kk=0;fprintf(fp, "\n");}
	}
	// GEOMETRY
    fprintf(fp, "GEOMETRY\n");
    fprintf(fp, "F=POINT\n");
    fprintf(fp, "CS=GRID\n");
    fprintf(fp, "X=0.00,Y=0.00,Z=0.00\n");
    fprintf(fp, "C=BLACK\n");
    fprintf(fp, "S=GLOBAL\n");
    fprintf(fp, "L=SOLID\n");
    fprintf(fp, "PL=4\n");
    fprintf(fp, "LT=0.1\n");
    fprintf(fp, "CLIPPING=CLIPTOVIEWPORT\n");
    fprintf(fp, "DRAWORDER=AFTERDATA\n");
    fprintf(fp, "MFC=\"\"\n");
    fprintf(fp, "T=RECTANGLE %g %g\n",xlength,ylength);
    for(i=0;i<CondNUM;i++){
        fprintf(fp, "GEOMETRY\n");
        fprintf(fp, "F=POINT\n");
        fprintf(fp, "CS=GRID\n");
        fprintf(fp, "X=%g,Y=%g,Z=0.00\n",CondX0[i]*dx,CondY0[i]*dy);
        fprintf(fp, "C=BLACK\n");
        fprintf(fp, "S=GLOBAL\n");
        fprintf(fp, "L=SOLID\n");
        fprintf(fp, "PL=4\n");
        fprintf(fp, "LT=0.1\n");
        fprintf(fp, "CLIPPING=CLIPTOVIEWPORT\n");
        fprintf(fp, "DRAWORDER=AFTERDATA\n");
        fprintf(fp, "MFC=\"\"\n");
        fprintf(fp, "T=RECTANGLE %g %g\n",CondX1[i]*dx-CondX0[i]*dx,CondY1[i]*dy-CondY0[i]*dy);
    }
    for(i=0;i<DielNUM;i++){
        fprintf(fp, "GEOMETRY\n");
        fprintf(fp, "F=POINT\n");
        fprintf(fp, "CS=GRID\n");
        fprintf(fp, "X=%g,Y=%g,Z=0.00\n",DielX0[i]*dx,DielY0[i]*dy);
        fprintf(fp, "C=BLACK\n");
        fprintf(fp, "S=GLOBAL\n");
        fprintf(fp, "L=SOLID\n");
        fprintf(fp, "PL=4\n");
        fprintf(fp, "LT=0.1\n");
        fprintf(fp, "CLIPPING=CLIPTOVIEWPORT\n");
        fprintf(fp, "DRAWORDER=AFTERDATA\n");
        fprintf(fp, "MFC=\"\"\n");
        fprintf(fp, "T=RECTANGLE %g %g\n",DielX1[i]*dx-DielX0[i]*dx,DielY1[i]*dy-DielY0[i]*dy);
    }
    fclose(fp);
}