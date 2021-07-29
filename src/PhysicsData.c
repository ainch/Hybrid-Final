#include "PhysicsData.h"

void PT_DATA(int MainG, Species* SP){

	if (spnum == 10) { // Electron
		strcpy(sp[isp].name,"E");
		sp[isp].m = 5.485799486e-4 * AMU;
		sp[isp].q = -1.0 * CQ;
	} else if (spnum == 13) { //Ar+   argon ion
		strcpy(sp[isp].name,"AR+");
		sp[isp].m = 39.9500 * AMU;
		sp[isp].q = 1.0* CQ;
	} else if (spnum == 21) { // O-  Negative ion
		strcpy(sp[isp].name,"O-");
		sp[isp].m = 16.0000 * AMU;
		sp[isp].q = -1.0* CQ;
	} else if (spnum == 22) { // O+ O ion
		strcpy(sp[isp].name,"O+");
		sp[isp].m = 16.0000 * AMU;
		sp[isp].q = 1.0* CQ;
	} else if (spnum == 31) { // O2+ ion
		strcpy(sp[isp].name,"O2+");
		sp[isp].m = 32.0000 * AMU;
		sp[isp].q = 1.0* CQ;
	} else {
		fprintf(stderr, "particle number  %d is incorrect\n", spnum);
		exit(1);
	}
}
void FD_DATA(int isp, int spnum, Fluid * FDs){
	if (spnum == 12) { // Argon metastable
		strcpy(fluid_sp[isp].name,"AR*");
		fluid_sp[isp].m = 39.9500 * AMU;
	} else if (spnum == 20) { // O atom
		strcpy(fluid_sp[isp].name,"OP");
		fluid_sp[isp].m = 16.0000 * AMU;
	} else if (spnum == 23) { // OD O metastable atom
		strcpy(fluid_sp[isp].name,"OD");
		fluid_sp[isp].m = 16.0000 * AMU;
	} else if (spnum == 32) { // O2A metastable
		strcpy(fluid_sp[isp].name,"O2A");
		fluid_sp[isp].m = 32.0000 * AMU;
	} else if (spnum == 33) { // O2B metastable
		strcpy(fluid_sp[isp].name,"O2B");
		fluid_sp[isp].m = 32.0000 * AMU;
	} else {
		fprintf(stderr, "Fluid species number  %d is incorrect\n", spnum);
		exit(1);
	}
}
void BG_DATA(int isp, int spnum, BackG *BGs){
	if (spnum == 11) { // Argon
		strcpy(BG_sp[isp].name,"AR");
		BG_sp[isp].m = 39.9500 * AMU;;
	} else if (spnum == 30) { // O2 Oxygen
		strcpy(BG_sp[isp].name,"O2");
		BG_sp[isp].m = 32.0000 * AMU;;
	} else {
		fprintf(stderr, "Background species number  %d is incorrect\n", spnum);
		exit(1);
	}
}