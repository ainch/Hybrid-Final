#include "main.h"

void display_title()
{
	puts("------------------------------------");
	puts("|    XY PIC Simulation using PCG   |");
	puts("|           2021. 03. 01           |");
	puts("|              G.PARK              |");
	puts("|           Hae June Lee           |");
	puts("|   Pusan Plasma Research Center   |");
	puts("|    Pusan National University     |");
	puts("------------------------------------");
}
void display_finish()
{
	fprintf(stderr,"---------------------------\n");
	fprintf(stderr,"|   Simulation finished   |\n");
	fprintf(stderr,"---------------------------\n");
}
int main(int argc, char *argv[])
{
	t = 0.0; // Time step
	tstep = 0; // INT Time step
	cstep = 0; // INT Cycle step
	DumpFlag = 0;
	if (argc == 5) {
		DumpFlag = 1;
		strcpy(InputFile, argv[2]);
		strcpy(DumpFile, argv[4]);
	} else if (argc == 3) {
		strcpy(InputFile, argv[2]);
		DumpFlag = 0;
	} else {
		puts(" 2D_Hybrid_PIC -i <InputFile>");
		puts(" \tor");
		puts(" 2D_Hybrid_PIC -i <InputFile> -d <DumpFile>");
		exit(1);
	}
	//seed = 1;
	seed = (long) getpid();
	display_title();       /* Display 2D Field Solver title */
	InputRead();	 	/* InputFile Parsing */
	Geometry_setting(); 	
	Source_setting();
	FieldSolverSetting();
	GasSetting();
	if(TecplotS_Gsize_Flag) 	Main_Variable_printorSave();
	if(TecplotS_CX_Flag) 		Cross_Section_data_Save();
	if(TecplotS_Particle_Flag) 	Initial_Particle_Save(TecplotS_Particle_Num,PtD);
	if(DumpFlag) LoadDumpFile();	/* DumpFile Read */
	main_cuda();			/* GPU Start */
 	display_finish();
 	return 0;
}