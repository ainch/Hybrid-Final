#include "main.h"

int main(int argc, char *argv[])
{
	t = 0.0;
	tstep = 0;
	cstep = 0;
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
	//display_title();       	/* Display 2D Field Solver title */
	InputRead(argc,argv);	/* InputFile Parsing */
	start();   				/* Initialized */

	DumpRead(argc,argv);	/* DumpFile Read */
	main_cuda();			/* GPU Start */
 	display_finish();
 	return 0;
}

/***************************************************************/
void display_title()
{
	puts("------------------------------------");
	puts("|    XY PIC Simulation using PCG   |");
	puts("|           2021. 03. 01           |");
	puts("|              G.PARK              |");
	puts("|           Hae June Lee           |");
	puts("|   Pusan Plasma Research Center   |");
	puts("|    Pusan National University     |");
	puts("------------------------------------\n");
}
/***************************************************************/

void display_finish()
{
	fprintf(stderr,"---------------------------\n");
	fprintf(stderr,"|   Simulation finished   |\n");
	fprintf(stderr,"---------------------------\n");
}
