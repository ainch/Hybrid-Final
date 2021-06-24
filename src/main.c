#include "main.h"

void display_title();
void display_finish();

int main(int argc, char *argv[])
{
	int Post_Processing_flag = 0;
	t = 0.0;
	tstep = 0;
	cstep = 0;
	if (argc == 6) {
		Post_Processing_flag = 1;
	} else if (argc == 5) {
		DumpFlag = 1;
		strcpy(InputFile, argv[2]);
		strcpy(DumpFile, argv[4]);
	} else if (argc == 3) {
		strcpy(InputFile, argv[2]);
		DumpFlag = 0;
	} else {
		puts("\t2D XY PIC -i <inputdeck>");
		exit(1);
	}
	display_title();       	/* Display 2D Field Solver title */
	InputRead(argc,argv);	/* InputFile Parsing */
	start();   				/* Initialized */
	DumpRead(argc,argv);	/* DumpFile Read */
	if(Post_Processing_flag){

	} else{
		//main_cuda();			/* GPU Start */
	}
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
