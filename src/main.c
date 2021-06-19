#include "main.h"
#include "parson.h"

void display_title();
void display_finish();

int main(int argc, char *argv[])
{
	//t = 0;
	display_title();       	/* Display 2D Field Solver title */
	//start(argc,argv);   	/* Allocate arrays and initialize */
	//main_cuda();			/* GPU Start */
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
