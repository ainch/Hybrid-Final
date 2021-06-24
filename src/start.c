#include "xypic.h"
#include "global.h"
#include <libgen.h>
#include <string.h>
#include <stdio.h>
#include "parson.h"

FILE *InputDeck;

void InputRead(int argc, char *argv[]) {
   int i;
   JSON_Value *InputValue;
   JSON_Object *MainObject;
   JSON_Object *SubObject1,*SubObject2,*SubObject3,*SubObject4,*SubObject5;
   JSON_Object *SystemObject;
   JSON_Object *BufObject;
   JSON_Array *BCArray,*CondArray,*DielArray;
   // Total read
   InputValue = json_parse_file(InputFile);
   MainObject = json_value_get_object(InputValue);
   //
   SubObject1 = json_object_get_object(MainObject,"GPU_Device_Number");
   printf("Device # = %d\n",(int)json_object_get_number(SubObject1,"Device"));
   //
   SubObject2 = json_object_get_object(MainObject,"Geometry");
   SystemObject = json_object_get_object(SubObject2,"SystemSpec");
   printf("X_length(m) = %d\n",(float)json_object_get_number(SystemObject,"X_length(m)"));
   printf("Y_length(m) = %d\n",(float)json_object_get_number(SystemObject,"Y_length(m)"));
   printf("Z_length(m) = %d\n",(float)json_object_get_number(SystemObject,"Z_length(m)"));
   printf("NumCellx = %d\n",(int)json_object_get_number(SystemObject,"NumGridx"));
   printf("NumCelly = %d\n",(int)json_object_get_number(SystemObject,"NumGridy"));
   BCArray = json_object_get_array(SubObject2,"BoundaryCondition");
   for (i=0;i<json_array_get_count(BCArray);i++){
      BufObject = json_array_get_object(BCArray,i);
      printf("X0 = %d,",(int)json_object_get_number(BufObject,"X0"));
      printf("X1 = %d,",(int)json_object_get_number(BufObject,"X1"));
      printf("Y0 = %d,",(int)json_object_get_number(BufObject,"Y0"));
      printf("Y1 = %d,",(int)json_object_get_number(BufObject,"Y1"));
      printf("B.C = %d,",(int)json_object_get_number(BufObject,"B.C"));
      printf("Temp = %g\n",(float)json_object_get_number(BufObject,"Temp"));
   }
   CondArray = json_object_get_array(SubObject2,"ConductorSpec");
   for (i=0;i<json_array_get_count(CondArray);i++){
      BufObject = json_array_get_object(CondArray,i);
      printf("X0 = %d,",(int)json_object_get_number(BufObject,"X0"));
      printf("X1 = %d,",(int)json_object_get_number(BufObject,"X1"));
      printf("Y0 = %d,",(int)json_object_get_number(BufObject,"Y0"));
      printf("Y1 = %d,",(int)json_object_get_number(BufObject,"Y1"));
      printf("Temp = %g",(float)json_object_get_number(BufObject,"Temp"));
      printf("DC = %g,",(float)json_object_get_number(BufObject,"DC(V)"));
      printf("Power1 = %g,",(float)json_object_get_number(BufObject,"Power1(W)"));
      printf("AC1 = %g,",(float)json_object_get_number(BufObject,"AC1(V)"));
      printf("Freq1 = %g,",(float)json_object_get_number(BufObject,"Freq1(1/s)"));
      printf("Phase1 = %g,",(float)json_object_get_number(BufObject,"Phase1(deg)"));
      printf("Power2 = %g,",(float)json_object_get_number(BufObject,"Power2(W)"));
      printf("AC2 = %g,",(float)json_object_get_number(BufObject,"AC2(V)"));
      printf("Freq2 = %g,",(float)json_object_get_number(BufObject,"Freq2(1/s)"));
      printf("Phase2 = %g,",(float)json_object_get_number(BufObject,"Phase2(deg)"));
      printf("R = %g,",(float)json_object_get_number(BufObject,"R(Ohm)"));
      printf("L = %g,",(float)json_object_get_number(BufObject,"L(H)"));
      printf("C = %g\n",(float)json_object_get_number(BufObject,"C(F)"));
   }
   DielArray = json_object_get_array(SubObject2,"DielectricSpec");
   for (i=0;i<json_array_get_count(DielArray);i++){
      BufObject = json_array_get_object(DielArray,i);
      printf("X0 = %d,",(int)json_object_get_number(BufObject,"X0"));
      printf("X1 = %d,",(int)json_object_get_number(BufObject,"X1"));
      printf("Y0 = %d,",(int)json_object_get_number(BufObject,"Y0"));
      printf("Y1 = %d,",(int)json_object_get_number(BufObject,"Y1"));
      printf("Epsilon = %g\n",(float)json_object_get_number(BufObject,"Epsilon"));
   }
   //
   SubObject3 = json_object_get_object(MainObject,"GasSpecies");
   //
   SubObject4 = json_object_get_object(MainObject,"SimulationMethod");
   printf("TimeStepForPIC = %d\n",(int)json_object_get_number(SubObject4,"TimeStepForPIC"));
   printf("TimeStepForConti = %d\n",(int)json_object_get_number(SubObject4,"TimeStepForConti"));
   SubObject5 = json_object_get_object(MainObject,"Diagnostics");


}
void start() {
   fprintf(stderr,"Tstrp = %d\n",tstep);
}
void DumpRead(int argc, char *argv[]) {
   fprintf(stderr,"Tstrp = %d\n",tstep);
}