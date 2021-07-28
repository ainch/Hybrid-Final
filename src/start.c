#include "xypic.h"
#include "global.h"
#include <libgen.h>
#include "parson.h"
#include "main.h"

FILE *InputDeck;

void InputRead(int argc, char *argv[]) {
   int i,j;
   int buf;
   JSON_Value *InputValue;
   JSON_Object *MainObject;
   JSON_Object *SubObject1,*SubObject2,*SubObject3;
   JSON_Object *SubObject4,*SubObject5,*SubObject6;
   JSON_Object *BufObject,*BufObject2;
   JSON_Array *BufArray,*CondArray,*DielArray;
   // Total read
   InputValue = json_parse_file(InputFile);
   MainObject = json_value_get_object(InputValue);
   //-------------------------//
   //----GPU_Device_Number----//
   //-------------------------//
   SubObject1 = json_object_get_object(MainObject,"GPU_Device_Number");
   device_num = (int)json_object_get_number(SubObject1,"GPUDevice");
   //-------------------------//
   //--------Geometry---------//
   //-------------------------//
   SubObject2 = json_object_get_object(MainObject,"Geometry");
   BufObject = json_object_get_object(SubObject2,"SystemSpec");
   xlength = (float)json_object_get_number(BufObject,"X_length(m)");
   ylength = (float)json_object_get_number(BufObject,"Y_length(m)");
   zlength = (float)json_object_get_number(BufObject,"Z_length(m)");
   ngx = (int)json_object_get_number(BufObject,"NumGridx");
   ngy = (int)json_object_get_number(BufObject,"NumGridy");
   if (ngx < 5 || ngy < 5){
      printf("\"ngx\",\"ngy\" is too small!\n");
      exit(1);
   }
   BufArray = json_object_get_array(SubObject2,"BoundaryCondition");
   BoundaryNUN = (int)json_array_get_count(BufArray);
   if (BoundaryNUN < 4){
      printf("\"BoundaryCondition\" is too small!(BoundaryNum>=4)\n");
      exit(1);
   }
   BoundaryX0 = VIMalloc(BoundaryNUN);
   BoundaryY0 = VIMalloc(BoundaryNUN);
   BoundaryX1 = VIMalloc(BoundaryNUN);
   BoundaryY1 = VIMalloc(BoundaryNUN);
   BoundaryBC = VIMalloc(BoundaryNUN);
   BoundaryTEMP = VFMalloc(BoundaryNUN);
   for (i=0;i<BoundaryNUN;i++){
      BufObject = json_array_get_object(BufArray,i);
      BoundaryX0[i]=(int)json_object_get_number(BufObject,"X0");
      BoundaryX1[i]=(int)json_object_get_number(BufObject,"X1");
      if(BoundaryX0[i] > BoundaryX1[i]){
         printf("Line %d X position is error in \"BoundaryCondition\"(\"X0\" <= \"X1\")\n",i+1);
         exit(1);
      } 
      BoundaryY0[i]=(int)json_object_get_number(BufObject,"Y0");
      BoundaryY1[i]=(int)json_object_get_number(BufObject,"Y1");
      if(BoundaryY0[i] > BoundaryY1[i]){
         printf("Line %d Y position is error in \"BoundaryCondition\"(\"Y0\" <= \"Y1\")\n",i+1);
         exit(1);
      } 
      BoundaryBC[i]=(int)json_object_get_number(BufObject,"B.C");
      if(BoundaryBC[i] > 4){
         printf("Line %d B.C is error in \"BoundaryCondition\"(\"B.C\" < 5)\n",i+1);
         exit(1);
      } 
      BoundaryTEMP[i]=(float)json_object_get_number(BufObject,"Temp");
   }
   BufArray = json_object_get_array(SubObject2,"ConductorSpec");
   CondNUN = (int)json_array_get_count(BufArray);
   if (CondNUN < 2){
      printf("\"Conductor\" is too small!(CondNum>=2)\n");
      exit(1);
   }
   CondM_ID = VIMalloc(CondNUN);
   CondX0 = VIMalloc(CondNUN);
   CondX1 = VIMalloc(CondNUN);
   CondY0 = VIMalloc(CondNUN);
   CondY1 = VIMalloc(CondNUN);
   CondTEMP = VFMalloc(CondNUN);
   for (i=0;i<CondNUN;i++){
      BufObject = json_array_get_object(BufArray,i);
      CondM_ID[i] = (int)json_object_get_number(BufObject,"M_ID");
      CondX0[i] = (int)json_object_get_number(BufObject,"X0");
      CondX1[i] = (int)json_object_get_number(BufObject,"X1");
      if(CondX0[i] > CondX1[i]){
         printf("Line %d X position is error in \"ConductorSpec\"(\"X0\" <= \"X1\")\n",i+1);
         exit(1);
      } 
      CondY0[i] = (int)json_object_get_number(BufObject,"Y0");
      CondY1[i] = (int)json_object_get_number(BufObject,"Y1");
      if(CondY0[i] > CondY1[i]){
         printf("Line %d Y position is error in \"ConductorSpec\"(\"Y0\" <= \"Y1\")\n",i+1);
         exit(1);
      } 
      CondTEMP[i] = (float)json_object_get_number(BufObject,"Temp");
   }
   BufArray = json_object_get_array(SubObject2,"Source");
   SrcNUN = (int)json_array_get_count(BufArray);
   SrcM_ID = VIMalloc(SrcNUN);
   SrcDC = VFMalloc(SrcNUN);
   SrcPOWER = VFMalloc(SrcNUN);
   SrcAC = VFMalloc(SrcNUN);
   SrcFREQ = VFMalloc(SrcNUN);
   SrcPHASE = VFMalloc(SrcNUN);
   SrcR = VFMalloc(SrcNUN);
   SrcL = VFMalloc(SrcNUN);
   SrcC = VFMalloc(SrcNUN);
   Min_FREQ  = 1e200;
   for (i=0;i<SrcNUN;i++){
      buf = 0;
      BufObject = json_array_get_object(BufArray,i);
      SrcM_ID[i] = (int)json_object_get_number(BufObject,"M_ID");
      for (j=0;j<CondNUN;j++){
         if(CondM_ID[j]==SrcM_ID[i])
            buf = 1;
      }
      if(buf == 0){
         printf("Line %d Source is no matching Conductor\n",i+1);
         exit(1);
      }
      SrcDC[i] = (float)json_object_get_number(BufObject,"DC(V)");
      SrcPOWER[i] = (float)json_object_get_number(BufObject,"Power(W)");
      SrcAC[i] = (float)json_object_get_number(BufObject,"AC(V)");
      SrcFREQ[i] = (float)json_object_get_number(BufObject,"Freq(1/s)");
      SrcPHASE[i] = (float)json_object_get_number(BufObject,"Phase(deg)");
      SrcR[i] = (float)json_object_get_number(BufObject,"R(Ohm)");
      SrcL[i] = (float)json_object_get_number(BufObject,"L(H)");
      SrcC[i] = (float)json_object_get_number(BufObject,"C(F)");
      if(SrcFREQ[i] != 0){
         if(SrcFREQ[i]<Min_FREQ)
            Min_FREQ  = SrcFREQ[i];
      }
   }
   BufArray = json_object_get_array(SubObject2,"DielectricSpec");
   DielNUN = (int)json_array_get_count(BufArray);
   DielM_ID = VIMalloc(DielNUN);
   DielX0 = VIMalloc(DielNUN);
   DielX1 = VIMalloc(DielNUN);
   DielY0 = VIMalloc(DielNUN);
   DielY1 = VIMalloc(DielNUN);
   DielEPS = VFMalloc(DielNUN);
   for (i=0;i<DielNUN;i++){
      BufObject = json_array_get_object(BufArray,i);
      DielM_ID[i] = (int)json_object_get_number(BufObject,"M_ID");
      for (j=0;j<CondNUN;j++){
         if(DielM_ID[i]==CondM_ID[j]){
            printf("Line %d M_ID is duplicated with conductor %d in \"DielectricSpec\"\n",i+1,CondM_ID[j]);
            exit(1);
         }
      }
      DielX0[i] = (int)json_object_get_number(BufObject,"X0");
      DielX1[i] = (int)json_object_get_number(BufObject,"X1");
      if(DielX0[i] > DielX1[i]){
         printf("Line %d X position is error in \"DielectricSpec\"(\"X0\" <= \"X1\")\n",i+1);
         exit(1);
      } 
      DielY0[i] = (int)json_object_get_number(BufObject,"Y0");
      DielY1[i] = (int)json_object_get_number(BufObject,"Y1");
      if(DielY0[i] > DielY1[i]){
         printf("Line %d Y position is error in \"DielelctricSpec\"(\"Y0\" <= \"Y1\")\n",i+1);
         exit(1);
      } 
      DielEPS[i] = (float)json_object_get_number(BufObject,"Epsilon");
   }
   exit(1);
   //-------------------------//
   //-------GasSpecies--------//
   //-------------------------//
   SubObject3 = json_object_get_object(MainObject,"GasSpecies");
   MainGas = (int)json_object_get_number(SubObject3,"Type(0:Ar,1:O2,2:Ar/O2)");
   //printf("MainGas = %d,\n",(int)json_object_get_number(SubObject3,"Type(0:Ar,1:O2,2:Ar/O2)"));   
   printf("TotalPres = %g,\n",(float)json_object_get_number(SubObject3,"TotalPres(Torr)"));
   if(MainGas==ARGON){
      BufObject = json_object_get_object(SubObject3,"Background");
      BufObject2 = json_object_get_object(BufObject,"Argon");
      printf("Ratio = %g,\n",(float)json_object_get_number(BufObject2,"Ratio"));
      printf("Temp = %g,\n",(float)json_object_get_number(BufObject2,"Temp(eV)"));
      BufObject = json_object_get_object(SubObject3,"NeutralSpecies");
      BufObject2 = json_object_get_object(BufObject,"Ar*");
      printf("Density = %g,\n",(float)json_object_get_number(BufObject2,"Density"));
      printf("Temp = %g,\n",(float)json_object_get_number(BufObject2,"Temp(eV)"));
      printf("Neutral Loadtype = %d,\n",(int)json_object_get_number(BufObject,"Loadtype"));
      BufObject2 = json_object_get_object(BufObject,"LoadPosition(m)");
      printf("X0 = %g,",(float)json_object_get_number(BufObject2,"X0"));
      printf("X1 = %g,",(float)json_object_get_number(BufObject2,"X1"));
      printf("Y0 = %g,",(float)json_object_get_number(BufObject2,"Y0"));
      printf("Y1 = %g,\n",(float)json_object_get_number(BufObject2,"Y1"));
      BufObject = json_object_get_object(SubObject3,"ChargeSpecies");
      BufObject2 = json_object_get_object(BufObject,"Electron");
      printf("S_ID = %d,\n",(int)json_object_get_number(BufObject2,"S_ID"));
      printf("Density = %g,\n",(float)json_object_get_number(BufObject2,"Density"));
      printf("Temp = %g,\n",(float)json_object_get_number(BufObject2,"Temp(eV)"));
      printf("np2c = %g,\n",(float)json_object_get_number(BufObject2,"np2c"));
      printf("Max_np = %d,\n",(int)json_object_get_number(BufObject2,"Max_np"));
      BufObject2 = json_object_get_object(BufObject,"Ar+");
      printf("S_ID = %d,\n",(int)json_object_get_number(BufObject2,"S_ID"));
      printf("Density = %g,\n",(float)json_object_get_number(BufObject2,"Density"));
      printf("Temp = %g,\n",(float)json_object_get_number(BufObject2,"Temp(eV)"));
      printf("np2c = %g,\n",(float)json_object_get_number(BufObject2,"np2c"));
      printf("Max_np = %d,\n",(int)json_object_get_number(BufObject2,"Max_np"));
      printf("Charged Loadtype = %d,\n",(int)json_object_get_number(BufObject,"Loadtype"));
      BufObject2 = json_object_get_object(BufObject,"LoadPosition(m)");
      printf("X0 = %g,",(float)json_object_get_number(BufObject2,"X0"));
      printf("X1 = %g,",(float)json_object_get_number(BufObject2,"X1"));
      printf("Y0 = %g,",(float)json_object_get_number(BufObject2,"Y0"));
      printf("Y1 = %g,\n",(float)json_object_get_number(BufObject2,"Y1"));
   }else if(MainGas==OXYGEN){
      BufObject = json_object_get_object(SubObject3,"Background");
      BufObject2 = json_object_get_object(BufObject,"Oxygen");
      printf("Ratio = %g,\n",(float)json_object_get_number(BufObject2,"Ratio"));
      printf("Temp = %g,\n",(float)json_object_get_number(BufObject2,"Temp(eV)"));
      BufObject = json_object_get_object(SubObject3,"NeutralSpecies");
      BufObject2 = json_object_get_object(BufObject,"OP");
      printf("Density = %g,\n",(float)json_object_get_number(BufObject2,"Density"));
      printf("Temp = %g,\n",(float)json_object_get_number(BufObject2,"Temp(eV)"));
      BufObject2 = json_object_get_object(BufObject,"OD");
      printf("Density = %g,\n",(float)json_object_get_number(BufObject2,"Density"));
      printf("Temp = %g,\n",(float)json_object_get_number(BufObject2,"Temp(eV)"));
      BufObject2 = json_object_get_object(BufObject,"O2A");
      printf("Density = %g,\n",(float)json_object_get_number(BufObject2,"Density"));
      printf("Temp = %g,\n",(float)json_object_get_number(BufObject2,"Temp(eV)"));
      BufObject2 = json_object_get_object(BufObject,"O2B");
      printf("Density = %g,\n",(float)json_object_get_number(BufObject2,"Density"));
      printf("Temp = %g,\n",(float)json_object_get_number(BufObject2,"Temp(eV)"));
      printf("Neutral Loadtype = %d,\n",(int)json_object_get_number(BufObject,"Loadtype"));
      BufObject2 = json_object_get_object(BufObject,"LoadPosition(m)");
      printf("X0 = %g,",(float)json_object_get_number(BufObject2,"X0"));
      printf("X1 = %g,",(float)json_object_get_number(BufObject2,"X1"));
      printf("Y0 = %g,",(float)json_object_get_number(BufObject2,"Y0"));
      printf("Y1 = %g,\n",(float)json_object_get_number(BufObject2,"Y1"));
      BufObject = json_object_get_object(SubObject3,"ChargeSpecies");
      BufObject2 = json_object_get_object(BufObject,"Electron");
      printf("S_ID = %d,\n",(int)json_object_get_number(BufObject2,"S_ID"));
      printf("Density = %g,\n",(float)json_object_get_number(BufObject2,"Density"));
      printf("Temp = %g,\n",(float)json_object_get_number(BufObject2,"Temp(eV)"));
      printf("np2c = %g,\n",(float)json_object_get_number(BufObject2,"np2c"));
      printf("Max_np = %d,\n",(int)json_object_get_number(BufObject2,"Max_np"));
      BufObject2 = json_object_get_object(BufObject,"O2+");
      printf("S_ID = %d,\n",(int)json_object_get_number(BufObject2,"S_ID"));
      printf("Density = %g,\n",(float)json_object_get_number(BufObject2,"Density"));
      printf("Temp = %g,\n",(float)json_object_get_number(BufObject2,"Temp(eV)"));
      printf("np2c = %g,\n",(float)json_object_get_number(BufObject2,"np2c"));
      printf("Max_np = %d,\n",(int)json_object_get_number(BufObject2,"Max_np"));
      BufObject2 = json_object_get_object(BufObject,"O+");
      printf("S_ID = %d,\n",(int)json_object_get_number(BufObject2,"S_ID"));
      printf("Density = %g,\n",(float)json_object_get_number(BufObject2,"Density"));
      printf("Temp = %g,\n",(float)json_object_get_number(BufObject2,"Temp(eV)"));
      printf("np2c = %g,\n",(float)json_object_get_number(BufObject2,"np2c"));
      printf("Max_np = %d,\n",(int)json_object_get_number(BufObject2,"Max_np"));
      BufObject2 = json_object_get_object(BufObject,"O-");
      printf("S_ID = %d,\n",(int)json_object_get_number(BufObject2,"S_ID"));
      printf("Density = %g,\n",(float)json_object_get_number(BufObject2,"Density"));
      printf("Temp = %g,\n",(float)json_object_get_number(BufObject2,"Temp(eV)"));
      printf("np2c = %g,\n",(float)json_object_get_number(BufObject2,"np2c"));
      printf("Max_np = %d,\n",(int)json_object_get_number(BufObject2,"Max_np"));
      printf("Charged Loadtype = %d,\n",(int)json_object_get_number(BufObject,"Loadtype"));
      BufObject2 = json_object_get_object(BufObject,"LoadPosition(m)");
      printf("X0 = %g,",(float)json_object_get_number(BufObject2,"X0"));
      printf("X1 = %g,",(float)json_object_get_number(BufObject2,"X1"));
      printf("Y0 = %g,",(float)json_object_get_number(BufObject2,"Y0"));
      printf("Y1 = %g,\n",(float)json_object_get_number(BufObject2,"Y1"));
   }else if(MainGas==ARO2){
      BufObject = json_object_get_object(SubObject3,"Background");
      BufObject2 = json_object_get_object(BufObject,"Oxygen");
      printf("Ratio = %g,\n",(float)json_object_get_number(BufObject2,"Ratio"));
      printf("Temp = %g,\n",(float)json_object_get_number(BufObject2,"Temp(eV)"));
      BufObject2 = json_object_get_object(BufObject,"Argon");
      printf("Ratio = %g,\n",(float)json_object_get_number(BufObject2,"Ratio"));
      printf("Temp = %g,\n",(float)json_object_get_number(BufObject2,"Temp(eV)"));
      BufObject = json_object_get_object(SubObject3,"NeutralSpecies");
      BufObject2 = json_object_get_object(BufObject,"Ar*");
      printf("Density = %g,\n",(float)json_object_get_number(BufObject2,"Density"));
      printf("Temp = %g,\n",(float)json_object_get_number(BufObject2,"Temp(eV)"));
      BufObject2 = json_object_get_object(BufObject,"OP");
      printf("Density = %g,\n",(float)json_object_get_number(BufObject2,"Density"));
      printf("Temp = %g,\n",(float)json_object_get_number(BufObject2,"Temp(eV)"));
      BufObject2 = json_object_get_object(BufObject,"OD");
      printf("Density = %g,\n",(float)json_object_get_number(BufObject2,"Density"));
      printf("Temp = %g,\n",(float)json_object_get_number(BufObject2,"Temp(eV)"));
      BufObject2 = json_object_get_object(BufObject,"O2A");
      printf("Density = %g,\n",(float)json_object_get_number(BufObject2,"Density"));
      printf("Temp = %g,\n",(float)json_object_get_number(BufObject2,"Temp(eV)"));
      BufObject2 = json_object_get_object(BufObject,"O2B");
      printf("Density = %g,\n",(float)json_object_get_number(BufObject2,"Density"));
      printf("Temp = %g,\n",(float)json_object_get_number(BufObject2,"Temp(eV)"));
      printf("Neutral Loadtype = %d,\n",(int)json_object_get_number(BufObject,"Loadtype"));
      BufObject2 = json_object_get_object(BufObject,"LoadPosition(m)");
      printf("X0 = %g,",(float)json_object_get_number(BufObject2,"X0"));
      printf("X1 = %g,",(float)json_object_get_number(BufObject2,"X1"));
      printf("Y0 = %g,",(float)json_object_get_number(BufObject2,"Y0"));
      printf("Y1 = %g,\n",(float)json_object_get_number(BufObject2,"Y1"));
      BufObject = json_object_get_object(SubObject3,"ChargeSpecies");
      BufObject2 = json_object_get_object(BufObject,"Electron");
      printf("S_ID = %d,\n",(int)json_object_get_number(BufObject2,"S_ID"));
      printf("Density = %g,\n",(float)json_object_get_number(BufObject2,"Density"));
      printf("Temp = %g,\n",(float)json_object_get_number(BufObject2,"Temp(eV)"));
      printf("np2c = %g,\n",(float)json_object_get_number(BufObject2,"np2c"));
      printf("Max_np = %d,\n",(int)json_object_get_number(BufObject2,"Max_np"));
      BufObject2 = json_object_get_object(BufObject,"Ar+");
      printf("S_ID = %d,\n",(int)json_object_get_number(BufObject2,"S_ID"));
      printf("Density = %g,\n",(float)json_object_get_number(BufObject2,"Density"));
      printf("Temp = %g,\n",(float)json_object_get_number(BufObject2,"Temp(eV)"));
      printf("np2c = %g,\n",(float)json_object_get_number(BufObject2,"np2c"));
      printf("Max_np = %d,\n",(int)json_object_get_number(BufObject2,"Max_np"));
      BufObject2 = json_object_get_object(BufObject,"O2+");
      printf("S_ID = %d,\n",(int)json_object_get_number(BufObject2,"S_ID"));
      printf("Density = %g,\n",(float)json_object_get_number(BufObject2,"Density"));
      printf("Temp = %g,\n",(float)json_object_get_number(BufObject2,"Temp(eV)"));
      printf("np2c = %g,\n",(float)json_object_get_number(BufObject2,"np2c"));
      printf("Max_np = %d,\n",(int)json_object_get_number(BufObject2,"Max_np"));
      BufObject2 = json_object_get_object(BufObject,"O+");
      printf("S_ID = %d,\n",(int)json_object_get_number(BufObject2,"S_ID"));
      printf("Density = %g,\n",(float)json_object_get_number(BufObject2,"Density"));
      printf("Temp = %g,\n",(float)json_object_get_number(BufObject2,"Temp(eV)"));
      printf("np2c = %g,\n",(float)json_object_get_number(BufObject2,"np2c"));
      printf("Max_np = %d,\n",(int)json_object_get_number(BufObject2,"Max_np"));
      BufObject2 = json_object_get_object(BufObject,"O-");
      printf("S_ID = %d,\n",(int)json_object_get_number(BufObject2,"S_ID"));
      printf("Density = %g,\n",(float)json_object_get_number(BufObject2,"Density"));
      printf("Temp = %g,\n",(float)json_object_get_number(BufObject2,"Temp(eV)"));
      printf("np2c = %g,\n",(float)json_object_get_number(BufObject2,"np2c"));
      printf("Max_np = %d,\n",(int)json_object_get_number(BufObject2,"Max_np"));
      printf("Charged Loadtype = %d,\n",(int)json_object_get_number(BufObject,"Loadtype"));
      BufObject2 = json_object_get_object(BufObject,"LoadPosition(m)");
      printf("X0 = %g,",(float)json_object_get_number(BufObject2,"X0"));
      printf("X1 = %g,",(float)json_object_get_number(BufObject2,"X1"));
      printf("Y0 = %g,",(float)json_object_get_number(BufObject2,"Y0"));
      printf("Y1 = %g,\n",(float)json_object_get_number(BufObject2,"Y1"));
   }else{
      exit(1);
   }
   //-------------------------//
   //----SimulationMethod-----//
   //-------------------------//
   SubObject4 = json_object_get_object(MainObject,"SimulationMethod");
   printf("TimeStepForPIC = %d\n",(int)json_object_get_number(SubObject4,"TimeStepForPIC"));
   printf("TimeStepForConti = %d\n",(int)json_object_get_number(SubObject4,"TimeStepForConti"));
   printf("PCGMarginOfError = %g\n",(float)json_object_get_number(SubObject4,"PCGMarginOfError"));
   printf("HistoryMax = %d\n",(int)json_object_get_number(SubObject4,"HistoryMax"));
   printf("HistoryDivide = %d\n",(int)json_object_get_number(SubObject4,"HistoryDivide"));
   printf("AverageIter = %d\n",(int)json_object_get_number(SubObject4,"AverageIter"));
   printf("ParticleLimit = %d\n",(int)json_object_get_number(SubObject4,"ParticleLimit"));
   printf("B_flag = %s\n",(char*)json_object_get_string(SubObject4,"B_flag"));
   BufObject = json_object_get_object(SubObject4,"Initialing_Int");
   printf("Initialing_Init_flag = %d\n",(int)json_object_get_number(BufObject,"Flag"));
   printf("Time_order = %d\n",(int)json_object_get_number(BufObject,"Time_order"));
   printf("Dump_order = %d\n",(int)json_object_get_number(BufObject,"Dump_order"));
   BufObject = json_object_get_object(SubObject4,"Addtion_IZ");
   printf("Add ionization_flag = %d\n",(int)json_object_get_number(BufObject,"Flag"));
   BufArray = json_object_get_array(BufObject,"AdditionIonization");
   for (i=0;i<json_array_get_count(BufArray);i++){
      BufObject2 = json_array_get_object(BufArray,i);
      printf("LoadType = %d,\n",(int)json_object_get_number(BufObject2,"LoadType"));
      printf("X0 = %g,",(float)json_object_get_number(BufObject2,"X0"));
      printf("X1 = %g,",(float)json_object_get_number(BufObject2,"X1"));
      printf("Y0 = %g,",(float)json_object_get_number(BufObject2,"Y0"));
      printf("Y1 = %g,",(float)json_object_get_number(BufObject2,"Y1"));
      printf("IonizationRate = %g\n",(float)json_object_get_number(BufObject2,"IonizationRate"));
   }
   BufObject = json_object_get_object(SubObject4,"ContinuityBreakPoint");
   printf("Margin(%) = %g\n",(float)json_object_get_number(BufObject,"Margin(%)"));
   BufObject = json_object_get_object(SubObject4,"PowerDrivenOption");
   printf("VoltageUpRatio = %g\n",(float)json_object_get_number(BufObject,"VoltageUpRatio"));
   BufObject = json_object_get_object(SubObject4,"SimulationStop");
   printf("EndTime(s) = %g\n",(float)json_object_get_number(BufObject,"EndTime(s)"));
   BufObject2 = json_object_get_object(BufObject,"AutoSaturation");
   printf("Flag = %d\n",(int)json_object_get_number(BufObject2,"Flag"));
   printf("AveMargin(%) = %g\n",(float)json_object_get_number(BufObject2,"AveMargin(%)"));
   printf("SatGoal = %d\n",(int)json_object_get_number(BufObject2,"SatGoal"));
   //-------------------------//
   //-------Diagnostics-------//
   //-------------------------//
   SubObject5 = json_object_get_object(MainObject,"Diagnostics");
   printf("BasicDiagnostic = %d\n",(int)json_object_get_number(SubObject5,"BasicDiagnostic"));
   BufObject = json_object_get_object(SubObject5,"ElectronRefection");
   printf("BD_Flag = %d\n",(int)json_object_get_number(BufObject,"Flag"));
   BufArray = json_object_get_array(BufObject,"Material_ID");
   for (i=0;i<json_array_get_count(BufArray);i++){
      printf("Material_ID = %d,\n",(int)json_array_get_number(BufArray,i));
   }
   BufArray = json_object_get_array(BufObject,"Coefficient");
   for (i=0;i<json_array_get_count(BufArray);i++){
      printf("Coefficient = %g,\n",(float)json_array_get_number(BufArray,i));
   }
   BufObject = json_object_get_object(SubObject5,"SecondaryElectronEmission");
   printf("SEE_Flag = %d\n",(int)json_object_get_number(BufObject,"Flag"));
   BufArray = json_object_get_array(BufObject,"Material_ID");
   for (i=0;i<json_array_get_count(BufArray);i++){
      printf("Material_ID = %d,\n",(int)json_array_get_number(BufArray,i));
   }
   BufArray = json_object_get_array(BufObject,"Coefficient");
   for (i=0;i<json_array_get_count(BufArray);i++){
      printf("Coefficient = %g,\n",(float)json_array_get_number(BufArray,i));
   }
   BufObject = json_object_get_object(SubObject5,"EnergyDistributionFunction");
   printf("SEE_Flag = %d\n",(int)json_object_get_number(BufObject,"Flag"));
   BufArray = json_object_get_array(BufObject,"RangeOfEnergyforElectron");
   for (i=0;i<json_array_get_count(BufArray);i++){
      printf("position = %g,\n",(float)json_array_get_number(BufArray,i));
   }
   BufArray = json_object_get_array(BufObject,"RegionforElectron");
   for (i=0;i<json_array_get_count(BufArray);i++){
      BufObject2 = json_array_get_object(BufArray,i);
      printf("X0 = %d,",(int)json_object_get_number(BufObject2,"X0"));
      printf("X1 = %d,",(int)json_object_get_number(BufObject2,"X1"));
      printf("Y0 = %d,",(int)json_object_get_number(BufObject2,"Y0"));
      printf("Y1 = %d,",(int)json_object_get_number(BufObject2,"Y1"));
   }
   BufArray = json_object_get_array(BufObject,"RangeOfEnergyforIon");
   for (i=0;i<json_array_get_count(BufArray);i++){
      printf("RangeOfEnergyforIon = %g,\n",(float)json_array_get_number(BufArray,i));
   }
   BufArray = json_object_get_array(BufObject,"RegionforIon");
   for (i=0;i<json_array_get_count(BufArray);i++){
      BufObject2 = json_array_get_object(BufArray,i);
      printf("S_ID = %d,",(int)json_object_get_number(BufObject2,"S_ID"));
      printf("X0 = %d,",(int)json_object_get_number(BufObject2,"X0"));
      printf("X1 = %d,",(int)json_object_get_number(BufObject2,"X1"));
      printf("Y0 = %d,",(int)json_object_get_number(BufObject2,"Y0"));
      printf("Y1 = %d,",(int)json_object_get_number(BufObject2,"Y1"));
   }
   printf("Tecplot_Save = %d\n",(int)json_object_get_number(BufObject,"Tecplot_Save"));
   BufObject2 = json_object_get_object(BufObject,"1CycleMovie");
   printf("CMovieFlag = %d\n",(int)json_object_get_number(BufObject2,"Flag"));
   printf("PhaseResolveNum = %d\n",(int)json_object_get_number(BufObject2,"PhaseResolveNum"));
   printf("AccumulationNUM = %d\n",(int)json_object_get_number(BufObject2,"AccumulationNUM"));
   BufObject = json_object_get_object(SubObject5,"IonEnergyAngularDistribution");
   BufArray = json_object_get_array(BufObject,"RangeOfEnergy");
   for (i=0;i<json_array_get_count(BufArray);i++){
      printf("RangeOfEnergy = %g,\n",(float)json_array_get_number(BufArray,i));
   }
   BufArray = json_object_get_array(BufObject,"RangeOfAngle");
   for (i=0;i<json_array_get_count(BufArray);i++){
      printf("RangeOfAngle = %g,\n",(float)json_array_get_number(BufArray,i));
   }
   BufArray = json_object_get_array(BufObject,"Position");
   for (i=0;i<json_array_get_count(BufArray);i++){
      BufObject2 = json_array_get_object(BufArray,i);
      printf("M_ID = %d,",(int)json_object_get_number(BufObject2,"M_ID"));
      printf("X0 = %d,",(int)json_object_get_number(BufObject2,"X0"));
      printf("X1 = %d,",(int)json_object_get_number(BufObject2,"X1"));
      printf("Y0 = %d,",(int)json_object_get_number(BufObject2,"Y0"));
      printf("Y1 = %d,",(int)json_object_get_number(BufObject2,"Y1"));
   }
   printf("CMovieFlag = %d\n",(int)json_object_get_number(BufObject,"Tecplot_Save"));
   BufObject2 = json_object_get_object(BufObject,"1CycleMovie");
   printf("CMovieFlag = %d\n",(int)json_object_get_number(BufObject2,"Flag"));
   printf("PhaseResolveNum = %d\n",(int)json_object_get_number(BufObject2,"PhaseResolveNum"));
   printf("AccumulationNUM = %d\n",(int)json_object_get_number(BufObject2,"AccumulationNUM"));
   BufObject = json_object_get_object(SubObject5,"Test_Particle");
   printf("Flag = %d\n",(int)json_object_get_number(BufObject,"Flag"));
   printf("MCC_Flag = %d\n",(int)json_object_get_number(BufObject,"MCC_Flag"));
   printf("NumberOfParticle = %d\n",(int)json_object_get_number(BufObject,"NumberOfParticle"));
   printf("Init_Position = %d\n",(int)json_object_get_number(BufObject,"Init_Position"));
   printf("Init_Velocity = %d\n",(int)json_object_get_number(BufObject,"Init_Velocity"));
   printf("Tecplot_Save = %d\n",(int)json_object_get_number(BufObject,"Tecplot_Save")); 
   BufObject = json_object_get_object(SubObject5,"TecplotSave");
   printf("Tecplot2D = %d\n",(int)json_object_get_number(BufObject,"Tecplot2D"));
   printf("Tec_Movie = %d\n",(int)json_object_get_number(BufObject,"Tec_Movie"));
   printf("Tec_Movie_FrameNum = %d\n",(int)json_object_get_number(BufObject,"Tec_Movie_FrameNum"));
   printf("Tec_Ave_Movie = %d\n",(int)json_object_get_number(BufObject,"Tec_Ave_Movie"));
   printf("Tec_Ave_Movie_Interval = %d\n",(int)json_object_get_number(BufObject,"Tec_Ave_Movie_Interval"));
   printf("Tec_Ave_Movie_num = %d\n",(int)json_object_get_number(BufObject,"Tec_Ave_Movie_num")); 
   BufObject = json_object_get_object(SubObject5,"DumpFileSave");
   BufArray = json_object_get_array(BufObject,"Cycle");
   for (i=0;i<json_array_get_count(BufArray);i++){
      BufObject2 = json_array_get_object(BufArray,i);
      printf("Start = %g,",(float)json_object_get_number(BufObject2,"Start"));
      printf("End = %g,",(float)json_object_get_number(BufObject2,"End"));
      printf("Term = %g\n",(float)json_object_get_number(BufObject2,"Term"));
   }
   //-------------------------//
   //-----CollisionSelect-----//
   //-------------------------//
   SubObject6 = json_object_get_object(MainObject,"CollisionSelect");
   if(MainGas==0){
BufObject = json_object_get_object(SubObject6,"ArgonCase");
   printf("R0=%d,",(int)json_object_get_number(BufObject,"0.e+Ar>e+Ar"));
   printf("R1=%d,",(int)json_object_get_number(BufObject,"1.e+Ar>e+Ar*"));
   printf("R2=%d,",(int)json_object_get_number(BufObject,"2.e+Ar>e+Ar*m"));
   printf("R3=%d,",(int)json_object_get_number(BufObject,"3.e+Ar>2e+Ar^"));
   printf("R4=%d,",(int)json_object_get_number(BufObject,"4.e+Ar*m>2e+Ar^"));
   printf("R5=%d,",(int)json_object_get_number(BufObject,"5.Ar+Ar^>Ar^+Ar"));
   printf("R6=%d,",(int)json_object_get_number(BufObject,"6.Ar+Ar^>Ar+Ar^"));
   printf("R7=%d,",(int)json_object_get_number(BufObject,"7.e+Ar*m>e+Ar"));
   printf("R8=%d,",(int)json_object_get_number(BufObject,"8.Ar*m+Ar>e+Ar^+Ar"));
   printf("R9=%d,",(int)json_object_get_number(BufObject,"9.Ar+Ar*m>Ar+Ar"));
   printf("\n");
   }else if(MainGas==1){
        BufObject = json_object_get_object(SubObject6,"OxygenCase");
printf("R0=%d,",(int)json_object_get_number(BufObject,"0.e+O2>e+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"1.e+O2>e+O2*"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"2.e+O2>e+O2*"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"3.e+O2>e+O2A"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"4.e+O2>e+O2B"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"5.e+O2>e+O2*"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"6.e+O2>OP+O-"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"7.e+O2>e+2OP"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"8.e+O2>e+OP+OD"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"9.e+O2>e+2OD"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"10.e+O2>2e+O2^"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"11.e+O2>e+OP+O*"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"12.e+O2>e+O^+O-"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"13.e+O2>2e+O^+OP"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"14.e+O2A>2e+O2+"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"15.e+O2A>OP+O-"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"16.e+O2A>e+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"17.e+O2A>e+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"18.e+O2A>e+2OP"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"19.e+O2A>e+OP+OD"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"20.e+O2A>e+2OD"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"21.e+O2A>2e+O^+OP"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"22.e+O2B>2e+O2^"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"23.e+O2B>OP+O-"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"24.e+O2B>e+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"25.e+O2B>e+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"26.e+O2B>e+2O"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"27.e+O2B>e+OP+OD"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"28.e+O2B>e+2OD"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"29.e+O2B>2e+O^+OP"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"30.e+O->2e+OP"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"31.e+O2^>OP+OD"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"32.e+OP>e+OP"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"33.e+OP>e+OD"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"34.e+OP>e+O*"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"35.e+OP>e+O*"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"36.e+OP>e+O*"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"37.e+OP>2e+O^"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"38.e+OP>e+O*"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"39.e+OD>2e+O^"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"40.e+OD>e+OP"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"41.O-+O2>O-+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"42.O-+O2>e+OP+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"43.O-+OP>e+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"44.O-+O2^>OP+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"45.O-+O^>2OP"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"46.O-+O2A>e+OP+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"47.O2^+OP>O2+O^"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"48.O2^+O2>O2+O2^"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"49.O2^+O2>O2^+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"50.O2^+O2>O^+OP+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"51.O2^+O2A>O2+O2^"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"52.O2^+O2B>O2+O2^"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"53.O^+O2>OP+O2^"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"54.O^+O2>O^+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"55.O^+OP>OP+O^"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"56.O^+O2A>O2^+OP"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"57.O^+O2B>O2^+OP"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"58.O-+O2B>e+OP+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"59.OP+OD>2OP"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"60.OD+O2>OP+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"61.OD+O2>OP+O2A"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"62.OD+O2>OP+O2B"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"63.O2A+OP>OP+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"64.O2A+O2>2O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"65.O2A+O2A>2O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"66.O2B+O2>2O2"));
 printf("\n");
   }else if(MainGas==2){
         BufObject = json_object_get_object(SubObject6,"Argon/OxygenCase");
printf("R0=%d,",(int)json_object_get_number(BufObject,"00.e+Ar>e+Ar"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"01.e+Ar>e+Ar*"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"02.e+Ar>e+Ar*m"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"03.e+Ar>2e+Ar+"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"04.e+Ar*m>2e+Ar^"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"05.e+O2>e+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"06.e+O2>e+O2*"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"07.e+O2>e+O2*"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"08.e+O2>e+O2A"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"09.e+O2>e+O2B"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"10.e+O2>e+O2*"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"11.e+O2>OP+O-"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"12.e+O2>e+2OP"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"13.e+O2>e+OP+OD"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"14.e+O2>e+2OD"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"15.e+O2>2e+O2+"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"16.e+O2>e+OP+O*"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"17.e+O2>e+O++O-"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"18.e+O2>2e+O^+OP"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"19.e+O2A>2e+O2^"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"20.e+O2A>OP+O-"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"21.e+O2A>e+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"22.e+O2A>e+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"23.e+O2A>e+2O")); 
printf("R0=%d,",(int)json_object_get_number(BufObject,"24.e+O2A>e+OP+OD"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"25.e+O2A>e+2OD"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"26.e+O2A>2e+O^+OP"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"27.e+O2B>2e+O2^"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"28.e+O2B>OP+O-"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"29.e+O2B>e+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"30.e+O2B>e+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"31.e+O2B>e+2O"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"32.e+O2B>e+OP+OD"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"33.e+O2B>e+2OD"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"34.e+O2B>2e+O++OP"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"35.e+O->2e+OP"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"36.e+O2+>OP+OD "));
printf("R0=%d,",(int)json_object_get_number(BufObject,"37.e+OP>e+OP"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"38.e+OP>e+OD"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"39.e+OP>e+O*"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"40.e+OP>e+O*"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"41.e+OP>e+O*"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"42.e+OP>2e+O^"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"43.e+OP>e+O*"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"44.e+OD>2e+O+"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"45.e+OD>e+O"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"46.O-+O2>O-+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"47.O-+O2>e+OP+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"48.O-+OP>e+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"49.O-+O2^>OP+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"50.O-+O^>2OP"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"51.O-+O2A>e+OP+O2"));  
printf("R0=%d,",(int)json_object_get_number(BufObject,"52.O2^+OP>O2+O^"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"53.O2^+O2>O2+O2^"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"54.O2^+O2>O2^+O2"));   
printf("R0=%d,",(int)json_object_get_number(BufObject,"55.O2^+O2>O^+OP+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"56.O2^+O2A>O2+O2^"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"57.O2^+O2B>O2+O2^"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"58.O2^+Ar>O2+Ar^"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"59.O2^+Ar>O2^+Ar^"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"60.O^+O2>OP+O2^"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"61.O^+O2>O^+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"62.O^+OP>OP+O^"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"63.O^+O2A>O2^+OP"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"64.O^+O2B>O2^+OP"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"65.Ar^+Ar>Ar+Ar^"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"66.Ar^+Ar>Ar++Ar"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"67.Ar^+O2>O2+Ar^"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"68.e+Ar*>e+Ar"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"69.O-+Ar^>OP+AR"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"70.O-+O2B>e+OP+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"71.OD+OP>2OP"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"72.OD+O2>OP+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"73.OD+O2>OP+O2A"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"74.OD+O2>OP+O2B"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"75.O2A+OP>OP+O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"76.O2A+O2>2O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"77.O2A+O2A>2O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"78.O2B+O2>2O2"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"79.Ar^+OP>Ar+O^"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"80.Ar^+O2>Ar+O2^"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"81.Ar*+Ar*>e+Ar+Ar^"));                
printf("R0=%d,",(int)json_object_get_number(BufObject,"82.Ar*+Ar>2Ar"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"83.Ar*+OP>OD+Ar"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"84.Ar*+OP>OP+Ar"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"85.Ar*+O2>2OP+Ar "));
printf("R0=%d,",(int)json_object_get_number(BufObject,"86.Ar*+O2>OP+OD+Ar"));
printf("R0=%d,",(int)json_object_get_number(BufObject,"87.Ar*+O2>O2+Ar"));
printf("\n");
   }else{
      exit(1);
   }   
}
void start() {
   int i;
   //-------------------------//
   //--------Geometry---------//
   //-------------------------//
   Gsize = ngx * ngy;
   ncx = ngx - 1;
   ncy = ngy - 1;
   Csize = ncx * ncy;
   dx=xlength/ncx;
   dy=ylength/ncy;
   idx=1/dx;
   idy=1/dy;
   dx2=dx*dx;
   dy2=dy*dy;
   dxdy2=dx2/dy2;
   hdx=0.5/dx;
   hdy=0.5/dy;
   r_eps0=1/EPS0;
   fncx=(float)ncx;
   fncy=(float)ncy;
   fngx=(float)ngx;
   fngy=(float)ngy;
   x_Garray=VFMalloc(ngx);
   for(i=0;i<ngx;i++) x_Garray[i]=i*dx;
   y_Garray=VFMalloc(ngy);
   for(i=0;i<ngy;i++) y_Garray[i]=i*dy;
}
void DumpRead(int argc, char *argv[]) {
   fprintf(stderr,"Tstrp = %d\n",tstep);
}
float **MFMalloc(int sizeX,int sizeY)
{
    int i;
    float **M;
    M=(float **)malloc(sizeX*sizeof(float *));
    for(i=0;i<sizeX;i++) M[i]=(float *)malloc(sizeY*sizeof(float));
    return M;
}
int ***TIMalloc(int sizeX,int sizeY,int sizeZ)
{
    int i,j;
    int ***T;
    T=(int ***)malloc(sizeX*sizeof(int **));
    for(i=0;i<sizeX;i++){
        T[i]=(int **)malloc(sizeY*sizeof(int *));
        for(j=0;j<sizeY;j++)	T[i][j]=(int *)malloc(sizeZ*sizeof(int));
    }
    return T;
}
float ***TFMalloc(int sizeX,int sizeY,int sizeZ)
{
    int i,j;
    float ***T;
    T=(float ***)malloc(sizeX*sizeof(float **));
    for(i=0;i<sizeX;i++){
        T[i]=(float **)malloc(sizeY*sizeof(float *));
        for(j=0;j<sizeY;j++)	T[i][j]=(float *)malloc(sizeZ*sizeof(float));
    }
    return T;
}
int **MIMalloc(int sizeX,int sizeY)
{
    int i;
    int **M;
    M=(int **)malloc(sizeX*sizeof(int *));
    for(i=0;i<sizeX;i++) M[i]=(int *)malloc(sizeY*sizeof(int));
    return M;
}
float *VFMalloc(int size)
{
    float *V;
  V=(float *)malloc(size*sizeof(float));
    return V;
}
int *VIMalloc(int size)
{
    int *V;
    V=(int *)malloc(size*sizeof(int));
    return V;
}
void TFFree(float ***T,int sizeX,int sizeY)
{
    int i,j;
    for(i=0;i<sizeX;i++){
        for(j=0;j<sizeY;j++) free(T[i][j]);
        free(T[i]);
    }
    free(T);
    return ;
}
void TIFree(int ***T,int sizeX,int sizeY)
{
    int i,j;
    for(i=0;i<sizeX;i++){
        for(j=0;j<sizeY;j++) free(T[i][j]);
        free(T[i]);
    }
    free(T);
    return ;
}
void MFFree(float **M,int sizeX)
{
    int i;
    for(i=0;i<sizeX;i++) free(M[i]);
    free(M);
    return ;
}
void MIFree(int **M,int sizeX)
{
    int i;
    for(i=0;i<sizeX;i++) free(M[i]);
    free(M);
    return ;
}
void TFInit(float ***T,float C,int sizeX,int sizeY,int sizeZ)
{
    int i,j,k;
    for(i=0;i<sizeX;i++)
        for(j=0;j<sizeY;j++)
            for(k=0;k<sizeZ;k++) T[i][j][k]=C;

    return;
}
void TIInit(int ***T,int C,int sizeX,int sizeY,int sizeZ)
{
    int i,j,k;
    for(i=0;i<sizeX;i++)
        for(j=0;j<sizeY;j++)
            for(k=0;k<sizeZ;k++) T[i][j][k]=C;

    return;
}
void MFInit(float **M,float C,int sizeX,int sizeY)
{
    int i,j;
    for(i=0;i<sizeX;i++)
        for(j=0;j<sizeY;j++) M[i][j]=C;

    return;
}
void MIInit(int **M,int C,int sizeX,int sizeY)
{
    int i,j;
    for(i=0;i<sizeX;i++)
        for(j=0;j<sizeY;j++) M[i][j]=C;

    return;
}
void VFInit(float *V,float C,int size)
{
    int i;
    for(i=0;i<size;i++) V[i]=C;

    return;
}
void VIInit(int *V,int C,int size)
{
    int i;
    for(i=0;i<size;i++) V[i]=C;

    return;
}
void MFCopy(float **M,float **C,int sizeX,int sizeY)
{
    int i,j;
    for(i=0;i<sizeX;i++){
        for(j=0;j<sizeY;j++){
         M[i][j]=C[i][j];
    }
}

    return;
}
void MICopy(int **M,int **C,int sizeX,int sizeY)
{
    int i,j;
    for(i=0;i<sizeX;i++)
        for(j=0;j<sizeY;j++) M[i][j]=C[i][j];

    return;
}
void VFCopy(float *V,float *C,int size)
{
    int i;
    for(i=0;i<size;i++) V[i]=C[i];

    return;
}
void VICopy(int *V,int *C,int size)
{
    int i;
    for(i=0;i<size;i++) V[i]=C[i];

    return;
}
