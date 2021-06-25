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
   printf("Device # = %d\n",(int)json_object_get_number(SubObject1,"Device"));
   //-------------------------//
   //--------Geometry---------//
   //-------------------------//
   SubObject2 = json_object_get_object(MainObject,"Geometry");
   BufObject = json_object_get_object(SubObject2,"SystemSpec");
   printf("X_length(m) = %g\n",(float)json_object_get_number(BufObject,"X_length(m)"));
   printf("Y_length(m) = %g\n",(float)json_object_get_number(BufObject,"Y_length(m)"));
   printf("Z_length(m) = %g\n",(float)json_object_get_number(BufObject,"Z_length(m)"));
   printf("NumCellx = %d\n",(int)json_object_get_number(BufObject,"NumGridx"));
   printf("NumCelly = %d\n",(int)json_object_get_number(BufObject,"NumGridy"));
   BufArray = json_object_get_array(SubObject2,"BoundaryCondition");
   for (i=0;i<json_array_get_count(BufArray);i++){
      BufObject = json_array_get_object(BufArray,i);
      printf("X0 = %d,",(int)json_object_get_number(BufObject,"X0"));
      printf("X1 = %d,",(int)json_object_get_number(BufObject,"X1"));
      printf("Y0 = %d,",(int)json_object_get_number(BufObject,"Y0"));
      printf("Y1 = %d,",(int)json_object_get_number(BufObject,"Y1"));
      printf("B.C = %d,",(int)json_object_get_number(BufObject,"B.C"));
      printf("Temp = %g\n",(float)json_object_get_number(BufObject,"Temp"));
   }
   BufArray = json_object_get_array(SubObject2,"ConductorSpec");
   for (i=0;i<json_array_get_count(BufArray);i++){
      BufObject = json_array_get_object(BufArray,i);
      printf("M_ID = %d,\n",(int)json_object_get_number(BufObject,"M_ID"));
      printf("X0 = %d,",(int)json_object_get_number(BufObject,"X0"));
      printf("X1 = %d,",(int)json_object_get_number(BufObject,"X1"));
      printf("Y0 = %d,",(int)json_object_get_number(BufObject,"Y0"));
      printf("Y1 = %d,",(int)json_object_get_number(BufObject,"Y1"));
      printf("Temp = %g",(float)json_object_get_number(BufObject,"Temp"));
      printf("DC = %g,\n",(float)json_object_get_number(BufObject,"DC(V)"));
      printf("Power1 = %g,",(float)json_object_get_number(BufObject,"Power1(W)"));
      printf("AC1 = %g,",(float)json_object_get_number(BufObject,"AC1(V)"));
      printf("Freq1 = %g,",(float)json_object_get_number(BufObject,"Freq1(1/s)"));
      printf("Phase1 = %g,\n",(float)json_object_get_number(BufObject,"Phase1(deg)"));
      printf("Power2 = %g,",(float)json_object_get_number(BufObject,"Power2(W)"));
      printf("AC2 = %g,",(float)json_object_get_number(BufObject,"AC2(V)"));
      printf("Freq2 = %g,",(float)json_object_get_number(BufObject,"Freq2(1/s)"));
      printf("Phase2 = %g,\n",(float)json_object_get_number(BufObject,"Phase2(deg)"));
      printf("R = %g,",(float)json_object_get_number(BufObject,"R(Ohm)"));
      printf("L = %g,",(float)json_object_get_number(BufObject,"L(H)"));
      printf("C = %g\n",(float)json_object_get_number(BufObject,"C(F)"));
   }
   BufArray = json_object_get_array(SubObject2,"DielectricSpec");
   for (i=0;i<json_array_get_count(BufArray);i++){
      BufObject = json_array_get_object(BufArray,i);
      printf("M_ID = %d,\n",(int)json_object_get_number(BufObject,"M_ID"));
      printf("X0 = %d,",(int)json_object_get_number(BufObject,"X0"));
      printf("X1 = %d,",(int)json_object_get_number(BufObject,"X1"));
      printf("Y0 = %d,",(int)json_object_get_number(BufObject,"Y0"));
      printf("Y1 = %d,",(int)json_object_get_number(BufObject,"Y1"));
      printf("Epsilon = %g\n",(float)json_object_get_number(BufObject,"Epsilon"));
   }
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
   fprintf(stderr,"Tstrp = %d\n",tstep);
}
void DumpRead(int argc, char *argv[]) {
   fprintf(stderr,"Tstrp = %d\n",tstep);
}