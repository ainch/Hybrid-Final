#include "start.h"

int IVnZC(char A[50],int Value){ // Int Value non Zero CHECK
   if(Value == 0){
      printf("\nERROR! - ");
      printf("\"%s\" = 0 \n",A);
      exit(1);
   }
   return Value;
}
float FVnZC(char A[50],float Value){ // float Value non Zero CHECK
   if(Value == 0.0){
      printf("\nERROR! - "); 
      printf("\"%s\" = 0 \n",A);
      exit(1);
   }
   return Value;
}
void InputRead() {
   int i,j;
   int buf;
   float fbuf1,fbuf2;
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
   printf("Read Geommetry\n"); 
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
   BoundaryNUM = (int)json_array_get_count(BufArray);
   printf("\tBoundary # = %d\n",BoundaryNUM); 
   if (BoundaryNUM < 4){
      printf("\"BoundaryCondition\" is too small!(BoundaryNum>=4)\n");
      exit(1);
   }
   BoundaryX0 = VIMalloc(BoundaryNUM);
   BoundaryY0 = VIMalloc(BoundaryNUM);
   BoundaryX1 = VIMalloc(BoundaryNUM);
   BoundaryY1 = VIMalloc(BoundaryNUM);
   BoundaryBC = VIMalloc(BoundaryNUM);
   BoundaryTEMP = VFMalloc(BoundaryNUM);
   for (i=0;i<BoundaryNUM;i++){
      BufObject = json_array_get_object(BufArray,i);
      BoundaryX0[i]=(int)json_object_get_number(BufObject,"X0");
      BoundaryX1[i]=(int)json_object_get_number(BufObject,"X1");
      if(BoundaryX0[i] > BoundaryX1[i]){
         printf("Line %d X position is error in \"BoundaryCondition\"(\"X0\" <= \"X1\")\n",i+1);
         exit(1);
      } 
      if(BoundaryX1[i]>=ngx){
         printf("Line %d X1 position is error in \"BoundaryCondition\"(\"X1\" < \"NumGridx\")\n",i+1);
         exit(1);
      }
      BoundaryY0[i]=(int)json_object_get_number(BufObject,"Y0");
      BoundaryY1[i]=(int)json_object_get_number(BufObject,"Y1");
      if(BoundaryY0[i] > BoundaryY1[i]){
         printf("Line %d Y position is error in \"BoundaryCondition\"(\"Y0\" <= \"Y1\")\n",i+1);
         exit(1);
      } 
      if(BoundaryY1[i]>=ngy){
         printf("Line %d Y1 position is error in \"BoundaryCondition\"(\"X1\" < \"NumGridy\")\n",i+1);
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
   CondNUM = (int)json_array_get_count(BufArray);
   printf("\tConductor # = %d\n",CondNUM); 
   if (CondNUM < 2){
      printf("\"Conductor\" is too small!(CondNum>=2)\n");
      exit(1);
   }
   CondM_ID = VIMalloc(CondNUM);
   CondX0 = VIMalloc(CondNUM);
   CondX1 = VIMalloc(CondNUM);
   CondY0 = VIMalloc(CondNUM);
   CondY1 = VIMalloc(CondNUM);
   CondTEMP = VFMalloc(CondNUM);
   for (i=0;i<CondNUM;i++){
      BufObject = json_array_get_object(BufArray,i);
      CondM_ID[i] = (int)json_object_get_number(BufObject,"M_ID");
      CondX0[i] = (int)json_object_get_number(BufObject,"X0");
      CondX1[i] = (int)json_object_get_number(BufObject,"X1");
      if(CondX0[i] > CondX1[i]){
         printf("Line %d X position is error in \"ConductorSpec\"(\"X0\" <= \"X1\")\n",i+1);
         exit(1);
      } 
      if(CondX1[i]>=ngx){
         printf("Line %d X1 position is error in \"ConductorSpec\"(\"X1\" < \"NumGridx\")\n",i+1);
         exit(1);
      }
      CondY0[i] = (int)json_object_get_number(BufObject,"Y0");
      CondY1[i] = (int)json_object_get_number(BufObject,"Y1");
      if(CondY0[i] > CondY1[i]){
         printf("Line %d Y position is error in \"ConductorSpec\"(\"Y0\" <= \"Y1\")\n",i+1);
         exit(1);
      } 
      if(CondY1[i]>=ngy){
         printf("Line %d Y1 position is error in \"ConductorSpec\"(\"Y1\" < \"NumGridy\")\n",i+1);
         exit(1);
      }
      CondTEMP[i] = (float)json_object_get_number(BufObject,"Temp");
   }
   BufArray = json_object_get_array(SubObject2,"Source");
   SrcNUM = (int)json_array_get_count(BufArray);
   printf("\tSource # = %d\n",SrcNUM); 
   SrcM_ID = VIMalloc(SrcNUM);
   SrcDC = VFMalloc(SrcNUM);
   SrcPOWER = VFMalloc(SrcNUM);
   SrcAC = VFMalloc(SrcNUM);
   SrcFREQ = VFMalloc(SrcNUM);
   SrcPHASE = VFMalloc(SrcNUM);
   SrcR = VFMalloc(SrcNUM);
   SrcL = VFMalloc(SrcNUM);
   SrcC = VFMalloc(SrcNUM);
   Min_FREQ = 1e200;
   Max_FREQ = 0.0;
   for (i=0;i<SrcNUM;i++){
      buf = 0;
      BufObject = json_array_get_object(BufArray,i);
      SrcM_ID[i] = (int)json_object_get_number(BufObject,"M_ID");
      for (j=0;j<CondNUM;j++){
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
         Max_FREQ = max(Max_FREQ,SrcFREQ[i]);
         if(SrcFREQ[i]<Min_FREQ)
            Min_FREQ  = SrcFREQ[i];
      }
      if(SrcDC[i] == 0 && SrcPOWER[i] == 0 && SrcAC[i] == 0){
         printf("Line %d Source is error(delete Line)\n",i+1);
         exit(1);
      }
   }
   BufArray = json_object_get_array(SubObject2,"DielectricSpec");
   DielNUM = (int)json_array_get_count(BufArray);
   printf("\tDielectric # = %d\n",DielNUM); 
   DielM_ID = VIMalloc(DielNUM);
   DielX0 = VIMalloc(DielNUM);
   DielX1 = VIMalloc(DielNUM);
   DielY0 = VIMalloc(DielNUM);
   DielY1 = VIMalloc(DielNUM);
   DielEPS = VFMalloc(DielNUM);
   for (i=0;i<DielNUM;i++){
      BufObject = json_array_get_object(BufArray,i);
      DielM_ID[i] = (int)json_object_get_number(BufObject,"M_ID");
      for (j=0;j<CondNUM;j++){
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
      if(DielX1[i]>=ngx){
         printf("Line %d X1 position is error in \"DielectricSpec\"(\"X1\" < \"NumGridx\")\n",i+1);
         exit(1);
      }
      DielY0[i] = (int)json_object_get_number(BufObject,"Y0");
      DielY1[i] = (int)json_object_get_number(BufObject,"Y1");
      if(DielY0[i] > DielY1[i]){
         printf("Line %d Y position is error in \"DielelctricSpec\"(\"Y0\" <= \"Y1\")\n",i+1);
         exit(1);
      } 
      if(DielY1[i]>=ngy){
         printf("Line %d Y1 position is error in \"DielectricSpec\"(\"Y1\" < \"NumGridy\")\n",i+1);
         exit(1);
      }
      DielEPS[i] = (float)json_object_get_number(BufObject,"Epsilon");
   }
   //-------------------------//
   //-------GasSpecies--------//
   //-------------------------//
   printf("Read GasSpecies\n");
   SubObject3 = json_object_get_object(MainObject,"GasSpecies");
   MainGas = (int)json_object_get_number(SubObject3,"Type(0:Ar,1:O2,2:Ar/O2)");
   switch (MainGas) {
	case ARGON: {	// ONLY ARGON
      printf("\tArgon gas\n"); 
		nsp = 2;		nfsp = 1;		nBG = 1;
		break;
	}
	case OXYGEN: { // Oxygen
      printf("\tOxygen gas\n"); 
		nsp = 4;		nfsp = 4;		nBG = 1;
		break;
	}
	case ARO2: { // Argon.Oxygen
      printf("\tArgon + Oxygen gas\n"); 
		nsp = 5;		nfsp = 5;		nBG = 2;
		break;
	}
   default: {
      printf("\t\"Type\" is error in GasSpecies.\n"); 
      exit(1);
   }
	}
   SP = (Species *) malloc(nsp * sizeof(Species));
	FG = (Fluid *) malloc(nfsp * sizeof(Fluid));
	BG = (BackG *) malloc(nBG * sizeof(BackG));
   Total_Pressure = (float)json_object_get_number(SubObject3,"TotalPres(Torr)");
   fbuf1 = 0.0;
   fbuf2 = 0.0;
   switch (MainGas) {
	case ARGON: {	// ONLY ARGON
      BufObject = json_object_get_object(SubObject3,"Background");
      BufObject2 = json_object_get_object(BufObject,"Argon");
      strcpy(BG[0].name,"Argon");
      BG[0].Ratio = 1.0;
      BG[0].Pres = Total_Pressure;
      BG[0].Temp = (float)json_object_get_number(BufObject2,"Temp(eV)");
      BG[0].mass = 39.950 * AMU;
      printf("\tAr Ratio = 100 %\n");
      BufObject = json_object_get_object(SubObject3,"NeutralSpecies");
      FG[0].Loadtype = (int)json_object_get_number(BufObject,"Loadtype");
      BufObject2 = json_object_get_object(BufObject,"LoadPosition(m)");
      FG[0].x_center = (float)json_object_get_number(BufObject2,"X0");
      FG[0].x_fall = (float)json_object_get_number(BufObject2,"X1");
      FG[0].y_center = (float)json_object_get_number(BufObject2,"Y0");
      FG[0].y_fall = (float)json_object_get_number(BufObject2,"Y1");
      if(FG[0].Loadtype == SMARTLOAD){
      }else if(FG[0].Loadtype == UNIFORM){
         if(FG[0].x_center>=FG[0].x_fall || FG[0].y_center>=FG[0].y_fall)
            exit(1);
      }else if(FG[0].Loadtype == EXPONETIAL){
         if(FG[0].x_center*FG[0].x_fall*FG[0].y_center*FG[0].y_fall < 0)
            exit(1);
      }else if(FG[0].Loadtype == COSINE){
         if(FG[0].x_center*FG[0].x_fall*FG[0].y_center*FG[0].y_fall < 0)
            exit(1);
      }else{
         printf("\t\"Loadtype\" is error in NeutralSpecies.\n"); 
         exit(1);
      }
      BufObject2 = json_object_get_object(BufObject,"Ar*");
      strcpy(FG[0].name,"Ar*");
      FG[0].InitDens =(float)json_object_get_number(BufObject2,"Density");
      FG[0].Temp = (float)json_object_get_number(BufObject2,"Temp(eV)");
      FG[0].mass = 39.950 * AMU;
      BufObject = json_object_get_object(SubObject3,"ChargeSpecies");
      SP[0].Loadtype = (int)json_object_get_number(BufObject,"Loadtype");
      BufObject2 = json_object_get_object(BufObject,"LoadPosition(m)");
      SP[0].x_center = (float)json_object_get_number(BufObject2,"X0");
      SP[0].x_fall = (float)json_object_get_number(BufObject2,"X1");
      SP[0].y_center = (float)json_object_get_number(BufObject2,"Y0");
      SP[0].y_fall = (float)json_object_get_number(BufObject2,"Y1");
      if(SP[0].Loadtype == SMARTLOAD){
      }else if(SP[0].Loadtype == UNIFORM){
         if(SP[0].x_center>=SP[0].x_fall || SP[0].y_center>=SP[0].y_fall)
            exit(1);
      }else if(SP[0].Loadtype == EXPONETIAL){
         if(SP[0].x_center*SP[0].x_fall*SP[0].y_center*SP[0].y_fall < 0)
            exit(1);
      }else if(FG[0].Loadtype == COSINE){
         if(SP[0].x_center*SP[0].x_fall*SP[0].y_center*SP[0].y_fall < 0)
            exit(1);
      }else{
         printf("\t\"Loadtype\" is error in ChargeSpecies.\n"); 
         exit(1);
      }
      BufObject2 = json_object_get_object(BufObject,"Electron");
      strcpy(SP[0].name,"Electron");
      SP[0].S_ID = (int)json_object_get_number(BufObject2,"S_ID");
      SP[0].InitDens = (float)json_object_get_number(BufObject2,"Density");
      SP[0].Temp = (float)json_object_get_number(BufObject2,"Temp(eV)");
      SP[0].np2c = (float)json_object_get_number(BufObject2,"np2c");
      SP[0].MAXNP = (int)json_object_get_number(BufObject2,"Max_np");
      SP[0].mass = 5.485799486e-4 * AMU;
      SP[0].q = -1.0 * CQ;
      BufObject2 = json_object_get_object(BufObject,"Ar+");
      strcpy(SP[1].name,"Ar+");
      SP[1].S_ID = (int)json_object_get_number(BufObject2,"S_ID");
      SP[1].InitDens = (float)json_object_get_number(BufObject2,"Density");
      SP[1].Temp = (float)json_object_get_number(BufObject2,"Temp(eV)");
      SP[1].np2c = (float)json_object_get_number(BufObject2,"np2c");
      SP[1].MAXNP = (int)json_object_get_number(BufObject2,"Max_np");
      SP[1].mass = 39.9500 * AMU;
      SP[1].q = 1.0 * CQ;
      break;
   }
	case OXYGEN: { // Oxygen
      BufObject = json_object_get_object(SubObject3,"Background");
      BufObject2 = json_object_get_object(BufObject,"Oxygen");
      strcpy(BG[0].name,"Oxygen");
      BG[0].Ratio = 1.0;
      BG[0].Pres = Total_Pressure;
      BG[0].Temp = (float)json_object_get_number(BufObject2,"Temp(eV)");
      BG[0].mass = 32.0000 * AMU;
      printf("\tO2 Ratio = 100 %\n");
      BufObject = json_object_get_object(SubObject3,"NeutralSpecies");
      FG[0].Loadtype = (int)json_object_get_number(BufObject,"Loadtype");
      BufObject2 = json_object_get_object(BufObject,"LoadPosition(m)");
      FG[0].x_center = (float)json_object_get_number(BufObject2,"X0");
      FG[0].x_fall = (float)json_object_get_number(BufObject2,"X1");
      FG[0].y_center = (float)json_object_get_number(BufObject2,"Y0");
      FG[0].y_fall = (float)json_object_get_number(BufObject2,"Y1");
      if(FG[0].Loadtype == SMARTLOAD){
      }else if(FG[0].Loadtype == UNIFORM){
         if(FG[0].x_center>=FG[0].x_fall || FG[0].y_center>=FG[0].y_fall)
            exit(1);
      }else if(FG[0].Loadtype == EXPONETIAL){
         if(FG[0].x_center*FG[0].x_fall*FG[0].y_center*FG[0].y_fall < 0)
            exit(1);
      }else if(FG[0].Loadtype == COSINE){
         if(FG[0].x_center*FG[0].x_fall*FG[0].y_center*FG[0].y_fall < 0)
            exit(1);
      }else{
         printf("\t\"Loadtype\" is error in NeutralSpecies.\n"); 
         exit(1);
      }
      BufObject2 = json_object_get_object(BufObject,"OP");
      strcpy(FG[0].name,"OP");
      FG[0].InitDens =(float)json_object_get_number(BufObject2,"Density");
      FG[0].Temp = (float)json_object_get_number(BufObject2,"Temp(eV)");
      FG[0].mass = 16.0000 * AMU;
      BufObject2 = json_object_get_object(BufObject,"OD");
      strcpy(FG[1].name,"OD");
      FG[1].InitDens =(float)json_object_get_number(BufObject2,"Density");
      FG[1].Temp = (float)json_object_get_number(BufObject2,"Temp(eV)");
      FG[1].mass = 16.0000 * AMU;
      BufObject2 = json_object_get_object(BufObject,"O2A");
      strcpy(FG[2].name,"O2A");
      FG[2].InitDens =(float)json_object_get_number(BufObject2,"Density");
      FG[2].Temp = (float)json_object_get_number(BufObject2,"Temp(eV)");
      FG[2].mass = 32.0000 * AMU;
      BufObject2 = json_object_get_object(BufObject,"O2B");
      strcpy(FG[3].name,"O2B");
      FG[3].InitDens =(float)json_object_get_number(BufObject2,"Density");
      FG[3].Temp = (float)json_object_get_number(BufObject2,"Temp(eV)");
      FG[3].mass = 32.0000 * AMU;
      BufObject = json_object_get_object(SubObject3,"ChargeSpecies");
      SP[0].Loadtype = (int)json_object_get_number(BufObject,"Loadtype");
      BufObject2 = json_object_get_object(BufObject,"LoadPosition(m)");
      SP[0].x_center = (float)json_object_get_number(BufObject2,"X0");
      SP[0].x_fall = (float)json_object_get_number(BufObject2,"X1");
      SP[0].y_center = (float)json_object_get_number(BufObject2,"Y0");
      SP[0].y_fall = (float)json_object_get_number(BufObject2,"Y1");
      if(SP[0].Loadtype == SMARTLOAD){
      }else if(SP[0].Loadtype == UNIFORM){
         if(SP[0].x_center>=SP[0].x_fall || SP[0].y_center>=SP[0].y_fall)
            exit(1);
      }else if(SP[0].Loadtype == EXPONETIAL){
         if(SP[0].x_center*SP[0].x_fall*SP[0].y_center*SP[0].y_fall < 0)
            exit(1);
      }else if(FG[0].Loadtype == COSINE){
         if(SP[0].x_center*SP[0].x_fall*SP[0].y_center*SP[0].y_fall < 0)
            exit(1);
      }else{
         printf("\t\"Loadtype\" is error in ChargeSpecies.\n"); 
         exit(1);
      }
      BufObject2 = json_object_get_object(BufObject,"Electron");
      strcpy(SP[0].name,"Electron");
      SP[0].S_ID = (int)json_object_get_number(BufObject2,"S_ID");
      SP[0].InitDens = (float)json_object_get_number(BufObject2,"Density");
      SP[0].Temp = (float)json_object_get_number(BufObject2,"Temp(eV)");
      SP[0].np2c = (float)json_object_get_number(BufObject2,"np2c");
      SP[0].MAXNP = (int)json_object_get_number(BufObject2,"Max_np");
      SP[0].mass = 5.485799486e-4 * AMU;
      SP[0].q = -1.0 * CQ;
      BufObject2 = json_object_get_object(BufObject,"O2+");
      strcpy(SP[1].name,"O2+");
      SP[1].S_ID = (int)json_object_get_number(BufObject2,"S_ID");
      SP[1].InitDens = (float)json_object_get_number(BufObject2,"Density");
      SP[1].Temp = (float)json_object_get_number(BufObject2,"Temp(eV)");
      SP[1].np2c = (float)json_object_get_number(BufObject2,"np2c");
      SP[1].MAXNP = (int)json_object_get_number(BufObject2,"Max_np");
      SP[1].mass = 32.0000 * AMU;
      SP[1].q = 1.0 * CQ;
      BufObject2 = json_object_get_object(BufObject,"O+");
      strcpy(SP[2].name,"O+");
      SP[2].S_ID = (int)json_object_get_number(BufObject2,"S_ID");
      SP[2].InitDens = (float)json_object_get_number(BufObject2,"Density");
      SP[2].Temp = (float)json_object_get_number(BufObject2,"Temp(eV)");
      SP[2].np2c = (float)json_object_get_number(BufObject2,"np2c");
      SP[2].MAXNP = (int)json_object_get_number(BufObject2,"Max_np");
      SP[2].mass = 16.000 * AMU;
      SP[2].q = 1.0 * CQ;
      BufObject2 = json_object_get_object(BufObject,"O-");
      strcpy(SP[3].name,"O-");
      SP[3].S_ID = (int)json_object_get_number(BufObject2,"S_ID");
      SP[3].InitDens = (float)json_object_get_number(BufObject2,"Density");
      SP[3].Temp = (float)json_object_get_number(BufObject2,"Temp(eV)");
      SP[3].np2c = (float)json_object_get_number(BufObject2,"np2c");
      SP[3].MAXNP = (int)json_object_get_number(BufObject2,"Max_np");
      SP[3].mass = 16.000 * AMU;
      SP[3].q = -1.0 * CQ;
      break;
   }
	case ARO2: { // Argon + Oxygen
      BufObject = json_object_get_object(SubObject3,"Background");
      BufObject2 = json_object_get_object(BufObject,"Argon");
      strcpy(BG[0].name,"Argon");
      fbuf1 = (float)json_object_get_number(BufObject2,"Ratio");
      BG[0].Temp = (float)json_object_get_number(BufObject2,"Temp(eV)");
      BG[0].mass = 39.950 * AMU;
      BufObject2 = json_object_get_object(BufObject,"Oxygen");
      strcpy(BG[1].name,"Oxygen");
      fbuf2 = (float)json_object_get_number(BufObject2,"Ratio");
      fbuf2 = fbuf2 / (fbuf1 + fbuf2);
      BG[1].Temp = (float)json_object_get_number(BufObject2,"Temp(eV)");
      BG[1].mass = 32.0000 * AMU;
      BG[0].Ratio = 1.0 - fbuf2;
      BG[0].Pres = Total_Pressure*BG[0].Ratio;
      BG[1].Ratio = fbuf2;
      BG[1].Pres = Total_Pressure*BG[1].Ratio;
      printf("\tAr Ratio = %3.3g, (%0.3g %)\n",(1-fbuf2)/(1-fbuf2),(1-fbuf2)*100);
      printf("\tO2 Ratio = %0.3g, (%0.3g %)\n",fbuf2/(1-fbuf2),(fbuf2)*100);
      fbuf1 = 0.0;
      fbuf2 = 0.0;
      BufObject = json_object_get_object(SubObject3,"NeutralSpecies");
      FG[0].Loadtype = (int)json_object_get_number(BufObject,"Loadtype");
      BufObject2 = json_object_get_object(BufObject,"LoadPosition(m)");
      FG[0].x_center = (float)json_object_get_number(BufObject2,"X0");
      FG[0].x_fall = (float)json_object_get_number(BufObject2,"X1");
      FG[0].y_center = (float)json_object_get_number(BufObject2,"Y0");
      FG[0].y_fall = (float)json_object_get_number(BufObject2,"Y1");
      if(FG[0].Loadtype == SMARTLOAD){
      }else if(FG[0].Loadtype == UNIFORM){
         if(FG[0].x_center>=FG[0].x_fall || FG[0].y_center>=FG[0].y_fall)
            exit(1);
      }else if(FG[0].Loadtype == EXPONETIAL){
         if(FG[0].x_center*FG[0].x_fall*FG[0].y_center*FG[0].y_fall < 0)
            exit(1);
      }else if(FG[0].Loadtype == COSINE){
         if(FG[0].x_center*FG[0].x_fall*FG[0].y_center*FG[0].y_fall < 0)
            exit(1);
      }else{
         printf("\t\"Loadtype\" is error in NeutralSpecies.\n"); 
         exit(1);
      }
      BufObject2 = json_object_get_object(BufObject,"Ar*");
      strcpy(FG[0].name,"Ar*");
      FG[0].InitDens =(float)json_object_get_number(BufObject2,"Density");
      FG[0].Temp = (float)json_object_get_number(BufObject2,"Temp(eV)");
      FG[0].mass = 39.950 * AMU;
      BufObject2 = json_object_get_object(BufObject,"OP");
      strcpy(FG[1].name,"OP");
      FG[1].InitDens =(float)json_object_get_number(BufObject2,"Density");
      FG[1].Temp = (float)json_object_get_number(BufObject2,"Temp(eV)");
      FG[1].mass = 16.0000 * AMU;
      BufObject2 = json_object_get_object(BufObject,"OD");
      strcpy(FG[2].name,"OD");
      FG[2].InitDens =(float)json_object_get_number(BufObject2,"Density");
      FG[2].Temp = (float)json_object_get_number(BufObject2,"Temp(eV)");
      FG[2].mass = 16.0000 * AMU;
      BufObject2 = json_object_get_object(BufObject,"O2A");
      strcpy(FG[3].name,"O2A");
      FG[3].InitDens =(float)json_object_get_number(BufObject2,"Density");
      FG[3].Temp = (float)json_object_get_number(BufObject2,"Temp(eV)");
      FG[3].mass = 32.0000 * AMU;
      BufObject2 = json_object_get_object(BufObject,"O2B");
      strcpy(FG[4].name,"O2B");
      FG[4].InitDens =(float)json_object_get_number(BufObject2,"Density");
      FG[4].Temp = (float)json_object_get_number(BufObject2,"Temp(eV)");
      FG[4].mass = 32.0000 * AMU;
      BufObject = json_object_get_object(SubObject3,"ChargeSpecies");
      SP[0].Loadtype = (int)json_object_get_number(BufObject,"Loadtype");
      BufObject2 = json_object_get_object(BufObject,"LoadPosition(m)");
      SP[0].x_center = (float)json_object_get_number(BufObject2,"X0");
      SP[0].x_fall = (float)json_object_get_number(BufObject2,"X1");
      SP[0].y_center = (float)json_object_get_number(BufObject2,"Y0");
      SP[0].y_fall = (float)json_object_get_number(BufObject2,"Y1");
      if(SP[0].Loadtype == SMARTLOAD){
      }else if(SP[0].Loadtype == UNIFORM){
         if(SP[0].x_center>=SP[0].x_fall || SP[0].y_center>=SP[0].y_fall)
            exit(1);
      }else if(SP[0].Loadtype == EXPONETIAL){
         if(SP[0].x_center*SP[0].x_fall*SP[0].y_center*SP[0].y_fall < 0)
            exit(1);
      }else if(FG[0].Loadtype == COSINE){
         if(SP[0].x_center*SP[0].x_fall*SP[0].y_center*SP[0].y_fall < 0)
            exit(1);
      }else{
         printf("\t\"Loadtype\" is error in ChargeSpecies.\n"); 
         exit(1);
      }
      BufObject2 = json_object_get_object(BufObject,"Electron");
      strcpy(SP[0].name,"Electron");
      SP[0].S_ID = (int)json_object_get_number(BufObject2,"S_ID");
      SP[0].InitDens = (float)json_object_get_number(BufObject2,"Density");
      SP[0].Temp = (float)json_object_get_number(BufObject2,"Temp(eV)");
      SP[0].np2c = (float)json_object_get_number(BufObject2,"np2c");
      SP[0].MAXNP = (int)json_object_get_number(BufObject2,"Max_np");
      SP[0].mass = 5.485799486e-4 * AMU;
      SP[0].q = -1.0 * CQ;
      BufObject2 = json_object_get_object(BufObject,"Ar+");
      strcpy(SP[1].name,"Ar+");
      SP[1].S_ID = (int)json_object_get_number(BufObject2,"S_ID");
      SP[1].InitDens = (float)json_object_get_number(BufObject2,"Density");
      SP[1].Temp = (float)json_object_get_number(BufObject2,"Temp(eV)");
      SP[1].np2c = (float)json_object_get_number(BufObject2,"np2c");
      SP[1].MAXNP = (int)json_object_get_number(BufObject2,"Max_np");
      SP[1].mass = 39.9500 * AMU;
      SP[1].q = 1.0 * CQ;
      BufObject2 = json_object_get_object(BufObject,"O2+");
      strcpy(SP[2].name,"O2+");
      SP[2].S_ID = (int)json_object_get_number(BufObject2,"S_ID");
      SP[2].InitDens = (float)json_object_get_number(BufObject2,"Density");
      SP[2].Temp = (float)json_object_get_number(BufObject2,"Temp(eV)");
      SP[2].np2c = (float)json_object_get_number(BufObject2,"np2c");
      SP[2].MAXNP = (int)json_object_get_number(BufObject2,"Max_np");
      SP[2].mass = 32.0000 * AMU;
      SP[2].q = 1.0 * CQ;
      BufObject2 = json_object_get_object(BufObject,"O+");
      strcpy(SP[3].name,"O+");
      SP[3].S_ID = (int)json_object_get_number(BufObject2,"S_ID");
      SP[3].InitDens = (float)json_object_get_number(BufObject2,"Density");
      SP[3].Temp = (float)json_object_get_number(BufObject2,"Temp(eV)");
      SP[3].np2c = (float)json_object_get_number(BufObject2,"np2c");
      SP[3].MAXNP = (int)json_object_get_number(BufObject2,"Max_np");
      SP[3].mass = 16.0000 * AMU;
      SP[3].q = 1.0 * CQ;
      BufObject2 = json_object_get_object(BufObject,"O-");
      strcpy(SP[4].name,"O-");
      SP[4].S_ID = (int)json_object_get_number(BufObject2,"S_ID");
      SP[4].InitDens = (float)json_object_get_number(BufObject2,"Density");
      SP[4].Temp = (float)json_object_get_number(BufObject2,"Temp(eV)");
      SP[4].np2c = (float)json_object_get_number(BufObject2,"np2c");
      SP[4].MAXNP = (int)json_object_get_number(BufObject2,"Max_np");
      SP[4].mass = 16.0000 * AMU;
      SP[4].q = -1.0 * CQ;
   }
   }
   //-------------------------//
   //----SimulationMethod-----//
   //-------------------------//
   printf("Read SimulationMethod\n");
   SubObject4 = json_object_get_object(MainObject,"SimulationMethod");   
   DT_PIC = IVnZC("TimeStep_PIC",(int)json_object_get_number(SubObject4,"TimeStep_PIC"));
   DT_CONTI = IVnZC("TimeStep_Conti",(int)json_object_get_number(SubObject4,"TimeStep_Conti"));
   dt = 1.0 / Max_FREQ / (float)DT_PIC;
   dtc = dt * (float)DT_CONTI;
   printf("\tPIC TimeStep = %g s\n",dt);
   printf("\tContinuity TimeStep = %g s\n",dtc);
   PCGtol = FVnZC("PCGMarginOfError",(float)json_object_get_number(SubObject4,"PCGMarginOfError"));
   HISTMAX = IVnZC("HistoryMax",(int)json_object_get_number(SubObject4,"HistoryMax"));
   dHIST = IVnZC("HistoryDivide",(int)json_object_get_number(SubObject4,"HistoryDivide"));
   N_ave = IVnZC("AverageIter",(int)json_object_get_number(SubObject4,"AverageIter"));
   N_smt = (int)json_object_get_number(SubObject4,"DensitySmoothing");
   np_lim = IVnZC("ParticleLimit",(int)json_object_get_number(SubObject4,"ParticleLimit"));
   ConstB_Flag = (int)json_object_get_number(SubObject4,"ConstB_Flag");
   if(ConstB_Flag){
      ConstBFile = (char*)json_object_get_string(SubObject4,"ConstB_File");
      printf("\tConstant magnetic DATA = %s \n",ConstBFile);
   }
   /* 
   BufObject = json_object_get_object(SubObject4,"Initialing_Int");
   printf("Initialing_Init_flag = %d\n",(int)json_object_get_number(BufObject,"Flag"));
   printf("Time_order = %d\n",(int)json_object_get_number(BufObject,"Time_order"));
   printf("Dump_order = %d\n",(int)json_object_get_number(BufObject,"Dump_order"));
   */
   /*
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
   */
   /*
   BufObject = json_object_get_object(SubObject4,"ContinuityBreakPoint");
   printf("Margin(%) = %g\n",(float)json_object_get_number(BufObject,"Margin(%)"));
   */
   /*
   BufObject = json_object_get_object(SubObject4,"PowerDrivenOption");
   printf("VoltageUpRatio = %g\n",(float)json_object_get_number(BufObject,"VoltageUpRatio"));
   */
   /*
   BufObject = json_object_get_object(SubObject4,"SimulationStop");
   printf("EndTime(s) = %g\n",(float)json_object_get_number(BufObject,"EndTime(s)"));
   BufObject2 = json_object_get_object(BufObject,"AutoSaturation");
   printf("Flag = %d\n",(int)json_object_get_number(BufObject2,"Flag"));
   printf("AveMargin(%) = %g\n",(float)json_object_get_number(BufObject2,"AveMargin(%)"));
   printf("SatGoal = %d\n",(int)json_object_get_number(BufObject2,"SatGoal"));
   */
   //-------------------------//
   //-------Diagnostics-------//
   //-------------------------//
   printf("Read Diagnostics\n");
   SubObject5 = json_object_get_object(MainObject,"Diagnostics");
   Basic_Flag = (int)json_object_get_number(SubObject5,"BasicDiagnostic");
   if(Basic_Flag){
      printf("\tBasic_Flag = %d",Basic_Flag);
      printf(", General Diagostics mode\n");
   }else{
      printf("\tMinimized Diagnostics mode\n");
   }
   /*
   BufObject = json_object_get_object(SubObject5,"TecplotSave");
   printf("Tecplot2D = %d\n",(int)json_object_get_number(BufObject,"Tecplot2D"));
   printf("Tec_Movie = %d\n",(int)json_object_get_number(BufObject,"Tec_Movie"));
   printf("Tec_Movie_FrameNum = %d\n",(int)json_object_get_number(BufObject,"Tec_Movie_FrameNum"));
   printf("Tec_Ave_Movie = %d\n",(int)json_object_get_number(BufObject,"Tec_Ave_Movie"));
   printf("Tec_Ave_Movie_Interval = %d\n",(int)json_object_get_number(BufObject,"Tec_Ave_Movie_Interval"));
   printf("Tec_Ave_Movie_num = %d\n",(int)json_object_get_number(BufObject,"Tec_Ave_Movie_num")); 
   */
   /*
   BufObject = json_object_get_object(SubObject5,"DumpFileSave");
   BufArray = json_object_get_array(BufObject,"Cycle");
   for (i=0;i<json_array_get_count(BufArray);i++){
      BufObject2 = json_array_get_object(BufArray,i);
      printf("Start = %g,",(float)json_object_get_number(BufObject2,"Start"));
      printf("End = %g,",(float)json_object_get_number(BufObject2,"End"));
      printf("Term = %g\n",(float)json_object_get_number(BufObject2,"Term"));
   }
   */
   /*
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
   */
   /*
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
   */
   /*
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
   */
   /*
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
   */
   /*
   BufObject = json_object_get_object(SubObject5,"Test_Particle");
   printf("Flag = %d\n",(int)json_object_get_number(BufObject,"Flag"));
   printf("MCC_Flag = %d\n",(int)json_object_get_number(BufObject,"MCC_Flag"));
   printf("NumberOfParticle = %d\n",(int)json_object_get_number(BufObject,"NumberOfParticle"));
   printf("Init_Position = %d\n",(int)json_object_get_number(BufObject,"Init_Position"));
   printf("Init_Velocity = %d\n",(int)json_object_get_number(BufObject,"Init_Velocity"));
   printf("Tecplot_Save = %d\n",(int)json_object_get_number(BufObject,"Tecplot_Save")); 
   */   
  
   //-------------------------//
   //-----Collision_On_Off-----//
   //-------------------------//
   printf("Read Collision_On_Off\n");
   SubObject6 = json_object_get_object(MainObject,"Collision_On_Off");
   if(MainGas==0){
      nRct_cx = 7;	// Number of reaction using Cross section data
		nRct_rc = 3;	// Number of reaction using Rate coefficient(constant)
      TnRct = nRct_cx + nRct_rc;
		mMnum = 1;
      Coll_Flag = (CollF *) malloc(TnRct * sizeof(CollF));
      for(i=0;i<TnRct;i++){
         Coll_Flag[i].Flag = 0;
	      Coll_Flag[i].mofM = 0.0; // target/projectile
	      Coll_Flag[i].Th_e = 0.0; //Threshold energy
	      Coll_Flag[i].RR = 0.0; // Reaction RATE
      }
      BufObject = json_object_get_object(SubObject6,"ArgonCase");
      CX_TEC_Flag = (int)json_object_get_number(BufObject,"TecplotSave");
      Coll_Flag[0].Flag = (float)json_object_get_number(BufObject,"0.e+Ar>e+Ar");
      Coll_Flag[0].mofM = SP[0].mass/BG[0].mass;
      Coll_Flag[1].Flag = (float)json_object_get_number(BufObject,"1.e+Ar>e+Ar*");
      Coll_Flag[1].mofM = SP[0].mass/BG[0].mass;
      Coll_Flag[2].Flag = (float)json_object_get_number(BufObject,"2.e+Ar>e+Ar*m");
      Coll_Flag[2].mofM = SP[0].mass/BG[0].mass;
      Coll_Flag[3].Flag = (float)json_object_get_number(BufObject,"3.e+Ar>2e+Ar^");
      Coll_Flag[3].mofM = SP[0].mass/BG[0].mass;
      Coll_Flag[4].Flag = (float)json_object_get_number(BufObject,"4.e+Ar*m>2e+Ar^");
      Coll_Flag[4].mofM = SP[0].mass/FG[0].mass;
      Coll_Flag[5].Flag = (float)json_object_get_number(BufObject,"5.Ar+Ar^>Ar^+Ar");
      Coll_Flag[5].mofM = BG[0].mass/SP[1].mass;
      Coll_Flag[6].Flag = (float)json_object_get_number(BufObject,"6.Ar+Ar^>Ar+Ar^");
      Coll_Flag[6].mofM = BG[0].mass/SP[1].mass;
      Coll_Flag[7].Flag = (float)json_object_get_number(BufObject,"7.e+Ar*m>e+Ar");
      Coll_Flag[7].mofM = SP[0].mass/FG[0].mass;
      Coll_Flag[8].Flag = (float)json_object_get_number(BufObject,"8.Ar*m+Ar>e+Ar^+Ar");
      Coll_Flag[8].mofM = FG[0].mass/BG[0].mass;
      Coll_Flag[9].Flag = (float)json_object_get_number(BufObject,"9.Ar+Ar*m>Ar+Ar");
      Coll_Flag[9].mofM = BG[0].mass/FG[0].mass;
      printf("\tReaction = %d\n",TnRct);
      printf("\tCross section = %d\n",nRct_cx);
      printf("\tReaction rate = %d\n",nRct_rc);
      printf("\tCheck reaction on or off...\n");
      for(i=0;i<TnRct;i++){
         if(Coll_Flag[i].Flag < 0)
            Coll_Flag[i].Flag = 0.0;
         if(Coll_Flag[i].Flag == 0)
            printf("\tReaction %d is off\n",i);
      }
      Argon_CrossSectionSET(Coll_Flag);
      printf("\tReaction setting Complete.\n");
   }else if(MainGas==1){
      nRct_cx = 58;	// Number of reaction using Cross section data
		nRct_rc = 9;	// Number of reaction using Rate coefficient(constant)
      TnRct = nRct_cx + nRct_rc;
		mMnum = 6;
      Coll_Flag = (CollF *) malloc(TnRct * sizeof(CollF));
      for(i=0;i<TnRct;i++){
         Coll_Flag[i].Flag = 0;
	      Coll_Flag[i].mofM = 0.0; // target/projectile
	      Coll_Flag[i].Th_e = 0.0; //Threshold energy
	      Coll_Flag[i].RR = 0.0; // Reaction RATE
      }
      BufObject = json_object_get_object(SubObject6,"OxygenCase");
      CX_TEC_Flag = (int)json_object_get_number(BufObject,"TecplotSave");
      Coll_Flag[0].Flag = (float)json_object_get_number(BufObject,"0.e+O2>e+O2");
      Coll_Flag[0].mofM = SP[0].mass/BG[0].mass;
      Coll_Flag[1].Flag = (float)json_object_get_number(BufObject,"1.e+O2>e+O2*");
      Coll_Flag[1].mofM = SP[0].mass/BG[0].mass;
      Coll_Flag[2].Flag = (float)json_object_get_number(BufObject,"2.e+O2>e+O2*");
      Coll_Flag[2].mofM = SP[0].mass/BG[0].mass;
      Coll_Flag[3].Flag = (float)json_object_get_number(BufObject,"3.e+O2>e+O2A");
      Coll_Flag[3].mofM = SP[0].mass/BG[0].mass;
      Coll_Flag[4].Flag = (float)json_object_get_number(BufObject,"4.e+O2>e+O2B");
      Coll_Flag[4].mofM = SP[0].mass/BG[0].mass;
      Coll_Flag[5].Flag = (float)json_object_get_number(BufObject,"5.e+O2>e+O2*");
      Coll_Flag[5].mofM = SP[0].mass/BG[0].mass;
      Coll_Flag[6].Flag = (float)json_object_get_number(BufObject,"6.e+O2>OP+O-");
      Coll_Flag[6].mofM = SP[0].mass/BG[0].mass;
      Coll_Flag[7].Flag = (float)json_object_get_number(BufObject,"7.e+O2>e+2OP");
      Coll_Flag[7].mofM = SP[0].mass/BG[0].mass;
      Coll_Flag[8].Flag = (float)json_object_get_number(BufObject,"8.e+O2>e+OP+OD");
      Coll_Flag[8].mofM = SP[0].mass/BG[0].mass;
      Coll_Flag[9].Flag = (float)json_object_get_number(BufObject,"9.e+O2>e+2OD");
      Coll_Flag[9].mofM = SP[0].mass/BG[0].mass;
      Coll_Flag[10].Flag = (float)json_object_get_number(BufObject,"10.e+O2>2e+O2^");
      Coll_Flag[10].mofM = SP[0].mass/BG[0].mass;
      Coll_Flag[11].Flag = (float)json_object_get_number(BufObject,"11.e+O2>e+OP+O*");
      Coll_Flag[11].mofM = SP[0].mass/BG[0].mass;
      Coll_Flag[12].Flag = (float)json_object_get_number(BufObject,"12.e+O2>e+O^+O-");
      Coll_Flag[12].mofM = SP[0].mass/BG[0].mass;
      Coll_Flag[13].Flag = (float)json_object_get_number(BufObject,"13.e+O2>2e+O^+OP");
      Coll_Flag[13].mofM = SP[0].mass/BG[0].mass;
      Coll_Flag[14].Flag = (float)json_object_get_number(BufObject,"14.e+O2A>2e+O2+");
      Coll_Flag[14].mofM = SP[0].mass/FG[0].mass;
      Coll_Flag[15].Flag = (float)json_object_get_number(BufObject,"15.e+O2A>OP+O-");
      Coll_Flag[15].mofM = SP[0].mass/FG[0].mass;
      Coll_Flag[16].Flag = (float)json_object_get_number(BufObject,"16.e+O2A>e+O2");
      Coll_Flag[16].mofM = SP[0].mass/FG[0].mass;
      Coll_Flag[17].Flag = (float)json_object_get_number(BufObject,"17.e+O2A>e+O2");
      Coll_Flag[17].mofM = SP[0].mass/FG[0].mass;
      Coll_Flag[18].Flag = (float)json_object_get_number(BufObject,"18.e+O2A>e+2OP");
      Coll_Flag[18].mofM = SP[0].mass/FG[0].mass;
      Coll_Flag[19].Flag = (float)json_object_get_number(BufObject,"19.e+O2A>e+OP+OD");
      Coll_Flag[19].mofM = SP[0].mass/FG[0].mass;
      Coll_Flag[20].Flag = (float)json_object_get_number(BufObject,"20.e+O2A>e+2OD");
      Coll_Flag[20].mofM = SP[0].mass/FG[0].mass;
      Coll_Flag[21].Flag = (float)json_object_get_number(BufObject,"21.e+O2A>2e+O^+OP");
      Coll_Flag[21].mofM = SP[0].mass/FG[0].mass;
      Coll_Flag[22].Flag = (float)json_object_get_number(BufObject,"22.e+O2B>2e+O2^");
      Coll_Flag[22].mofM = SP[0].mass/FG[1].mass;
      Coll_Flag[23].Flag = (float)json_object_get_number(BufObject,"23.e+O2B>OP+O-");
      Coll_Flag[23].mofM = SP[0].mass/FG[1].mass;
      Coll_Flag[24].Flag = (float)json_object_get_number(BufObject,"24.e+O2B>e+O2");
      Coll_Flag[24].mofM = SP[0].mass/FG[1].mass;
      Coll_Flag[25].Flag = (float)json_object_get_number(BufObject,"25.e+O2B>e+O2");
      Coll_Flag[25].mofM = SP[0].mass/FG[1].mass;
      Coll_Flag[26].Flag = (float)json_object_get_number(BufObject,"26.e+O2B>e+2O");
      Coll_Flag[26].mofM = SP[0].mass/FG[1].mass;
      Coll_Flag[27].Flag = (float)json_object_get_number(BufObject,"27.e+O2B>e+OP+OD");
      Coll_Flag[27].mofM = SP[0].mass/FG[1].mass;
      Coll_Flag[28].Flag = (float)json_object_get_number(BufObject,"28.e+O2B>e+2OD");
      Coll_Flag[28].mofM = SP[0].mass/FG[1].mass;
      Coll_Flag[29].Flag = (float)json_object_get_number(BufObject,"29.e+O2B>2e+O^+OP");
      Coll_Flag[29].mofM = SP[0].mass/FG[1].mass;
      Coll_Flag[30].Flag = (float)json_object_get_number(BufObject,"30.e+O->2e+OP");
      Coll_Flag[30].mofM = SP[0].mass/SP[3].mass;
      Coll_Flag[31].Flag = (float)json_object_get_number(BufObject,"31.e+O2^>OP+OD");
      Coll_Flag[31].mofM = SP[0].mass/SP[1].mass;
      Coll_Flag[32].Flag = (float)json_object_get_number(BufObject,"32.e+OP>e+OP");
      Coll_Flag[32].mofM = SP[0].mass/FG[2].mass;
      Coll_Flag[33].Flag = (float)json_object_get_number(BufObject,"33.e+OP>e+OD");
      Coll_Flag[33].mofM = SP[0].mass/FG[2].mass;
      Coll_Flag[34].Flag = (float)json_object_get_number(BufObject,"34.e+OP>e+O*");
      Coll_Flag[34].mofM = SP[0].mass/FG[2].mass;
      Coll_Flag[35].Flag = (float)json_object_get_number(BufObject,"35.e+OP>e+O*");
      Coll_Flag[35].mofM = SP[0].mass/FG[2].mass;
      Coll_Flag[36].Flag = (float)json_object_get_number(BufObject,"36.e+OP>e+O*");
      Coll_Flag[36].mofM = SP[0].mass/FG[2].mass;
      Coll_Flag[37].Flag = (float)json_object_get_number(BufObject,"37.e+OP>2e+O^");
      Coll_Flag[37].mofM = SP[0].mass/FG[2].mass;
      Coll_Flag[38].Flag = (float)json_object_get_number(BufObject,"38.e+OP>e+O*");
      Coll_Flag[38].mofM = SP[0].mass/FG[2].mass;
      Coll_Flag[39].Flag = (float)json_object_get_number(BufObject,"39.e+OD>2e+O^");
      Coll_Flag[39].mofM = SP[0].mass/FG[3].mass;
      Coll_Flag[40].Flag = (float)json_object_get_number(BufObject,"40.e+OD>e+OP");
      Coll_Flag[40].mofM = SP[0].mass/FG[3].mass;
      Coll_Flag[41].Flag = (float)json_object_get_number(BufObject,"41.O-+O2>O-+O2");
      Coll_Flag[41].mofM = SP[3].mass/BG[0].mass;
      Coll_Flag[42].Flag = (float)json_object_get_number(BufObject,"42.O-+O2>e+OP+O2");
      Coll_Flag[42].mofM = SP[3].mass/BG[0].mass;
      Coll_Flag[43].Flag = (float)json_object_get_number(BufObject,"43.O-+OP>e+O2");
      Coll_Flag[43].mofM = SP[3].mass/FG[2].mass;
      Coll_Flag[44].Flag = (float)json_object_get_number(BufObject,"44.O-+O2^>OP+O2");
      Coll_Flag[44].mofM = SP[3].mass/SP[1].mass;
      Coll_Flag[45].Flag = (float)json_object_get_number(BufObject,"45.O-+O^>2OP");
      Coll_Flag[45].mofM = SP[3].mass/SP[2].mass;
      Coll_Flag[46].Flag = (float)json_object_get_number(BufObject,"46.O-+O2A>e+OP+O2");
      Coll_Flag[46].mofM = SP[3].mass/FG[0].mass;
      Coll_Flag[47].Flag = (float)json_object_get_number(BufObject,"47.O2^+OP>O2+O^");
      Coll_Flag[47].mofM = SP[1].mass/FG[2].mass;
      Coll_Flag[48].Flag = (float)json_object_get_number(BufObject,"48.O2^+O2>O2+O2^");
      Coll_Flag[48].mofM = SP[1].mass/BG[0].mass;
      Coll_Flag[49].Flag = (float)json_object_get_number(BufObject,"49.O2^+O2>O2^+O2");
      Coll_Flag[49].mofM = SP[1].mass/BG[0].mass;
      Coll_Flag[50].Flag = (float)json_object_get_number(BufObject,"50.O2^+O2>O^+OP+O2");
      Coll_Flag[50].mofM = SP[1].mass/BG[0].mass;
      Coll_Flag[51].Flag = (float)json_object_get_number(BufObject,"51.O2^+O2A>O2+O2^");
      Coll_Flag[51].mofM = SP[1].mass/FG[0].mass;
      Coll_Flag[52].Flag = (float)json_object_get_number(BufObject,"52.O2^+O2B>O2+O2^");
      Coll_Flag[52].mofM = SP[1].mass/FG[1].mass;
      Coll_Flag[53].Flag = (float)json_object_get_number(BufObject,"53.O^+O2>OP+O2^");
      Coll_Flag[53].mofM = SP[2].mass/BG[0].mass;
      Coll_Flag[54].Flag = (float)json_object_get_number(BufObject,"54.O^+O2>O^+O2");
      Coll_Flag[54].mofM = SP[2].mass/BG[0].mass;
      Coll_Flag[55].Flag = (float)json_object_get_number(BufObject,"55.O^+OP>OP+O^");
      Coll_Flag[55].mofM = SP[2].mass/FG[2].mass;
      Coll_Flag[56].Flag = (float)json_object_get_number(BufObject,"56.O^+O2A>O2^+OP");
      Coll_Flag[56].mofM = SP[2].mass/FG[0].mass;
      Coll_Flag[57].Flag = (float)json_object_get_number(BufObject,"57.O^+O2B>O2^+OP");
      Coll_Flag[57].mofM = SP[2].mass/FG[1].mass;
      Coll_Flag[58].Flag = (float)json_object_get_number(BufObject,"58.O-+O2B>e+OP+O2");
      Coll_Flag[58].mofM = SP[3].mass/FG[1].mass;
      Coll_Flag[59].Flag = (float)json_object_get_number(BufObject,"59.OP+OD>2OP");
      Coll_Flag[59].mofM = FG[2].mass/FG[3].mass;
      Coll_Flag[60].Flag = (float)json_object_get_number(BufObject,"60.OD+O2>OP+O2");
      Coll_Flag[60].mofM = FG[3].mass/BG[0].mass;
      Coll_Flag[61].Flag = (float)json_object_get_number(BufObject,"61.OD+O2>OP+O2A");
      Coll_Flag[61].mofM = FG[3].mass/BG[0].mass;
      Coll_Flag[62].Flag = (float)json_object_get_number(BufObject,"62.OD+O2>OP+O2B");
      Coll_Flag[62].mofM = FG[3].mass/BG[0].mass;
      Coll_Flag[63].Flag = (float)json_object_get_number(BufObject,"63.O2A+OP>OP+O2");
      Coll_Flag[63].mofM = FG[1].mass/FG[2].mass;
      Coll_Flag[64].Flag = (float)json_object_get_number(BufObject,"64.O2A+O2>2O2");
      Coll_Flag[64].mofM = FG[0].mass/BG[0].mass;
      Coll_Flag[65].Flag = (float)json_object_get_number(BufObject,"65.O2A+O2A>2O2");
      Coll_Flag[65].mofM = FG[0].mass/FG[0].mass;
      Coll_Flag[66].Flag = (float)json_object_get_number(BufObject,"66.O2B+O2>2O2");
      Coll_Flag[66].mofM = FG[1].mass/BG[0].mass;
      printf("\tReaction = %d\n",TnRct);
      printf("\tCross section = %d\n",nRct_cx);
      printf("\tReaction rate = %d\n",nRct_rc);
      printf("\tCheck reaction on or off...\n");
      for(i=0;i<TnRct;i++){
         if(Coll_Flag[i].Flag < 0)
            Coll_Flag[i].Flag = 0.0;
         if(Coll_Flag[i].Flag == 0)
            printf("\tReaction %d is off\n",i);
      }
      Oxygen_CrossSectionSET(Coll_Flag);     
      printf("\tReaction setting Complete.\n");
   }else if(MainGas==2){
      nRct_cx = 68;	// Number of reaction using Cross section data
		nRct_rc = 20;	// Number of reaction using Rate coefficient(constant)
      TnRct = nRct_cx + nRct_rc;
		mMnum = 10;
      Coll_Flag = (CollF *) malloc(TnRct * sizeof(CollF));
      for(i=0;i<TnRct;i++){
         Coll_Flag[i].Flag = 0;
	      Coll_Flag[i].mofM = 0.0; // target/projectile
	      Coll_Flag[i].Th_e = 0.0; //Threshold energy
	      Coll_Flag[i].RR = 0.0; // Reaction RATE
      }
      BufObject = json_object_get_object(SubObject6,"Argon/OxygenCase");
      CX_TEC_Flag = (int)json_object_get_number(BufObject,"TecplotSave");
      Coll_Flag[0].Flag = (float)json_object_get_number(BufObject,"0.e+Ar>e+Ar");
      Coll_Flag[0].mofM = SP[0].mass/BG[0].mass;
      Coll_Flag[1].Flag = (float)json_object_get_number(BufObject,"1.e+Ar>e+Ar*");
      Coll_Flag[1].mofM = SP[0].mass/BG[0].mass;
      Coll_Flag[2].Flag = (float)json_object_get_number(BufObject,"2.e+Ar>e+Ar*m");
      Coll_Flag[2].mofM = SP[0].mass/BG[0].mass;
      Coll_Flag[3].Flag = (float)json_object_get_number(BufObject,"3.e+Ar>2e+Ar^");
      Coll_Flag[3].mofM = SP[0].mass/BG[0].mass;
      Coll_Flag[4].Flag = (float)json_object_get_number(BufObject,"4.e+Ar*m>2e+Ar^");
      Coll_Flag[4].mofM = SP[0].mass/FG[0].mass;
      Coll_Flag[5].Flag = (float)json_object_get_number(BufObject,"5.e+O2>e+O2");
      Coll_Flag[5].mofM = SP[0].mass/BG[1].mass;
      Coll_Flag[6].Flag = (float)json_object_get_number(BufObject,"6.e+O2>e+O2*");
      Coll_Flag[6].mofM = SP[0].mass/BG[1].mass;
      Coll_Flag[7].Flag = (float)json_object_get_number(BufObject,"7.e+O2>e+O2*");
      Coll_Flag[7].mofM = SP[0].mass/BG[1].mass;
      Coll_Flag[8].Flag = (float)json_object_get_number(BufObject,"8.e+O2>e+O2A");
      Coll_Flag[8].mofM = SP[0].mass/BG[1].mass;
      Coll_Flag[9].Flag = (float)json_object_get_number(BufObject,"9.e+O2>e+O2B");
      Coll_Flag[9].mofM = SP[0].mass/BG[1].mass;
      Coll_Flag[10].Flag = (float)json_object_get_number(BufObject,"10.e+O2>e+O2*");
      Coll_Flag[10].mofM = SP[0].mass/BG[1].mass;
      Coll_Flag[11].Flag = (float)json_object_get_number(BufObject,"11.e+O2>OP+O-");
      Coll_Flag[11].mofM = SP[0].mass/BG[1].mass;
      Coll_Flag[12].Flag = (float)json_object_get_number(BufObject,"12.e+O2>e+2OP");
      Coll_Flag[12].mofM = SP[0].mass/BG[1].mass;
      Coll_Flag[13].Flag = (float)json_object_get_number(BufObject,"13.e+O2>e+OP+OD");
      Coll_Flag[13].mofM = SP[0].mass/BG[1].mass;
      Coll_Flag[14].Flag = (float)json_object_get_number(BufObject,"14.e+O2>e+2OD");
      Coll_Flag[14].mofM = SP[0].mass/BG[1].mass;
      Coll_Flag[15].Flag = (float)json_object_get_number(BufObject,"15.e+O2>2e+O2^");
      Coll_Flag[15].mofM = SP[0].mass/BG[1].mass;
      Coll_Flag[16].Flag = (float)json_object_get_number(BufObject,"16.e+O2>e+OP+O*");
      Coll_Flag[16].mofM = SP[0].mass/BG[1].mass;
      Coll_Flag[17].Flag = (float)json_object_get_number(BufObject,"17.e+O2>e+O^+O-");
      Coll_Flag[17].mofM = SP[0].mass/BG[1].mass;
      Coll_Flag[18].Flag = (float)json_object_get_number(BufObject,"18.e+O2>2e+O^+OP");
      Coll_Flag[18].mofM = SP[0].mass/BG[1].mass;
      Coll_Flag[19].Flag = (float)json_object_get_number(BufObject,"19.e+O2A>2e+O2^");
      Coll_Flag[19].mofM = SP[0].mass/FG[1].mass;
      Coll_Flag[20].Flag = (float)json_object_get_number(BufObject,"20.e+O2A>OP+O-");
      Coll_Flag[20].mofM = SP[0].mass/FG[1].mass;
      Coll_Flag[21].Flag = (float)json_object_get_number(BufObject,"21.e+O2A>e+O2");
      Coll_Flag[21].mofM = SP[0].mass/FG[1].mass;
      Coll_Flag[22].Flag = (float)json_object_get_number(BufObject,"22.e+O2A>e+O2");
      Coll_Flag[22].mofM = SP[0].mass/FG[1].mass;
      Coll_Flag[23].Flag = (float)json_object_get_number(BufObject,"23.e+O2A>e+2OP");
      Coll_Flag[23].mofM = SP[0].mass/FG[1].mass;
      Coll_Flag[24].Flag = (float)json_object_get_number(BufObject,"24.e+O2A>e+OP+OD");
      Coll_Flag[24].mofM = SP[0].mass/FG[1].mass;
      Coll_Flag[25].Flag = (float)json_object_get_number(BufObject,"25.e+O2A>e+2OD");
      Coll_Flag[25].mofM = SP[0].mass/FG[1].mass;
      Coll_Flag[26].Flag = (float)json_object_get_number(BufObject,"26.e+O2A>2e+O^+OP");
      Coll_Flag[26].mofM = SP[0].mass/FG[1].mass;
      Coll_Flag[27].Flag = (float)json_object_get_number(BufObject,"27.e+O2B>2e+O2^");
      Coll_Flag[27].mofM = SP[0].mass/FG[2].mass;
      Coll_Flag[28].Flag = (float)json_object_get_number(BufObject,"28.e+O2B>OP+O-");
      Coll_Flag[28].mofM = SP[0].mass/FG[2].mass;
      Coll_Flag[29].Flag = (float)json_object_get_number(BufObject,"29.e+O2B>e+O2");
      Coll_Flag[29].mofM = SP[0].mass/FG[2].mass;
      Coll_Flag[30].Flag = (float)json_object_get_number(BufObject,"30.e+O2B>e+O2");
      Coll_Flag[30].mofM = SP[0].mass/FG[2].mass;
      Coll_Flag[31].Flag = (float)json_object_get_number(BufObject,"31.e+O2B>e+2O");
      Coll_Flag[31].mofM = SP[0].mass/FG[2].mass;
      Coll_Flag[32].Flag = (float)json_object_get_number(BufObject,"32.e+O2B>e+OP+OD");
      Coll_Flag[32].mofM = SP[0].mass/FG[2].mass;
      Coll_Flag[33].Flag = (float)json_object_get_number(BufObject,"33.e+O2B>e+2OD");
      Coll_Flag[33].mofM = SP[0].mass/FG[2].mass;
      Coll_Flag[34].Flag = (float)json_object_get_number(BufObject,"34.e+O2B>2e+O^+OP");
      Coll_Flag[34].mofM = SP[0].mass/FG[2].mass;
      Coll_Flag[35].Flag = (float)json_object_get_number(BufObject,"35.e+O->2e+OP");
      Coll_Flag[35].mofM = SP[0].mass/SP[4].mass;
      Coll_Flag[36].Flag = (float)json_object_get_number(BufObject,"36.e+O2^>OP+OD");
      Coll_Flag[36].mofM = SP[0].mass/SP[2].mass;
      Coll_Flag[37].Flag = (float)json_object_get_number(BufObject,"37.e+OP>e+OP");
      Coll_Flag[37].mofM = SP[0].mass/FG[3].mass;
      Coll_Flag[38].Flag = (float)json_object_get_number(BufObject,"38.e+OP>e+OD");
      Coll_Flag[38].mofM = SP[0].mass/FG[3].mass;
      Coll_Flag[39].Flag = (float)json_object_get_number(BufObject,"39.e+OP>e+O*");
      Coll_Flag[39].mofM = SP[0].mass/FG[3].mass;
      Coll_Flag[40].Flag = (float)json_object_get_number(BufObject,"40.e+OP>e+O*");
      Coll_Flag[40].mofM = SP[0].mass/FG[3].mass;
      Coll_Flag[41].Flag = (float)json_object_get_number(BufObject,"41.e+OP>e+O*");
      Coll_Flag[41].mofM = SP[0].mass/FG[3].mass;
      Coll_Flag[42].Flag = (float)json_object_get_number(BufObject,"42.e+OP>2e+O^");
      Coll_Flag[42].mofM = SP[0].mass/FG[3].mass;
      Coll_Flag[43].Flag = (float)json_object_get_number(BufObject,"43.e+OP>e+O*");
      Coll_Flag[43].mofM = SP[0].mass/FG[3].mass;
      Coll_Flag[44].Flag = (float)json_object_get_number(BufObject,"44.e+OD>2e+O^");
      Coll_Flag[44].mofM = SP[0].mass/FG[4].mass;
      Coll_Flag[45].Flag = (float)json_object_get_number(BufObject,"45.e+OD>e+OP");
      Coll_Flag[45].mofM = SP[0].mass/FG[4].mass;
      Coll_Flag[46].Flag = (float)json_object_get_number(BufObject,"46.O-+O2>O-+O2");
      Coll_Flag[46].mofM = SP[4].mass/BG[1].mass;
      Coll_Flag[47].Flag = (float)json_object_get_number(BufObject,"47.O-+O2>e+OP+O2");
      Coll_Flag[47].mofM = SP[4].mass/BG[1].mass;
      Coll_Flag[48].Flag = (float)json_object_get_number(BufObject,"48.O-+OP>e+O2");
      Coll_Flag[48].mofM = SP[4].mass/FG[3].mass;
      Coll_Flag[49].Flag = (float)json_object_get_number(BufObject,"49.O-+O2^>OP+O2");
      Coll_Flag[49].mofM = SP[4].mass/SP[2].mass;
      Coll_Flag[50].Flag = (float)json_object_get_number(BufObject,"50.O-+O^>2OP");
      Coll_Flag[50].mofM = SP[4].mass/SP[3].mass;
      Coll_Flag[51].Flag = (float)json_object_get_number(BufObject,"51.O-+O2A>e+OP+O2");
      Coll_Flag[51].mofM = SP[4].mass/FG[1].mass;
      Coll_Flag[52].Flag = (float)json_object_get_number(BufObject,"52.O2^+OP>O2+O^");
      Coll_Flag[52].mofM = SP[2].mass/FG[3].mass;
      Coll_Flag[53].Flag = (float)json_object_get_number(BufObject,"53.O2^+O2>O2+O2^");
      Coll_Flag[53].mofM = SP[2].mass/BG[1].mass;
      Coll_Flag[54].Flag = (float)json_object_get_number(BufObject,"54.O2^+O2>O2^+O2");
      Coll_Flag[54].mofM = SP[2].mass/BG[1].mass;
      Coll_Flag[55].Flag = (float)json_object_get_number(BufObject,"55.O2^+O2>O^+OP+O2");
      Coll_Flag[55].mofM = SP[2].mass/BG[1].mass;
      Coll_Flag[56].Flag = (float)json_object_get_number(BufObject,"56.O2^+O2A>O2+O2^");
      Coll_Flag[56].mofM = SP[2].mass/FG[1].mass;
      Coll_Flag[57].Flag = (float)json_object_get_number(BufObject,"57.O2^+O2B>O2+O2^");
      Coll_Flag[57].mofM = SP[2].mass/FG[2].mass;
      Coll_Flag[58].Flag = (float)json_object_get_number(BufObject,"58.O2^+Ar>O2+Ar^");
      Coll_Flag[58].mofM = SP[2].mass/BG[0].mass;
      Coll_Flag[59].Flag = (float)json_object_get_number(BufObject,"59.O2^+Ar>O2^+Ar^");
      Coll_Flag[59].mofM = SP[2].mass/BG[0].mass;
      Coll_Flag[60].Flag = (float)json_object_get_number(BufObject,"60.O^+O2>OP+O2^");
      Coll_Flag[60].mofM = SP[3].mass/BG[1].mass;
      Coll_Flag[61].Flag = (float)json_object_get_number(BufObject,"61.O^+O2>O^+O2");
      Coll_Flag[61].mofM = SP[3].mass/BG[1].mass;
      Coll_Flag[62].Flag = (float)json_object_get_number(BufObject,"62.O^+OP>OP+O^");
      Coll_Flag[62].mofM = SP[3].mass/FG[3].mass;
      Coll_Flag[63].Flag = (float)json_object_get_number(BufObject,"63.O^+O2A>O2^+OP");
      Coll_Flag[63].mofM = SP[3].mass/FG[1].mass;
      Coll_Flag[64].Flag = (float)json_object_get_number(BufObject,"64.O^+O2B>O2^+OP");
      Coll_Flag[64].mofM = SP[3].mass/FG[2].mass;
      Coll_Flag[65].Flag = (float)json_object_get_number(BufObject,"65.Ar+Ar^>Ar^+Ar");
      Coll_Flag[65].mofM = BG[0].mass/SP[1].mass;
      Coll_Flag[66].Flag = (float)json_object_get_number(BufObject,"66.Ar^+Ar>Ar^+Ar");
      Coll_Flag[66].mofM = SP[1].mass/BG[0].mass;
      Coll_Flag[67].Flag = (float)json_object_get_number(BufObject,"67.Ar^+O2>O2+Ar^");
      Coll_Flag[67].mofM = SP[1].mass/BG[1].mass;
      Coll_Flag[68].Flag = (float)json_object_get_number(BufObject,"68.e+Ar*>e+Ar");
      Coll_Flag[68].mofM = SP[0].mass/FG[0].mass;
      Coll_Flag[69].Flag = (float)json_object_get_number(BufObject,"69.O-+Ar^>OP+AR");
      Coll_Flag[69].mofM = SP[4].mass/SP[1].mass;
      Coll_Flag[70].Flag = (float)json_object_get_number(BufObject,"70.O-+O2B>e+OP+O2");
      Coll_Flag[70].mofM = SP[4].mass/FG[2].mass;
      Coll_Flag[71].Flag = (float)json_object_get_number(BufObject,"71.OP+OD>2OP");
      Coll_Flag[71].mofM = FG[3].mass/FG[4].mass;
      Coll_Flag[72].Flag = (float)json_object_get_number(BufObject,"72.OD+O2>OP+O2");
      Coll_Flag[72].mofM = FG[4].mass/BG[1].mass;
      Coll_Flag[73].Flag = (float)json_object_get_number(BufObject,"73.OD+O2>OP+O2A");
      Coll_Flag[73].mofM = FG[4].mass/BG[1].mass;
      Coll_Flag[74].Flag = (float)json_object_get_number(BufObject,"74.OD+O2>OP+O2B");
      Coll_Flag[74].mofM = FG[4].mass/BG[1].mass;
      Coll_Flag[75].Flag = (float)json_object_get_number(BufObject,"75.O2A+OP>OP+O2");
      Coll_Flag[75].mofM = FG[1].mass/FG[3].mass;
      Coll_Flag[76].Flag = (float)json_object_get_number(BufObject,"76.O2A+O2>2O2");
      Coll_Flag[76].mofM = FG[1].mass/BG[1].mass;
      Coll_Flag[77].Flag = (float)json_object_get_number(BufObject,"77.O2A+O2A>2O2");
      Coll_Flag[77].mofM = FG[1].mass/FG[1].mass;
      Coll_Flag[78].Flag = (float)json_object_get_number(BufObject,"78.O2B+O2>2O2");
      Coll_Flag[78].mofM = FG[2].mass/BG[1].mass;
      Coll_Flag[79].Flag = (float)json_object_get_number(BufObject,"79.Ar^+OP>Ar+O^");
      Coll_Flag[79].mofM = SP[1].mass/FG[3].mass;
      Coll_Flag[80].Flag = (float)json_object_get_number(BufObject,"80.Ar^+O2>Ar+O2^");
      Coll_Flag[80].mofM = SP[1].mass/BG[1].mass;
      Coll_Flag[81].Flag = (float)json_object_get_number(BufObject,"81.Ar*+Ar*>e+Ar+Ar^");
      Coll_Flag[81].mofM = FG[0].mass/FG[0].mass;
      Coll_Flag[82].Flag = (float)json_object_get_number(BufObject,"82.Ar*+Ar>2Ar");
      Coll_Flag[82].mofM = FG[0].mass/BG[0].mass;
      Coll_Flag[83].Flag = (float)json_object_get_number(BufObject,"83.Ar*+OP>OD+Ar");
      Coll_Flag[83].mofM = FG[0].mass/FG[3].mass;
      Coll_Flag[84].Flag = (float)json_object_get_number(BufObject,"84.Ar*+OP>OP+Ar");
      Coll_Flag[84].mofM = FG[0].mass/FG[3].mass;
      Coll_Flag[85].Flag = (float)json_object_get_number(BufObject,"85.Ar*+O2>2OP+Ar");
      Coll_Flag[85].mofM = FG[0].mass/BG[1].mass;
      Coll_Flag[86].Flag = (float)json_object_get_number(BufObject,"86.Ar*+O2>OP+OD+Ar");
      Coll_Flag[86].mofM = FG[0].mass/BG[1].mass;
      Coll_Flag[87].Flag = (float)json_object_get_number(BufObject,"87.Ar*+O2>O2+Ar");
      Coll_Flag[87].mofM = FG[0].mass/BG[1].mass;
      printf("\tReaction = %d\n",TnRct);
      printf("\tCross section = %d\n",nRct_cx);
      printf("\tReaction rate = %d\n",nRct_rc);
      printf("\tCheck reaction on or off...\n");
      for(i=0;i<TnRct;i++){
         if(Coll_Flag[i].Flag < 0)
            Coll_Flag[i].Flag = 0.0;
         if(Coll_Flag[i].Flag == 0)
            printf("\tReaction %d is off\n",i);
      }
      ArO2_CrossSectionSET(Coll_Flag);
      printf("\tReaction setting Complete.\n");
   }else{
      exit(1);
   }
   printf("Finish parsing \"%s\"\n",InputFile); 
}
void Geometry_setting() {
   int i,j,k,ID,CID,SID;
   printf("Initial-setting start\n"); 
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
   //
   vec_G = (HGA *) malloc(Gsize * sizeof(HGA));
   for(i=0;i<Gsize;i++){
      vec_G[i].Boundary = (int) 0;
      vec_G[i].CondID = (int) 0;
      vec_G[i].Temp = (float) 0.0;
      vec_G[i].BackDens = (float) 0.0;
      vec_G[i].face = NO_FACE;
      vec_G[i].area = (float) dy*zlength;
   }
   vec_C = (HCA *) malloc(Csize * sizeof(HCA));
   for(i=0;i<Csize;i++){
      vec_C[i].PlasmaRegion = (int) 1;
      vec_C[i].eps_r = (float) 1.0;
   }
   for (k=0;k<BoundaryNUM;k++){
      for (i=BoundaryX0[k];i<=BoundaryX1[k];i++){
         for (j=BoundaryY0[k];j<=BoundaryY1[k];j++){
            ID = i*ngy+j;
            vec_G[ID].Temp = BoundaryTEMP[k];
            if(vec_G[ID].Boundary != DIRICHLET)
               vec_G[ID].Boundary = BoundaryBC[k];
         }
      }
   }
   for (k=0;k<DielNUM;k++){
      for (i=DielX0[k];i<=DielX1[k];i++){
         for (j=DielY0[k];j<=DielY1[k];j++){
            ID = i*ngy+j;
            if(vec_G[ID].Boundary == 0)
               vec_G[ID].Boundary = DIELECTRIC;
            if(i < DielX1[k] && j <DielY1[k]){
               CID = i*ncy + j;
               vec_C[CID].PlasmaRegion = (int) 0;
               vec_C[CID].eps_r = DielEPS[k];
            }
         }
      }
   }
   for (k=0;k<CondNUM;k++){
      for (i=CondX0[k];i<=CondX1[k];i++){
         for (j=CondY0[k];j<=CondY1[k];j++){
            ID = i*ngy + j;
            vec_G[ID].CondID = CondM_ID[k];
            if(vec_G[ID].Boundary == 0 || vec_G[ID].Boundary == DIELECTRIC)
               vec_G[ID].Boundary = CONDUCTOR;
            if(i < CondX1[k] && j <CondY1[k]){
               CID = i*ncy + j;
               vec_C[CID].PlasmaRegion = (int) 0;
               vec_C[CID].eps_r = (float) 0.0;
            }
         }
      }
   }
   int P1,P2,P3,P4;
   for(i=0+1; i<ngx-1; i++){
      for(j=0+1; j<ngy-1; j++){
         ID = i*ngy+j;
         CID = i*ncy+j;
         // P2 | P1
         // ---+----
         // P4 | P3
         P1 = vec_C[CID].eps_r;
         P2 = vec_C[CID-ncy].eps_r;
         P3 = vec_C[CID-1].eps_r;
         P4 = vec_C[CID-ncy-1].eps_r;
         if(vec_G[ID].Boundary == DIELECTRIC){
            if((P2==1.0 && P1==1.0 && P4==1.0) || (P2!=0.0 && P1!=0.0 && P4!=0.0 && P3==0.0)) {
					//          //
					//          //
					//     -----//
					//     |----//
					//     |----//               
               vec_G[ID].face=UL_CORN;
               vec_G[ID].area= 0.5*(dx+dy)*zlength;
            }else if((P2==1.0 && P1==1.0 && P3==1.0) || (P2!=0.0 && P1!=0.0 && P3!=0.0 && P4==0.0)) {
					//          //
					//          //
					//------    //
					//-----|    //
					//-----|    //               
               vec_G[ID].face=UR_CORN;
               vec_G[ID].area= 0.5*(dx+dy)*zlength;
            }else if((P2==1.0 && P1==1.0 && P4!=1.0 && P3!=1.0) || (P2!=0.0 && P1!=0.0 && P4==0.0 && P3==0.0)) {
					//          //
					//          //
					//----------//
					//----------//
					//----------//
               vec_G[ID].face=UP;
               vec_G[ID].area= dx*zlength;
            }else if((P4==1.0 && P3==1.0 && P2==1.0) || (P2!=0.0 && P4!=0.0 && P3!=0.0 && P1==0.0)) {
					//     |----//
					//     |----//
					//     |----//
					//          //
					//          //               
               vec_G[ID].face=LL_CORN;
               vec_G[ID].area= 0.5*(dx+dy)*zlength;
            }else if((P4==1.0 && P3==1.0 && P1==1.0) || (P1!=0.0 && P4!=0.0 && P3!=0.0 && P2==0.0)) {
					//----|     //
					//----|     //
					//----|     //
					//          //
					//          //               
               vec_G[ID].face=LR_CORN;
               vec_G[ID].area= 0.5*(dx+dy)*zlength;
            }else if((P4==1.0 && P3==1.0 && P1!=1.0 && P2!=1.0) || (P4!=0.0 && P3!=0.0 && P1==0.0 && P2==0.0)) {
					//----------//
					//----------//
					//----------//
					//          //
					//          //               
               vec_G[ID].face=DOWN;
               vec_G[ID].area= dx*zlength;
            }else if((P4==1.0 && P2==1.0 && P1!=1.0 && P3!=1.0) || (P4!=0.0 && P2!=0.0 && P1==0.0 && P3==0.0)) {
					//     |----//
					//     |----//
					//     |----//
					//     |----//
					//     |----//               
               vec_G[ID].face=LEFT;
               vec_G[ID].area= dy*zlength;
            }else if((P4!=1.0 && P2!=1.0 && P1==1.0 && P3==1.0) || (P4==0.0 && P2==0.0 && P1!=0.0 && P3!=0.0)) {
					//----|     //
					//----|     //
					//----|     //
					//----|     //
					//----|     //               
               vec_G[ID].face=RIGHT;
               vec_G[ID].area= dy*zlength;
            }  
         }
         if(vec_G[ID].Boundary == CONDUCTOR){
            if((P2==1.0 && P1==1.0 && P4==1.0)||(P2!=0.0 && P1!=0.0 && P4!=0.0 && P3==0.0)) {
               vec_G[ID].face=UL_CORN;
               vec_G[ID].area= 0.5*(dx+dy)*zlength;
            }
            else if((P2==1.0 && P1==1.0 && P3==1.0)||(P2!=0.0 && P1!=0.0 && P3!=0.0 && P4==0.0)) {
               vec_G[ID].face=UR_CORN;
               vec_G[ID].area= 0.5*(dx+dy)*zlength;
            }
            else if((P2==1.0 && P1==1.0 && P4!=1.0 && P3!=1.0)||(P2!=0.0 && P1!=0.0 && P4==0.0 && P3==0.0)) {
               vec_G[ID].face=UP;
               vec_G[ID].area= dx*zlength;
            }
            else if((P4==1.0 && P3==1.0 && P2==1.0)||(P2!=0.0 && P4!=0.0 && P3!=0.0 && P1==0.0)) {
               vec_G[ID].face=LL_CORN;
               vec_G[ID].area= 0.5*(dx+dy)*zlength;
            }
            else if((P4==1.0 && P3==1.0 && P1==1.0)||(P1!=0.0 && P4!=0.0 && P3!=0.0 && P2==0.0)) {
               vec_G[ID].face=LR_CORN;
               vec_G[ID].area= 0.5*(dx+dy)*zlength;
            }
            else if((P4==1.0 && P3==1.0 && P1!=1.0 && P2!=1.0)||(P4!=0.0 && P3!=0.0 && P1==0.0 && P2==0.0)) {
               vec_G[ID].face=DOWN;
               vec_G[ID].area= dx*zlength;
            }
            else if((P4==1.0 && P2==1.0 && P1!=1.0 && P3!=1.0)||(P4!=0.0 && P2!=0.0 && P1==0.0 && P3==0.0)) {
               vec_G[ID].face=LEFT;
               vec_G[ID].area= dy*zlength;
            }
            else if((P4!=1.0 && P2!=1.0 && P1==1.0 && P3==1.0)||(P4==0.0 && P2==0.0 && P1!=0.0 && P3!=0.0)) {
               vec_G[ID].face=RIGHT;
               vec_G[ID].area= dy*zlength;
            }
         }
      }
   }
   for(j=1;j<ncy;j++) {
      ID = j;
      if(vec_G[ID].Boundary==DIRICHLET) {
         vec_G[ID].face=RIGHT;
         vec_G[ID].area=dy*zlength;
         if(vec_G[ID+1].CondID==0 && vec_G[ID-1].CondID!=0 && vec_G[ID].CondID!=0) {
            vec_G[ID].face=UP;
            vec_G[ID].area=dx*zlength;
         }else if(vec_G[ID+1].CondID!=0 && vec_G[ID-1].CondID==0 && vec_G[ID].CondID!=0) {
            vec_G[ID].face=DOWN;
            vec_G[ID].area=dx*zlength;
         }else if(vec_G[ID+1].CondID!=0 && vec_G[ID-1].CondID!=0 && vec_G[ID].CondID!=0) {
            vec_G[ID].face=NO_FACE;
            vec_G[ID].area=dy*zlength;
         }else if(vec_G[ID+1].CondID==0 && vec_G[ID-1].CondID==0 && vec_G[ID].CondID!=0) {
            vec_G[ID].face=NO_FACE;
            vec_G[ID].area=dy*zlength;
         }
      }else if(vec_G[ID].Boundary==NEUMANN) {
         vec_G[ID].Boundary=NO_FACE;
         vec_G[ID].area=dy*zlength;
         if(vec_G[ID+1].CondID==0 && vec_G[ID-1].CondID!=0 && vec_G[ID].CondID!=0) {
            vec_G[ID].face=UP;
            vec_G[ID].area=dx*zlength;
         }else if(vec_G[ID+1].CondID!=0 && vec_G[ID-1].CondID==0 && vec_G[ID].CondID!=0) {
            vec_G[ID].face=DOWN;
            vec_G[ID].area=dx*zlength;
         }else if(vec_G[ID+1].CondID!=0 && vec_G[ID-1].CondID!=0 && vec_G[ID].CondID==0) {
            vec_G[ID].face=RIGHT;
            vec_G[ID].area=dy*zlength;
         }
      }
      ID = ncx*ngy+j;
      if(vec_G[ID].Boundary==DIRICHLET) {
         vec_G[ID].face=LEFT;
         vec_G[ID].area=dy*zlength;
         if(vec_G[ID+1].CondID==0 && vec_G[ID-1].CondID!=0 && vec_G[ID-ngy].CondID!=0) {
            vec_G[ID].face=UP;
            vec_G[ID].area=dx*zlength;
         }else if(vec_G[ID+1].CondID!=0 && vec_G[ID-1].CondID==0 && vec_G[ID-ngy].CondID!=0) {
            vec_G[ID].face=DOWN;
            vec_G[ID].area=dx*zlength;
         }else if(vec_G[ID+1].CondID!=0 && vec_G[ID-1].CondID!=0 && vec_G[ID-ngy].CondID!=0) {
            vec_G[ID].face=NO_FACE;
            vec_G[ID].area=dy*zlength;
         }else if(vec_G[ID+1].CondID==0 && vec_G[ID-1].CondID==0 && vec_G[ID-ngy].CondID!=0) {
            vec_G[ID].face=NO_FACE;
            vec_G[ID].area=dy*zlength;
         }
      }else if(vec_G[ID].Boundary==NEUMANN) {
         vec_G[ID].face=NO_FACE;
         vec_G[ID].area=dy*zlength;
         if(vec_G[ID+1].CondID==0 && vec_G[ID-1].CondID!=0 && vec_G[ID-ngy].CondID!=0) {
            vec_G[ID].face=UP;
            vec_G[ID].area=dx*zlength;
         }else if(vec_G[ID+1].CondID!=0 && vec_G[ID-1].CondID==0 && vec_G[ID-ngy].CondID!=0) {
            vec_G[ID].face=DOWN;
            vec_G[ID].area=dx*zlength;
         }else if(vec_G[ID+1].CondID!=0 && vec_G[ID-1].CondID!=0 && vec_G[ID-ngy].CondID==0) {
            vec_G[ID].face=LEFT;
            vec_G[ID].area=dy*zlength;
         }
      }
   }
   for(i=1;i<ncx;i++) {
      ID = i*ngy;
      if(vec_G[ID].Boundary==DIRICHLET) {
         vec_G[ID].face=UP;
         vec_G[ID].area=dx*zlength;
         if(vec_G[ID+ngy].CondID==0 && vec_G[ID-ngy].CondID!=0 && vec_G[ID].CondID!=0) {
            vec_G[ID].face=RIGHT;
            vec_G[ID].area=dy*zlength;
         }else if(vec_G[ID+ngy].CondID!=0 && vec_G[ID-ngy].CondID==0 && vec_G[ID].CondID!=0) {
            vec_G[ID].face=LEFT;
            vec_G[ID].area=dy*zlength;
         }else if(vec_G[ID+ngy].CondID!=0 && vec_G[ID-ngy].CondID!=0 && vec_G[ID].CondID!=0) {
            vec_G[ID].face=NO_FACE;
            vec_G[ID].area=dy*zlength;
         }else if(vec_G[ID+ngy].CondID==0 && vec_G[ID-ngy].CondID==0 && vec_G[ID].CondID!=0) {
            vec_G[ID].face=NO_FACE;
            vec_G[ID].area=dy*zlength;
         }
      }else if(vec_G[ID].Boundary==NEUMANN) {
         vec_G[ID].face=NO_FACE;
         vec_G[ID].area=dy*zlength;
         if(vec_G[ID+ngy].CondID==0 && vec_G[ID-ngy].CondID!=0 && vec_G[ID].CondID!=0) {
            vec_G[ID].face=RIGHT;
            vec_G[ID].area=dy*zlength;
         }else if(vec_G[ID+ngy].CondID!=0 && vec_G[ID-ngy].CondID==0 && vec_G[ID].CondID!=0) {
            vec_G[ID].face=LEFT;
            vec_G[ID].area=dy*zlength;
         }else if(vec_G[ID+ngy].CondID!=0 && vec_G[ID-ngy].CondID!=0 && vec_G[ID].CondID==0) {
            vec_G[ID].face=UP;
            vec_G[ID].area=dx*zlength;
         }
      }
      ID = i*ngy+ncy;
      if(vec_G[ID].Boundary==DIRICHLET) {
         vec_G[ID].face=DOWN;
         vec_G[ID].area=dx*zlength;
         if(vec_G[ID+ngy].CondID==0 && vec_G[ID-ngy].CondID!=0 && vec_G[ID-1].CondID!=0) {
            vec_G[ID].face=RIGHT;
            vec_G[ID].area=dy*zlength;
         }else if(vec_G[ID+ngy].CondID!=0 && vec_G[ID-ngy].CondID==0 && vec_G[ID-1].CondID!=0) {
            vec_G[ID].face=LEFT;
            vec_G[ID].area=dy*zlength;
         }else if(vec_G[ID+ngy].CondID!=0 && vec_G[ID-ngy].CondID!=0 && vec_G[ID-1].CondID!=0) {
            vec_G[ID].face=NO_FACE;
            vec_G[ID].area=dy*zlength;
         }else if(vec_G[ID+ngy].CondID==0 && vec_G[ID-ngy].CondID==0 && vec_G[ID-1].CondID!=0) {
            vec_G[ID].face=NO_FACE;
            vec_G[ID].area=dy*zlength;
         }
      }else if(vec_G[ID].Boundary==NEUMANN) {
         vec_G[ID].face=NO_FACE;
         vec_G[ID].area=dy*zlength;
         if(vec_G[ID+ngy].CondID==0 && vec_G[ID-ngy].CondID!=0 && vec_G[ID-1].CondID!=0) {
            vec_G[ID].face=RIGHT;
            vec_G[ID].area=dy*zlength;
         }else if(vec_G[ID+ngy].CondID!=0 && vec_G[ID-ngy].CondID==0 && vec_G[ID-1].CondID!=0) {
            vec_G[ID].face=LEFT;
            vec_G[ID].area=dy*zlength;
         }else if(vec_G[ID+ngy].CondID!=0 && vec_G[ID-ngy].CondID!=0 && vec_G[ID-1].CondID==0) {
            vec_G[ID].face=DOWN;
            vec_G[ID].area=dx*zlength;
         }
      }
   }
   ID = 0;
   if(vec_G[ID].Boundary==DIRICHLET  && vec_G[ID+1].CondID==0) {
      vec_G[ID].face=UR_CORN;
      vec_G[ID].area= 0.5*(dx+dy)*zlength;
   }
   ID = ncx*ngy;
   if(vec_G[ID].Boundary==DIRICHLET  && vec_G[ID+1].CondID==0) {
      vec_G[ID].face=UL_CORN;
      vec_G[ID].area= 0.5*(dx+dy)*zlength;
   }
   ID = ncy;
   if(vec_G[ID].Boundary==DIRICHLET  && vec_G[ID+1].CondID==0) {
      vec_G[ID].face=LR_CORN;
      vec_G[ID].area= 0.5*(dx+dy)*zlength;
   }
   ID = ncx*ngy+ncy;
   if(vec_G[ID].Boundary==DIRICHLET  && vec_G[ID+1].CondID==0) {
      vec_G[ID].face=LL_CORN;
      vec_G[ID].area= 0.5*(dx+dy)*zlength;
   }
   // Set Structure for GPU
    // make a Set Structure Index
	// 1~99:dielectric 0:Plasma 100~:conductor //-1~:boundary
	StructureIndex = MIMalloc((ncx + 2),(ncy + 2));
   MIInit(StructureIndex, 0, ncx + 2, ncy + 2);
	vec_StructureIndex = VIMalloc((ncx + 2) * (ncy + 2));
	VIInit(vec_StructureIndex, 0, (ncx + 2) * (ncy + 2));
   //Material
	for (i = 0; i < ncx; i++) {
		for (j = 0; j < ncy; j++) {
			StructureIndex[i + 1][j + 1] = 0;
			for (k = 0; k < DielNUM; k++)
				if (i >= DielX0[k] && i < DielX1[k] && j >= DielY0[k] && j < DielY1[k])
					StructureIndex[i + 1][j + 1] = DielM_ID[k];
			for (k = 0; k < CondNUM; k++)
				if (i >= CondX0[k] && i < CondX1[k] && j >= CondY0[k] && j < CondY1[k])
					StructureIndex[i + 1][j + 1] = 100 + CondM_ID[k];
		}
	}
   // BOUNDARY check excluding edges
	// -2, -4, -6, -8 : Neumann
	// -12, -16, -20, -24 : Neumann edge
	for (i = 1; i < ngx; i++) {
      // Bottom boundary check, excluding edges
		if ((vec_G[(i-1)*ngy].Boundary == NEUMANN || vec_G[i*ngy].Boundary == NEUMANN)&& (StructureIndex[i][1] == 0)) {
         StructureIndex[i][0] = -vec_G[(i-1)*ngy].Boundary - 2;
		} else {
			if (StructureIndex[i][1]) { // material
				StructureIndex[i][0] = StructureIndex[i][1];
			} else
            // Dirichlet Boundary
				StructureIndex[i][0] = -1;
		}
      // Top boundary check, excluding edges
		if ((vec_G[(i-1)*ngy+ncy].Boundary == NEUMANN|| vec_G[i*ngy+ncy].Boundary == NEUMANN)&& (StructureIndex[i][ngy - 1] == 0)) {
			StructureIndex[i][ngy] = -vec_G[(i-1)*ngy+ncy].Boundary - 6;
		} else {
			if (StructureIndex[i][ngy - 1]) {
				StructureIndex[i][ngy] = StructureIndex[i][ngy - 1];
			} else
				StructureIndex[i][ngy] = -1;
		}
	}
	for (j = 1; j < ngy; j++) {
      // Left boundary check, excluding edges
		if ((vec_G[j-1].Boundary == NEUMANN
				|| vec_G[j].Boundary == NEUMANN)
				&& (StructureIndex[1][j] == 0)) {
			StructureIndex[0][j] = -vec_G[j-1].Boundary;
		} else {
			if (StructureIndex[1][j]) {
				StructureIndex[0][j] = StructureIndex[1][j];
			} else
				StructureIndex[0][j] = -1;
		}
      // Right boundary check, excluding edges
		if ((vec_G[ncx*ngy+j-1].Boundary == NEUMANN
				|| vec_G[ncx*ngy+j].Boundary == NEUMANN)
				&& (StructureIndex[ngx - 1][j] == 0)) {
			StructureIndex[ngx][j] = -vec_G[ncx*ngy+j].Boundary - 4;
		} else {
			if (StructureIndex[ngx - 1][j]) {
				StructureIndex[ngx][j] = StructureIndex[ngx - 1][j];
			} else
				StructureIndex[ngx][j] = -1;
		}
	}
   // Boundary edge check
	StructureIndex[0][0] = -vec_G[ngy].Boundary - vec_G[1].Boundary - 8;
	StructureIndex[ncx + 1][0] = -vec_G[(ncx-1)*ngy].Boundary - vec_G[ncx*ngy+1].Boundary - 12;
	StructureIndex[0][ncy + 1] = -vec_G[ngy+ncy].Boundary	- vec_G[ncy-1].Boundary - 16;
	StructureIndex[ncx + 1][ncy + 1] = -vec_G[(ncy-1)*ngy+ncy].Boundary - vec_G[(ncx)*ngy+ncy-1].Boundary- 20;
	for (i = 0; i < ncx + 2; i++)
		for (j = 0; j < ncy + 2; j++)
			vec_StructureIndex[i * (ncy + 2) + j] = StructureIndex[i][j];
   /*
   for(j=ncy-1;j>=0;j--){
      for(i=0;i<ncx;i++){
         CID = i*ncy + j;
         printf("%d",vec_C[CID].PlasmaRegion);
      }
      printf("\n");
   }
   exit(1);
   */
  /*
   for(j=ngy-1;j>=0;j--){
      for(i=0;i<ngx;i++){
         ID = i*ngy + j;
         printf("%d",vec_G[ID].Boundary);
      }
      printf("\n");
   }printf("\n");
   exit(1);
   */
   /*
   for(j=ngy;j>=0;j--){
      for(i=0;i<ngx+1;i++){
         printf("%d",vec_StructureIndex[i * (ncy + 2) + j]);
      }
      printf("\n");
   }printf("\n");
   exit(1);
   */
   //exit(1); // devolep point

}
void FieldSolverSetting(){
   exit(1);
}
void GasSetting(){

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

