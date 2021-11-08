#include "Efield.h"

void Source_setting(){
   int i,j;
   int buf,buf1,TScnt,S_ID;
   printf("Initial Source setting\n"); 
   dt = 1.0 / Max_FREQ / (float)DT_PIC;
   dtc = dt * (float)DT_CONTI;
   CYCLE_NUM =  DT_PIC;
   DT_MCCn = 0;
   dt_mcc = 0.0;
   dt_dx = dt/dx;
   dt_dy = dt/dy;
   printf("\tPIC TimeStep = %g (s)\n",dt);
   printf("\tContinuity TimeStep = %g (s)\n",dtc);
   printf("\tMinimum Freq Cycle step # = %d (#)\n",CYCLE_NUM);
   if(CYCLE_NUM == 0){
      printf("Error : Minimum Freq Cycle step.\n");
      exit(1);
   }
   Efield_Flag = VIMalloc(CondNUMR); // 0:ground 1~ : Source ID
   VIInit(Efield_Flag,0,CondNUMR);
   Cond_Source_num = VIMalloc(CondNUMR);  // Number of Fix power frequency 
   VIInit(Cond_Source_num,0,CondNUMR);
   Cond_count = VIMalloc(CondNUMR);
   VIInit(Cond_count,0,CondNUMR);
   Cond_Power = VFMalloc(CondNUMR);
   VFInit(Cond_Power,0.0f,CondNUMR);
   Cond_Power_ID = MIMalloc(CondNUMR,2);
   MIInit(Cond_Power_ID,0,CondNUMR,2);
   //Find Source num each of Conductor
   TScnt = 0;
   for (i=0;i<CondNUMR;i++){
      buf = 0;
      buf1 = 0;
      for (j=0;j<SrcNUM;j++){
         if(i+1==SrcM_ID[j]){
            buf++;
            if(SrcPOWER[j] !=0 ){
               if(buf1 == 2){
                  printf("Error! : There should be no more than 3 different Power each conductor.\n"); 
                  exit(1);
               }
               Cond_Power_ID[i][buf1] = j;
               buf1++;
            }
         }
      }
      if(buf == 0){
         printf("\tConductor ID[%d] : Grounded Electrode",i+1);
         Efield_Flag[i] = 0;
      }else if(buf < 4){
         printf("\tConductor ID[%d] : Powered Electrode",i+1);
         TScnt++;
         Efield_Flag[i] = 1;
      }
      if(buf1 == 0){
         printf(", Voltage Driven \n");
      }else if(buf1 == 1){
         printf(", Single Power Driven \n");
      }else if(buf1 == 2){
         printf(", Dual Power Driven \n");
      }else{
         printf("Error! : There should be no more than 3 different frequencies each conductor.\n"); 
         exit(1);
      }
      Cond_Source_num[i] = buf1;
   }
   printf("\tTotal Number of Source set : [%d]\n",TScnt);
   if(External_Flag){
      CC_a  = MFMalloc(CondNUMR,5);
      for (i=0;i<CondNUM;i++){
         CC_a[CondM_ID[i]-1][0] = 2.25*CondL[i]/dt/dt + 1.5*CondR[i]/dt + 1/CondC[i];
         CC_a[CondM_ID[i]-1][1] = -6*CondL[i]/dt/dt - 2*CondR[i]/dt;
         CC_a[CondM_ID[i]-1][2] = 5.5*CondL[i]/dt/dt + .5*CondR[i]/dt;
         CC_a[CondM_ID[i]-1][3] = -2*CondL[i]/dt/dt;
         CC_a[CondM_ID[i]-1][4] = .25*CondL[i]/dt/dt;
      }
   }
}