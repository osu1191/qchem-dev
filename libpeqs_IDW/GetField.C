#include "peqs_class.h"

void PEQSclass::Get_Electric_Field_on_ACG(double *PAv)
{
  int i,j,k,ijk,l,m;
  double sum1=0.0,sum2=0.0,sum3=0.0,sumC1=0.0,sumC2=0.0;
  INTEGER NBasis = bSetMgr.crntShlsStats(STAT_NBASIS);
  INTEGER N2 = NBasis*NBasis;
  double *kPts2 = QAllocDouble(count*3);  // for ACG testing
  double *kWts = QAllocDouble(count);   // for ACG testing
  FileMan(FM_READ,FILE_PEQS_GRID,FM_DP,count*3,0,FM_BEG,kPts2);  // for ACG testing
  FileMan(FM_READ,FILE_ATM_BATCH_WTS,FM_DP,count,0,FM_BEG,kWts);   // for ACG testing

  double *jDMA = QAllocDouble(N2);
  double *jDMB = QAllocDouble(N2);
  FileMan(FM_READ,FILE_DENSITY_MATRIX,FM_DP,N2,0,FM_BEG,jDMA);
  FileMan(FM_READ,FILE_DENSITY_MATRIX,FM_DP,N2,N2,FM_BEG,jDMB);
  double *jChiFree, *aChi;
  aChi = QAllocDoubleWithInit(NBasis*count);
  EvlBasis(aChi,NBasis,kPts2,count);

  double *RhoA = QAllocDoubleWithInit(count);
  double *RhoB = QAllocDoubleWithInit(count);
  double *RhoACG = QAllocDoubleWithInit(count);
  double *unwtRhoACG = QAllocDoubleWithInit(count);
  VRload(RhoACG,count,0.0);
  VRload(unwtRhoACG,count,0.0);

  for(j = 0; j < count; j++){
     for (l = 0; l < NBasis; l++) {
       for (m = 0; m < NBasis; m++) {
       RhoA[j] += jDMA[l*NBasis+m]*aChi[l*count+j]*aChi[m*count+j];
       RhoB[j] += jDMB[l*NBasis+m]*aChi[l*count+j]*aChi[m*count+j];
       }
     }
  }
  for(j = 0; j < count; j++){
       sum2 += kWts[j]*RhoA[j];   // Integration on ACG
       sum2 += kWts[j]*RhoB[j];   // Integration on ACG
       RhoACG[j] = kWts[j]*(RhoA[j] + RhoB[j]);
       unwtRhoACG[j] = (RhoA[j] + RhoB[j]);
  }

  qtime_t shackle08;
  double shackle09[3];
  shackle08 = QTimerOn();
  cout << "Will calculate potential now" << endl;
  int cycles=MaxPPer;
  cout << "cycles in elec = " << cycles << endl;
  int acgPPer=cycles, offset=0;
  double *irregPts=QAllocDouble(3*cycles); //,*acg_phi=QAllocDouble(cycles);
  double *AField=QAllocDouble(3*cycles);
  for(l=0;l<count;l+=cycles){
  if (cycles <= (count-l)) acgPPer = cycles;
  else acgPPer = (count-l);
    VRload(AField,3*acgPPer,0.0);
    VRcopy(irregPts,&kPts2[(3*l)],3*acgPPer);   // THIS IS THE DOUBTFUL PART... you have to copy all the ACG points using this line!!!
    AOints(AField,NULL,NULL,NULL,PAv,irregPts,NULL,NULL,&acgPPer,112);   // EField on ACG pts
    FileMan(FM_WRITE,FILE_PEQS_GRID,FM_DP,3*acgPPer,count*3+offset,FM_BEG,AField);
    offset += 3*acgPPer;
  }
  QFree(AField);
  cout << "Electric field have been calculated" << endl;

  QTimerOff(shackle09,shackle08);
  printf(" Field calculation in ACG points completed in %2.2f (cpu) %2.2f (wall) seconds. \n",shackle09[0],shackle09[2]);

  QFree(jDMA), QFree(jDMB);
  QFree(RhoA), QFree(RhoB), QFree(RhoACG);
  QFree(aChi);
  QFree(kPts2); QFree(kWts);
  QFree(unwtRhoACG);
}

void PEQSclass::GetDensityFromEfield()
{
   int ijk,i,j,k,Nx=NPoints[0],Ny=NPoints[1],Nz=NPoints[2];
   double dx,dy,dz;
   dx = delX[0];
   dy = delX[1];
   dz = delX[2];

   VRscale(valEField,3*NTotPts,-1.0/(4.0*Pi));
   VRscale(valGf,3*NTotPts,-1.0/(Pi));
   VRload(rho_from_field,NTotPts,0.0);

   for(ijk=0; ijk<NTotPts; ijk++){
     get_3ind_from_multInd(i,j,k,ijk);
     
     if ((i-4)>=0 && i<=(Nx-5)){
        rho_from_field[ijk] -= (1.0/280.0)* valEField[3*((i+4)*Ny*Nz+j*Nz+k)]/dx;
        rho_from_field[ijk] += (1.0/280.0)* valEField[3*((i-4)*Ny*Nz+j*Nz+k)]/dx;
        rho_from_field[ijk] += (4.0/105.0)* valEField[3*((i+3)*Ny*Nz+j*Nz+k)]/dx;
        rho_from_field[ijk] -= (4.0/105.0)* valEField[3*((i-3)*Ny*Nz+j*Nz+k)]/dx;
        rho_from_field[ijk] -= (1.0/5.0)* valEField[3*((i+2)*Ny*Nz+j*Nz+k)]/dx;
        rho_from_field[ijk] += (1.0/5.0)* valEField[3*((i-2)*Ny*Nz+j*Nz+k)]/dx;
        rho_from_field[ijk] += (4.0/5.0)* valEField[3*((i+1)*Ny*Nz+j*Nz+k)]/dx;
        rho_from_field[ijk] -= (4.0/5.0)* valEField[3*((i-1)*Ny*Nz+j*Nz+k)]/dx;
     }
     if ((j-4)>=0 && j<=(Ny-5)){
        rho_from_field[ijk] -= (1.0/280.0)* valEField[3*(i*Ny*Nz+(j+4)*Nz+k)+1]/dy;
        rho_from_field[ijk] += (1.0/280.0)* valEField[3*(i*Ny*Nz+(j-4)*Nz+k)+1]/dy;
        rho_from_field[ijk] += (4.0/105.0)* valEField[3*(i*Ny*Nz+(j+3)*Nz+k)+1]/dy;
        rho_from_field[ijk] -= (4.0/105.0)* valEField[3*(i*Ny*Nz+(j-3)*Nz+k)+1]/dy;
        rho_from_field[ijk] -= (1.0/5.0)* valEField[3*(i*Ny*Nz+(j+2)*Nz+k)+1]/dy;
        rho_from_field[ijk] += (1.0/5.0)* valEField[3*(i*Ny*Nz+(j-2)*Nz+k)+1]/dy;
        rho_from_field[ijk] += (4.0/5.0)* valEField[3*(i*Ny*Nz+(j+1)*Nz+k)+1]/dy;
        rho_from_field[ijk] -= (4.0/5.0)* valEField[3*(i*Ny*Nz+(j-1)*Nz+k)+1]/dy;
     
     }
     if ((k-4)>=0 && k<=(Nz-5)){
        rho_from_field[ijk] -= (1.0/280.0)* valEField[3*(i*Ny*Nz+j*Nz+(k+4))+2]/dz;
        rho_from_field[ijk] += (1.0/280.0)* valEField[3*(i*Ny*Nz+j*Nz+(k-4))+2]/dz;
        rho_from_field[ijk] += (4.0/105.0)* valEField[3*(i*Ny*Nz+j*Nz+(k+3))+2]/dz;
        rho_from_field[ijk] -= (4.0/105.0)* valEField[3*(i*Ny*Nz+j*Nz+(k-3))+2]/dz;
        rho_from_field[ijk] -= (1.0/5.0)* valEField[3*(i*Ny*Nz+j*Nz+(k+2))+2]/dz;
        rho_from_field[ijk] += (1.0/5.0)* valEField[3*(i*Ny*Nz+j*Nz+(k-2))+2]/dz;
        rho_from_field[ijk] += (4.0/5.0)* valEField[3*(i*Ny*Nz+j*Nz+(k+1))+2]/dz;
        rho_from_field[ijk] -= (4.0/5.0)* valEField[3*(i*Ny*Nz+j*Nz+(k-1))+2]/dz;
     }
   }
   
   VRscale(rho_from_field,NTotPts,-1.0/(4.0*Pi));
   double total_Rho_from_Field = Integrator(rho_from_field);
   cout << "totRho density from field = " << total_Rho_from_Field << endl;
}



void PEQSclass::SwitchZone()
{
   double *xcg,*ycg,*zcg;
   INTEGER *MidCG=QAllocINTEGER(NTotPts);
   int Nx=NPoints[0],Ny=NPoints[1],Nz=NPoints[2];
   int ijk,ibox,jbox,kbox,i,j,k,l,m,n=0,sum_occ=0;
   VRload(valSf,3*NTotPts,0.0);
   VRcopy(valSf,valEField,3*NTotPts);
   xcg=QAllocDouble(Nx);
   ycg=QAllocDouble(Ny);
   zcg=QAllocDouble(Nz);
   VRload(MidCG,NTotPts,0);
   VRload(xcg,Nx,0.0);
   VRload(ycg,Ny,0.0);
   VRload(zcg,Nz,0.0);

//==========Automating Switching zones=========================//

   double xmin=Carts[0],ymin=Carts[1],zmin=Carts[2];
   double xmax=xmin, ymax=ymin, zmax=zmin;
   double xcar,ycar,zcar;
   for(l=0;l<NAtoms;l++){
     xcar=Carts[3*l+0];
     ycar=Carts[3*l+1];
     zcar=Carts[3*l+2];
     xmin=min(xcar,xmin);
     ymin=min(ycar,ymin);
     zmin=min(zcar,zmin);
     xmax=max(xcar,xmax);
     ymax=max(ycar,ymax);
     zmax=max(zcar,zmax);
   }
   cout << "shift1 = " << shift1 << endl;
   cout << "shift3 = " << shift3 << endl;
   double InXn=xmin-shift1,InYn=ymin-shift1,InZn=zmin-shift1; // limits for Inside[]
   double InXp=xmax+shift1,InYp=ymax+shift1,InZp=zmax+shift1; // limits for Inside[]
   double OutXn=xmin-shift3,OutYn=ymin-shift3,OutZn=zmin-shift3; // limits for Outside[]
   double OutXp=xmax+shift3,OutYp=ymax+shift3,OutZp=zmax+shift3; // limits for Outside[]

   cout << "InXn = " << InXn/ang2bohr << '\t' << "InXp = " << InXp/ang2bohr << endl;
   cout << "InYn = " << InYn/ang2bohr << '\t' << "InYp = " << InYp/ang2bohr << endl;
   cout << "InZn = " << InZn/ang2bohr << '\t' << "InZp = " << InZp/ang2bohr << endl;

   cout << "OutXn = " << OutXn/ang2bohr << '\t' << "OutXp = " << OutXp/ang2bohr << endl;
   cout << "OutYn = " << OutYn/ang2bohr << '\t' << "OutYp = " << OutYp/ang2bohr << endl;
   cout << "OutZn = " << OutZn/ang2bohr << '\t' << "OutZp = " << OutZp/ang2bohr << endl;
  
//==== Collecting the ijk-points in the switching zone =====================//   

   int c1=0,switch_count;
   for(ijk=0;ijk<NTotPts;ijk++){
     if(((InXn > XYZ[3*ijk] && XYZ[3*ijk] > OutXn) || (InXp < XYZ[3*ijk] && XYZ[3*ijk] < OutXp)) &&
        ((InYn > XYZ[3*ijk+1] && XYZ[3*ijk+1] > OutYn) || (InYp < XYZ[3*ijk] && XYZ[3*ijk] < OutYp))  &&
        ((InZn > XYZ[3*ijk+2] && XYZ[3*ijk+2] > OutZn) || (InZp < XYZ[3*ijk] && XYZ[3*ijk] < OutZp))){

        MidCG[c1]=ijk;
        c1++;
     }
   }
   switch_count=c1;

   double ax,ay,az,b=abs(shift3-shift1),Lamb;   
   for(ijk=0;ijk<switch_count;ijk++){
     get_3ind_from_multInd(i,j,k,MidCG[ijk]);
     int prevX=(i-1)*Ny*Nz+j*Nz+k;
     int nextX=(i+1)*Ny*Nz+j*Nz+k;
     int prevY=i*Ny*Nz+(j-1)*Nz+k;
     int nextY=i*Ny*Nz+(j+1)*Nz+k;
     int prevZ=i*Ny*Nz+j*Nz+(k-1);
     int nextZ=i*Ny*Nz+j*Nz+(k+1);

     double x=valueOfX(i);
     double y=valueOfY(j);
     double z=valueOfZ(k);

     if(x<0.0)ax=(OutXn-InXn)/2;
     else ax=(OutXp-InXp)/2;
     if(y<0.0)ay=(OutYn-InYn)/2;
     else ay=(OutYp-InYp)/2;
     if(z<0.0)az=(OutZn-InZn)/2;
     else az=(OutZp-InZp)/2;

     Lamb=hypertan(x,ax,b)*hypertan(y,ay,b)*hypertan(z,az,b);

     valSf[3*MidCG[ijk]]=Lamb*valEField[3*prevX]+(1-Lamb)*valEField[3*nextX];
     valSf[3*MidCG[ijk]+1]=Lamb*valEField[3*prevY+1]+(1-Lamb)*valEField[3*nextY+1];
     valSf[3*MidCG[ijk]+2]=Lamb*valEField[3*prevZ+2]+(1-Lamb)*valEField[3*nextZ+2];
   }
   VRscale(valSf,3*NTotPts,4.0*Pi);
}



