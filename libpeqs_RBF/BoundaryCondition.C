#include "peqs_class.h"

void PEQSclass::BoundaryCondition()
{
  int i,j,k,l,Nx,Ny,Nz,index1=0,index2=0,index3=0;
  for (l=0;l<mg_order;l++){
    Nx=NPoints[3*l+0];Ny=NPoints[3*l+1];Nz=NPoints[3*l+2];
    for (i=0;i<Nx;i++){
      BCX_mg[index1+9*i+0] = (i-4)*Ny*Nz; 
      BCX_mg[index1+9*i+1] = (i-3)*Ny*Nz; 
      BCX_mg[index1+9*i+2] = (i-2)*Ny*Nz; 
      BCX_mg[index1+9*i+3] = (i-1)*Ny*Nz;
      BCX_mg[index1+9*i+4] = i*Ny*Nz;
      if (i<Nx-1) BCX_mg[index1+9*i+5] = (i+1)*Ny*Nz; 
      if (i<Nx-2) BCX_mg[index1+9*i+6] = (i+2)*Ny*Nz;
      if (i<Nx-3) BCX_mg[index1+9*i+7] = (i+3)*Ny*Nz;
      if (i<Nx-4) BCX_mg[index1+9*i+8] = (i+4)*Ny*Nz;
    }
    for (j=0;j<Ny;j++){
      BCY_mg[index2+9*j+0] = (j-4)*Nz; 
      BCY_mg[index2+9*j+1] = (j-3)*Nz;
      BCY_mg[index2+9*j+2] = (j-2)*Nz; 
      BCY_mg[index2+9*j+3] = (j-1)*Nz; 
      BCY_mg[index2+9*j+4] = j*Nz;
      if (j<Ny-1) BCY_mg[index2+9*j+5] = (j+1)*Nz; 
      if (j<Ny-2) BCY_mg[index2+9*j+6] = (j+2)*Nz;
      if (j<Ny-3) BCY_mg[index2+9*j+7] = (j+3)*Nz;
      if (j<Ny-4) BCY_mg[index2+9*j+8] = (j+4)*Nz;
    }
    for (k=0;k<Nz;k++){
      BCZ_mg[index3+9*k+0] = k-4;
      BCZ_mg[index3+9*k+1] = k-3; 
      BCZ_mg[index3+9*k+2] = k-2; 
      BCZ_mg[index3+9*k+3] = k-1;
      BCZ_mg[index3+9*k+4] = k;
      if (k<Nz-1) BCZ_mg[index3+9*k+5] = k+1; 
      if (k<Nz-2) BCZ_mg[index3+9*k+6] = k+2; 
      if (k<Nz-3) BCZ_mg[index3+9*k+7] = k+3; 
      if (k<Nz-4) BCZ_mg[index3+9*k+8] = k+4;
    }
    index1 += 9*Nx; index2 += 9*Ny; index3 += 9*Nz;
  }
  VIcopy(BCX,&BCX_mg[0],9*NPoints[0]);
  VIcopy(BCY,&BCY_mg[0],9*NPoints[1]);
  VIcopy(BCZ,&BCZ_mg[0],9*NPoints[2]);
}
