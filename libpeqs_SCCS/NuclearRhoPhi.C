#include "peqs_class.h"

// 
// Nuclear part of the electrostatic potential is computed here.
// For the electronic part, see ElectronicRho()
//
// MPC 2017
//

void PEQSclass::NuclearRhoPhi()
{
  int *EffAtNo = QAllocINTEGER(NAtoms);
  VIcopy(EffAtNo,AtNo,NAtoms);

  if (rem_read(REM_PSEUDOPOTENTIAL) > 0)
  {  // JMH(4/2020): for ECPs, need to modify the atomic numbers  

     int* atomic_numbers0; // actual atno
     get_carts(NULL, NULL, &atomic_numbers0, NULL);
     int *atomic_numbers = QAllocINTEGER(NAtoms); // effective atno
     VIcopy(atomic_numbers, atomic_numbers0, NAtoms);

     INTEGER pseudodata[104];
     GetPseudoData(pseudodata);
     for(int i=0; i < NAtoms; i++) { 
        atomic_numbers[i] -= pseudodata[atomic_numbers[i]];
     } 
     VIcopy(EffAtNo,atomic_numbers,NAtoms);
     QFree(atomic_numbers);
  }

  double charge=1.0,dist=0.0;
  int i,j,k,ijk,l;
  int Nx=NPoints[0],Ny=NPoints[1],Nz=NPoints[2],x0,y0,z0;
  double dx,dx2,dy,dy2,dz,dz2;
  dx = delX[0]; dy = delX[1]; dz = delX[2];
  dx2 = dx*dx; dy2 = dy*dy; dz2 = dz*dz;
  for (l=0;l<NAtoms;l++){
    charge = -1.0*(double)EffAtNo[l];
    if (GaussBlur){ //gaussian blur nuclear charge 
      double sigma2=sigma*sigma;
      double cnst=1.0/(pow(2.0*Pi,1.5)*sigma*sigma*sigma);
      for (ijk=0;ijk<NTotPts;ijk++){
        get_3ind_from_multInd(i,j,k,ijk);
        double x = valueOfX(i) - Carts[3*l+0],x2=x*x;
        double y = valueOfY(j) - Carts[3*l+1],y2=y*y;
        double z = valueOfZ(k) - Carts[3*l+2],z2=z*z;
        double expX = exp(-0.5*x2/sigma2);
        double expY = exp(-0.5*y2/sigma2);
        double expZ = exp(-0.5*z2/sigma2);
        dist = sqrt(x2+y2+z2);
        if (dist != 0.0) phi_nuke[ijk] += (charge/dist)*erf(dist/(sqrt(2.0)*sigma));
        else phi_nuke[ijk] += sqrt(2)*charge/(sigma*sqrt(Pi));
      }
    }
    else{ //point charge
      double *rX=QAllocDouble(3),tolerance=0.001;
      for (ijk=0;ijk<NTotPts;ijk++){
        VRsub(rX,&XYZ[3*ijk],&Carts[3*l],3);
        dist = sqrt(VRdot(rX,rX,3));
        if (dist > tolerance) phi_nuke[ijk] += charge/dist;
      }
      QFree(rX);
    }
  }
  VRscale(phi_nuke,NTotPts,1.0/(4.0*Pi));
  QFree(EffAtNo);

  #pragma omp parallel for private(i,j,k,ijk) firstprivate(dx2,dy2,dz2)
  for (ijk=0;ijk<NTotPts;ijk++)
  {
     get_3ind_from_multInd(i,j,k,ijk,false);
     if (fdif_order == 2) 
        rho_nuke[ijk] = (2.0)*(1.0/dx2 + 1.0/dy2 + 1.0/dz2)*phi_nuke[ijk];
     else if (fdif_order == 4) 
        rho_nuke[ijk] = (5.0/2.0)*(1.0/dx2 + 1.0/dy2 + 1.0/dz2)*phi_nuke[ijk];
    else if (fdif_order == 6) 
        rho_nuke[ijk] = (49.0/18.0)*(1.0/dx2 + 1.0/dy2 + 1.0/dz2)*phi_nuke[ijk];
    else if (fdif_order == 8) 
        rho_nuke[ijk] = (205.0/72.0)*(1.0/dx2 + 1.0/dy2 + 1.0/dz2)*phi_nuke[ijk];

    int x0=BCX[9*i+4]; int y0=BCY[9*j+4]; int z0=BCZ[9*k+4];
    if (fdif_order == 2)
    {
       if (BCX[9*i+3] >= 0) rho_nuke[ijk] -= phi_nuke[BCX[9*i+3]+y0+z0]/dx2;
       if (BCX[9*i+5] >= 0) rho_nuke[ijk] -= phi_nuke[BCX[9*i+5]+y0+z0]/dx2;
       if (BCY[9*j+3] >= 0) rho_nuke[ijk] -= phi_nuke[x0+BCY[9*j+3]+z0]/dy2;
       if (BCY[9*j+5] >= 0) rho_nuke[ijk] -= phi_nuke[x0+BCY[9*j+5]+z0]/dy2;
       if (BCZ[9*k+3] >= 0) rho_nuke[ijk] -= phi_nuke[x0+y0+BCZ[9*k+3]]/dz2;
       if (BCZ[9*k+5] >= 0) rho_nuke[ijk] -= phi_nuke[x0+y0+BCZ[9*k+5]]/dz2;
       if (i == 0 || i == (Nx-1) || j == 0 || j == (Ny-1) || k == 0 || k == (Nz-1))
          rho_nuke[ijk] = 0.0;
    }
    else if (fdif_order == 4)
    {
       if (BCX[9*i+2] >= 0) rho_nuke[ijk] += (1.0/12.0)*phi_nuke[BCX[9*i+2]+y0+z0]/dx2;
       if (BCX[9*i+3] >= 0) rho_nuke[ijk] -= (4.0/3.0)*phi_nuke[BCX[9*i+3]+y0+z0]/dx2;
       if (BCX[9*i+5] >= 0) rho_nuke[ijk] -= (4.0/3.0)*phi_nuke[BCX[9*i+5]+y0+z0]/dx2;
       if (BCX[9*i+6] >= 0) rho_nuke[ijk] += (1.0/12.0)*phi_nuke[BCX[9*i+6]+y0+z0]/dx2;
       if (BCY[9*j+2] >= 0) rho_nuke[ijk] += (1.0/12.0)*phi_nuke[x0+BCY[9*j+2]+z0]/dy2;
       if (BCY[9*j+3] >= 0) rho_nuke[ijk] -= (4.0/3.0)*phi_nuke[x0+BCY[9*j+3]+z0]/dy2;
       if (BCY[9*j+5] >= 0) rho_nuke[ijk] -= (4.0/3.0)*phi_nuke[x0+BCY[9*j+5]+z0]/dy2;
       if (BCY[9*j+6] >= 0) rho_nuke[ijk] += (1.0/12.0)*phi_nuke[x0+BCY[9*j+6]+z0]/dy2;
       if (BCZ[9*k+2] >= 0) rho_nuke[ijk] += (1.0/12.0)*phi_nuke[x0+y0+BCZ[9*k+2]]/dz2;
       if (BCZ[9*k+3] >= 0) rho_nuke[ijk] -= (4.0/3.0)*phi_nuke[x0+y0+BCZ[9*k+3]]/dz2;
       if (BCZ[9*k+5] >= 0) rho_nuke[ijk] -= (4.0/3.0)*phi_nuke[x0+y0+BCZ[9*k+5]]/dz2;
       if (BCZ[9*k+6] >= 0) rho_nuke[ijk] += (1.0/12.0)*phi_nuke[x0+y0+BCZ[9*k+6]]/dz2;
       if (i == 0 || i == 1 || i == (Nx-1) || i == (Nx-2) ||
           j == 0 || j == 1 || j == (Ny-1) || j == (Ny-2) ||
           k == 0 || k == 1 || k == (Nz-1) || k == (Nz-2) 
          )rho_nuke[ijk] = 0.0;
    }
    else if (fdif_order == 6)
    {
       if (BCX[9*i+1] >= 0) rho_nuke[ijk] -= (1.0/90.0)*phi_nuke[BCX[9*i+1]+y0+z0]/dx2;
       if (BCX[9*i+2] >= 0) rho_nuke[ijk] += (3.0/20.0)*phi_nuke[BCX[9*i+2]+y0+z0]/dx2;
       if (BCX[9*i+3] >= 0) rho_nuke[ijk] -= (3.0/2.0)*phi_nuke[BCX[9*i+3]+y0+z0]/dx2;
       if (BCX[9*i+5] >= 0) rho_nuke[ijk] -= (3.0/2.0)*phi_nuke[BCX[9*i+5]+y0+z0]/dx2;
       if (BCX[9*i+6] >= 0) rho_nuke[ijk] += (3.0/20.0)*phi_nuke[BCX[9*i+6]+y0+z0]/dx2;
       if (BCX[9*i+7] >= 0) rho_nuke[ijk] -= (1.0/90.0)*phi_nuke[BCX[9*i+7]+y0+z0]/dx2;
       if (BCY[9*j+1] >= 0) rho_nuke[ijk] -= (1.0/90.0)*phi_nuke[x0+BCY[9*j+1]+z0]/dy2;
       if (BCY[9*j+2] >= 0) rho_nuke[ijk] += (3.0/20.0)*phi_nuke[x0+BCY[9*j+2]+z0]/dy2;
       if (BCY[9*j+3] >= 0) rho_nuke[ijk] -= (3.0/2.0)*phi_nuke[x0+BCY[9*j+3]+z0]/dy2;
       if (BCY[9*j+5] >= 0) rho_nuke[ijk] -= (3.0/2.0)*phi_nuke[x0+BCY[9*j+5]+z0]/dy2;
       if (BCY[9*j+6] >= 0) rho_nuke[ijk] += (3.0/20.0)*phi_nuke[x0+BCY[9*j+6]+z0]/dy2;
       if (BCY[9*j+7] >= 0) rho_nuke[ijk] -= (1.0/90.0)*phi_nuke[x0+BCY[9*j+7]+z0]/dy2;
       if (BCZ[9*k+1] >= 0) rho_nuke[ijk] -= (1.0/90.0)*phi_nuke[x0+y0+BCZ[9*k+1]]/dz2;
       if (BCZ[9*k+2] >= 0) rho_nuke[ijk] += (3.0/20.0)*phi_nuke[x0+y0+BCZ[9*k+2]]/dz2;
       if (BCZ[9*k+3] >= 0) rho_nuke[ijk] -= (3.0/2.0)*phi_nuke[x0+y0+BCZ[9*k+3]]/dz2;
       if (BCZ[9*k+5] >= 0) rho_nuke[ijk] -= (3.0/2.0)*phi_nuke[x0+y0+BCZ[9*k+5]]/dz2;
       if (BCZ[9*k+6] >= 0) rho_nuke[ijk] += (3.0/20.0)*phi_nuke[x0+y0+BCZ[9*k+6]]/dz2;
       if (BCZ[9*k+7] >= 0) rho_nuke[ijk] -= (1.0/90.0)*phi_nuke[x0+y0+BCZ[9*k+7]]/dz2;
       if (i == 0 || i == 1 || i == 2 || i == (Nx-1) || i == (Nx-2) || i == (Nx-3) || 
           j == 0 || j == 1 || j == 2 || j == (Ny-1) || j == (Ny-2) || j == (Ny-3) ||
           k == 0 || k == 1 || k == 2 || k == (Nz-1) || k == (Nz-2) || k == (Nz-3) 
          )rho_nuke[ijk] = 0.0;
    }
    else if (fdif_order == 8)
    {
       if (BCX[9*i+0] >= 0) rho_nuke[ijk] += (1.0/560.0)*phi_nuke[BCX[9*i+0]+y0+z0]/dx2;
       if (BCX[9*i+1] >= 0) rho_nuke[ijk] -= (8.0/315.0)*phi_nuke[BCX[9*i+1]+y0+z0]/dx2;
       if (BCX[9*i+2] >= 0) rho_nuke[ijk] += (1.0/5.0)*phi_nuke[BCX[9*i+2]+y0+z0]/dx2;
       if (BCX[9*i+3] >= 0) rho_nuke[ijk] -= (8.0/5.0)*phi_nuke[BCX[9*i+3]+y0+z0]/dx2;
       if (BCX[9*i+5] >= 0) rho_nuke[ijk] -= (8.0/5.0)*phi_nuke[BCX[9*i+5]+y0+z0]/dx2;
       if (BCX[9*i+6] >= 0) rho_nuke[ijk] += (1.0/5.0)*phi_nuke[BCX[9*i+6]+y0+z0]/dx2;
       if (BCX[9*i+7] >= 0) rho_nuke[ijk] -= (8.0/315.0)*phi_nuke[BCX[9*i+7]+y0+z0]/dx2;
       if (BCX[9*i+8] >= 0) rho_nuke[ijk] += (1.0/560.0)*phi_nuke[BCX[9*i+8]+y0+z0]/dx2;
       if (BCY[9*j+0] >= 0) rho_nuke[ijk] += (1.0/560.0)*phi_nuke[x0+BCY[9*j+0]+z0]/dy2;
       if (BCY[9*j+1] >= 0) rho_nuke[ijk] -= (8.0/315.0)*phi_nuke[x0+BCY[9*j+1]+z0]/dy2;
       if (BCY[9*j+2] >= 0) rho_nuke[ijk] += (1.0/5.0)*phi_nuke[x0+BCY[9*j+2]+z0]/dy2;
       if (BCY[9*j+3] >= 0) rho_nuke[ijk] -= (8.0/5.0)*phi_nuke[x0+BCY[9*j+3]+z0]/dy2;
       if (BCY[9*j+5] >= 0) rho_nuke[ijk] -= (8.0/5.0)*phi_nuke[x0+BCY[9*j+5]+z0]/dy2;
       if (BCY[9*j+6] >= 0) rho_nuke[ijk] += (1.0/5.0)*phi_nuke[x0+BCY[9*j+6]+z0]/dy2;
       if (BCY[9*j+7] >= 0) rho_nuke[ijk] -= (8.0/315.0)*phi_nuke[x0+BCY[9*j+7]+z0]/dy2;
       if (BCY[9*j+8] >= 0) rho_nuke[ijk] += (1.0/560.0)*phi_nuke[x0+BCY[9*j+8]+z0]/dy2;
       if (BCZ[9*k+0] >= 0) rho_nuke[ijk] += (1.0/560.0)*phi_nuke[x0+y0+BCZ[9*k+0]]/dz2;
       if (BCZ[9*k+1] >= 0) rho_nuke[ijk] -= (8.0/315.0)*phi_nuke[x0+y0+BCZ[9*k+1]]/dz2;
       if (BCZ[9*k+2] >= 0) rho_nuke[ijk] += (1.0/5.0)*phi_nuke[x0+y0+BCZ[9*k+2]]/dz2;
       if (BCZ[9*k+3] >= 0) rho_nuke[ijk] -= (8.0/5.0)*phi_nuke[x0+y0+BCZ[9*k+3]]/dz2;
       if (BCZ[9*k+5] >= 0) rho_nuke[ijk] -= (8.0/5.0)*phi_nuke[x0+y0+BCZ[9*k+5]]/dz2;
       if (BCZ[9*k+6] >= 0) rho_nuke[ijk] += (1.0/5.0)*phi_nuke[x0+y0+BCZ[9*k+6]]/dz2;
       if (BCZ[9*k+7] >= 0) rho_nuke[ijk] -= (8.0/315.0)*phi_nuke[x0+y0+BCZ[9*k+7]]/dz2;
       if (BCZ[9*k+8] >= 0) rho_nuke[ijk] += (1.0/560.0)*phi_nuke[x0+y0+BCZ[9*k+8]]/dz2;
       if (i == 0 || i == 1 || i == 2 || i == 3 ||
           i == (Nx-1) || i == (Nx-2) || i == (Nx-3) || i == (Nx-4) ||
           j == 0 || j == 1 || j == 2 || j == 3 ||
           j == (Ny-1) || j == (Ny-2) || j == (Ny-3) || j == (Ny-4) ||
           k == 0 || k == 1 || k == 2 || k == 3 ||
           k == (Nz-1) || k == (Nz-2) || k == (Nz-3) || k == (Nz-4)
          )rho_nuke[ijk] = 0.0;
    }      
  }
}
