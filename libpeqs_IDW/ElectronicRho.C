#include "peqs_class.h"

//
// Electronic part of the electrostatic potential is computed here
// For the nuclear part, see NuclearRhoPhi()
//
// MPC 2017
//

void PEQSclass::ElectronicRho(double *PAv)
{
  int i,j,k,ijk,NPPer=MaxPPer,Nx,Ny,Nz;
  double dx,dx2,dy,dy2,dz,dz2,*BatchPts=QAllocDouble(3*MaxPPer),*ESP=QAllocDouble(MaxPPer);
  Nx = NPoints[0]; Ny = NPoints[1]; Nz = NPoints[2];
  dx = delX[0]; dx2 = dx*dx;
  dy = delX[1]; dy2 = dy*dy;
  dz = delX[2]; dz2 = dz*dz;
  VRload(rho_elec,NTotPts,0.0);
  VRload(phi_elec,NTotPts,0.0);

  qtime_t shackle02;
  double shackle03[3];
  shackle02 = QTimerOn();
  double *data=QAllocDouble(NTotPts);

  VRscale(phi_elec,NTotPts,1.0/(4.0*Pi));
  #pragma omp parallel for private(i,j,k,ijk) firstprivate(dx2,dy2,dz2)
  for (ijk=0;ijk<NTotPts;ijk++){
    get_3ind_from_multInd(i,j,k,ijk,false);
    int x0=BCX[9*i+4]; int y0=BCY[9*j+4]; int z0=BCZ[9*k+4];
    if (fdif_order == 2)
    {
       rho_elec[ijk] = (2.0)*(1.0/dx2 + 1.0/dy2 + 1.0/dz2)*phi_elec[ijk];
       if (BCX[9*i+3] >= 0) rho_elec[ijk] -= phi_elec[BCX[9*i+3]+y0+z0]/dx2;
       if (BCX[9*i+5] >= 0) rho_elec[ijk] -= phi_elec[BCX[9*i+5]+y0+z0]/dx2;
       if (BCY[9*j+3] >= 0) rho_elec[ijk] -= phi_elec[x0+BCY[9*j+3]+z0]/dy2;
       if (BCY[9*j+5] >= 0) rho_elec[ijk] -= phi_elec[x0+BCY[9*j+5]+z0]/dy2;
       if (BCZ[9*k+3] >= 0) rho_elec[ijk] -= phi_elec[x0+y0+BCZ[9*k+3]]/dz2;
       if (BCZ[9*k+5] >= 0) rho_elec[ijk] -= phi_elec[x0+y0+BCZ[9*k+5]]/dz2;
       if (i == 0 || i == (Nx-1) || j == 0 || j == (Ny-1) || k == 0 || k == (Nz-1))
           rho_elec[ijk] = 0.0;
    }
    else if (fdif_order == 4)
    {
       rho_elec[ijk] = (5.0/2.0)*(1.0/dx2 + 1.0/dy2 + 1.0/dz2)*phi_elec[ijk];
       if (BCX[9*i+2] >= 0) rho_elec[ijk] += (1.0/12.0)*phi_elec[BCX[9*i+2]+y0+z0]/dx2;
       if (BCX[9*i+3] >= 0) rho_elec[ijk] -= (4.0/3.0)*phi_elec[BCX[9*i+3]+y0+z0]/dx2;
       if (BCX[9*i+5] >= 0) rho_elec[ijk] -= (4.0/3.0)*phi_elec[BCX[9*i+5]+y0+z0]/dx2;
       if (BCX[9*i+6] >= 0) rho_elec[ijk] += (1.0/12.0)*phi_elec[BCX[9*i+6]+y0+z0]/dx2;
       if (BCY[9*j+2] >= 0) rho_elec[ijk] += (1.0/12.0)*phi_elec[x0+BCY[9*j+2]+z0]/dy2;
       if (BCY[9*j+3] >= 0) rho_elec[ijk] -= (4.0/3.0)*phi_elec[x0+BCY[9*j+3]+z0]/dy2;
       if (BCY[9*j+5] >= 0) rho_elec[ijk] -= (4.0/3.0)*phi_elec[x0+BCY[9*j+5]+z0]/dy2;
       if (BCY[9*j+6] >= 0) rho_elec[ijk] += (1.0/12.0)*phi_elec[x0+BCY[9*j+6]+z0]/dy2;
       if (BCZ[9*k+2] >= 0) rho_elec[ijk] += (1.0/12.0)*phi_elec[x0+y0+BCZ[9*k+2]]/dz2;
       if (BCZ[9*k+3] >= 0) rho_elec[ijk] -= (4.0/3.0)*phi_elec[x0+y0+BCZ[9*k+3]]/dz2;
       if (BCZ[9*k+5] >= 0) rho_elec[ijk] -= (4.0/3.0)*phi_elec[x0+y0+BCZ[9*k+5]]/dz2;
       if (BCZ[9*k+6] >= 0) rho_elec[ijk] += (1.0/12.0)*phi_elec[x0+y0+BCZ[9*k+6]]/dz2;
       if (i == 0 || i == 1 || i == (Nx-1) || i == (Nx-2) ||
           j == 0 || j == 1 || j == (Ny-1) || j == (Ny-2) ||
           k == 0 || k == 1 || k == (Nz-1) || k == (Nz-2) 
          )rho_elec[ijk] = 0.0;
    }
    else if (fdif_order == 6)
    {
       rho_elec[ijk] = (49.0/18.0)*(1.0/dx2 + 1.0/dy2 + 1.0/dz2)*phi_elec[ijk];
       if (BCX[9*i+1] >= 0) rho_elec[ijk] -= (1.0/90.0)*phi_elec[BCX[9*i+1]+y0+z0]/dx2;
       if (BCX[9*i+2] >= 0) rho_elec[ijk] += (3.0/20.0)*phi_elec[BCX[9*i+2]+y0+z0]/dx2;
       if (BCX[9*i+3] >= 0) rho_elec[ijk] -= (3.0/2.0)*phi_elec[BCX[9*i+3]+y0+z0]/dx2;
       if (BCX[9*i+5] >= 0) rho_elec[ijk] -= (3.0/2.0)*phi_elec[BCX[9*i+5]+y0+z0]/dx2;
       if (BCX[9*i+6] >= 0) rho_elec[ijk] += (3.0/20.0)*phi_elec[BCX[9*i+6]+y0+z0]/dx2;
       if (BCX[9*i+7] >= 0) rho_elec[ijk] -= (1.0/90.0)*phi_elec[BCX[9*i+7]+y0+z0]/dx2;
       if (BCY[9*j+1] >= 0) rho_elec[ijk] -= (1.0/90.0)*phi_elec[x0+BCY[9*j+1]+z0]/dy2;
       if (BCY[9*j+2] >= 0) rho_elec[ijk] += (3.0/20.0)*phi_elec[x0+BCY[9*j+2]+z0]/dy2;
       if (BCY[9*j+3] >= 0) rho_elec[ijk] -= (3.0/2.0)*phi_elec[x0+BCY[9*j+3]+z0]/dy2;
       if (BCY[9*j+5] >= 0) rho_elec[ijk] -= (3.0/2.0)*phi_elec[x0+BCY[9*j+5]+z0]/dy2;
       if (BCY[9*j+6] >= 0) rho_elec[ijk] += (3.0/20.0)*phi_elec[x0+BCY[9*j+6]+z0]/dy2;
       if (BCY[9*j+7] >= 0) rho_elec[ijk] -= (1.0/90.0)*phi_elec[x0+BCY[9*j+7]+z0]/dy2;
       if (BCZ[9*k+1] >= 0) rho_elec[ijk] -= (1.0/90.0)*phi_elec[x0+y0+BCZ[9*k+1]]/dz2;
       if (BCZ[9*k+2] >= 0) rho_elec[ijk] += (3.0/20.0)*phi_elec[x0+y0+BCZ[9*k+2]]/dz2;
       if (BCZ[9*k+3] >= 0) rho_elec[ijk] -= (3.0/2.0)*phi_elec[x0+y0+BCZ[9*k+3]]/dz2;
       if (BCZ[9*k+5] >= 0) rho_elec[ijk] -= (3.0/2.0)*phi_elec[x0+y0+BCZ[9*k+5]]/dz2;
       if (BCZ[9*k+6] >= 0) rho_elec[ijk] += (3.0/20.0)*phi_elec[x0+y0+BCZ[9*k+6]]/dz2;
       if (BCZ[9*k+7] >= 0) rho_elec[ijk] -= (1.0/90.0)*phi_elec[x0+y0+BCZ[9*k+7]]/dz2;
       if (i == 0 || i == 1 || i == 2 || i == (Nx-1) || i == (Nx-2) || i == (Nx-3) ||
           j == 0 || j == 1 || j == 2 || j == (Ny-1) || j == (Ny-2) || j == (Ny-3) ||
           k == 0 || k == 1 || k == 2 || k == (Nz-1) || k == (Nz-2) || k == (Nz-3) 
          )rho_elec[ijk] = 0.0;
    }
    else if (fdif_order == 8)
    {
       rho_elec[ijk] = (205.0/72.0)*(1.0/dx2 + 1.0/dy2 + 1.0/dz2)*phi_elec[ijk];
       if (BCX[9*i+0] >= 0) rho_elec[ijk] += (1.0/560.0)*phi_elec[BCX[9*i+0]+y0+z0]/dx2;
       if (BCX[9*i+1] >= 0) rho_elec[ijk] -= (8.0/315.0)*phi_elec[BCX[9*i+1]+y0+z0]/dx2;
       if (BCX[9*i+2] >= 0) rho_elec[ijk] += (1.0/5.0)*phi_elec[BCX[9*i+2]+y0+z0]/dx2;
       if (BCX[9*i+3] >= 0) rho_elec[ijk] -= (8.0/5.0)*phi_elec[BCX[9*i+3]+y0+z0]/dx2;
       if (BCX[9*i+5] >= 0) rho_elec[ijk] -= (8.0/5.0)*phi_elec[BCX[9*i+5]+y0+z0]/dx2;
       if (BCX[9*i+6] >= 0) rho_elec[ijk] += (1.0/5.0)*phi_elec[BCX[9*i+6]+y0+z0]/dx2;
       if (BCX[9*i+7] >= 0) rho_elec[ijk] -= (8.0/315.0)*phi_elec[BCX[9*i+7]+y0+z0]/dx2;
       if (BCX[9*i+8] >= 0) rho_elec[ijk] += (1.0/560.0)*phi_elec[BCX[9*i+8]+y0+z0]/dx2;
       if (BCY[9*j+0] >= 0) rho_elec[ijk] += (1.0/560.0)*phi_elec[x0+BCY[9*j+0]+z0]/dy2;
       if (BCY[9*j+1] >= 0) rho_elec[ijk] -= (8.0/315.0)*phi_elec[x0+BCY[9*j+1]+z0]/dy2;
       if (BCY[9*j+2] >= 0) rho_elec[ijk] += (1.0/5.0)*phi_elec[x0+BCY[9*j+2]+z0]/dy2;
       if (BCY[9*j+3] >= 0) rho_elec[ijk] -= (8.0/5.0)*phi_elec[x0+BCY[9*j+3]+z0]/dy2;
       if (BCY[9*j+5] >= 0) rho_elec[ijk] -= (8.0/5.0)*phi_elec[x0+BCY[9*j+5]+z0]/dy2;
       if (BCY[9*j+6] >= 0) rho_elec[ijk] += (1.0/5.0)*phi_elec[x0+BCY[9*j+6]+z0]/dy2;
       if (BCY[9*j+7] >= 0) rho_elec[ijk] -= (8.0/315.0)*phi_elec[x0+BCY[9*j+7]+z0]/dy2;
       if (BCY[9*j+8] >= 0) rho_elec[ijk] += (1.0/560.0)*phi_elec[x0+BCY[9*j+8]+z0]/dy2;
       if (BCZ[9*k+0] >= 0) rho_elec[ijk] += (1.0/560.0)*phi_elec[x0+y0+BCZ[9*k+0]]/dz2;
       if (BCZ[9*k+1] >= 0) rho_elec[ijk] -= (8.0/315.0)*phi_elec[x0+y0+BCZ[9*k+1]]/dz2;
       if (BCZ[9*k+2] >= 0) rho_elec[ijk] += (1.0/5.0)*phi_elec[x0+y0+BCZ[9*k+2]]/dz2;
       if (BCZ[9*k+3] >= 0) rho_elec[ijk] -= (8.0/5.0)*phi_elec[x0+y0+BCZ[9*k+3]]/dz2;
       if (BCZ[9*k+5] >= 0) rho_elec[ijk] -= (8.0/5.0)*phi_elec[x0+y0+BCZ[9*k+5]]/dz2;
       if (BCZ[9*k+6] >= 0) rho_elec[ijk] += (1.0/5.0)*phi_elec[x0+y0+BCZ[9*k+6]]/dz2;
       if (BCZ[9*k+7] >= 0) rho_elec[ijk] -= (8.0/315.0)*phi_elec[x0+y0+BCZ[9*k+7]]/dz2;
       if (BCZ[9*k+8] >= 0) rho_elec[ijk] += (1.0/560.0)*phi_elec[x0+y0+BCZ[9*k+8]]/dz2;
       if (i == 0 || i == 1 || i == 2 || i == 3 ||
           i == (Nx-1) || i == (Nx-2) || i == (Nx-3) || i == (Nx-4) ||
           j == 0 || j == 1 || j == 2 || j == 3 ||
           j == (Ny-1) || j == (Ny-2) || j == (Ny-3) || j == (Ny-4) ||
           k == 0 || k == 1 || k == 2 || k == 3 ||
           k == (Nz-1) || k == (Nz-2) || k == (Nz-3) || k == (Nz-4)
          )rho_elec[ijk] = 0.0;
    }
  }

  QTimerOff(shackle03,shackle02);
  if (peqs_print > 0)
     printf("\tPEqS:  Computed the e- density in %2.2f (cpu) %2.2f (wall) sec. \n",
           shackle03[0],shackle03[2]);

}
