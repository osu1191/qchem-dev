#include "peqs_class.h"

void PEQSclass::ComputeCharge(double *Q, double *vector, bool restriction)
{
  //note that this is -1.0 * laplacian
  double dx,dx2,dy,dy2,dz,dz2;
  int i,j,k,ijk,Nx,Ny,Nz,index=mg_order-mg_level,index1=0,index2=0,index3=0;
  if (restriction) index -= 1;
  Nx = NPoints[3*index+0]; Ny = NPoints[3*index+1]; Nz = NPoints[3*index+2];
  dx = delX[3*index+0]; dy = delX[3*index+1]; dz = delX[3*index+2];
  dx2 = dx*dx; dy2 = dy*dy; dz2 = dz*dz;
  for (ijk=0;ijk<index;ijk++){
    index1 += 9*NPoints[3*ijk+0];
    index2 += 9*NPoints[3*ijk+1];
    index3 += 9*NPoints[3*ijk+2];
  }
  #pragma omp parallel \
  private(i,j,k,ijk) \
  firstprivate(index,index1,index2,index3,restriction,dx2,dy2,dz2) 
  {
  #pragma omp for
  for (ijk=0;ijk<NTotPts_mg[index];ijk++){
    get_3ind_from_multInd(i,j,k,ijk,restriction);  
    int x0=BCX_mg[index1+9*i+4]; int y0=BCY_mg[index2+9*j+4]; int z0=BCZ_mg[index3+9*k+4];
    if (fdif_order == 2){
    Q[ijk] =  (2.0)*(1.0/dx2 + 1.0/dy2 + 1.0/dz2)*vector[ijk];
    if (BCX_mg[index1+9*i+3] >= 0) Q[ijk] -= vector[BCX_mg[index1+9*i+3]+y0+z0]/dx2;
    if (BCX_mg[index1+9*i+5] >= 0) Q[ijk] -= vector[BCX_mg[index1+9*i+5]+y0+z0]/dx2;
    if (BCY_mg[index2+9*j+3] >= 0) Q[ijk] -= vector[x0+BCY_mg[index2+9*j+3]+z0]/dy2;
    if (BCY_mg[index2+9*j+5] >= 0) Q[ijk] -= vector[x0+BCY_mg[index2+9*j+5]+z0]/dy2;
    if (BCZ_mg[index3+9*k+3] >= 0) Q[ijk] -= vector[x0+y0+BCZ_mg[index3+9*k+3]]/dz2;
    if (BCZ_mg[index3+9*k+5] >= 0) Q[ijk] -= vector[x0+y0+BCZ_mg[index3+9*k+5]]/dz2;
    }
    else if (fdif_order == 4){
    Q[ijk] =  (5.0/2.0)*(1.0/dx2 + 1.0/dy2 + 1.0/dz2)*vector[ijk];
    if (BCX_mg[index1+9*i+2] >= 0) Q[ijk] += (1.0/12.0)*vector[BCX_mg[index1+9*i+2]+y0+z0]/dx2;
    if (BCX_mg[index1+9*i+3] >= 0) Q[ijk] -= (4.0/3.0)*vector[BCX_mg[index1+9*i+3]+y0+z0]/dx2;
    if (BCX_mg[index1+9*i+5] >= 0) Q[ijk] -= (4.0/3.0)*vector[BCX_mg[index1+9*i+5]+y0+z0]/dx2;
    if (BCX_mg[index1+9*i+6] >= 0) Q[ijk] += (1.0/12.0)*vector[BCX_mg[index1+9*i+6]+y0+z0]/dx2;
    if (BCY_mg[index2+9*j+2] >= 0) Q[ijk] += (1.0/12.0)*vector[x0+BCY_mg[index2+9*j+2]+z0]/dy2;
    if (BCY_mg[index2+9*j+3] >= 0) Q[ijk] -= (4.0/3.0)*vector[x0+BCY_mg[index2+9*j+3]+z0]/dy2;
    if (BCY_mg[index2+9*j+5] >= 0) Q[ijk] -= (4.0/3.0)*vector[x0+BCY_mg[index2+9*j+5]+z0]/dy2;
    if (BCY_mg[index2+9*j+6] >= 0) Q[ijk] += (1.0/12.0)*vector[x0+BCY_mg[index2+9*j+6]+z0]/dy2;
    if (BCZ_mg[index3+9*k+2] >= 0) Q[ijk] += (1.0/12.0)*vector[x0+y0+BCZ_mg[index3+9*k+2]]/dz2;
    if (BCZ_mg[index3+9*k+3] >= 0) Q[ijk] -= (4.0/3.0)*vector[x0+y0+BCZ_mg[index3+9*k+3]]/dz2;
    if (BCZ_mg[index3+9*k+5] >= 0) Q[ijk] -= (4.0/3.0)*vector[x0+y0+BCZ_mg[index3+9*k+5]]/dz2;
    if (BCZ_mg[index3+9*k+6] >= 0) Q[ijk] += (1.0/12.0)*vector[x0+y0+BCZ_mg[index3+9*k+6]]/dz2;
    }
    else if (fdif_order == 6){
    Q[ijk] =  (49.0/18.0)*(1.0/dx2 + 1.0/dy2 + 1.0/dz2)*vector[ijk];
    if (BCX_mg[index1+9*i+1] >= 0) Q[ijk] -= (1.0/90.0)*vector[BCX_mg[index1+9*i+1]+y0+z0]/dx2;
    if (BCX_mg[index1+9*i+2] >= 0) Q[ijk] += (3.0/20.0)*vector[BCX_mg[index1+9*i+2]+y0+z0]/dx2;
    if (BCX_mg[index1+9*i+3] >= 0) Q[ijk] -= (3.0/2.0)*vector[BCX_mg[index1+9*i+3]+y0+z0]/dx2;
    if (BCX_mg[index1+9*i+5] >= 0) Q[ijk] -= (3.0/2.0)*vector[BCX_mg[index1+9*i+5]+y0+z0]/dx2;
    if (BCX_mg[index1+9*i+6] >= 0) Q[ijk] += (3.0/20.0)*vector[BCX_mg[index1+9*i+6]+y0+z0]/dx2;
    if (BCX_mg[index1+9*i+7] >= 0) Q[ijk] -= (1.0/90.0)*vector[BCX_mg[index1+9*i+7]+y0+z0]/dx2;
    if (BCY_mg[index2+9*j+1] >= 0) Q[ijk] -= (1.0/90.0)*vector[x0+BCY_mg[index2+9*j+1]+z0]/dy2;
    if (BCY_mg[index2+9*j+2] >= 0) Q[ijk] += (3.0/20.0)*vector[x0+BCY_mg[index2+9*j+2]+z0]/dy2;
    if (BCY_mg[index2+9*j+3] >= 0) Q[ijk] -= (3.0/2.0)*vector[x0+BCY_mg[index2+9*j+3]+z0]/dy2;
    if (BCY_mg[index2+9*j+5] >= 0) Q[ijk] -= (3.0/2.0)*vector[x0+BCY_mg[index2+9*j+5]+z0]/dy2;
    if (BCY_mg[index2+9*j+6] >= 0) Q[ijk] += (3.0/20.0)*vector[x0+BCY_mg[index2+9*j+6]+z0]/dy2;
    if (BCY_mg[index2+9*j+7] >= 0) Q[ijk] -= (1.0/90.0)*vector[x0+BCY_mg[index2+9*j+7]+z0]/dy2;
    if (BCZ_mg[index3+9*k+1] >= 0) Q[ijk] -= (1.0/90.0)*vector[x0+y0+BCZ_mg[index3+9*k+1]]/dz2;
    if (BCZ_mg[index3+9*k+2] >= 0) Q[ijk] += (3.0/20.0)*vector[x0+y0+BCZ_mg[index3+9*k+2]]/dz2;
    if (BCZ_mg[index3+9*k+3] >= 0) Q[ijk] -= (3.0/2.0)*vector[x0+y0+BCZ_mg[index3+9*k+3]]/dz2;
    if (BCZ_mg[index3+9*k+5] >= 0) Q[ijk] -= (3.0/2.0)*vector[x0+y0+BCZ_mg[index3+9*k+5]]/dz2;
    if (BCZ_mg[index3+9*k+6] >= 0) Q[ijk] += (3.0/20.0)*vector[x0+y0+BCZ_mg[index3+9*k+6]]/dz2;
    if (BCZ_mg[index3+9*k+7] >= 0) Q[ijk] -= (1.0/90.0)*vector[x0+y0+BCZ_mg[index3+9*k+7]]/dz2;
    }
    else if (fdif_order == 8){ 
    Q[ijk] =  (205.0/72.0)*(1.0/dx2 + 1.0/dy2 + 1.0/dz2)*vector[ijk];
    if (BCX_mg[index1+9*i+0] >= 0) Q[ijk] += (1.0/560.0)*vector[BCX_mg[index1+9*i+0]+y0+z0]/dx2;
    if (BCX_mg[index1+9*i+1] >= 0) Q[ijk] -= (8.0/315.0)*vector[BCX_mg[index1+9*i+1]+y0+z0]/dx2;
    if (BCX_mg[index1+9*i+2] >= 0) Q[ijk] += (1.0/5.0)*vector[BCX_mg[index1+9*i+2]+y0+z0]/dx2;
    if (BCX_mg[index1+9*i+3] >= 0) Q[ijk] -= (8.0/5.0)*vector[BCX_mg[index1+9*i+3]+y0+z0]/dx2;
    if (BCX_mg[index1+9*i+5] >= 0) Q[ijk] -= (8.0/5.0)*vector[BCX_mg[index1+9*i+5]+y0+z0]/dx2;
    if (BCX_mg[index1+9*i+6] >= 0) Q[ijk] += (1.0/5.0)*vector[BCX_mg[index1+9*i+6]+y0+z0]/dx2;
    if (BCX_mg[index1+9*i+7] >= 0) Q[ijk] -= (8.0/315.0)*vector[BCX_mg[index1+9*i+7]+y0+z0]/dx2;
    if (BCX_mg[index1+9*i+8] >= 0) Q[ijk] += (1.0/560.0)*vector[BCX_mg[index1+9*i+8]+y0+z0]/dx2;
    if (BCY_mg[index2+9*j+0] >= 0) Q[ijk] += (1.0/560.0)*vector[x0+BCY_mg[index2+9*j+0]+z0]/dy2;
    if (BCY_mg[index2+9*j+1] >= 0) Q[ijk] -= (8.0/315.0)*vector[x0+BCY_mg[index2+9*j+1]+z0]/dy2;
    if (BCY_mg[index2+9*j+2] >= 0) Q[ijk] += (1.0/5.0)*vector[x0+BCY_mg[index2+9*j+2]+z0]/dy2;
    if (BCY_mg[index2+9*j+3] >= 0) Q[ijk] -= (8.0/5.0)*vector[x0+BCY_mg[index2+9*j+3]+z0]/dy2;
    if (BCY_mg[index2+9*j+5] >= 0) Q[ijk] -= (8.0/5.0)*vector[x0+BCY_mg[index2+9*j+5]+z0]/dy2;
    if (BCY_mg[index2+9*j+6] >= 0) Q[ijk] += (1.0/5.0)*vector[x0+BCY_mg[index2+9*j+6]+z0]/dy2;
    if (BCY_mg[index2+9*j+7] >= 0) Q[ijk] -= (8.0/315.0)*vector[x0+BCY_mg[index2+9*j+7]+z0]/dy2;
    if (BCY_mg[index2+9*j+8] >= 0) Q[ijk] += (1.0/560.0)*vector[x0+BCY_mg[index2+9*j+8]+z0]/dy2;
    if (BCZ_mg[index3+9*k+0] >= 0) Q[ijk] += (1.0/560.0)*vector[x0+y0+BCZ_mg[index3+9*k+0]]/dz2;
    if (BCZ_mg[index3+9*k+1] >= 0) Q[ijk] -= (8.0/315.0)*vector[x0+y0+BCZ_mg[index3+9*k+1]]/dz2;
    if (BCZ_mg[index3+9*k+2] >= 0) Q[ijk] += (1.0/5.0)*vector[x0+y0+BCZ_mg[index3+9*k+2]]/dz2;
    if (BCZ_mg[index3+9*k+3] >= 0) Q[ijk] -= (8.0/5.0)*vector[x0+y0+BCZ_mg[index3+9*k+3]]/dz2;
    if (BCZ_mg[index3+9*k+5] >= 0) Q[ijk] -= (8.0/5.0)*vector[x0+y0+BCZ_mg[index3+9*k+5]]/dz2;
    if (BCZ_mg[index3+9*k+6] >= 0) Q[ijk] += (1.0/5.0)*vector[x0+y0+BCZ_mg[index3+9*k+6]]/dz2;
    if (BCZ_mg[index3+9*k+7] >= 0) Q[ijk] -= (8.0/315.0)*vector[x0+y0+BCZ_mg[index3+9*k+7]]/dz2;
    if (BCZ_mg[index3+9*k+8] >= 0) Q[ijk] += (1.0/560.0)*vector[x0+y0+BCZ_mg[index3+9*k+8]]/dz2;
    }
  }
  }
}

void PEQSclass::FinDif1stDeriv(double *func, double *deriv)
{
  int i,j,k,ijk,Nx=NPoints[0],Ny=NPoints[1],Nz=NPoints[2];
  double dx=delX[0],dy=delX[1],dz=delX[2];
  VRload(deriv,3*NTotPts,0.0);
  #pragma omp parallel \
  private(i,j,k,ijk) \
  firstprivate(dx,dy,dz) 
  {
  #pragma omp for
  for (ijk=0;ijk<NTotPts;ijk++){
    get_3ind_from_multInd(i,j,k,ijk);
    int x0=BCX[9*i+4]; int y0=BCY[9*j+4]; int z0=BCZ[9*k+4];
    if (fdif_order == 2){
    if (BCX[9*i+3] >= 0) deriv[3*ijk+0] -= (1.0/2.0)*func[BCX[9*i+3]+y0+z0]/dx;
    if (BCX[9*i+5] >= 0) deriv[3*ijk+0] += (1.0/2.0)*func[BCX[9*i+5]+y0+z0]/dx;
    if (BCY[9*j+3] >= 0) deriv[3*ijk+1] -= (1.0/2.0)*func[x0+BCY[9*j+3]+z0]/dy;
    if (BCY[9*j+5] >= 0) deriv[3*ijk+1] += (1.0/2.0)*func[x0+BCY[9*j+5]+z0]/dy;
    if (BCZ[9*k+3] >= 0) deriv[3*ijk+2] -= (1.0/2.0)*func[x0+y0+BCZ[9*k+3]]/dz;
    if (BCZ[9*k+5] >= 0) deriv[3*ijk+2] += (1.0/2.0)*func[x0+y0+BCZ[9*k+5]]/dz;
    if (i == 0 || i == (Nx-1))
         deriv[3*ijk+0] = 0.0;
    if (j == 0 || j == (Ny-1))
         deriv[3*ijk+1] = 0.0;
    if (k == 0 || k == (Nz-1))
         deriv[3*ijk+2] = 0.0;
    }
    else if (fdif_order == 4){
    if (BCX[9*i+2] >= 0) deriv[3*ijk+0] += (1.0/12.0)*func[BCX[9*i+2]+y0+z0]/dx;
    if (BCX[9*i+3] >= 0) deriv[3*ijk+0] -= (2.0/3.0)*func[BCX[9*i+3]+y0+z0]/dx;
    if (BCX[9*i+5] >= 0) deriv[3*ijk+0] += (2.0/3.0)*func[BCX[9*i+5]+y0+z0]/dx;
    if (BCX[9*i+6] >= 0) deriv[3*ijk+0] -= (1.0/12.0)*func[BCX[9*i+6]+y0+z0]/dx;
    if (BCY[9*j+2] >= 0) deriv[3*ijk+1] += (1.0/12.0)*func[x0+BCY[9*j+2]+z0]/dy;
    if (BCY[9*j+3] >= 0) deriv[3*ijk+1] -= (2.0/3.0)*func[x0+BCY[9*j+3]+z0]/dy;
    if (BCY[9*j+5] >= 0) deriv[3*ijk+1] += (2.0/3.0)*func[x0+BCY[9*j+5]+z0]/dy;
    if (BCY[9*j+6] >= 0) deriv[3*ijk+1] -= (1.0/12.0)*func[x0+BCY[9*j+6]+z0]/dy;
    if (BCZ[9*k+2] >= 0) deriv[3*ijk+2] += (1.0/12.0)*func[x0+y0+BCZ[9*k+2]]/dz;
    if (BCZ[9*k+3] >= 0) deriv[3*ijk+2] -= (2.0/3.0)*func[x0+y0+BCZ[9*k+3]]/dz;
    if (BCZ[9*k+5] >= 0) deriv[3*ijk+2] += (2.0/3.0)*func[x0+y0+BCZ[9*k+5]]/dz;
    if (BCZ[9*k+6] >= 0) deriv[3*ijk+2] -= (1.0/12.0)*func[x0+y0+BCZ[9*k+6]]/dz;
    if (i == 0 || i == 1 || i == (Nx-1) || i == (Nx-2))
         deriv[3*ijk+0] = 0.0;
    if (j == 0 || j == 1 || j == (Ny-1) || j == (Ny-2))
         deriv[3*ijk+1] = 0.0;
    if (k == 0 || k == 1 || k == (Nz-1) || k == (Nz-2))
         deriv[3*ijk+2] = 0.0;
    }
    else if (fdif_order == 6){
    if (BCX[9*i+1] >= 0) deriv[3*ijk+0] -= (1.0/60.0)*func[BCX[9*i+1]+y0+z0]/dx;
    if (BCX[9*i+2] >= 0) deriv[3*ijk+0] += (3.0/20.0)*func[BCX[9*i+2]+y0+z0]/dx;
    if (BCX[9*i+3] >= 0) deriv[3*ijk+0] -= (3.0/4.0)*func[BCX[9*i+3]+y0+z0]/dx;
    if (BCX[9*i+5] >= 0) deriv[3*ijk+0] += (3.0/4.0)*func[BCX[9*i+5]+y0+z0]/dx;
    if (BCX[9*i+6] >= 0) deriv[3*ijk+0] -= (3.0/20.0)*func[BCX[9*i+6]+y0+z0]/dx;
    if (BCX[9*i+7] >= 0) deriv[3*ijk+0] += (1.0/60.0)*func[BCX[9*i+7]+y0+z0]/dx;
    if (BCY[9*j+1] >= 0) deriv[3*ijk+1] -= (1.0/60.0)*func[x0+BCY[9*j+1]+z0]/dy;
    if (BCY[9*j+2] >= 0) deriv[3*ijk+1] += (3.0/20.0)*func[x0+BCY[9*j+2]+z0]/dy;
    if (BCY[9*j+3] >= 0) deriv[3*ijk+1] -= (3.0/4.0)*func[x0+BCY[9*j+3]+z0]/dy;
    if (BCY[9*j+5] >= 0) deriv[3*ijk+1] += (3.0/4.0)*func[x0+BCY[9*j+5]+z0]/dy;
    if (BCY[9*j+6] >= 0) deriv[3*ijk+1] -= (3.0/20.0)*func[x0+BCY[9*j+6]+z0]/dy;
    if (BCY[9*j+7] >= 0) deriv[3*ijk+1] += (1.0/60.0)*func[x0+BCY[9*j+7]+z0]/dy;
    if (BCZ[9*k+1] >= 0) deriv[3*ijk+2] -= (1.0/60.0)*func[x0+y0+BCZ[9*k+1]]/dz;
    if (BCZ[9*k+2] >= 0) deriv[3*ijk+2] += (3.0/20.0)*func[x0+y0+BCZ[9*k+2]]/dz;
    if (BCZ[9*k+3] >= 0) deriv[3*ijk+2] -= (3.0/4.0)*func[x0+y0+BCZ[9*k+3]]/dz;
    if (BCZ[9*k+5] >= 0) deriv[3*ijk+2] += (3.0/4.0)*func[x0+y0+BCZ[9*k+5]]/dz;
    if (BCZ[9*k+6] >= 0) deriv[3*ijk+2] -= (3.0/20.0)*func[x0+y0+BCZ[9*k+6]]/dz;
    if (BCZ[9*k+7] >= 0) deriv[3*ijk+2] += (1.0/60.0)*func[x0+y0+BCZ[9*k+7]]/dz;
    if (i == 0 || i == 1 || i == 2 || i == (Nx-1) || i == (Nx-2) || i == (Nx-3))
         deriv[3*ijk+0] = 0.0;
    if (j == 0 || j == 1 || j == 2 || j == (Ny-1) || j == (Ny-2) || j == (Ny-3))
         deriv[3*ijk+1] = 0.0;
    if (k == 0 || k == 1 || k == 2 || k == (Nz-1) || k == (Nz-2) || k == (Nz-3))
         deriv[3*ijk+2] = 0.0;
    }
    else if (fdif_order == 8){
    if (BCX[9*i+0] >= 0) deriv[3*ijk+0] += (1.0/280.0)*func[BCX[9*i+0]+y0+z0]/dx;
    if (BCX[9*i+1] >= 0) deriv[3*ijk+0] -= (4.0/105.0)*func[BCX[9*i+1]+y0+z0]/dx;
    if (BCX[9*i+2] >= 0) deriv[3*ijk+0] += (1.0/5.0)*func[BCX[9*i+2]+y0+z0]/dx;
    if (BCX[9*i+3] >= 0) deriv[3*ijk+0] -= (4.0/5.0)*func[BCX[9*i+3]+y0+z0]/dx;
    if (BCX[9*i+5] >= 0) deriv[3*ijk+0] += (4.0/5.0)*func[BCX[9*i+5]+y0+z0]/dx;
    if (BCX[9*i+6] >= 0) deriv[3*ijk+0] -= (1.0/5.0)*func[BCX[9*i+6]+y0+z0]/dx;
    if (BCX[9*i+7] >= 0) deriv[3*ijk+0] += (4.0/105.0)*func[BCX[9*i+7]+y0+z0]/dx;
    if (BCX[9*i+8] >= 0) deriv[3*ijk+0] -= (1.0/280.0)*func[BCX[9*i+8]+y0+z0]/dx;
    if (BCY[9*j+0] >= 0) deriv[3*ijk+1] += (1.0/280.0)*func[x0+BCY[9*j+0]+z0]/dy;
    if (BCY[9*j+1] >= 0) deriv[3*ijk+1] -= (4.0/105.0)*func[x0+BCY[9*j+1]+z0]/dy;
    if (BCY[9*j+2] >= 0) deriv[3*ijk+1] += (1.0/5.0)*func[x0+BCY[9*j+2]+z0]/dy;
    if (BCY[9*j+3] >= 0) deriv[3*ijk+1] -= (4.0/5.0)*func[x0+BCY[9*j+3]+z0]/dy;
    if (BCY[9*j+5] >= 0) deriv[3*ijk+1] += (4.0/5.0)*func[x0+BCY[9*j+5]+z0]/dy;
    if (BCY[9*j+6] >= 0) deriv[3*ijk+1] -= (1.0/5.0)*func[x0+BCY[9*j+6]+z0]/dy;
    if (BCY[9*j+7] >= 0) deriv[3*ijk+1] += (4.0/105.0)*func[x0+BCY[9*j+7]+z0]/dy;
    if (BCY[9*j+8] >= 0) deriv[3*ijk+1] -= (1.0/280.0)*func[x0+BCY[9*j+8]+z0]/dy;
    if (BCZ[9*k+0] >= 0) deriv[3*ijk+2] += (1.0/280.0)*func[x0+y0+BCZ[9*k+0]]/dz;
    if (BCZ[9*k+1] >= 0) deriv[3*ijk+2] -= (4.0/105.0)*func[x0+y0+BCZ[9*k+1]]/dz;
    if (BCZ[9*k+2] >= 0) deriv[3*ijk+2] += (1.0/5.0)*func[x0+y0+BCZ[9*k+2]]/dz;
    if (BCZ[9*k+3] >= 0) deriv[3*ijk+2] -= (4.0/5.0)*func[x0+y0+BCZ[9*k+3]]/dz;
    if (BCZ[9*k+5] >= 0) deriv[3*ijk+2] += (4.0/5.0)*func[x0+y0+BCZ[9*k+5]]/dz;
    if (BCZ[9*k+6] >= 0) deriv[3*ijk+2] -= (1.0/5.0)*func[x0+y0+BCZ[9*k+6]]/dz;
    if (BCZ[9*k+7] >= 0) deriv[3*ijk+2] += (4.0/105.0)*func[x0+y0+BCZ[9*k+7]]/dz;
    if (BCZ[9*k+8] >= 0) deriv[3*ijk+2] -= (1.0/280.0)*func[x0+y0+BCZ[9*k+8]]/dz;
    if (i == 0 || i == 1 || i == 2 || i == 3 ||
        i == (Nx-1) || i == (Nx-2) || i == (Nx-3) || i == (Nx-4))
         deriv[3*ijk+0] = 0.0;
    if (j == 0 || j == 1 || j == 2 || j == 3 ||
        j == (Ny-1) || j == (Ny-2) || j == (Ny-3) || j == (Ny-4))
         deriv[3*ijk+1] = 0.0;
    if (k == 0 || k == 1 || k == 2 || k == 3 ||
        k == (Nz-1) || k == (Nz-2) || k == (Nz-3) || k == (Nz-4))
         deriv[3*ijk+2] = 0.0;
    }
  }
  }
}
