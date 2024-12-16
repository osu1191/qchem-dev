#include "peqs_class.h"

bool PEQSclass::MultiGridControl(int system)
{
  double resid = 1.0;
  if ((multigrid && PBcallNo == 0 && (system == 1 || system == 2))
       || (multigrid && (system == 3 || system == 4 ))){
    for (int i=0; i<sweep1; i++){
      resid = ConjGrad(1,system);
      if (converge_cg) break;
    }
    ResetCGIter();SetMGLevel(3);
    for (int i=0; i<sweep1; i++) resid = ConjGrad(1,system);
    ResetCGIter();SetMGLevel(2);
    for (int i=0; i<sweep1; i++) resid = ConjGrad(1,system);
    ResetCGIter();SetMGLevel(1);
    for (int i=0; i<max_iter; i++){
      resid = ConjGrad(1,system);
      if (converge_cg) break;
    }
    if (mg_cycle == 2){
      ResetCGIter();SetMGLevel(2);
      for (int i=0; i<sweep3; i++) resid = ConjGrad(2,system);
      ResetCGIter();SetMGLevel(1);
      for (int i=0; i<max_iter; i++){
        resid = ConjGrad(1,system);
        if (converge_cg) break;
      }
      ResetCGIter();SetMGLevel(2);
      for (int i=0; i<sweep2; i++) resid = ConjGrad(2,system);
      ResetCGIter();SetMGLevel(3);
      for (int i=0; i<sweep3; i++) resid = ConjGrad(2,system);
      ResetCGIter();SetMGLevel(2);
      for (int i=0; i<sweep1; i++) resid = ConjGrad(1,system);
      ResetCGIter();SetMGLevel(1);
      for (int i=0; i<max_iter; i++){
      resid = ConjGrad(1,system);
        if (converge_cg) break;
      }
      ResetCGIter();SetMGLevel(2);
      for (int i=0; i<sweep3; i++) resid = ConjGrad(2,system);
      ResetCGIter();SetMGLevel(1);
      for (int i=0; i<max_iter; i++){
        resid = ConjGrad(1,system);
        if (converge_cg) break;
      }
    }
    ResetCGIter();SetMGLevel(2);
    for (int i=0; i<sweep2; i++) resid = ConjGrad(2,system);
    ResetCGIter();SetMGLevel(3);
    for (int i=0; i<sweep2; i++) resid = ConjGrad(2,system);
    ResetCGIter();SetMGLevel(4);
    for (int i=0; i<sweep2; i++){
      resid = ConjGrad(2,system);
      if (converge_cg) break;
    }
  }
  else{
    ResetCGIter();
    for (int i=0; i<max_iter; i++){
      resid = ConjGrad(0,system);
      if (converge_cg) break;
    }
  }
  return converge_cg;
}

void PEQSclass::RestrictError(double *rstrct_e, double *rho, double *phi)
{
  int i,j,k,ijk,i_2,j_2,k_2,ijk_2,Nx,Ny,Nz,index1=0,index2=0,index3=0,mg_index=mg_order-mg_level,rstrct_index=mg_index-1;
  double *error=QAllocDouble(NTotPts_mg[rstrct_index]),*temp_del2phi=QAllocDouble(NTotPts_mg[rstrct_index]);
  ComputeCharge(temp_del2phi,phi,true);
  VRsub(error,rho,temp_del2phi,NTotPts_mg[rstrct_index]);
  Nx = NPoints[3*rstrct_index+0]; Ny = NPoints[3*rstrct_index+1]; Nz = NPoints[3*rstrct_index+2];
  for (ijk=0;ijk<rstrct_index;ijk++){
    index1 += 9*NPoints[3*ijk+0];
    index2 += 9*NPoints[3*ijk+1];
    index3 += 9*NPoints[3*ijk+2];
  }
  #pragma omp parallel for private(i,j,k,i_2,j_2,k_2,ijk,ijk_2) firstprivate(index1,index2,index3)
  for (i=0;i<Nx;i+=2){
    for (j=0;j<Ny;j+=2){
      for (k=0;k<Nz;k+=2){
        ijk = get_multInd(i,j,k,true,false);
        i_2 = i/2; j_2 = j/2; k_2 = k/2;
        ijk_2 = get_multInd(i_2,j_2,k_2,false,false);
        int x0=BCX_mg[index1+9*i+4]; int y0=BCY_mg[index2+9*j+4]; int z0=BCZ_mg[index3+9*k+4];
        rstrct_e[ijk_2] = 0.125*error[ijk];
        if (BCX_mg[index1+9*i+3] >= 0)
          rstrct_e[ijk_2] += 0.0625*error[BCX_mg[index1+9*i+3]+y0+z0];
        if (BCX_mg[index1+9*i+5] >= 0)
          rstrct_e[ijk_2] += 0.0625*error[BCX_mg[index1+9*i+5]+y0+z0];
        if (BCY_mg[index2+9*j+3] >= 0)
          rstrct_e[ijk_2] += 0.0625*error[x0+BCY_mg[index2+9*j+3]+z0];
        if (BCY_mg[index2+9*j+5] >= 0)
          rstrct_e[ijk_2] += 0.0625*error[x0+BCY_mg[index2+9*j+5]+z0];
        if (BCZ_mg[index3+9*k+3] >= 0)
          rstrct_e[ijk_2] += 0.0625*error[x0+y0+BCZ_mg[index3+9*k+3]];
        if (BCZ_mg[index3+9*k+5] >= 0)
          rstrct_e[ijk_2] += 0.0625*error[x0+y0+BCZ_mg[index3+9*k+5]];
        if (BCX_mg[index1+9*i+3] >= 0 && BCY_mg[index2+9*j+3] >= 0)
          rstrct_e[ijk_2] += 0.03125*error[BCX_mg[index1+9*i+3]+BCY_mg[index2+9*j+3]+z0];
        if (BCX_mg[index1+9*i+3] >= 0 && BCY_mg[index2+9*j+5] >= 0)
          rstrct_e[ijk_2] += 0.03125*error[BCX_mg[index1+9*i+3]+BCY_mg[index2+9*j+5]+z0];
        if (BCX_mg[index1+9*i+5] >= 0 && BCY_mg[index2+9*j+3] >= 0)
          rstrct_e[ijk_2] += 0.03125*error[BCX_mg[index1+9*i+5]+BCY_mg[index2+9*j+3]+z0];
        if (BCX_mg[index1+9*i+5] >= 0 && BCY_mg[index2+9*j+5] >= 0)
          rstrct_e[ijk_2] += 0.03125*error[BCX_mg[index1+9*i+5]+BCY_mg[index2+9*j+5]+z0];
        if (BCX_mg[index1+9*i+3] >= 0 && BCZ_mg[index3+9*k+3] >= 0)
          rstrct_e[ijk_2] += 0.03125*error[BCX_mg[index1+9*i+3]+y0+BCZ_mg[index3+9*k+3]];
        if (BCX_mg[index1+9*i+3] >= 0 && BCZ_mg[index3+9*k+5] >= 0)
          rstrct_e[ijk_2] += 0.03125*error[BCX_mg[index1+9*i+3]+y0+BCZ_mg[index3+9*k+5]];
        if (BCX_mg[index1+9*i+5] >= 0 && BCZ_mg[index3+9*k+3] >= 0)
          rstrct_e[ijk_2] += 0.03125*error[BCX_mg[index1+9*i+5]+y0+BCZ_mg[index3+9*k+3]];
        if (BCX_mg[index1+9*i+5] >= 0 && BCZ_mg[index3+9*k+5] >= 0)
          rstrct_e[ijk_2] += 0.03125*error[BCX_mg[index1+9*i+5]+y0+BCZ_mg[index3+9*k+5]];
        if (BCY_mg[index2+9*j+3] >= 0 && BCZ_mg[index3+9*k+3] >= 0)
          rstrct_e[ijk_2] += 0.03125*error[x0+BCY_mg[index2+9*j+3]+BCZ_mg[index3+9*k+3]];
        if (BCY_mg[index2+9*j+3] >= 0 && BCZ_mg[index3+9*k+5] >= 0)
          rstrct_e[ijk_2] += 0.03125*error[x0+BCY_mg[index2+9*j+3]+BCZ_mg[index3+9*k+5]];
        if (BCY_mg[index2+9*j+5] >= 0 && BCZ_mg[index3+9*k+3] >= 0)
          rstrct_e[ijk_2] += 0.03125*error[x0+BCY_mg[index2+9*j+5]+BCZ_mg[index3+9*k+3]];
        if (BCY_mg[index2+9*j+5] >= 0 && BCZ_mg[index3+9*k+5] >= 0)
          rstrct_e[ijk_2] += 0.03125*error[x0+BCY_mg[index2+9*j+5]+BCZ_mg[index3+9*k+5]];
        if (BCX_mg[index1+9*i+3] >= 0 && BCY_mg[index2+9*j+3] >= 0 && BCZ_mg[index3+9*k+3] >= 0)
          rstrct_e[ijk_2] += 0.015625*error[BCX_mg[index1+9*i+3]+BCY_mg[index2+9*j+3]+BCZ_mg[index3+9*k+3]];
        if (BCX_mg[index1+9*i+3] >= 0 && BCY_mg[index2+9*j+3] >= 0 && BCZ_mg[index3+9*k+5] >= 0)
          rstrct_e[ijk_2] += 0.015625*error[BCX_mg[index1+9*i+3]+BCY_mg[index2+9*j+3]+BCZ_mg[index3+9*k+5]];
        if (BCX_mg[index1+9*i+3] >= 0 && BCY_mg[index2+9*j+5] >= 0 && BCZ_mg[index3+9*k+3] >= 0)
          rstrct_e[ijk_2] += 0.015625*error[BCX_mg[index1+9*i+3]+BCY_mg[index2+9*j+5]+BCZ_mg[index3+9*k+3]];
        if (BCX_mg[index1+9*i+3] >= 0 && BCY_mg[index2+9*j+5] >= 0 && BCZ_mg[index3+9*k+5] >= 0)
          rstrct_e[ijk_2] += 0.015625*error[BCX_mg[index1+9*i+3]+BCY_mg[index2+9*j+5]+BCZ_mg[index3+9*k+5]];
        if (BCX_mg[index1+9*i+5] >= 0 && BCY_mg[index2+9*j+3] >= 0 && BCZ_mg[index3+9*k+3] >= 0)
          rstrct_e[ijk_2] += 0.015625*error[BCX_mg[index1+9*i+5]+BCY_mg[index2+9*j+3]+BCZ_mg[index3+9*k+3]];
        if (BCX_mg[index1+9*i+5] >= 0 && BCY_mg[index2+9*j+3] >= 0 && BCZ_mg[index3+9*k+5] >= 0)
          rstrct_e[ijk_2] += 0.015625*error[BCX_mg[index1+9*i+5]+BCY_mg[index2+9*j+3]+BCZ_mg[index3+9*k+5]];
        if (BCX_mg[index1+9*i+5] >= 0 && BCY_mg[index2+9*j+5] >= 0 && BCZ_mg[index3+9*k+3] >= 0)
          rstrct_e[ijk_2] += 0.015625*error[BCX_mg[index1+9*i+5]+BCY_mg[index2+9*j+5]+BCZ_mg[index3+9*k+3]];
        if (BCX_mg[index1+9*i+5] >= 0 && BCY_mg[index2+9*j+5] >= 0 && BCZ_mg[index3+9*k+5] >= 0)
          rstrct_e[ijk_2] += 0.015625*error[BCX_mg[index1+9*i+5]+BCY_mg[index2+9*j+5]+BCZ_mg[index3+9*k+5]];
      }
    }
  }
  QFree(error);QFree(temp_del2phi);
}

void PEQSclass::InterpolError(double *phi, double *phi_corr)
{
  int i,j,k,ijk,i2,j2,k2,ijk2,Nx,Ny,Nz,index1=0,index2=0,index3=0,index1_2=0,index2_2=0,index3_2=0,mg_index=mg_order-mg_level,interpol_index=mg_index+1;
  double *interpolation=QAllocDouble(NTotPts_mg[mg_index]);
  Nx = NPoints[3*interpol_index+0]; Ny = NPoints[3*interpol_index+1]; Nz = NPoints[3*interpol_index+2];
  for (ijk=0;ijk<interpol_index;ijk++){
    index1 += 9*NPoints[3*ijk+0];
    index2 += 9*NPoints[3*ijk+1];
    index3 += 9*NPoints[3*ijk+2];
  }
  for (ijk=0;ijk<mg_index;ijk++){
    index1_2 += 9*NPoints[3*ijk+0];
    index2_2 += 9*NPoints[3*ijk+1];
    index3_2 += 9*NPoints[3*ijk+2];
  }
  #pragma omp parallel for private(i,j,k,i2,j2,k2,ijk,ijk2) firstprivate(index1,index2,index3)
  for (i=0;i<Nx;i++){
    for (j=0;j<Ny;j++){
      for (k=0;k<Nz;k++){
        ijk = get_multInd(i,j,k,false,true);
        i2 = 2*i; j2 = 2*j; k2 = 2*k;
        ijk2 = get_multInd(i2,j2,k2,false,false);
        int x0=BCX_mg[index1+9*i+4]; int y0=BCY_mg[index2+9*j+4]; int z0=BCZ_mg[index3+9*k+4];
        int x0_2=BCX_mg[index1_2+9*i2+4]; int y0_2=BCY_mg[index2_2+9*j2+4]; int z0_2=BCZ_mg[index3_2+9*k2+4];
        interpolation[ijk2] = phi_corr[ijk];
        if (BCX_mg[index1+9*i+5] >= 0 && BCX_mg[index1_2+9*i2+5] >= 0){
          interpolation[BCX_mg[index1_2+9*i2+5]+y0_2+z0_2] = 0.5*(phi_corr[BCX_mg[index1+9*i+5]+y0+z0]+phi_corr[ijk]);
        }
        if (BCY_mg[index2+9*j+5] >= 0 && BCY_mg[index2_2+9*j2+5] >= 0){
          interpolation[x0_2+BCY_mg[index2_2+9*j2+5]+z0_2] = 0.5*(phi_corr[x0+BCY_mg[index2+9*j+5]+z0]+phi_corr[ijk]);
        }
        if (BCZ_mg[index3+9*k+5] >= 0 && BCZ_mg[index3_2+9*k2+5] >= 0){
          interpolation[x0_2+y0_2+BCZ_mg[index3_2+9*k2+5]] = 0.5*(phi_corr[x0+y0+BCZ_mg[index3+9*k+5]]+phi_corr[ijk]);
        }
        if (BCX_mg[index1+9*i+5] >= 0 && BCY_mg[index2+9*j+5] >= 0 &&
            BCX_mg[index1_2+9*i2+5] >= 0 && BCY_mg[index2_2+9*j2+5] >= 0){
          interpolation[BCX_mg[index1_2+9*i2+5]+BCY_mg[index2_2+9*j2+5]+z0_2] =
            0.25*(phi_corr[BCX_mg[index1+9*i+5]+BCY_mg[index2+9*j+5]+z0]+
            phi_corr[BCX_mg[index1+9*i+5]+y0+z0]+phi_corr[x0+BCY_mg[index2+9*j+5]+z0]+phi_corr[ijk]);
        }
        if (BCX_mg[index1+9*i+5] >= 0 && BCZ_mg[index3+9*k+5] >= 0 &&
            BCX_mg[index1_2+9*i2+5] >= 0 && BCZ_mg[index3_2+9*k2+5] >= 0){
          interpolation[BCX_mg[index1_2+9*i2+5]+y0_2+BCZ_mg[index3_2+9*k2+5]] =
            0.25*(phi_corr[BCX_mg[index1+9*i+5]+y0+BCZ_mg[index3+9*k+5]]+
            phi_corr[BCX_mg[index1+9*i+5]+y0+z0]+phi_corr[x0+y0+BCZ_mg[index3+9*k+5]]+phi_corr[ijk]);
        }
        if (BCY_mg[index2+9*j+5] >= 0 && BCZ_mg[index3+9*k+5] >= 0 &&
            BCY_mg[index2_2+9*j2+5] >= 0 && BCZ_mg[index3_2+9*k2+5] >= 0){
          interpolation[x0_2+BCY_mg[index2_2+9*j2+5]+BCZ_mg[index3_2+9*k2+5]] =
            0.25*(phi_corr[x0+BCY_mg[index2+9*j+5]+BCZ_mg[index3+9*k+5]]+
            phi_corr[x0+BCY_mg[index2+9*j+5]+z0]+phi_corr[x0+y0+BCZ_mg[index3+9*k+5]]+phi_corr[ijk]);
        }
        if (BCX_mg[index1+9*i+5] >= 0 && BCY_mg[index2+9*j+5] >= 0 && BCZ_mg[index3+9*k+5] >= 0 &&
            BCX_mg[index1_2+9*i2+5] >= 0 && BCY_mg[index2_2+9*j2+5] >= 0 && BCZ_mg[index3_2+9*k2+5] >= 0){
          interpolation[BCX_mg[index1_2+9*i2+5]+BCY_mg[index2_2+9*j2+5]+BCZ_mg[index3_2+9*k2+5]] =
            0.125*(phi_corr[BCX_mg[index1+9*i+5]+BCY_mg[index2+9*j+5]+BCZ_mg[index3+9*k+5]]+
            phi_corr[BCX_mg[index1+9*i+5]+BCY_mg[index2+9*j+5]+z0]+phi_corr[BCX_mg[index1+9*i+5]+y0+BCZ_mg[index3+9*k+5]]+
            phi_corr[x0+BCY_mg[index2+9*j+5]+BCZ_mg[index3+9*k+5]]+phi_corr[x0+y0+BCZ_mg[index3+9*k+5]]+
            phi_corr[x0+BCY_mg[index2+9*j+5]+z0]+phi_corr[BCX_mg[index1+9*i+5]+y0+z0]+phi_corr[ijk]);
        }
      }
    }
  }
  VRadd(phi,phi,interpolation,NTotPts_mg[mg_index]);
  QFree(interpolation);
}
