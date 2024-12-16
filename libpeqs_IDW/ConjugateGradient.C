#include "peqs_class.h"

double PEQSclass::ConjGrad(int operation, int system)
{
  int ijk,index=mg_order-mg_level;
  double normR=0.0,normQ=0.0,ratio=0.0,top=0.0,bottom=0.0,beta=0.0,diagonal=1.0;
  double *APHI,*PHI,*RESID,*RHO,*ZVEC,*BASIS,*AB;
  if (multigrid){
    if (mg_level == 4){
      if (iter_cg == 0 && operation == 2){
        if (system == 1) InterpolError(phi_elec,phi_corr3);
        else if (system == 2) InterpolError(phi_nuke,phi_corr3);
        else if (system == 3) InterpolError(phi_total,phi_corr3);
        else if (system == 4) InterpolError(phitot_fast,phi_corr3);
      }
      if (system == 1){PHI=phi_elec; RHO=rho_elec;}
      else if (system == 2){PHI=phi_nuke; RHO=rho_nuke;}
      else if (system == 3){PHI=phi_total; RHO=rho_total;}
      else if (system == 4){PHI=phitot_fast; RHO=rhotot_fast;}
      APHI=Aphi; RESID=resid; ZVEC=zvec; BASIS=basis; AB=Abasis;
    }
    else if (mg_level == 3){
      if (iter_cg == 0 && operation == 1){
        if (system == 1) RestrictError(rstrct_3,rho_elec,phi_elec);
        else if (system == 2) RestrictError(rstrct_3,rho_nuke,phi_nuke);
        else if (system == 3) RestrictError(rstrct_3,rho_total,phi_total);
        else if (system == 4) RestrictError(rstrct_3,rhotot_fast,phitot_fast);
      }
      else if (iter_cg == 0 && operation == 2) InterpolError(phi_corr3,phi_corr2);
      APHI=Aphi_3; PHI=phi_corr3; RHO=rstrct_3; RESID=resid_3; ZVEC=zvec_3; BASIS=basis_3; AB=Abasis_3;
    }
    else if (mg_level == 2){
      if (iter_cg == 0 && operation == 1) RestrictError(rstrct_2,rstrct_3,phi_corr3);
      else if (iter_cg == 0 && operation == 2) InterpolError(phi_corr2,phi_corr1);
      APHI=Aphi_2; PHI=phi_corr2; RHO=rstrct_2; RESID=resid_2; ZVEC=zvec_2; BASIS=basis_2; AB=Abasis_2;
    }
    else if (mg_level == 1 && multigrid){
      if (iter_cg == 0 && operation == 1) RestrictError(rstrct_1,rstrct_2,phi_corr2);
      APHI=Aphi_1; PHI=phi_corr1; RHO=rstrct_1; RESID=resid_1; ZVEC=zvec_1; BASIS=basis_1; AB=Abasis_1;
    }
  }
  if (mg_level == 1 && !multigrid){
    if (system == 1){PHI=phi_elec; RHO=rho_elec;}
    else if (system == 2){PHI=phi_nuke; RHO=rho_nuke;}
    else if (system == 3){PHI=phi_total; RHO=rho_total;}
    else if (system == 4){PHI=phitot_fast; RHO=rhotot_fast;}
    APHI=Aphi; RESID=resid; ZVEC=zvec; BASIS=basis; AB=Abasis;
  }

  VRload(APHI,NTotPts_mg[index],0.0);

  ComputeCharge(APHI,PHI,false);
  if (iter_cg == 0){
  VRload(RESID,NTotPts_mg[index],0.0);
  VRload(AB,NTotPts_mg[index],0.0);
    #pragma omp parallel \
    private(ijk) \
    firstprivate(index,diagonal) 
    {
      #pragma omp for
      for (ijk=0;ijk<NTotPts_mg[index];ijk++){
        RESID[ijk] = RHO[ijk] - APHI[ijk];
        ZVEC[ijk] = RESID[ijk]/diagonal;
        BASIS[ijk] = ZVEC[ijk];
      }
    }
  }
  else{
    #pragma omp parallel \
    private(ijk) \
    firstprivate(index,diagonal) \
    reduction(+:bottom,top)
    {
      #pragma omp for
      for (ijk=0;ijk<NTotPts_mg[index];ijk++){
        bottom += RESID[ijk]*ZVEC[ijk];
        RESID[ijk] = RESID[ijk] - alpha_cg*AB[ijk];
        ZVEC[ijk] = RESID[ijk]/diagonal;
        top += RESID[ijk]*ZVEC[ijk];
      }
    }
  }
  #pragma omp parallel \
  private(ijk) \
  firstprivate(index) \
  reduction(+:normR,normQ)
  {
    #pragma omp for
    for (ijk=0;ijk<NTotPts_mg[index];ijk++){
      normR += RESID[ijk]*RESID[ijk];
      normQ += RHO[ijk]*RHO[ijk];
    }
  }
  normR = pow(normR,0.5);
  normQ = pow(normQ,0.5);
  if (normR == 0.0 || normQ == 0.0) normQ = 1.0;
  ratio = normR/normQ;
  if ((ratio - cg_thresh) < 0.0) converge_cg = true;
  else converge_cg = false;
  if (iter_cg > 0 && !converge_cg){
    beta = top/bottom;
    #pragma omp parallel for private(ijk) firstprivate(beta)
    for (ijk=0; ijk<NTotPts_mg[index]; ijk++) BASIS[ijk] = ZVEC[ijk] + beta*BASIS[ijk];
  }
  iter_cg++;
  if (!converge_cg){
    if (iter_cg > max_iter)
      QCrash(" Maximum conjugate gradient iterations reached.");
    top=0.0,bottom=0.0;
    ComputeCharge(AB,BASIS,false);
    #pragma omp parallel \
    private(ijk) \
    firstprivate(index)\
    reduction(+:top,bottom)
    {
      #pragma omp for
      for (ijk=0;ijk<NTotPts_mg[index];ijk++){
        top += RESID[ijk]*ZVEC[ijk];
        bottom += BASIS[ijk]*AB[ijk];
      }
    }
    alpha_cg = top/bottom;
    #pragma omp parallel for private(ijk) 
    for (ijk=0;ijk<NTotPts_mg[index];ijk++) PHI[ijk] = PHI[ijk] + alpha_cg*BASIS[ijk];
  }
  return ratio;
}

