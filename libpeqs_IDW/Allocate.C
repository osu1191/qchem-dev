#include "peqs_class.h"

void PEQSclass::Allocate()
{
  Aphi = QAllocDouble(NTotPts);
  Abasis = QAllocDouble(NTotPts);
  basis = QAllocDouble(NTotPts);
  resid = QAllocDouble(NTotPts);
  zvec = QAllocDouble(NTotPts);
  BCX = QAllocINTEGER(9*NPoints[0]);
  BCY = QAllocINTEGER(9*NPoints[1]);
  BCZ = QAllocINTEGER(9*NPoints[2]);
  vdWScale = QAllocDouble(NTotPts);
  epsilon = QAllocDouble(NTotPts);
  log_epsilon = QAllocDouble(NTotPts);
  delLogEps = QAllocDouble(3*NTotPts);
  delPhi = QAllocDouble(3*NTotPts);
  rho_total = QAllocDouble(NTotPts);
  phi_total = QAllocDouble(NTotPts);
  rho_elec = QAllocDouble(NTotPts);
  phi_elec = QAllocDouble(NTotPts);
  rho_nuke = QAllocDouble(NTotPts);
  phi_nuke = QAllocDouble(NTotPts);
  valCart = QAllocDouble(NTotPts);
  rho_solute = QAllocDouble(NTotPts);
  phi_solute = QAllocDouble(NTotPts);
  rho_iter = QAllocDouble(NTotPts);
  rho_polar = QAllocDouble(NTotPts);
  phi_polar = QAllocDouble(NTotPts);

//============= Interpolation =====================//

  val3 = QAllocDouble(NTotPts); 
  mono = QAllocDouble(NTotPts); 
  dip = QAllocDouble(NTotPts);  
  quad = QAllocDouble(NTotPts); 
  octa = QAllocDouble(NTotPts); 
  hexa = QAllocDouble(NTotPts); 
  fifth = QAllocDouble(NTotPts);
  sixth = QAllocDouble(NTotPts); 
  seventh = QAllocDouble(NTotPts);
  eight = QAllocDouble(NTotPts);  

//============================================//

//CS: added this for Poisson-Boltzmann code
  lambda = QAllocDouble(NTotPts);
  rho_ions = QAllocDouble(NTotPts);
  phi_ions = QAllocDouble(NTotPts);
  if (doNeqJob == 1 && noneq_partition == 0){
    rhopol_ref_slow = QAllocDouble(NTotPts);
    rhopol_ref_fast = QAllocDouble(NTotPts);
    phipol_ref_slow = QAllocDouble(NTotPts);
  }
  else if (doNeqJob == 1 && noneq_partition == 1){
    if (noneq_state == 0){
      delPhiFast = QAllocDouble(3*NTotPts);
      rhotot_fast = QAllocDouble(NTotPts);
      rhopol_ref_fast = QAllocDouble(NTotPts);
      phitot_fast = QAllocDouble(NTotPts);
      phipol_ref_fast = QAllocDouble(NTotPts);
    }
    rhopol_ref_slow = QAllocDouble(NTotPts);
    phipol_ref_slow = QAllocDouble(NTotPts);
  }
}
