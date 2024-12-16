#include "peqs_class.h"

void PEQSclass::WriteRhoPhi()
{
  if (doNeqJob == 0 || (doNeqJob == 1 && noneq_state == 0)){
    VRscale(phi_nuke,NTotPts,1.0/(4.0*Pi));
    FileMan(FM_WRITE,FILE_REF_PHI,FM_DP,NTotPts,0*NTotPts,FM_BEG,phi_nuke);
    FileMan(FM_WRITE,FILE_REF_RHO,FM_DP,NTotPts,0,FM_BEG,rho_nuke);
  }
  else if (doNeqJob == 1 && noneq_state == 1){
    VRscale(phi_nuke,NTotPts,1.0/(4.0*Pi));
    FileMan(FM_WRITE,FILE_REF_PHI,FM_DP,NTotPts,1*NTotPts,FM_BEG,phi_nuke);
    FileMan(FM_WRITE,FILE_REF_RHO,FM_DP,NTotPts,0*NTotPts,FM_BEG,rho_nuke);
  }
}

void PEQSclass::ReadRhoPhi()
{
  if ((doNeqJob == 0 || (doNeqJob == 1 && noneq_state == 0)) && PBcallNo > 0){
    FileMan(FM_READ,FILE_REF_PHI,FM_DP,NTotPts,0*NTotPts,FM_BEG,phi_nuke);
    FileMan(FM_READ,FILE_REF_RHO,FM_DP,NTotPts,0*NTotPts,FM_BEG,rho_nuke);
  }
  else if ((doNeqJob == 1 && noneq_state == 1) && PBcallNo == 0){
    FileMan(FM_READ,FILE_REF_PHI,FM_DP,NTotPts,0*NTotPts,FM_BEG,phi_nuke);
    FileMan(FM_READ,FILE_REF_RHO,FM_DP,NTotPts,0*NTotPts,FM_BEG,rho_nuke);
  }
  else if ((doNeqJob == 1 && noneq_state == 1) && PBcallNo > 0){
    FileMan(FM_READ,FILE_REF_PHI,FM_DP,NTotPts,1*NTotPts,FM_BEG,phi_nuke);
    FileMan(FM_READ,FILE_REF_RHO,FM_DP,NTotPts,0*NTotPts,FM_BEG,rho_nuke);
  }
}
