#include "peqs_class.h"
#include "qchem.h"
#include "rem_values.h"
#include "BasisSet.hh"
#include "OneEMtrx.hh"
#include "hirshfeld.h"

extern "C" void cnvtmm(double*,int*);

void PEQSclass::PrepareSolver()
{
  cout << "in PrepSolver()" << endl;
//  octave_testing(); // testing octave in a lame way
  int factorX = (135*NPoints[0]+153)/8;
  int factorY = (135*NPoints[1]+153)/8;
  int factorZ = (135*NPoints[2]+153)/8;
  iter_polar=0,iter_cg=0;
  alpha_cg=0.0;
  VRload(Aphi,NTotPts,0.0);
  VRload(Abasis,NTotPts,0.0);
  VRload(basis,NTotPts,0.0);
  VRload(resid,NTotPts,0.0);
  VRload(zvec,NTotPts,0.0);
  VRload(BCX,9*NPoints[0],-1);
  VRload(BCY,9*NPoints[1],-1);
  VRload(BCZ,9*NPoints[2],-1);
  BCX_mg = QAllocINTEGER(factorX);VRload(BCX_mg,factorX,-1);
  BCY_mg = QAllocINTEGER(factorY);VRload(BCY_mg,factorY,-1);
  BCZ_mg = QAllocINTEGER(factorZ);VRload(BCZ_mg,factorZ,-1);
  if(mg_order > 1){
    Abasis_3 = QAllocDouble(NTotPts_mg[1]);VRload(Abasis_3,NTotPts_mg[1],0.0);
    Aphi_3 = QAllocDouble(NTotPts_mg[1]);VRload(Aphi_3,NTotPts_mg[1],0.0);
    resid_3 = QAllocDouble(NTotPts_mg[1]);VRload(resid_3,NTotPts_mg[1],0.0);
    basis_3 = QAllocDouble(NTotPts_mg[1]);VRload(basis_3,NTotPts_mg[1],0.0);
    zvec_3 = QAllocDouble(NTotPts_mg[1]);VRload(zvec_3,NTotPts_mg[1],0.0);
    rstrct_3 = QAllocDouble(NTotPts_mg[1]);VRload(rstrct_3,NTotPts_mg[1],0.0);
    phi_corr3 = QAllocDouble(NTotPts_mg[1]);VRload(phi_corr3,NTotPts_mg[1],0.0);
  }
  if(mg_order > 2){
    Abasis_2 = QAllocDouble(NTotPts_mg[2]);VRload(Abasis_2,NTotPts_mg[2],0.0);
    Aphi_2 = QAllocDouble(NTotPts_mg[2]);VRload(Aphi_2,NTotPts_mg[2],0.0);
    resid_2 = QAllocDouble(NTotPts_mg[2]);VRload(resid_2,NTotPts_mg[2],0.0);
    basis_2 = QAllocDouble(NTotPts_mg[2]);VRload(basis_2,NTotPts_mg[2],0.0);
    zvec_2 = QAllocDouble(NTotPts_mg[2]);VRload(zvec_2,NTotPts_mg[2],0.0);
    rstrct_2 = QAllocDouble(NTotPts_mg[2]);VRload(rstrct_2,NTotPts_mg[2],0.0);
    phi_corr2 = QAllocDouble(NTotPts_mg[2]);VRload(phi_corr2,NTotPts_mg[2],0.0);
  }
  if(mg_order > 3){
    Abasis_1 = QAllocDouble(NTotPts_mg[3]);VRload(Abasis_1,NTotPts_mg[3],0.0);
    Aphi_1 = QAllocDouble(NTotPts_mg[3]);VRload(Aphi_1,NTotPts_mg[3],0.0);
    resid_1 = QAllocDouble(NTotPts_mg[3]);VRload(resid_1,NTotPts_mg[3],0.0);
    basis_1 = QAllocDouble(NTotPts_mg[3]);VRload(basis_1,NTotPts_mg[3],0.0);
    zvec_1 = QAllocDouble(NTotPts_mg[3]);VRload(zvec_1,NTotPts_mg[3],0.0);
    rstrct_1 = QAllocDouble(NTotPts_mg[3]);VRload(rstrct_1,NTotPts_mg[3],0.0);
    phi_corr1 = QAllocDouble(NTotPts_mg[3]);VRload(phi_corr1,NTotPts_mg[3],0.0);
  }
}

