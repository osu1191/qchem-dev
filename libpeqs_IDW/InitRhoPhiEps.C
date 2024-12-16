#include "peqs_class.h"

void PEQSclass::InitRhoPhiEps(double *PAv)
{
  double Qnuke=0.0,Qelec=0.0;
  double Velec=0.0;

  if(potgrid == 2){
     Make_Atom_Centered_Grid();
     Get_Electronic_Density_on_ACG(PAv);
     if(interp == 1) Expanded_IDW();
     if(interp == 1) Combine_IDW_and_MulPot(PAv);
  }

  VRload(rho_elec,NTotPts,0.0); // Gygi
  VRload(phi_elec,NTotPts,0.0); // Gygi

  ElectronicRho(PAv); // for Gygi

  for(int ijk=0; ijk<NTotPts; ijk++){
    cout << phi_elec[ijk] << endl;
  }

  ComputeDielectric(false);
  Qelec = Integrator(rho_elec);
  if (peqs_print > 0){
     printf("\tPEqS:  Total elec   charge on grid = %+8.3f\n",-Qelec);
  }

  VRload(rho_nuke,NTotPts,0.0);
  VRload(phi_nuke,NTotPts,0.0);
  NuclearRhoPhi();
  Qnuke = Integrator(rho_nuke);
  if (peqs_print > 0)
  {
    printf("\tPEqS:  Total nucl   charge on grid = %+8.3f\n",-Qnuke);
    printf("\tPEqS:  Total solute charge on grid = %+8.3f\n",-(Qelec+Qnuke));
  }

  VRload(rho_iter,NTotPts,0.0);
  VRload(rho_total,NTotPts,0.0);
  VRload(rho_polar,NTotPts,0.0);
  VRload(rho_solute,NTotPts,0.0);
  VRload(phi_solute,NTotPts,0.0);
  VRload(phi_total,NTotPts,0.0);
  VRload(phi_polar,NTotPts,0.0);
}
