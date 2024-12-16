#include "peqs_class.h"

PEQSclass::PEQSclass()
{
  int i,j,k;
  NumDim = 3;
  NTotPts_mg = NULL;
  NPoints = NULL;
  XMin = NULL;
  XMax = NULL;
  Center = NULL;
  delX = NULL;
  BCX = NULL;
  BCY = NULL;
  BCZ = NULL;
  vdWScale = NULL;
  epsilon = NULL;
  log_epsilon = NULL;
  delLogEps = NULL;
  delPhi = NULL;
  delPhiE = NULL;
  delPhiN = NULL;
  delPhiFast = NULL;
  delRhoElec = NULL;
  delCart = NULL;
  rho_total = NULL;
  rho_elec = NULL;
  rho_nuke = NULL;
  rho_solute = NULL;
  rho_iter = NULL;
  rho_polar = NULL;
  phi_total = NULL;
  phi_elec = NULL;
  phi_nuke = NULL;
  phi_solute = NULL;
  phi_polar = NULL;
  rhotot_fast = NULL;
  rhopol_ref_slow = NULL;
  rhopol_ref_fast = NULL;
  phitot_fast = NULL;
  phipol_ref_slow = NULL;
  phipol_ref_fast = NULL;
  Aphi = NULL;
  Abasis = NULL;
  resid = NULL;
  basis = NULL;
  zvec = NULL;
  BCX_mg = NULL;
  BCY_mg = NULL;
  BCZ_mg = NULL;
  Abasis_3 = NULL;
  Aphi_3 = NULL;
  resid_3 = NULL;
  basis_3 = NULL;
  zvec_3 = NULL;
  rstrct_3 = NULL;
  phi_corr3 = NULL;
  Abasis_2 = NULL;
  Aphi_2 = NULL;
  resid_2 = NULL;
  basis_2 = NULL;
  zvec_2 = NULL;
  rstrct_2 = NULL;
  phi_corr2 = NULL;
  Abasis_1 = NULL;
  Aphi_1 = NULL;
  resid_1 = NULL;
  basis_1 = NULL;
  zvec_1 = NULL;
  rstrct_1 = NULL;
  phi_corr1 = NULL;
  XYZ = NULL;

//====== Interpolation==========//

  mono = NULL;
  dip = NULL;
  quad = NULL;
  octa = NULL;
  hexa = NULL;
  fifth = NULL;
  sixth = NULL;
  seventh = NULL;
  eight = NULL;
  RhoCart = NULL;
  val_elec = NULL;
  valCart = NULL;
  valEField = NULL;
  valACG = NULL;
  ACGField = NULL;
  phiACG = NULL;
  phi_irr = NULL;
  Field = NULL;
  val3f = NULL;
  valSf = NULL;
  delVal = NULL;
  rho_from_field = NULL;
  val2 = NULL;
  val3 = NULL;
  xacg = NULL;
  yacg = NULL;
  zacg = NULL;
  occ = NULL;
  gauss = NULL;
  InSide = NULL;
  OutSide = NULL;
  MidZone = NULL;

  efx = NULL;
  efy = NULL;
  efz = NULL;
  valGf = NULL;

  power = 2.0; 
  shift1 = 1.0;
  shift2 = 2.0;
  shift3 = 3.0;
  filter = 0.2;
  inhouse = 1; 
  outhouse = 3;
  MulOrder = 2;
  potgrid = 1; 
  count = 0;   
  nGrid = 0;   
  GrdPat = 1;  
  interp = 1;  
  cutoff = 10; 

  interp = 2;
  max_der = 0;
  order = 5;
  orderp1 = order + 1;
  halfordp1 = orderp1/2;
  orderp1t2 = 2*orderp1;
  norder = orderp1 - max_der;
  nhalf = norder/2;
  nth = norder - 1;

//==============================//

  fdif_order = 8;
  peqs_print = 0;

  ang2bohr = ConvFac(ANGSTROMS_TO_BOHRS);
  au2kcal = ConvFac(HARTREES_TO_KCAL_MOL);
  Pi = ConvFac(VALUE_OF_PI); 

  PBcallNo = rem_read(REM_PEQS_CALLS);

  CavityType = PEQS_CAVTYPE_VDW;
  GaussBlur = false; 
  sigma = 0.45;
  Interface = false;
  vdwType = PEQS_VDW_UNSCALED;
  vdwScale = 1.0;
  vdwShift = 0.0;
  erf_delta = 0.50;
  
  r_cavity = INVALID_REAL;
  FindCenter = 3;
  AutoRadius = false;
  r_hybrid = 0.0;
  tanh_cut1 = 0.0; 
  tanh_cut2 = 1.0*ang2bohr;
  tanh_scale = 0.0;
  r_interface = 2.75*ang2bohr;
  interface_scale = 0.0;
  gibbsDS = 0.0;
  z_direction = 0;

  multigrid = true; 
  mg_order = 4;
  mg_cycle = 2;
  sweep1 = 2;
  sweep2 = 3;
  sweep3 = sweep1 + sweep2;
  max_iter = 500;
  solver_thresh = 0.0001;
  cg_thresh = 0.00001;
  eta = 0.6;
  MaxPPer = 0;

  doNeqJob = false;
  noneq_state = INVALID_INT;
  noneq_partition = INVALID_INT; 

  electron = 1;
  nuclear = 2;
  total = 3;
  total_fast = 4;

  converge_total = false;
  converge_rho = false;
  converge_peq = false;
  converge_cg = false;
  converge_polar = false;

  eps_vacuum = 1.00; 
  eps_solvent = 78.39; 
  eps_infi = 1.7778;
  dielectric = eps_vacuum;
  opt_dielec = eps_vacuum;

 //CS: added some variables for Poisson-Boltzmann code
  poisson_boltzmann = false;
  rho_ions_thresh = 1.0e-4;
  conv_temp_au = 3.1668154225e-6;
  temperature = 298.00*conv_temp_au;
  concentration = 0.0;
  eta_pb = 0.2;
  stern_thickness = 0.0;
  ion_radius = 0.0;
  probe_radius = 0.0;
  converge_rho_ions = false;
  converge_phi_ions = false;
  converge_phi_polar = false;
  linear = true;
  size_modification = false;
  lambda = NULL;
  rho_ions = NULL;
  phi_ions = NULL;
  phi_polar = NULL;
 
  Read_PEQS();
  Read_PEQSgrid();
  Read_EpsilonGrid();

  if (CavityType == PEQS_CAVTYPE_SPHERE)
  {
    if (tanh_cut2 != 0.0)
       tanh_scale = 4.0/tanh_cut2;
    else
       QCrash(" Need to specify the spherical cavity interpolation length!");
    if (AutoRadius == false && r_cavity == INVALID_REAL)
      QCrash(" Need to specify the spherical cavity radius!");
  }
  if (Interface == true){
    if (r_interface != 0.0)
       interface_scale = 4.0/r_interface;
    else
       QCrash(" Need to specify the interface interpolation length!");
    if (z_direction == 0 )
       QCrash (" Need to specify the z-direction for interface dielectric interpolation");
  }

  if (!multigrid) mg_order = 1;
  mg_level = mg_order;

  if (doNeqJob && noneq_partition == INVALID_INT)
    noneq_partition = 0;
  if (doNeqJob && noneq_partition == 0){
    marcus_slow = (eps_solvent-eps_infi)/(eps_solvent-1.0);
    marcus_fast = 1.0-marcus_slow;
    if (peqs_print > 2)
    {
      printf("\tPEqS:  Marcus partition slow scale factor = %.4f\n",marcus_slow);
      printf("\t                        fast scale factor = %.4f\n",marcus_fast);
    }
  }
  else{
    marcus_slow = 1.0;
    marcus_fast = 0.0;
  }

  Allocate();
  Coordinates();
}

PEQSclass::~PEQSclass()
{
//========Interpolation================//

  if (mono != NULL) QFree(mono); //ACG
  if (dip != NULL) QFree(dip);  //ACG
  if (quad != NULL) QFree(quad);  //ACG
  if (octa != NULL) QFree(octa);  //ACG
  if (hexa != NULL) QFree(hexa);  //ACG
  if (fifth != NULL) QFree(fifth);  //ACG
  if (sixth != NULL) QFree(sixth);  //ACG
  if (seventh != NULL) QFree(seventh);  //ACG
  if (eight != NULL) QFree(eight);  //ACG
  if (RhoCart != NULL) QFree(RhoCart);  //ACG
  if (phi_irr != NULL) QFree(phi_irr);  //ACG
  if (phiACG != NULL) QFree(phiACG);  //ACG
  if (valCart != NULL) QFree(valCart);  //ACG
  if (val_elec != NULL) QFree(val_elec);  //ACG
  if (valEField != NULL) QFree(valEField);  //ACG
  if (valACG != NULL) QFree(valACG);  //ACG
  if (ACGField != NULL) QFree(ACGField);  //ACG
  if (Field != NULL) QFree(Field);  // EField
  if (val3f != NULL) QFree(val3f);  // EField
  if (valSf != NULL) QFree(valSf);  // EField
  if (delVal != NULL) QFree(delVal);  // EField
  if (rho_from_field != NULL) QFree(rho_from_field);  // EField
  if (val2 != NULL) QFree(val2);  //ACG
  if (val3 != NULL) QFree(val3);  //ACG
  if (occ != NULL) QFree(occ);  //ACG
  if (gauss != NULL) QFree(gauss);  //ACG
  if (xacg != NULL) QFree(xacg);  //ACG
  if (yacg != NULL) QFree(yacg);  //ACG
  if (zacg != NULL) QFree(zacg);  //ACG  
  if (InSide != NULL) QFree(InSide);  //ACG
  if (OutSide != NULL) QFree(OutSide);  //ACG
  if (MidZone != NULL) QFree(MidZone);  //Do_IDW3
 
  if (efx != NULL) QFree(efx);
  if (efy != NULL) QFree(efy);
  if (efz != NULL) QFree(efz);
  if (valGf != NULL) QFree(valGf);
//=====================================//


   if (NTotPts_mg != NULL) QFree(NTotPts_mg);
   if (NPoints != NULL) QFree(NPoints);
   if (XMin != NULL) QFree(XMin);
   if (XMax != NULL) QFree(XMax);
   if (Center != NULL) QFree(Center);
   if (delX != NULL) QFree(delX);
   if (BCX != NULL) QFree(BCX);
   if (BCY != NULL) QFree(BCY);
   if (BCZ != NULL) QFree(BCZ);
   if (vdWScale != NULL) QFree(vdWScale);
   if (epsilon != NULL) QFree(epsilon);
   if (log_epsilon != NULL) QFree(log_epsilon);
   if (delLogEps != NULL) QFree(delLogEps);
   if (delPhi != NULL) QFree(delPhi);
   if (delPhiFast != NULL) QFree(delPhiFast);
   if (delPhiE != NULL) QFree(delPhiE);
   if (delPhiN != NULL) QFree(delPhiN);
   if (delRhoElec != NULL) QFree(delRhoElec);
   if (delCart != NULL) QFree(delCart);
   if (rho_total != NULL) QFree(rho_total);
   if (rho_elec != NULL) QFree(rho_elec);
   if (rho_nuke != NULL) QFree(rho_nuke);
   if (rho_solute != NULL) QFree(rho_solute);
   if (rho_iter != NULL) QFree(rho_iter);
   if (rho_polar != NULL) QFree(rho_polar);
   if (phi_total != NULL) QFree(phi_total);
   if (phi_elec != NULL) QFree(phi_elec);
   if (phi_nuke != NULL) QFree(phi_nuke);
   if (phi_solute != NULL) QFree(phi_solute);
   if (phi_polar != NULL) QFree(phi_polar);
   if (rhotot_fast != NULL) QFree(rhotot_fast);
   if (rhopol_ref_slow != NULL) QFree(rhopol_ref_slow);
   if (rhopol_ref_fast != NULL) QFree(rhopol_ref_fast);
   if (phitot_fast != NULL) QFree(phitot_fast);
   if (phipol_ref_slow != NULL) QFree(phipol_ref_slow);
   if (phipol_ref_fast != NULL) QFree(phipol_ref_fast);
   if (Aphi != NULL) QFree(Aphi);
   if (Abasis != NULL) QFree(Abasis);
   if (resid != NULL) QFree(resid);
   if (basis != NULL) QFree(basis);
   if (zvec != NULL) QFree(zvec);
   if (BCX_mg != NULL) QFree(BCX_mg);
   if (BCY_mg != NULL) QFree(BCY_mg);
   if (BCZ_mg != NULL) QFree(BCZ_mg);
   if (Abasis_3 != NULL) QFree(Abasis_3);
   if (Aphi_3 != NULL) QFree(Aphi_3);
   if (resid_3 != NULL) QFree(resid_3);
   if (basis_3 != NULL) QFree(basis_3);
   if (zvec_3 != NULL) QFree(zvec_3);
   if (rstrct_3 != NULL) QFree(rstrct_3);
   if (phi_corr3 != NULL) QFree(phi_corr3);
   if (Abasis_2 != NULL) QFree(Abasis_2);
   if (Aphi_2 != NULL) QFree(Aphi_2);
   if (resid_2 != NULL) QFree(resid_2);
   if (basis_2 != NULL) QFree(basis_2);
   if (zvec_2 != NULL) QFree(zvec_2);
   if (rstrct_2 != NULL) QFree(rstrct_2);
   if (phi_corr2 != NULL) QFree(phi_corr2);
   if (Abasis_1 != NULL) QFree(Abasis_1);
   if (Aphi_1 != NULL) QFree(Aphi_1);
   if (resid_1 != NULL) QFree(resid_1);
   if (basis_1 != NULL) QFree(basis_1);
   if (zvec_1 != NULL) QFree(zvec_1);
   if (rstrct_1 != NULL) QFree(rstrct_1);
   if (phi_corr1 != NULL) QFree(phi_corr1);
   if (XYZ != NULL){QFree(XYZ);}
   //CS: free Poisson-Boltzmann variables
   if (lambda != NULL) QFree(lambda);
   if (rho_ions != NULL) QFree(rho_ions);
   if (phi_ions != NULL) QFree(phi_ions);
   if (phi_polar != NULL) QFree(phi_polar);
}
