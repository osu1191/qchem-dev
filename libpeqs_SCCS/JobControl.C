#include "peqs_class.h"
#include "../libgen/evlbasis.h"

//CS: added this jobtype from ReferenceJob blueprint
void PEQSclass::PBJob()
{
  converge_total = false;
  converge_rho = false;
  converge_rho_ions = false;   

  qtime_t shackle00;
  double shackle01[3];
  shackle00 = QTimerOn();
  ComputeLambda();
  ReferenceRhoPhi();
  ResetPolarIter();
    
  //set up electrolytes: for now NaCl
  Ion anion(-1,ion_radius);
  Ion cation(1,ion_radius);  
  std::vector<Ion> Ions;
  Ions.push_back(anion);
  Ions.push_back(cation);
  el = Electrolyte(Ions,concentration);

  //initialize ion charges on grid
  std::fill_n(rho_ions, NTotPts, 0.0);
  std::fill_n(rho_iter, NTotPts, 0.0);
 
  do{
     do {
       do {
         converge_total = MultiGridControl(total);
         if (converge_total && iter_polar == 0 && peqs_print > 0){
            printf("\tPEqS:  Computing total reference ESP in diel medium");
         }
         else if (converge_total && iter_polar > 0 && peqs_print > 0){
            cout << ".";
         }
       }while(!converge_total);
       converge_rho = UpdateRhoIter();
       if (converge_rho && converge_total){
         if (peqs_print > 0)
            printf("\n\tPEqS:  Converged total reference ESP and polarization charge!\n");
       }
     }while(!converge_rho);
     converge_rho_ions = UpdateRhoIons(el);
     iter_polar =0;
     if (converge_rho && converge_total && converge_rho_ions){
       if (peqs_print > 0)
          printf("\n\tPEqS:  Converged total reference ESP, polarization and electrolyte ion charge!\n");
       SetConvPEQ(true);
     }
  }while(!converge_rho_ions);
  if (doNeqJob == 0 || (doNeqJob == 1 && noneq_partition == 0))
  {
    VRscale(phi_elec, NTotPts,4.0*Pi);
    VRscale(phi_nuke, NTotPts,4.0*Pi);
    VRscale(phi_total,NTotPts,4.0*Pi);
  }
  QTimerOff(shackle01,shackle00);
  if (peqs_print > 0)
     printf("\tPEqS:  Converged the reference potential in %2.2f (cpu) %2.2f (wall) seconds. \n",
        shackle01[0],shackle01[2]);
}

void PEQSclass::ReferenceJob()
{
  converge_total = false;
  converge_rho = false;

  qtime_t shackle00;
  double shackle01[3];
  shackle00 = QTimerOn();
  ReferenceRhoPhi();
  ResetPolarIter();
  do {
    do {
      converge_total = MultiGridControl(total);
      if (converge_total && iter_polar == 0 && peqs_print > 0)
         printf("\tPEqS:  Computing total reference ESP in diel medium");
      else if (converge_total && iter_polar > 0 && peqs_print > 0)
	cout << ".";
    }while(!converge_total);
    converge_rho = UpdateRhoIter();
    if (converge_rho && converge_total){
      if (peqs_print > 0)
         printf("\n\tPEqS:  Converged total reference ESP and polarization charge!\n");
      SetConvPEQ(true);
    }
  }while(!converge_rho);
  //Note: Poisson's equation is solved without the 4.0*Pi factor so incorporate it here.
  //Note: Do this for Pekar partition later on.
  if (doNeqJob == 0 || (doNeqJob == 1 && noneq_partition == 0))
  {
    VRscale(phi_elec, NTotPts,4.0*Pi);
    VRscale(phi_nuke, NTotPts,4.0*Pi);
    VRscale(phi_total,NTotPts,4.0*Pi);
    VRscale(phi_polar,NTotPts,4.0*Pi);
  }
  QTimerOff(shackle01,shackle00);
  if (peqs_print > 0)
     printf("\tPEqS:  Coverged the reference potential in %2.2f (cpu) %2.2f (wall) seconds. \n",
        shackle01[0],shackle01[2]);

//  PrintEpsilon();
//  Print();
}

void PEQSclass::ReferenceRhoPhi()
{
  VRadd(phi_solute,phi_elec,phi_nuke,NTotPts);
  VRadd(rho_solute,rho_elec,rho_nuke,NTotPts);
  VRcopy(phi_total,phi_solute,NTotPts);
  VRcopy(rho_total,rho_solute,NTotPts);
  VRload(rho_iter,NTotPts,0.0);
  VRload(phi_polar,NTotPts,0.0);
}

//CS: add class to update ionic charge density for Poisson-Boltzmann
bool PEQSclass::UpdateRhoIons(Electrolyte el){
  double residual=0.0;
  double rho_ions_unscaled = 0.0;
  double phi_local = 0.0;
  double sum=0.0;
  int ijk;
  //compute concentration from dielectric function and electrostatic potential
  if(linear){
    if(size_modification==false){
       #pragma omp for
       for (ijk=0;ijk<NTotPts;ijk++){
          // need 4*pi factor here
          phi_local = 4.0*Pi*phi_total[ijk];
          rho_ions_unscaled = -2.0*lambda[ijk]*el.get_concentration_au()*phi_local/temperature;
          residual += (rho_ions_unscaled-rho_ions[ijk])*(rho_ions_unscaled-rho_ions[ijk]);
          rho_ions[ijk] = eta_pb*rho_ions_unscaled + (1.0-eta_pb)*rho_ions[ijk];
       }
    }
    else{//with size modification
       double c = el.get_concentration_au()/el.get_c_comb();
       #pragma omp for
       for (ijk=0;ijk<NTotPts;ijk++){
          // need 4*pi factor here
          phi_local = 4.0*Pi*phi_total[ijk];
          rho_ions_unscaled = -2.0*lambda[ijk]*el.get_concentration_au()*phi_local/temperature
                              /(1.0-c*(1-lambda[ijk]));
          residual += (rho_ions_unscaled-rho_ions[ijk])*(rho_ions_unscaled-rho_ions[ijk]);
          rho_ions[ijk] = eta_pb*rho_ions_unscaled + (1.0-eta_pb)*rho_ions[ijk];
       }
    }
    residual = pow(residual,0.5);
    //std::cout << "Residual: " << residual << std::endl;
  }
  else{ //non-linear (full) Poisson Boltzmann
    double rho_ions_unscaled = 0.0;
    double rho_ions_scaled;
    double phi_local = 0.0;
    if(size_modification==false){
       #pragma omp for
       for (ijk=0;ijk<NTotPts;ijk++){
          if(lambda[ijk]>1e-15){
             phi_local = 4.0*Pi*phi_total[ijk];
             rho_ions_unscaled = -2.0*lambda[ijk]*el.get_concentration_au()*sinh(phi_local/temperature);
          }else{
             rho_ions_unscaled = 0;
          }
          rho_ions_scaled = eta_pb*rho_ions_unscaled + (1.0-eta_pb)*rho_ions[ijk];
          residual += (rho_ions_scaled-rho_ions[ijk])*(rho_ions_scaled-rho_ions[ijk]);
          rho_ions[ijk] = rho_ions_scaled;
       }
    }
    else{//non-linear with size modification
       double c = el.get_concentration_au()/el.get_c_comb();
       #pragma omp for
       for (ijk=0;ijk<NTotPts;ijk++){
          if(lambda[ijk]>1e-15){
             phi_local = 4.0*Pi*phi_total[ijk];
             rho_ions_unscaled = -2.0*lambda[ijk]*el.get_concentration_au()*sinh(phi_local/temperature)
                  /(1.0+c*2.0*(lambda[ijk]*cosh(phi_local/temperature)-1.0));
          }else{
             rho_ions_unscaled = 0.0;
          }
          rho_ions_scaled = eta_pb*rho_ions_unscaled + (1.0-eta_pb)*rho_ions[ijk];
          residual += (rho_ions_scaled-rho_ions[ijk])*(rho_ions_scaled-rho_ions[ijk]);
          rho_ions[ijk] = rho_ions_scaled;
       }
    }
    residual = pow(residual,0.5);
    //std::cout << "Residual: " << residual << std::endl;
  }
  if ((rho_ions_thresh-residual) > 0.0) {
    converge_rho_ions = true;
  }
return converge_rho_ions;
}

bool PEQSclass::UpdateRhoIter()
{
  int ijk,l;
  double residual=0.0;
//  cout << "eta = " << eta << endl;

  FinDif1stDeriv(phi_total,delPhi);
  #pragma omp parallel \
  private(ijk,l) \
  reduction(+:residual)
  {
    #pragma omp for
    for (ijk=0;ijk<NTotPts;ijk++)
    {
      double dot = 0.0;
      for (l=0;l<3;l++) 
         dot += delLogEps[3*ijk+l]*delPhi[3*ijk+l];

      // JMH (2/19dd) -- removed 4*pi factors
      //residual += eda*eda*(dot/(4.0*Pi)-rho_iter[ijk])*(dot/(4.0*Pi)-rho_iter[ijk]);
      //rho_iter[ijk] = (eda*dot/(4.0*Pi)) + (1.0-eda)*rho_iter[ijk];
      residual += eta*eta*(dot-rho_iter[ijk])*(dot-rho_iter[ijk]);
      rho_iter[ijk] = eta*dot + (1.0-eta)*rho_iter[ijk];
      if(poisson_boltzmann){
        rho_polar[ijk] = rho_iter[ijk] + ((1.0-epsilon[ijk])/epsilon[ijk])*(rho_solute[ijk]+rho_ions[ijk]);
        rho_total[ijk] = rho_solute[ijk] + rho_ions[ijk] + rho_polar[ijk];
      }
      else{
        rho_polar[ijk] = rho_iter[ijk] + ((1.0-epsilon[ijk])/epsilon[ijk])*rho_solute[ijk];
        rho_total[ijk] = rho_solute[ijk] + rho_polar[ijk];
        phi_polar[ijk] = phi_total[ijk] - phi_solute[ijk];
      }
    }
  }
  iter_polar++;
  residual = pow(residual,0.5);
//  cout << "resid = " << residual << endl;

  double total_pol;
  total_pol = Integrator(rho_polar);
  cout << "Total RhoPolar = " << total_pol << endl;  

  if ((residual - solver_thresh) < 0.0) converge_polar = true;
  if (!converge_polar && iter_polar > max_iter)
    QCrash(" Unable to converge the induced polarization charge density.");
  return converge_polar;
}


void PEQSclass::MarcusRefJob()
{
  VRcopy(phipol_ref_slow,phi_polar,NTotPts);
  VRcopy(rhopol_ref_slow,rho_polar,NTotPts);
  int ijk;

  VRscale(phipol_ref_slow,NTotPts,marcus_slow);
  VRscale(rhopol_ref_slow,NTotPts,marcus_slow);

  VRsub(rhopol_ref_fast,rho_polar,rhopol_ref_slow,NTotPts);

  FileMan(FM_WRITE,FILE_NONEQ_PHI,FM_DP,NTotPts,0,FM_BEG,phipol_ref_slow);
  FileMan(FM_WRITE,FILE_NONEQ_RHO,FM_DP,NTotPts,0,FM_BEG,rhopol_ref_slow);
  FileMan(FM_WRITE,FILE_NONEQ_RHO,FM_DP,NTotPts,NTotPts,FM_BEG,rhopol_ref_fast);
}

void PEQSclass::MarcusNeqJob()
{
  converge_total = false;
  converge_rho = false;
  MarcusNeqRhoPhi();
  ResetPolarIter();
//  cout << "Inside Marcus Neq Job()" << endl;
  do {
    do {
      converge_total = MultiGridControl(total);
      if (converge_total && iter_polar == 0 && peqs_print > 0)
         printf("\tPEqS:  Computing ionized ESP in optical diel medium");
      else if (converge_total && iter_polar > 0 && peqs_print > 0)
	cout << ".";
    }while(!converge_total);
    converge_rho = MarcusNeqRhoIter();
    if (converge_rho && converge_total){
      if (peqs_print > 0)
         printf("\n\tPEqS:  Converged total ionized ESP and polarization charge!\n");
      SetConvPEQ(true);
    }
  }while(!converge_rho);
  VRscale(phi_elec,NTotPts,4.0*Pi);
  VRscale(phi_nuke,NTotPts,4.0*Pi);
  VRscale(phi_total,NTotPts,4.0*Pi);
  VRscale(phi_polar,NTotPts,4.0*Pi);
  VRscale(phipol_ref_slow,NTotPts,4.0*Pi);
}

void PEQSclass::MarcusNeqRhoPhi()
{
  VRadd(phi_solute,phi_elec,phi_nuke,NTotPts);
  VRadd(rho_solute,rho_elec,rho_nuke,NTotPts);

  FileMan(FM_READ,FILE_NONEQ_PHI,FM_DP,NTotPts,0,FM_BEG,phipol_ref_slow);
  FileMan(FM_READ,FILE_NONEQ_RHO,FM_DP,NTotPts,0,FM_BEG,rhopol_ref_slow);
  FileMan(FM_READ,FILE_NONEQ_RHO,FM_DP,NTotPts,NTotPts,FM_BEG,rhopol_ref_fast);

  VRscale(phipol_ref_slow,NTotPts,1.0/(4.0*Pi));

  VRadd(phi_total,phi_solute,phipol_ref_slow,NTotPts);
  VRadd(rho_total,rho_solute,rhopol_ref_slow,NTotPts);

  VRload(rho_iter,NTotPts,0.0);
  VRload(phi_polar,NTotPts,0.0);
}

bool PEQSclass::MarcusNeqRhoIter()
{
  int ijk,l;
  double residual=0.0;
  FinDif1stDeriv(phi_total,delPhi);
//  cout << "Doing Marcus Neq Iter job" << endl;
  #pragma omp parallel \
  private(ijk,l) \
  reduction(+:residual)
  {
    #pragma omp for
    for (ijk=0;ijk<NTotPts;ijk++)
    {
      double dot = 0.0;
      for (l=0;l<3;l++) 
         dot += delLogEps[3*ijk+l]*delPhi[3*ijk+l];

      // JMH (2/19) -- removed 4*pi factors
      //residual += eda*eda*(dot/(4.0*Pi)-rho_iter[ijk])*(dot/(4.0*Pi)-rho_iter[ijk]);
      //rho_iter[ijk] = (eda*dot/(4.0*Pi)) + (1.0-eda)*rho_iter[ijk];
      residual += eta*eta*(dot-rho_iter[ijk])*(dot-rho_iter[ijk]);
      rho_iter[ijk] = eta*dot + (1.0-eta)*rho_iter[ijk];

      rho_polar[ijk] = rho_iter[ijk]
                     + ((1.0-epsilon[ijk])/epsilon[ijk])*(rho_solute[ijk]+rhopol_ref_slow[ijk]);
      rho_total[ijk] = rho_solute[ijk] + rhopol_ref_slow[ijk] + rho_polar[ijk];
      phi_polar[ijk] = phi_total[ijk] - phi_solute[ijk] - phipol_ref_slow[ijk];
    }
  }
  iter_polar++;
  residual = pow(residual,0.5);
  if ((residual - solver_thresh) < 0.0) 
     converge_polar = true;
  if (!converge_polar && iter_polar > max_iter)
    QCrash(" Unable to converge the induced polarization charge density.");
  return converge_polar;
}

void PEQSclass::PekarRefJob()
{
  int ijk;
  converge_total = false;
  converge_rho = false;

  if (!converge_peq) QCrash(" The reference state calculation has failed, should not be here!");
  else ResetConvPEQ();
 
  PekarNeqRhoPhi();
  ResetPolarIter();
  do {
    do {
      converge_total = MultiGridControl(total_fast);
      if (converge_total && iter_polar == 0 && peqs_print > 0)
         printf("\tPEqS:  Computing total reference ESP in optical diel medium");
      else if (converge_total && iter_polar > 0 && peqs_print > 0)
        cout << ".";
    }while(!converge_total);
      converge_rho = PekarNeqRhoIter();
    if (converge_rho && converge_total)
    {
      if (peqs_print > 0)
         printf("\n\tPEqS:  Converged total reference ESP and polarization charge!\n");
      SetConvPEQ(true);
      VRscale(phi_elec, NTotPts,4.0*Pi);
      VRscale(phi_nuke, NTotPts,4.0*Pi);
      VRscale(phi_total,NTotPts,4.0*Pi);
      VRscale(phi_polar,NTotPts,4.0*Pi);
      VRscale(phipol_ref_fast,NTotPts,4.0*Pi);
      VRsub(rhopol_ref_slow,rho_polar,rhopol_ref_fast,NTotPts);
      VRsub(phipol_ref_slow,phi_polar,phipol_ref_fast,NTotPts);
      double *qS1ESPT1 = QAllocDouble(NTotPts);
      #pragma omp parallel for private(ijk) 
      for(ijk=0;ijk<NTotPts;ijk++) 
        qS1ESPT1[ijk] = (phi_elec[ijk]+phi_nuke[ijk])*rhopol_ref_slow[ijk];
      double qS1_ESPT1 = Integrator(qS1ESPT1);
      FileMan(FM_WRITE,FILE_REF_SOLV_ENERGY,FM_DP,1,0,FM_BEG,&qS1_ESPT1);
      FileMan(FM_WRITE,FILE_NONEQ_RHO,FM_DP,NTotPts,0,FM_BEG,rhopol_ref_slow);       
      FileMan(FM_WRITE,FILE_NONEQ_PHI,FM_DP,NTotPts,0,FM_BEG,phipol_ref_slow);
      QFree(qS1ESPT1);
    }
  }while(!converge_rho);
}

void PEQSclass::PekarNeqJob()
{
  converge_total = false;
  converge_rho = false;
  PekarNeqRhoPhi();
  ResetPolarIter();
  do {
    do {
      converge_total = MultiGridControl(total);
      if (converge_total && iter_polar == 0 && peqs_print > 0)
         printf("\tPEqS:  Computing total ionized ESP in optical diel medium");
      else if (converge_total && iter_polar > 0 && peqs_print > 0)
        cout << ".";
    }while(!converge_total);
      converge_rho = PekarNeqRhoIter();
    if (converge_rho && converge_total){
      if (peqs_print > 0)
         printf("\n\tPEqS:  Converged total ionized ESP and polarization charge!\n");
      SetConvPEQ(true);
    }
  }while(!converge_rho);
  VRscale(phi_elec, NTotPts,4.0*Pi);
  VRscale(phi_nuke, NTotPts,4.0*Pi);
  VRscale(phi_total,NTotPts,4.0*Pi);
  VRscale(phi_polar,NTotPts,4.0*Pi);
}

void PEQSclass::PekarNeqRhoPhi()
{
  if (noneq_state == 0)
  {
    ComputeDielectric(true);
    VRadd(phitot_fast,phi_elec,phi_nuke,NTotPts);
    VRadd(rhotot_fast,rho_elec,rho_nuke,NTotPts);
    VRload(rho_iter,NTotPts,0.0);
    VRload(phipol_ref_fast,NTotPts,0.0);
  }
  else{
    FileMan(FM_READ,FILE_NONEQ_RHO,FM_DP,NTotPts,0,FM_BEG,rhopol_ref_slow);
    FileMan(FM_READ,FILE_NONEQ_PHI,FM_DP,NTotPts,0,FM_BEG,phipol_ref_slow);
    VRadd(phi_solute,phi_elec,phi_nuke,NTotPts);
    VRadd(rho_solute,rho_elec,rho_nuke,NTotPts);
    VRcopy(phi_total,phi_solute,NTotPts);
    VRcopy(rho_total,rho_solute,NTotPts);
    VRload(rho_iter,NTotPts,0.0);
    VRload(phi_polar,NTotPts,0.0);
  }
}

bool PEQSclass::PekarNeqRhoIter()
{
  int ijk,l;
  double residual=0.0;
  if (noneq_state == 0){
    FinDif1stDeriv(phitot_fast,delPhiFast);
    #pragma omp parallel \
    private(ijk,l) \
    reduction(+:residual)
    {
      #pragma omp for
      for (ijk=0;ijk<NTotPts;ijk++)
      {
        double dot = 0.0;
        for (l=0;l<3;l++) 
           dot += delLogEps[3*ijk+l]*delPhiFast[3*ijk+l];

        // JMH (2/19) -- removed 4*pi factors
        //residual += eda*eda*(dot/(4.0*Pi)-rho_iter[ijk])*(dot/(4.0*Pi)-rho_iter[ijk]);
        //rho_iter[ijk] = (eda*dot/(4.0*Pi)) + (1.0-eda)*rho_iter[ijk];
        residual += eta*eta*(dot-rho_iter[ijk])*(dot-rho_iter[ijk]);
        rho_iter[ijk] = eta*dot + (1.0-eta)*rho_iter[ijk];


        rhopol_ref_fast[ijk] = rho_iter[ijk] + ((1.0-epsilon[ijk])/epsilon[ijk])*rho_solute[ijk];
        rhotot_fast[ijk] = rho_solute[ijk] + rhopol_ref_fast[ijk];
        phipol_ref_fast[ijk] = phitot_fast[ijk] - phi_solute[ijk];
      }
    }
    iter_polar++;
    residual = pow(residual,0.5);
    if ((residual - solver_thresh) < 0.0) 
       converge_polar = true;
    if (!converge_polar && iter_polar > max_iter)
      QCrash(" Unable to converge the induced polarization charge density.");
  }
  else{
    FinDif1stDeriv(phi_total,delPhi);
    #pragma omp parallel \
    private(ijk,l) \
    reduction(+:residual)
    {
      #pragma omp for
      for (ijk=0;ijk<NTotPts;ijk++)
      {
        double dot = 0.0;
        for (l=0;l<3;l++) 
           dot += delLogEps[3*ijk+l]*delPhi[3*ijk+l];

        // JMH (2/19) -- removed 4*pi factors
        //residual += eda*eda*(dot/(4.0*Pi)-rho_iter[ijk])*(dot/(4.0*Pi)-rho_iter[ijk]);
        //rho_iter[ijk] = (eda*dot/(4.0*Pi)) + (1.0-eda)*rho_iter[ijk];
        residual += eta*eta*(dot-rho_iter[ijk])*(dot-rho_iter[ijk]);
        rho_iter[ijk] = eta*dot + (1.0-eta)*rho_iter[ijk];


        rho_polar[ijk] = rho_iter[ijk] + ((1.0-epsilon[ijk])/epsilon[ijk])*rho_solute[ijk];
        rho_total[ijk] = rho_solute[ijk] + rho_polar[ijk];
        phi_polar[ijk] = phi_total[ijk] - phi_solute[ijk];
      }
    }
    iter_polar++;
    residual = pow(residual,0.5);
    if ((residual - solver_thresh) < 0.0) 
       converge_polar = true;
    if (!converge_polar && iter_polar > max_iter)
      QCrash(" Unable to converge the induced polarization charge density.");
  }
  return converge_polar;
}
