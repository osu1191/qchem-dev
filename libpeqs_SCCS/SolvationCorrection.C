#include "peqs_class.h"
#include <math.h>
#include <armadillo>

double PEQSclass::SolvEnergyPB(Electrolyte el){
  int ijk;
  double G_total=0.0;
  if (!converge_peq) QCrash(" Must converge the electrostatic potentials first!");

    //calculate total charge of mobile ions
    double Qions=0.0,Qpol=0.0;
    Qions = Integrator(rho_ions);
    Qpol = Integrator(rho_polar);

    if (peqs_print > 0){
       printf("\tPEqS:  Total ionic   charge on grid = %+8.6f\n",-Qions);
       printf("\tPEqS:  Total polar   charge on grid = %+8.6f\n",-Qpol);
       printf("\tPEqS:  Sum of       charges on grid = %+8.6f\n",-Qions-Qpol);
    }
    double *qESPN=QAllocDouble(NTotPts),*qESPE=QAllocDouble(NTotPts),*qESPI=QAllocDouble(NTotPts),*qESPEOP=QAllocDouble(NTotPts);
    double q_ESPN=0.0, q_ESPE=0.0, q_ESPI=0.0,q_ESPEOP=0.0;

    double c = el.get_concentration_au()/el.get_c_comb(); 

    #pragma omp parallel for private(ijk)
    for (ijk=0;ijk<NTotPts;ijk++){
      qESPN[ijk] = phi_nuke[ijk]*(rho_polar[ijk]+rho_ions[ijk]);
      qESPE[ijk] = phi_elec[ijk]*(rho_polar[ijk]+rho_ions[ijk]);
      qESPI[ijk] = rho_ions[ijk]*phi_total[ijk];

      if(linear){
         if(size_modification==false){
            if(lambda[ijk]>1e-15){
               qESPEOP[ijk] = -lambda[ijk]*el.get_concentration_au()*std::pow(phi_total[ijk], 2.0)/temperature;
            }else{
               qESPEOP[ijk] = 0.0;
            }
         }else{//with size modification
            if(lambda[ijk]>1e-15){
               qESPEOP[ijk] = -lambda[ijk]*el.get_concentration_au()*phi_total[ijk]/temperature
                               /(1.0-c*(1-lambda[ijk]));
            }else{
               qESPEOP[ijk] = 0.0;
            }
         }  
      }else{ //nonlinear

         if(size_modification==false){
              if(lambda[ijk]>1e-15){
                 qESPEOP[ijk] = 2.0*el.get_concentration_au()*temperature*(1.0-lambda[ijk]*cosh(phi_total[ijk]/temperature));
              }else{
                 qESPEOP[ijk] = 2.0*el.get_concentration_au()*temperature;
              }
         }else{//with size modification
              if(lambda[ijk]>1e-15){
                qESPEOP[ijk] = -el.get_c_comb()*temperature*log(1.0+c*2.0*(lambda[ijk]*cosh(phi_total[ijk]/temperature)-1.0));
              }else{
                qESPEOP[ijk] = -el.get_c_comb()*temperature*log(1.0-c*2.0);
              }
         }      
      }
    }

    q_ESPN = Integrator(qESPN);
    q_ESPE = Integrator(qESPE);
    q_ESPI = Integrator(qESPI);
    q_ESPEOP = Integrator(qESPEOP);
 
    if(linear){
       G_total = 0.5*q_ESPN+0.5*q_ESPE;
    }
    else{
       G_total = 0.5*q_ESPN+0.5*q_ESPE-0.5*q_ESPI+q_ESPEOP;
    }
    printf("\tReference G(solv) = %14.8f a.u. = %8.5f kcal/mol\n",G_total,G_total*au2kcal);
    if (peqs_print > 0 ){
      printf("\tPEqS:  q_ESPN      = %14.8f a.u. = %8.5f kcal/mol\n",q_ESPN,q_ESPN*au2kcal);
      printf("\tPEqS:  q_ESPE      = %14.8f a.u. = %8.5f kcal/mol\n",q_ESPE,q_ESPE*au2kcal);
      printf("\tPEqS:  q_ESPI      = %14.8f a.u. = %8.5f kcal/mol\n",q_ESPI,q_ESPI*au2kcal);
      printf("\tPEqS:  q_ESPEOP    = %14.8f a.u. = %8.5f kcal/mol\n",q_ESPEOP,q_ESPEOP*au2kcal);
    }

    //CS:fix some stuff for use with 1e Integrals in STVman
    //CS: taken over from original Poisson code below
    if(linear==true){
       G_total = 0.5*q_ESPN - 0.5*q_ESPE; //in STVman, a call to MakeV adds q_ESPE to E1, need to account for that here.
    }else{
       G_total = 0.5*q_ESPN - 0.5*q_ESPE-0.5*q_ESPI+q_ESPEOP;
    }
    QFree(qESPN);QFree(qESPE);
    rem_write(1,REM_PEQ_CORR);

    FileMan(FM_WRITE,FILE_PEQS_DATA,FM_INT,1,0,FM_BEG,&NTotPts);
    //add rho_ions to rho_polar just for the 1e Integral part
    VRadd(rho_polar,rho_polar,rho_ions,NTotPts);
    VRscale(rho_polar,NTotPts,-1.0*volumeElmnt()); //turn densities into charges; account for -1 charge of e- in stvman
    FileMan(FM_WRITE,FILE_PEQS_DATA,FM_DP,NTotPts,0,FM_CUR,rho_polar);
    FileMan(FM_WRITE,FILE_PEQS_DATA,FM_DP,3*NTotPts,0,FM_CUR,XYZ);

    //free matrices and write stuff to file for use in next iteration
    QFree(qESPN);QFree(qESPE);QFree(qESPI),QFree(qESPEOP);

  PBcallNo++;
  rem_write(PBcallNo,REM_PEQ_CALLS);

  if (peqs_print > 1) {
    Print();
  }

  return G_total;
}


double PEQSclass::SolvEnergy(double EMax)
{
  int ijk;
  double G_total=0.0;
  if (!converge_peq) QCrash(" Must converge the electrostatic potentials first!");
  if (!doNeqJob || noneq_state == 0){
    double *qESPN=QAllocDouble(NTotPts),*qESPE=QAllocDouble(NTotPts);
    double q_ESPN=0.0,q_ESPE=0.0;
    #pragma omp parallel for private(ijk)
    for (ijk=0;ijk<NTotPts;ijk++){
      qESPN[ijk] = phi_nuke[ijk]*rho_polar[ijk];
      qESPE[ijk] = phi_elec[ijk]*rho_polar[ijk];
    }
    q_ESPN = Integrator(qESPN);
    q_ESPE = Integrator(qESPE);
    G_total = q_ESPN + q_ESPE;
    G_total *= 0.5; //this is the actual reference solvation free energy
    printf("\tPEqS:  Reference G(solv) = %14.8f a.u. = %5.2f kcal/mol\n",G_total,G_total*au2kcal);
    if (peqs_print > 0 ){
      printf("\tPEqS:  q_ESPN    = %8.4f\n",q_ESPN);
      printf("\tPEqS:  q_ESPE    = %8.4f\n",q_ESPE);
    }
    //REF_SOLV_ENERGY was written elsewhere for pekar partition.
    if (noneq_partition == 0) FileMan(FM_WRITE,FILE_REF_SOLV_ENERGY,FM_DP,1,0,FM_BEG,&G_total);
    G_total = q_ESPN - q_ESPE; //in STVman, a call to MakeV adds q_ESPE to E1, need to account for that here.
    G_total *= 0.5; //return this value to SCFman
    QFree(qESPN);QFree(qESPE);
    rem_write(1,REM_PEQ_CORR);
    FileMan(FM_WRITE,FILE_PEQS_DATA,FM_INT,1,0,FM_BEG,&NTotPts);
    FileMan(FM_WRITE,FILE_PEQS_DATA,FM_DP,1,0,FM_CUR,&EMax); // Solvent aware
    VRscale(rho_polar,NTotPts,-1.0*volumeElmnt()); //turn densities into charges; account for -1 charge of e- in stvman
    FileMan(FM_WRITE,FILE_PEQS_DATA,FM_DP,NTotPts,0,FM_CUR,rho_polar);
    FileMan(FM_WRITE,FILE_PEQS_DATA,FM_DP,3*NTotPts,0,FM_CUR,XYZ);

    if (CavityType != PEQS_CAVTYPE_GYGI){
	Veps=QAllocDouble(NBasis*NBasis);	
	VRload(Veps,NBasis*NBasis,0.0);
    }
    FileMan(FM_WRITE,FILE_PEQS_DATA,FM_DP,NBasis*NBasis,0,FM_CUR,Veps); // Solvent aware
  }
  else if (doNeqJob && noneq_state == 1){
    double V_ref=0.0;
    double *qF2ESPN2=QAllocDouble(NTotPts);
    double *qF2ESPE2=QAllocDouble(NTotPts);
    double *qF2ESPS1=QAllocDouble(NTotPts);
    double *qS1ESPN2=QAllocDouble(NTotPts);
    double *qS1ESPE2=QAllocDouble(NTotPts);
    double *qF1ESPS1=QAllocDouble(NTotPts);
    for (ijk=0;ijk<NTotPts;ijk++){
      qF2ESPN2[ijk] = phi_nuke[ijk]*rho_polar[ijk];
      qF2ESPE2[ijk] = phi_elec[ijk]*rho_polar[ijk];
      qS1ESPN2[ijk] = phi_nuke[ijk]*rhopol_ref_slow[ijk];
      qS1ESPE2[ijk] = phi_elec[ijk]*rhopol_ref_slow[ijk];
      if (noneq_partition == 0){
        qF2ESPS1[ijk] = rho_polar[ijk]*phipol_ref_slow[ijk];
        qF1ESPS1[ijk] = rhopol_ref_fast[ijk]*phipol_ref_slow[ijk];
      }
    }
    FileMan(FM_READ,FILE_REF_SOLV_ENERGY,FM_DP,1,0,FM_BEG,&V_ref);
    double qF2_ESPN2 = Integrator(qF2ESPN2);
    double qF2_ESPE2 = Integrator(qF2ESPE2);
    double qF2_ESPS1 = Integrator(qF2ESPS1);
    double qS1_ESPN2 = Integrator(qS1ESPN2);
    double qS1_ESPE2 = Integrator(qS1ESPE2);
    double qS1_ESPT1;
    if (noneq_partition == 0) qS1_ESPT1 = 2.0*V_ref;
    else if (noneq_partition == 1) qS1_ESPT1 = V_ref; 
    double qF1_ESPS1 = Integrator(qF1ESPS1);

    if (noneq_partition == 1) {qF2_ESPS1 = 0.0; qF1_ESPS1 = 0.0;}

    if (noneq_partition == 0)
      G_total = qF2_ESPN2 + qF2_ESPE2 + 2.0*qS1_ESPN2 + 2.0*qS1_ESPE2 - qS1_ESPT1 + qF2_ESPS1 - qF1_ESPS1;
    else if (noneq_partition == 1)
      G_total = qF2_ESPN2 + qF2_ESPE2 + 2.0*qS1_ESPN2 + 2.0*qS1_ESPE2 - qS1_ESPT1;

    G_total *= 0.5; //this is the actual solvation free energy
    printf("\tPEqS:  Ionized   G(solv) = %14.8f a.u. = %5.2f kcal/mol\n",G_total,G_total*au2kcal);

    if (peqs_print > 0){
      printf("\tPEqS:  qF2_ESPN2 = %8.4f\n",qF2_ESPN2);
      printf("\tPEqS:  qF2_ESPE2 = %8.4f\n",qF2_ESPE2);
      printf("\tPEqS:  qF2_ESPS1 = %8.4f\n",qF2_ESPS1);
      printf("\tPEqS:  qS1_ESPN2 = %8.4f\n",qS1_ESPN2);
      printf("\tPEqS:  qS1_ESPN2 = %8.4f\n",qS1_ESPN2);
      printf("\tPEqS:  qS1_ESPE2 = %8.4f\n",qS1_ESPE2);
      printf("\tPEqS:  qS1_ESPT1 = %8.4f\n",qS1_ESPT1);
      printf("\tPEqS:  qF1_ESPS1 = %8.4f\n",qF1_ESPS1);
    }

    //in STVman, MakeV adds qF2_ESPE2+qS1_ESPE2 to E1, need to account for that here.
    if (noneq_partition == 0)
      G_total = qF2_ESPN2 - qF2_ESPE2 + 2.0*qS1_ESPN2 - qS1_ESPT1 + qF2_ESPS1 - qF1_ESPS1;
    else if (noneq_partition == 1)
      G_total = qF2_ESPN2 - qF2_ESPE2 + 2.0*qS1_ESPN2 - qS1_ESPT1;
    G_total *= 0.5;

    rem_write(1,REM_PEQ_CORR);
    FileMan(FM_WRITE,FILE_PEQS_DATA,FM_INT,1,0,FM_BEG,&NTotPts);
    FileMan(FM_WRITE,FILE_PEQS_DATA,FM_DP,1,0,FM_CUR,&EMax); // Solvent aware
    //for non-eq state, total polarization response arises from rho^{fast}_{pol,ion} + rho^{slow}_{pol,ground}
    VRadd(rho_polar,rho_polar,rhopol_ref_slow,NTotPts);
    VRscale(rho_polar,NTotPts,-1.0*volumeElmnt());//turn densities into charges; account for -1 charge of e- in stvman
    FileMan(FM_WRITE,FILE_PEQS_DATA,FM_DP,NTotPts,0,FM_CUR,rho_polar);
    FileMan(FM_WRITE,FILE_PEQS_DATA,FM_DP,3*NTotPts,0,FM_CUR,XYZ);
    if (CavityType != PEQS_CAVTYPE_GYGI){
        Veps=QAllocDouble(NBasis*NBasis);
        VRload(Veps,NBasis*NBasis,0.0);
    }
    FileMan(FM_WRITE,FILE_PEQS_DATA,FM_DP,NBasis*NBasis,0,FM_CUR,Veps); // Solvent aware
    QFree(qF2ESPN2);QFree(qF2ESPE2);QFree(qF2ESPS1);
    QFree(qS1ESPN2);QFree(qS1ESPE2);QFree(qF1ESPS1);
  }
  PBcallNo++;
  rem_write(PBcallNo,REM_PEQ_CALLS);
  if (peqs_print > 1) {
    Print();
  }

  return G_total;
}
