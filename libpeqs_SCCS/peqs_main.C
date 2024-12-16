#include "peqs_class.h"

/******************************************************************
 *                                                                *
 *  PEQS solves the Poisson equation on a grid                    *
 *                                                                *
 *  MPC (04/15)                                                   *
 *                                                                *
 ******************************************************************/

double PEQS_Main(double *PAv, double Ediis)
{
  static bool firstCall=true;

  //cout << " --------------------------------------- " << endl;
  if (firstCall){
     printf("\tPEqS:  Begin Poisson Eqn. Solver\n");
  }

  PEQSclass PEQS;
  PEQS.PrepareSolver();
  PEQS.BoundaryCondition();
  PEQS.InitRhoPhiEps(PAv);
  if (firstCall){
     if(PEQS.poisson_boltzmann){
        std::cout <<"\tPEqs:  This is a Poisson-Boltzmann calculation." << std::endl;
        if(PEQS.linear){
          std::cout << "\tPEqs:  Using the linearized Poisson-Boltzmann equation." << std::endl;
        }
        else{
          std::cout << "\tPEqs:  Using the full (nonlinear) Poisson-Boltzmann equation." << std::endl;
        }
        std::cout << std::scientific;
        if(PEQS.size_modification){
          std::cout << "\tPEqs:  Calculations are solved for finite ion size with radius R = " << PEQS.ion_radius << " a.u." << std::endl;
        }
        else{
          std::cout << "\tPEqs:  Assuming point-like ions." << std::endl;
        }
        std::cout << "\tPEqs:  The Stern-layer thickness is set to a = " << PEQS.stern_thickness  << " a.u." << std::endl;
        std::cout << "\tPEqs:  The Boltzmann factor assumes the following temperature: " << PEQS.temperature/PEQS.conv_temp_au << " K" << std::endl;
        std::cout << "\tPEqs:  The electrolyte concentration is set to: " << PEQS.concentration << " mol/L." << std::endl;
     }
  }
  rem_write(0,REM_PEQ_CORR);

  //equil reference solvation calculation or nonequil reference state calculation 
  //if marcus job, create and store the slow polarization component
  if (PEQS.DoNeqJob() == 0 || (PEQS.DoNeqJob() == 1 && PEQS.NeqState() == 0)){
    //CS: added this for Poisson-Boltzmann code
    if (PEQS.poisson_boltzmann) PEQS.PBJob();
    else PEQS.ReferenceJob();
    if (PEQS.NeqPart() == 0) PEQS.MarcusRefJob();
  }

  //compute the reference slow/fast polarization component for pekar partition
  if (PEQS.DoNeqJob() == 1 && PEQS.NeqState() == 0 && PEQS.NeqPart() == 1)
    PEQS.PekarRefJob();

  //compute the ionized fast polarization component for marcus partition...
  if (PEQS.DoNeqJob() == 1 && PEQS.NeqState() == 1 && PEQS.NeqPart() == 0)
    PEQS.MarcusNeqJob();

  //...or for the pekar patition
  if (PEQS.DoNeqJob() == 1 && PEQS.NeqState() == 1 && PEQS.NeqPart() == 1)
    PEQS.PekarNeqJob();

  firstCall=false;
  if (PEQS.ConvPEQ()){
    if(PEQS.poisson_boltzmann){
      double G_solvation = PEQS.SolvEnergyPB(PEQS.el);
      return G_solvation;
    }
    else{
      double G_solvation = PEQS.SolvEnergy(Ediis);
    //cout << " --------------------------------------- " << endl;
      return G_solvation;
    }
  }
  else QCrash(" Sorry, the Poisson equation solver was unable to find a solution.");
}
