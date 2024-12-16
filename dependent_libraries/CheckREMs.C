#include "qchem.h"
#include "InputSection.h"
#include "TokenList.h"
#include "rem_values.h"
#include "functionals.h"
#include "dftcodes.h"
#include "JobTypes.h"
#include "sapt.h"

/*
 *
 * Here, we check the REMs for consistency and conflicts.  Mutually-exclusive combinations
 * of REMs should be added to this routine on a case-by-base basis.  The idea is to bail
 * out early if there are obvious signs of trouble in the input file.
 *
 * There is not a lot here at the moment (and I know there are more possible conflicts
 * than this), so please add to this routine...
 *
 * JMH (7/2014)
 *
*/

bool IsJobTypeOkay(int);

void ByeBye()
{
  printf("\n\n");

  // don't want to use QCrash here because that prompts user to submit a crash report
  printf("\n\n\tThe Q-Chem input file has failed to pass inspection\n\n");
  QKill();

}

void CheckREMs()
{ printf("\nChecking the input file for inconsistencies... ");

/* 
  // do this in PostProcessRemDefaults, else it doesn't trap some errors
  if (!IsJobTypeOkay(rem_read(REM_JOBTYP))){
     cout << "Urecognized job type\n"; ByeBye();
  }
*/

  ////////////// General ////////////////
  int NAtoms = rem_read(REM_NATOMS);
  int JobType= rem_read(REM_JOBTYP);
  if (NAtoms == 1)
  {
     if (JobType == JOBTYPE_GEOM_OPT || JobType == JOBTYPE_TS_OPT){
        cout << "\nGeometry optimization for one atom does not make sense.\n";
        ByeBye();
     }

/*
// JMH - need this for some relativistic jobs
     else if (JobType == JOBTYPE_FREQUENCY){
        cout << "Frequency calculation for one atom does not make sense.\n";
        ByeBye();
     }
*/

  }

  if (rem_read(REM_KONSCF) <= 0){
     cout << "\nSCF_CONVERGENCE set to <= 0 in the $rem section.\n"
          << "Please choose a positive value or use the default.\n"; 
     ByeBye();
  }

  //////////////// AIMD ////////////////////

  int IThermo = rem_read(REM_AIMD_THERMOSTAT);
  int Integr  = rem_read(REM_AIMD_INTEGRATION);
  if (IThermo == LANGEVIN && Integr != VVERLET){
     cout << "\nLangevin thermostat is available only with velocity Verlet integration.\n";
     ByeBye();
  }

  //////////////// SAPT ////////////////////
  if (rem_read(REM_SAPT_AO) == 1)
  {  // some features not available in AO-based SAPT code; need MO-based code instead
     if (rem_read(REM_SAPT_DISPERSION) == sapt_disp_SAPT0)
     {
        printf("\nYou have requested 2nd-order (SAPT0-style dispersion), but this is not available\n");
        printf("in the AO-based SAPT code.  Either choose an empirical dispersion correction or set\n");
        printf("ALGORITHM = MO in $sapt section to use 2nd-order dispersion.\n");
        ByeBye();
     }
     if (rem_read(REM_SAPT_COUPLE_IND) == 1){
        rem_write(0,REM_SAPT_COUPLE_IND);
        QWarn("Turning off 3-body induction couplings, which are not available for AO-based SAPT algorithm");
     }
  }
  if (rem_read(REM_SAPT_DISPERSION) == sapt_disp_SAPT0)
  { // SAPT empirical dispersion - check version number
     int iEmpDisp = rem_read(REM_SAPT_EMPIRICAL_DISP);
     if (iEmpDisp <= 0 || iEmpDisp > 3){
        printf("\n Requested SAPT + aiD empirical dispersion potential version %d\n",iEmpDisp);
        printf(" Only valid version numbers are 1, 2, and 3\n");
        ByeBye();
     }
  }

  ////////////// PCM Stuff ////////////////
  //
  // Solvation models use keyword input sections, and unrecognized keywords are ignored
  // by the PCM code.  Therefore let's check them here in case of misspellings, etc.
  // Right now we only check the keyword itself (TL[0]), not its value (TL[1], etc.).
  // Could check both at the expense of more lines of code...
  //
  //////////////////////////////////////////
  InputSection is("$pcm");
  if(is){
     TokenList TL;
     LOGICAL Done;
     do{
        is >> TL; 
        Done = TL && TL[0]  == "$end";
        if(!Done && TL.NSignifTokens()){ // list all allowed $pcm keywords here
        //if(!Done) // list all allowed $pcm keywords here
	   
           if (TL[0] != "EqS_Conv" && TL[0] != "EqSolv" && TL[0] != "EqState" && TL[0] != "Theory" && 
               TL[0] != "EqS_Ref" && TL[0] != "NonEquilibrium" && TL[0] != "Solver" && TL[0] != "SwitchThresh" && 
	       TL[0] != "Method" && TL[0] != "GridMethod" && TL[0] != "HPoints" && 
               TL[0] != "HeavyPoints" && TL[0] != "MMHPoints" &&
               TL[0] != "MMHeavyPoints" && TL[0] != "Radii" && TL[0] != "CavityRadius" && 
               TL[0] != "ProbeRadius" && TL[0] != "TorusRadius" && TL[0] != "alpha_smooth" &&
               TL[0] != "TorusGrid" && TL[0] != "TorusOnly" && TL[0] != "TorusZeta" &&
               TL[0] != "TPrint" && 
               TL[0] != "CavityCenter" && TL[0] != "vdwScale" && TL[0] != "SASradius" && 
               TL[0] != "PrintLevel" && TL[0] != "FirstIteReq" &&
               TL[0] != "K" && TL[0] != "Ksym" && TL[0] != "Ktrans" && TL[0] != "DComp" && 
               TL[0] != "PreCond" &&
               TL[0] != "doPara" && TL[0] != "NoBaseGuess" && TL[0] != "NoFastSSVPE" && 
               TL[0] != "NoFastCPCM" &&
               TL[0] != "AltSwitch" && TL[0] != "SumRule" && TL[0] != "SumSwitch" && TL[0] != "FixPVA2" &&
               TL[0] != "MidRange" && TL[0] != "FarRange" && TL[0] != "MulOrder" && TL[0] != "UseMultipole" &&
               TL[0] != "FineLength" && TL[0] != "ChargeSeparation" && TL[0] != "StateSpecific" &&
               TL[0] != "TdNonEq" &&
               TL[0] != "CGthresh" && TL[0] != "NoMatrix" && TL[0] != "NoneqIterMode" && 
               TL[0] != "NoneqRelax" &&
	       TL[0] != "SurfaceType" && TL[0] != "Resolution" && TL[0] != "DensityThresh" && 
               TL[0] != "RefineLevel" &&
	       TL[0] != "IsoElecConv" && TL[0] != "IsoElecMaxIter"  && TL[0] != "QKnorm"
              ){
                 cout << "\n\tAn unrecognized keyword \"" << TL[0] << "\" was found in the $pcm input section\n";
                 ByeBye();
              }
        }
// JMH - should check the various surface types and other options
     }while (!Done);
  }

  is.assign("$pcm_solvent");
  if(is){ 
      cout << "\nUse of $pcm_solvent section has been deprecated starting in Q-Chem v. 4.2.1.\n"
           << "Please use $solvent instead.  See the Q-Chem manual for more instructions.\n";
      ByeBye();
  }

  is.assign("$solvent");
  if (is){
     int iMethod = rem_read(REM_SOLVENT_METHOD);
     if (iMethod != PCM && iMethod != ONSAGER)
        QWarn("Found a $solvent section, which is only used for PCM and Kirkwood-Onsager solvation models");

     TokenList TL;
     LOGICAL Done;
     do{
        is >> TL; 
        Done = TL && TL[0]  == "$end";
        if(!Done)
        {  // list all allowed $solvent keywords here
           // they are different for different solvent models
           if (iMethod == PCM)
           {
              if (TL[0] != "Cav_water"      && TL[0] != "Cav_n-hexane"      && TL[0] != "Cav_benzene"  &&
                  TL[0] != "Cav_cholorform" && TL[0] != "Cav_cyclohexane"   && TL[0] != "Cav_methanol" &&
                  TL[0] != "Cav_ethanol"    && TL[0] != "Cav_toluene"       && TL[0] != "NonEls"       &&
                  TL[0] != "Dielectric"     && TL[0] != "OpticalDielectric" && TL[0] != "Temperature"  &&
                  TL[0] != "Pressure"       && TL[0] != "SolventRho"        && TL[0] != "SolventRadius"&& 
                  TL[0] != "ScreenLength"   && TL[0] != "NSolventAtoms"     && TL[0] != "SolventAtom"  &&
		  TL[0] != "Dielectric_Infi"
                 ){
                    cout << "\n\tAn unrecognized keyword \"" << TL[0] << "\" was found in the $solvent input section\n";
                    ByeBye();
              }
           }else if (iMethod == ONSAGER){
              if (TL[0] != "Dielectric" && TL[0] != "CavityRadius" && TL[0] != "MultipoleOrder")
              {
                 cout << "\n\tAn unrecognized keyword \"" << TL[0] << "\" was found in the $solvent input section\n";
                 ByeBye();
              }
           }
        }
     }while(!Done);
  } 

  is.assign("$chem_sol");
  if(is){
     int iMethod = rem_read(REM_SOLVENT_METHOD);
     if (iMethod != CHEM_SOL)
        QWarn("Found a $chem_sol input section but the Langevin Dipoles solvation model was not requested");

     TokenList TL;
     LOGICAL Done;
     do{
        is >> TL; 
        Done = TL && TL[0]  == "$end";
        if(!Done){  
           if (TL[0] != "EField" && TL[0] != "NGrids" && TL[0] != "Print" && TL[0] != "ReadRadii")
           {
              cout << "\n\tAn unrecognized keyword \"" << TL[0] << "\" was found in the $chem_sol input section\n";
              ByeBye();
           }
        }
     }while (!Done);
  }

  is.assign("$smx");
  if(is){
     int iMethod = rem_read(REM_SOLVENT_METHOD);
     if (iMethod != SM8 && iMethod != SM12 && iMethod != SMD)
        QWarn("Found a $smx input section but no SMx solvation model was requested");

     TokenList TL;
     LOGICAL Done;
     do{
        is >> TL;
        Done = TL && TL[0]  == "$end";
        if(!Done){  
           if (TL[0] == "Charges"){
             if(TL[1] == "MK")     rem_write(3,REM_SMX_SOLVATION);
             if(TL[1] == "CHELPG") rem_write(4,REM_SMX_SOLVATION);
           }
           if (TL[0] != "Solvent" && TL[0] != "Print" && TL[0] != "Charges")
           {
              cout << "\n\tAn unrecognized keyword \"" << TL[0] << "\" was found in the $smx input section\n";
              ByeBye();
           }
        }
     }while (!Done);
  }

///////////////////////////////////////////////
//ZQ- handle methods not supporting analytical gradient/Hessian
//Tickets- 107, 1524
// JMH - why do we not set IDERIV automatically, like every other Q-Chem job?
  const int isolvent = rem_read(REM_SOLVENT_METHOD);
  const int ibasis  = rem_read(REM_IBASIS);
  const int levcor = rem_read(REM_LEVCOR);
  const int eomcor = rem_read(REM_EOM_CORR);
  const int ideriv = rem_read(REM_IDERIV);
  int job_need_2nd = JobType == JOBTYPE_FREQUENCY;
  int job_need_1st = JobType == JOBTYPE_GEOM_OPT || JobType == JOBTYPE_TS_OPT ||
                     JobType == JOBTYPE_FORCE || JobType == JOBTYPE_REACTION_PATH ||
                     job_need_2nd == 1;
  if(ideriv > 0 && job_need_1st){
    if(isolvent == SM8 && (ibasis == 361010 || ibasis == 361011)){ // 6-31+G* and 6-31+G**
      printf("\n\n");
      printf("\tDesired analytical derivatives not available.\n");
      printf("\tPlease set IDERIV = 0 to run finite difference calculation.\n");
      ByeBye();
    }
  }

  if (isolvent > 0){
     // NMR is not supported by PCM
     if (JobType == NMR){
        cout << "\nSolvation models are not supported for use with NMR calculations.\n";
        ByeBye();
     }else if (rem_read(REM_D_SCF) > 0){
        cout << "\nThe D-SCF CPSCF code (for properties) does not support solvation models.\n";
        ByeBye();
     }
  }


  /////// End of PCM Stuff ///////////

///////////////////////
// $peqs
///////////////////////
  is.assign("$peqs");
  if(is){
     TokenList TL;
     LOGICAL Done;
     do{
        is >> TL;
        Done = TL && TL[0]  == "$end";
        if(!Done){  
           if (TL[0] != "SoluteCavity"    && TL[0] != "VDW_Type"       && TL[0] != "VDW_Scale"      &&
               TL[0] != "VDW_Shift"       && TL[0] != "RigidScale"     && TL[0] != "RHybrid"        &&
               TL[0] != "SphereRadius"    && TL[0] != "CavityCutoff"   && 
               TL[0] != "InterpolLength"  && TL[0] != "InterpolScale"  &&  
               TL[0] != "SphereCenter"    && TL[0] != "CenterType"     && 
               TL[0] != "AutoRadius"      && TL[0] != "Interface"      && TL[0] != "GibbsDS"        &&
               TL[0] != "InterfaceLength" && TL[0] != "InterfaceScale" && 
               TL[0] != "InterfaceDirection" &&
               TL[0] != "FDiffOrder"      && TL[0] != "Multigrid"      && TL[0] != "MG_Order"       &&
               TL[0] != "MG_Cycle"        && TL[0] != "MG_Sweeps"      && TL[0] != "MaxIter"        &&
               TL[0] != "CG_Thresh"       && TL[0] != "SolverThresh"   && TL[0] != "PolarIterScale" &&
               TL[0] != "NonEquilJob"     && TL[0] != "NonEquilState"  && TL[0] != "NonEquilPartition"  &&
               TL[0] != "SolventDiel"     && TL[0] != "OpticalDiel"    && TL[0] != "GaussBlur"     &&
               TL[0] != "GaussWidth"      && TL[0] != "BatchSize"      && TL[0] != "Print"         &&
               //CS: added some Poisson-Boltzmann keywords here
               TL[0] != "PoissonBoltzmann"&& TL[0] != "PB_Type"        && TL[0] != "ElectrolyteConcentration"  &&
               TL[0] != "Eta_PB"          && TL[0] != "Eta"            && TL[0] != "probe_radius"  &&
               TL[0] != "stern_thickness" && TL[0] != "rho_ions_thresh"&& TL[0] != "size_modification" &&
               TL[0] != "ion_radius"      && TL[0] != "MulOrder"       && TL[0] != "MaxThresh"       &&
               TL[0] != "MinThresh"       && TL[0] != "SolventRadius"  && TL[0] != "FillZero"       &&
               TL[0] != "Scale"           && TL[0] != "zhi_switch"     && TL[0] != "eta_switch"     &&
               TL[0] != "PotentialGrid"   && TL[0] != "Idw_power"      && TL[0] != "InterpType"    &&
               TL[0] != "ConvoLim"        && TL[0] != "InterpOrder"    && TL[0] != "MaxDer_Interp"   &&
               TL[0] != "LowLim"          && TL[0] != "HighLim"  
              )
           {
              cout << "\n\tAn unrecognized keyword \"" << TL[0] << "\" was found in the $peqs input section\n";
              ByeBye();
           }
        }
     }while (!Done);
  }
///////////////////////
// $xpol
///////////////////////
  is.assign("$xpol");
  if(is){
     TokenList TL;
     LOGICAL Done;
     do{
        is >> TL;
        Done = TL && TL[0]  == "$end";
        if(!Done){  
           if (TL[0] != "embed" && TL[0] != "charges" && TL[0] != "dft-lrc" && TL[0] != "print")
           {
              cout << "\n\tAn unrecognized keyword \"" << TL[0] << "\" was found in the $xpol input section\n";
              ByeBye();
           }
        }
     }while (!Done);
  }
///////////////////////
// $sapt
///////////////////////
  is.assign("$sapt");
  if(is){
     TokenList TL;
     LOGICAL Done;
     do{
        is >> TL;
        Done = TL && TL[0]  == "$end";
        if(!Done)
        {  
           if (TL[0] != "print"         && TL[0] != "order"     && TL[0] != "basis"    && 
               TL[0] != "dscf"          && TL[0] != "algorithm" && TL[0] != "cdft-eda" &&
               TL[0] != "EmpiricalDisp" && TL[0] != "exchange"  && TL[0] != "cphf"     &&
               TL[0] != "3b-ind"
              )
           {
              cout << "\n\tAn unrecognized keyword \"" << TL[0] << "\" was found in the $sapt input section\n";
              ByeBye();
           }
        }
     }while (!Done);
  }
///////////////////////
// $mbe 
///////////////////////
  is.assign("$mbe");
  if(is){
     TokenList TL;
     LOGICAL Done;
     do{
        is >> TL;
        Done = TL && TL[0] == "$end"; 
        if (!Done){
           if (TL[0] != "order"      && TL[0] != "BSSE"       && TL[0] != "embed" &&
               TL[0] != "bsse_order" && TL[0] != "bsse_type"
              )
           {
              cout << "\n\tAn unrecognized keyword \"" << TL[0] << "\" was found in the $mbe input section\n";
              ByeBye();
           }
        }
     }while (!Done);
  }

/////////////////////////////////////////////////////////////////////////

  int unrestricted = rem_read(REM_JUSTAL);
  if (rem_read(REM_FRACTIONAL_ELECTRON) > 0 && !unrestricted){
     printf("\nPlease use an unrestricted formalism for calculations with fractional electrons\n");
     ByeBye();
  }

  // There's some problem with reading the energy and/or dipole files
  // that isn't worth fixing.
  if (((JobType == JOBTYPE_DIPOLE && ideriv == 0) ||
       (JobType == JOBTYPE_POLARIZABILITY)) &&
      (useCCMAN(levcor) == true) &&
      (rem_read(REM_CCMAN2) == 0)) {
    cout << "\nCan't do FD dipole or FD polarizability with ccman1, use ccman2 instead." << endl;
    ByeBye();
  }

  const bool use_libresponse = rem_read(REM_RESPONSE) > 0;
  if (use_libresponse &&
      ((levcor > 0) || (eomcor > 0))) {
    cout << "\nCan't do general linear response with anything other than HF or DFT." << endl;
    ByeBye();
  }

  const bool use_polman = (JobType == JOBTYPE_POLARIZABILITY) && (rem_read(REM_RESPONSE_POLAR) > -1);
  const bool use_responseman = use_libresponse || use_polman;
  const int iNLC = rem_read(REM_NL_CORRELATION);
  if (use_responseman && (iNLC > 0)) {
      cout << "\nDensity functionals with NLC currently not supported for responseman-based routines." << endl;
      ByeBye();
  }

  ////////////// FDEman Stuff ////////////////
  InputSection InpSec("$fde");
  if(InpSec){
    TokenList TL;
    LOGICAL Done;
    bool isXC = false;
    bool isXC_B = false;
    bool isX = false;
    bool isX_B = false;
    do{
      InpSec >> TL;
      Done = TL && TL[0]  == "$end";
      if(!Done && TL.NSignifTokens()){
        //check for invalid keywords
        if (TL[0] != "T_func" && TL[0] != "X_func" && TL[0] != "C_func" && TL[0] != "XC_func" &&
            TL[0] != "expansion" && TL[0] != "rhoB_method" && TL[0] != "rhoB_basis" && TL[0] != "PrintLevel" && TL[0] != "Debug" &&
            TL[0] != "prepol" && TL[0] != "prepolarization"	&& TL[0] == "prepol_type" && TL[0] == "prepolarization_type" &&
            TL[0] != "polarization" && TL[0] != "pol_B" && TL[0] != "polarization_B" &&
            TL[0] != "polarization_type" && TL[0] != "pol_B_type" && TL[0] != "polarization_B_type" &&
            TL[0] != "X_func_B" && TL[0] != "C_func_B" && TL[0] != "XC_func_B"
           ){
          cout << "\n\t An unrecognized keyword \"" << TL[0] << "\" was found in the $fde input section\n";
          ByeBye();
        }
        // check whether TL[1] is a valid X or XC functional
        if (TL[0] == "XC_func" ) {
          isXC = true;
          if (TL[1] != "PBE" && TL[1] != "BLYP" && TL[1] != "BP86" && TL[1] != "PW91"){
            int xccode = XCode(TL[1]);
            if (xccode >= 0){
              cout << "\n\t XCcode = " << xccode << ". Has to be < 0 ! Wrong functional for XC_Func." << endl;
              ByeBye();
            }
          }
        }

        if (TL[0] == "XC_Func_B" ) {
          isXC_B = true;
          if (TL[1] != "PBE" && TL[1] != "BLYP" && TL[1] != "BP86" && TL[1] != "PW91"){
            int xccode = XCode(TL[1]);
            if (xccode >= 0){
              cout << "\n\t XCcode = " << xccode << ". Has to be < 0 ! Wrong functional for XC_Func_B." << endl;
              ByeBye();
            }
          }
        }


        if (TL[0] == "X_func" ){
          int xcode = XCode(TL[1]);
          if (xcode < 0) {
            cout << "\n\t Cannot use a XC functional as X functional for X_Func." << endl;
            ByeBye();
          }
        }

        if (TL[0] == "X_Func_B" ){
          int xcode = XCode(TL[1]);
          if (xcode < 0) {
            cout << "\n\t Cannot use a XC functional as X functional for X_Func_B" << endl;
            ByeBye();
          }
        }

        // check whether XC or X+C
        if (TL[0] == "X_func" || TL[0] == "C_func") {isX = true;}
        if (isXC && isX) {
          cout << "\n\t Multiple definitions of XC functional in $fde section! Use either XC_func or X_func + C_func!" << endl;
          ByeBye();
        }

        if (TL[0] == "X_func_B" || TL[0] == "C_func_B") {isX_B = true;}
        if (isXC_B && isX_B) {
          cout << "\n\t Multiple definitions of XC functional in $fde section! Use either XC_func_B or X_func_B + C_func_B!" << endl;
          ByeBye();
        }

      }
    }while (!Done);
  }

  ///////// NMR ///////////
  if (JobType == JOBTYPE_NMR  || 
      JobType == JOBTYPE_ISSC || 
      JobType == JOBTYPE_DYNPOLAR ||
      JobType == JOBTYPE_HYPERPOLAR ||
      rem_read(REM_MOPROP) > 0)
  {
    //printf(" NMR JobType = %d\n",JobType);
     //
     // check that the functional is actually supported for NMR (for electric jobs in MOPROPMAN also)
     // JMH (3/18): meta-GGAs, range-separated hybrids, and nonlocal correlation are definitely not 
     // supported but I'm not sure this list of unsupported functionals is complete
     //
     XCFunctional XCFunc = rem_read(REM_LEVEXC) != -1 ? read_xcfunc(rem_read(REM_LEVEXC),
         (rem_read(REM_LEVCOR) <= XCFUNC_MAX) ? rem_read(REM_LEVCOR):0) : 0;

     bool shouldExit=false;
     if (XCFunc.IsCoulombAtten() || rem_read(REM_LRC_DFT) > 0){
        cout << "\n\n\tRange-separated hybrid functionals are ";
        shouldExit=true;
     }else if (XCFunc.HasTau() || XCFunc.HasLap() || XCFunc.IsMeta()){
        cout << "\n\n\tMeta-GGA functionals are ";
        shouldExit=true;

     }else if (!XCFunc.CanUse3rdDer()){  // this will pick up VV10
        if (JobType == JOBTYPE_NMR || JobType == JOBTYPE_ISSC) {
           // NMR and ISSC should be fine.
        }else {
           cout << "\n\n\tThis functional is ";
           shouldExit=true;
        }
     }

     if (shouldExit){
        cout << "not supported for NMR or electric properties calculations\n";
        ByeBye();
     }
  }

  ///////////////////////////// All done  ////////////////////////////////////////////////
  printf("\t...done.\n");
}


/////////////////////////////////////

bool IsJobTypeOkay(int j)
{
// There is nothing else to check that user didn't try JOBTYPE = xxx where 'xxx'
// is an acceptable option for some *other* REM....

  bool okay=false; 
  if (j == JOBTYPE_SINGLE_POINT   ||
      j == JOBTYPE_GEOM_OPT       ||
      j == JOBTYPE_TS_OPT         ||
      j == JOBTYPE_POLARIZABILITY ||
      j == JOBTYPE_FREQUENCY      ||
      j == JOBTYPE_REACTION_PATH  ||
      j == JOBTYPE_FORCE          ||
      j == JOBTYPE_NMR            ||
      j == JOBTYPE_AIMD           ||
      j == JOBTYPE_BSSE           ||
      j == JOBTYPE_EDA            ||
      j == JOBTYPE_PIMC           ||
      j == JOBTYPE_PIMD           ||
      j == JOBTYPE_FSM            ||
      j == JOBTYPE_GSM            ||
      j == JOBTYPE_PES_SCAN       ||
      j == JOBTYPE_BH             ||
      j == JOBTYPE_RAND           ||
      j == JOBTYPE_DIPOLE         ||
      j == JOBTYPE_STATPOLAR      ||
      j == JOBTYPE_DYNPOLAR       ||
      j == JOBTYPE_HYPERPOLAR     ||
      j == JOBTYPE_ISSC           ||
//      j == JOBTYPE_XPOL           ||
      j == JOBTYPE_XSAPT          ||
      j == JOBTYPE_RRAMAN         ||
      j == JOBTYPE_SCF_GUESS      ||
      j == JOBTYPE_NUCGRAD_ECP    ||
      j == JOBTYPE_SIMPER         ||
      j == JOBTYPE_NUCENERGY_ECP  ||
      j == JOBTYPE_NUCGRAD_ECP_2 
  )okay=true;
  //cout << "JobType = " << j << ".  Okay? " << okay << endl;
  return okay; 
}
