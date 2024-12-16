#include "peqs_class.h"
#include "rem_values.h"
#include "InputSection.h"
#include "TokenList.h"

void PEQSclass::Read_PEQS()
{
  InputSection is("$peqs");
  if(is){
    TokenList TL;
    LOGICAL Done;
    do{
       is >> TL;
       Done = TL && TL[0]  == "$end";
       if(!Done){
         if(TL[0] == "SoluteCavity")
         {
           if (TL[1] == "None")
             CavityType = PEQS_CAVTYPE_NONE;
           else if(TL[1] == "RigidVDW")
             CavityType = PEQS_CAVTYPE_VDW;
           else if(TL[1] == "Spherical")
             CavityType = PEQS_CAVTYPE_SPHERE;
	   else if(TL[1] == "Arbitrary") 
             CavityType = PEQS_CAVTYPE_ARBI; 
           else QCrash (" Requested solute cavity not found.");
         }
         if(TL[0] == "VDW_Type")
         {
           if (TL[1] == "Unscaled")
             vdwType = PEQS_VDW_UNSCALED;
           else if (TL[1] == "Scaled")
             vdwType = PEQS_VDW_SCALED;
           else if (TL[1] == "Shifted")
             vdwType = PEQS_VDW_SHIFTED;
           else if (TL[1] == "Hybrid")
             vdwType = PEQS_VDW_HYBRID;
           else QCrash (" Requested VDWType not found.");
         }
         if(TL[0] == "VDW_Scale")
           vdwScale = TL[1].GetDouble();
         if(TL[0] == "VDW_Shift")
           vdwShift = TL[1].GetDouble();
         if(TL[0] == "RigidScale")
           erf_delta = TL[1].GetDouble()*ang2bohr;
         if(TL[0] == "RHybrid")
           r_hybrid = TL[1].GetDouble()*ang2bohr;        
         if(TL[0] == "SphereRadius")
           r_cavity = TL[1].GetDouble()*ang2bohr;
         if(TL[0] == "CavityCutoff")
           tanh_cut1 = TL[1].GetDouble()*ang2bohr;
         if(TL[0] == "InterpolLength")
           tanh_cut2 = TL[1].GetDouble()*ang2bohr;
         if(TL[0] == "InterpolScale")
           tanh_scale = TL[1].GetDouble()/ang2bohr;
         if(TL[0] == "SphereCenter"){
           CavityType = PEQS_CAVTYPE_SPHERE;
           FindCenter = 1;
           Center = QAllocDouble(3);
           Center[0] = TL[1].GetDouble()*ang2bohr;
           Center[1] = TL[2].GetDouble()*ang2bohr;
           Center[2] = TL[3].GetDouble()*ang2bohr;
         }
         if(TL[0] == "CenterType"){
           if (TL[1] == "Auto")
             FindCenter = 2;
           else if (TL[1] == "COM")
             FindCenter = 3;
           else QCrash (" Requested method for finding the spherical cavity center not found.");
         }
         if(TL[0] == "AutoRadius"){
           if (TL[1] == "true") 
             AutoRadius = true;
           else if (TL[1] == "false") 
             AutoRadius = false;
         }
         if(TL[0] == "Interface"){
           if (TL[1] == "true") 
             Interface = true;
           else if (TL[1] == "false") 
             Interface = false;
         }
         if(TL[0] == "GibbsDS")
           gibbsDS = TL[1].GetDouble()*ang2bohr;
         if(TL[0] == "InterfaceLength")
           r_interface = TL[1].GetDouble()*ang2bohr; 
         if(TL[0] == "InterfaceScale")
           interface_scale = TL[1].GetDouble()/ang2bohr;
         if(TL[0] == "InterfaceDirection")
         {
           if (TL[1] == "Positive") 
             z_direction = 1;
           else if (TL[1] == "Negative") 
             z_direction = -1;
         } 
         if(TL[0] == "FDiffOrder")
           fdif_order = TL[1].GetINTEGER(); 
         if(TL[0] == "Multigrid")
         {
           if (TL[1] == "true") 
             multigrid = true;
           else if (TL[1] == "false") 
             multigrid = false;
         }  
         if(TL[0] == "MG_Order")
           mg_order = TL[1].GetINTEGER();
         if(TL[0] == "MG_Cycle")
         {
           if (TL[1] == "VCycle")
             mg_cycle = 1;
           else if(TL[1] == "WCycle")
             mg_cycle = 2;
         }
         if(TL[0] == "MG_Sweeps"){
           sweep1 = TL[1].GetINTEGER();
           sweep2 = TL[2].GetINTEGER();
           sweep3 = sweep1 + sweep2;
         }
         if(TL[0] == "MaxIter")
           max_iter = TL[1].GetINTEGER();

         if(TL[0] == "CG_Thresh")
           cg_thresh = TL[1].GetDouble();
         if(TL[0] == "SolverThresh")
           solver_thresh = TL[1].GetDouble();
         if(TL[0] == "PolarIterScale")
           eta = TL[1].GetDouble();
         if(TL[0] == "NonequilJob")
           doNeqJob = true; 
         if(TL[0] == "NonequilState")
         {
           if(TL[1] == "Reference")
             noneq_state = 0;
           else if(TL[1] == "Ionized")
             noneq_state = 1;
           else QCrash (" Invalid state of the non-equilibrium system.");
         }
         if(TL[0] == "NonequilPartition"){
           if (TL[1] == "Marcus")
             noneq_partition = 0;
           else if (TL[1] == "Pekar")
             noneq_partition = 1;
           else QCrash (" Requested non-equilibrium partitioning scheme not found.");
         }

//================Interpolarion================//

         if(TL[0] == "PotentialGrid"){
           if(TL[1] == "Cartesian") potgrid = 1;
           if(TL[1] == "ACG") potgrid = 2;
         }
         if(TL[0] == "Occ-CutOff")
           cutoff = TL[1].GetINTEGER();
         if(TL[0] == "MulOrder")
           MulOrder = TL[1].GetINTEGER();
         if(TL[0] == "InterpType"){
           if(TL[1] == "IDW") interp = 1;
           if(TL[1] == "Spline") interp = 2;
	   if(TL[1] == "DivDiff") interp = 3;
         }
	 if(TL[0] == "InterpOrder"){
           order = TL[1].GetINTEGER();
         }
	 if(TL[0] == "MaxDer_Interp"){
           max_der = TL[1].GetINTEGER();
         }
         if(TL[0] == "Idw_power"){
           if(TL[1] == "1") power = 1;
           if(TL[1] == "2") power = 2;
           if(TL[1] == "3") power = 3;
         }
	 if(TL[0] == "InStart"){
           shift1 = TL[1].GetDouble()*ang2bohr;
         }
         if(TL[0] == "OutStart"){
           shift2 = TL[1].GetDouble()*ang2bohr;
         }
	 if(TL[0] == "PotStart"){
           shift3 = TL[1].GetDouble()*ang2bohr;
         }
	 if(TL[0] == "FilterWidth"){
           filter = TL[1].GetDouble();
         }
	 if(TL[0] == "InHouse"){
           inhouse = TL[1].GetINTEGER();
         }
         if(TL[0] == "OutHouse"){
           outhouse = TL[1].GetINTEGER();
         }

//=============================================//
          
         if(TL[0] == "SolventDiel")
           eps_solvent = TL[1].GetDouble();
         if(TL[0] == "OpticalDiel")
           eps_infi = TL[1].GetDouble();
         if(TL[0] == "GaussBlur"){
           if (TL[1] == "true") 
              GaussBlur = true;
           else if (TL[1] == "false") 
              GaussBlur = false;
         }
         if(TL[0] == "GaussWidth")
           sigma = TL[1].GetDouble()*ang2bohr;
         if(TL[0] == "BatchSize")
           MaxPPer = TL[1].GetINTEGER();
         if(TL[0] == "print")
           peqs_print = TL[1].GetINTEGER();
         //CS:add some input values for Poisson-Boltzmann
         if(TL[0] == "PoissonBoltzmann"){
           if (TL[1] == "true")
             poisson_boltzmann = true;
           else if (TL[1] == "false")
             poisson_boltzmann = false;
         }
         if(TL[0] == "PB_TYPE"){
           if (TL[1] == "linear")
              linear = true;
           else if (TL[1] == "nonlinear")
              linear = false;
           else QCrash ("Requested Poisson-Boltzmann type not implemented.");
         }
         if(TL[0] == "size_modification"){ 
           if (TL[1] == "true")
              size_modification = true;
           else if (TL[1] == "false")
              size_modification = false;
           else QCrash ("Size Modification can only be true or false.");
         }
         if(TL[0] == "ElectrolyteConcentration")
           concentration = TL[1].GetDouble();
         if(TL[0] == "Eta_PB")
           eta_pb = TL[1].GetDouble();
         if(TL[0] == "Eta")
           eta = TL[1].GetDouble();
         if(TL[0] == "rho_ions_thresh")
           rho_ions_thresh = TL[1].GetDouble();
         if(TL[0] == "stern_thickness")
           stern_thickness = TL[1].GetDouble()*ang2bohr;
         if(TL[0] == "ion_radius")
           ion_radius = TL[1].GetDouble()*ang2bohr;
         if(TL[0] == "probe_radius")
           probe_radius = TL[1].GetDouble()*ang2bohr;
       }
    }while(!Done);
  }
}

void PEQSclass::Read_PEQSgrid()
{
  if (NPoints == NULL) NPoints = QAllocINTEGER(mg_order*NumDim); 
  if (NTotPts_mg == NULL) NTotPts_mg = QAllocINTEGER(mg_order);
  if (XMin == NULL) XMin = QAllocDouble(NumDim);
  if (XMax == NULL) XMax = QAllocDouble(NumDim);
  if (delX == NULL) delX = QAllocDouble(mg_order*NumDim);
  if (CavityType == PEQS_CAVTYPE_SPHERE && FindCenter > 1) Center = QAllocDouble(3); 
  
  InputSection is("$peqs_grid");
  if(is){
    TokenList TL;
    LOGICAL Done;
    do{
       is >> TL;
       Done = TL && TL[0]  == "$end";
       if(!Done){
         if (TL[0] == "DimX"){
           NPoints[0] = TL[1].GetINTEGER();
           XMin[0] = TL[2].GetDouble()*ang2bohr;
           XMax[0] = TL[3].GetDouble()*ang2bohr;
         } 
         if (TL[0] == "DimY"){
           NPoints[1] = TL[1].GetINTEGER();
           XMin[1] = TL[2].GetDouble()*ang2bohr;
           XMax[1] = TL[3].GetDouble()*ang2bohr;
         }

         if (TL[0] == "DimZ"){
           NPoints[2] = TL[1].GetINTEGER();
           XMin[2] = TL[2].GetDouble()*ang2bohr;
           XMax[2] = TL[3].GetDouble()*ang2bohr;
         }
       }
    }while(!Done);
  }

  if(CavityType == PEQS_CAVTYPE_ARBI) Read_EpsilonGrid();

  NTotPts = NPoints[0]*NPoints[1]*NPoints[2];
  if (MaxPPer == 0){
    if (NTotPts % 3 == 0) MaxPPer = NTotPts/3;
    else MaxPPer = NTotPts/NPoints[0];
  }
}

void PEQSclass::Read_EpsilonGrid()
{
  double *pts = QAllocDouble(4);
  int NGridPoints = 0;
  InputSection is("$epsilon");
  if(is){
      TokenList TL;
      LOGICAL Done;
      do{
        is >> TL; 
        Done = TL && TL[0] == "$end";
        if(!Done){
        if (TL.NTokens() != 4) QCrash("Problem reading $epsilon section");
        else{ 
          for (int i=0; i<4; i++) pts[i]=TL[i].GetDouble();
        FileMan(FM_WRITE,FILE_PEQS_GRID,FM_DP,4,0,FM_CUR,pts);
        NGridPoints++;
        }
       }
      }while(!Done);
      QFree(pts);
   }
}

