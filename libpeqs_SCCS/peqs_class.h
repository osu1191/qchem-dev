#ifndef _PEQS_CLASS_H
#define _PEQS_CLASS_H

#include <stdio.h>
#include "qccfg.h"
#include "qcmath2.h"
#include "qcmath.h"
#include "qchem.h"
#include <cstdlib>
#include <iostream>
#include "convfac.h"
#include "DFTBasisSet.h"
#include "BasisSet.hh"
#include "BSetMgr.hh"
#include "ShlPrs.hh"
#include "Electrolyte.h"
#include <armadillo>

#define INVALID_INT  -99999
#define INVALID_REAL -99999.9

class PEQSclass {
private:

   enum PEQS_CavType {
      PEQS_CAVTYPE_NONE  =0,
      PEQS_CAVTYPE_VDW   =1,
      PEQS_CAVTYPE_SPHERE=2,
      PEQS_CAVTYPE_SCCS  =3,
      PEQS_CAVTYPE_ARBI  =4
   };
   enum PEQS_VDW {
      PEQS_VDW_UNSCALED=1,
      PEQS_VDW_SCALED  =2,
      PEQS_VDW_SHIFTED =3,
      PEQS_VDW_HYBRID  =4
   };

   int NumDim,NTotPts,*NPoints,mg_order,mg_level,fdif_order,peqs_print;
   int *NTotPts_mg;
   int CavityType,GaussBlur,z_direction,iter_polar,iter_cg,max_iter,PBcallNo;
   int mg_cycle,sweep1,sweep2,sweep3,noneq_state,noneq_partition,MaxPPer;
   int *BCX,*BCY,*BCZ,*BCX_mg,*BCY_mg,*BCZ_mg,*AtNo,NAtoms;
   int vdwType,electron,nuclear,total,total_fast,FindCenter;
   bool AllocateMem,multigrid,Interface,doNeqJob,converge_peq,converge_cg,converge_polar;
   bool converge_total,converge_rho,AutoRadius;
   double sigma,solver_thresh,cg_thresh;
   double *XMin,*XMax,*delX,*XYZ,*Carts,*Center;
   double *vdWScale,*epsilon,*log_epsilon,*delLogEps,*delPhi,*delPhiFast,*delRhoElec;
   double *rho_total,*rho_elec,*rho_nuke,*rho_iter,*rho_solute,*rho_polar;
   double *rhotot_fast,*rhopol_ref_slow,*rhopol_ref_fast;
   double *phi_total,*phi_elec,*phi_nuke,*phi_polar,*phi_solute;
   double *phitot_fast,*phipol_ref_slow,*phipol_ref_fast;
   double Pi,eps_solvent,eps_vacuum,eps_infi,dielectric,opt_dielec,marcus_slow,marcus_fast;
   double vdwScale,vdwShift,erf_delta,tanh_scale,tanh_cut1,tanh_cut2;
   double r_cavity,r_hybrid,r_interface,interface_scale,gibbsDS;
   double alpha_cg,eta,residual,ang2bohr,au2kcal;
   double *Abasis,*Aphi,*resid,*basis,*zvec;
   double *Abasis_3,*Abasis_2,*Abasis_1;
   double *Aphi_3,*Aphi_2,*Aphi_1;
   double *resid_3,*resid_2,*resid_1;
   double *basis_3,*basis_2,*basis_1;
   double *zvec_3,*zvec_2,*zvec_1;
   double *rstrct_3,*rstrct_2,*rstrct_1;
   double *phi_corr3,*phi_corr2,*phi_corr1;
   double *rho_total_mg,*rho_total_temp,*rho_elec_mg,*rho_elec_temp,*rho_nuke_mg,*rho_nuke_temp;

//===== Solvent aware surface caviy ==============================//
   double *phi_sccs, *Veps;
   double rho_max, rho_min;
   INTEGER NBasis;
   int limit;
   double sccs_switch;
   double Rsol,zhi_delta,eta_delta,FillZero,alfa;  
   double *npro,*t_eta,*FillFrac,*mod_sr, *tst_sr; 
   double *conv,*norm; 
   double *eps, *ifunc;
//===============================================//

//CS: added some variables for Poisson-Boltzmann code
   bool converge_rho_ions,converge_phi_ions,converge_phi_polar;
   double rho_ions_thresh;
   double eta_pb;
   double probe_radius;
   double *lambda,*rho_ions,*phi_ions;

public:
   PEQSclass();
  ~PEQSclass();

   int NZ()            const{return NPoints[2];}
   double switch1()      const{return sccs_switch;}
   int NPts()          const{return NTotPts;}
   int Batch()         const{return MaxPPer;}
   int NeqState()      const {return noneq_state;}
   int NeqPart()      const{return noneq_partition;}
   bool DoNeqJob()    const{return doNeqJob;}
   void SetMGLevel(int i) {mg_level = i;}
   void SetConvPEQ(bool conv) {converge_peq = conv;}
   void ResetConvPEQ() {converge_peq = false;}
   void ResetCGIter() {iter_cg=0;converge_cg=false;}
   void ResetPolarIter() {iter_polar=0;converge_polar=false;}
   bool ConvPEQ() const {return converge_peq;}
   double xMin(int i)  const {return XMin[i];}
   double xMax(int i)  const {return XMax[i];}
   double dX(int i)    const {return delX[i];}
   double Length(int i) const {return XMax[i]-XMin[i];}
   double valueOfX(int i) const {return XMin[0]+delX[0]*(double)i;}
   double valueOfY(int i) const {return XMin[1]+delX[1]*(double)i;}
   double valueOfZ(int i) const {return XMin[2]+delX[2]*(double)i;}

//CS: added some variables for Poisson-Boltzmann code
   bool poisson_boltzmann;
   bool size_modification;
   bool linear;
   double temperature;
   double stern_thickness;
   double ion_radius;
   double concentration;
   double conv_temp_au;
   Electrolyte el;



   int get_multInd(int i, int j=0, int k=0, bool restriction=false, bool interpolation=false) const {
     int l,mg_index=0;
     if (multigrid) mg_index = mg_order-mg_level;
     l = mg_index;
     if (restriction) l -= 1;
     else if (interpolation) l += 1;
     int NyNz=NPoints[3*l+1]*NPoints[3*l+2],Nz=NPoints[3*l+2];
     return i*NyNz+j*Nz+k;
   }

   void get_3ind_from_multInd(int& i, int& j, int& k, int ijk, bool restriction=false) const{
      int l,mg_index=0;
      if (NumDim != 3) QCrash(" NumDim must be three!!!");
      if (multigrid) mg_index = mg_order-mg_level;
      l = mg_index;
      if (restriction) l -= 1;
      int NyNz=NPoints[3*l+1]*NPoints[3*l+2],Nz=NPoints[3*l+2];
      k = ijk % Nz;
      j = ( (ijk % NyNz) - k )/Nz;
      i = ( ijk - j*Nz - k )/NyNz;
   }

   double volumeElmnt()  const
   {
      double volElt=1.0;
      for (int i=0; i < NumDim; i++) volElt *= delX[i];
      return volElt;
   }
  
   double Integrator(double *integrand) const{
      int i;
      double integral=0.0,dV=volumeElmnt();
      #pragma omp parallel \
      private(i) \
      firstprivate(dV) \
      reduction(+:integral)
      {
        #pragma omp for
        for (i=0; i < NTotPts; i++) integral += integrand[i]*volumeElmnt();
      }
      return integral;

   }
   
   void Print(void) {
     int i,j,k,ijk,Nx=NPoints[0],Ny=NPoints[1],Nz=NPoints[2];
     ofstream dataFileZ("PEQS.data");
     ofstream dataFileA("Debug_x.data");
     ofstream dataFileB("Debug_y.data");
     ofstream dataFileC("Debug_z.data");
     ofstream dataFileD("multipoleX.data");
     ofstream dataFileE("multipoleY.data");
     ofstream dataFileF("multipoleZ.data");
     ofstream dataFileG("epsALL.data");
     for (i=0;i<Nx;i++){
       for (j=0;j<Ny;j++){
         for (k=0;k<Nz;k++){
           ijk = get_multInd(i,j,k);
           double x = valueOfX(i);
           double y = valueOfY(j);
           double z = valueOfZ(k);
	   dataFileG << x/ang2bohr  << '\t' << y/ang2bohr  << '\t' << z/ang2bohr  << '\t' << epsilon[ijk] << endl;
           if (z==0.0) dataFileZ << x << '\t' << y << '\t' 
               << rho_total[ijk] << '\t' << rho_elec[ijk]  << '\t' << rho_nuke[ijk]   << '\t' << '\t' 
               << rho_polar[ijk] << '\t' << phi_total[ijk] << '\t' << phi_elec[ijk]   << '\t' 
               << phi_nuke[ijk]  << '\t' << phi_polar[ijk] << '\t' << delPhi[3*ijk+2] << '\t' 
               << epsilon[ijk]   << '\t' << delLogEps[3*ijk+2] << '\t' << rho_ions[ijk] << '\t' 
               << '\t' << rho_iter[ijk]  << '\t' << lambda[ijk] << '\t' << valCart[ijk] << endl; 
	   if(z==0.0 && y==0.0){
                dataFileA << x/ang2bohr  << '\t' << delLogEps[3*ijk+0] << '\t' << rho_polar[ijk] << endl; 
                }
           if(x==0.0 && z==0.0){
                dataFileB << y/ang2bohr  << '\t' << delLogEps[3*ijk+1] << '\t' << rho_polar[ijk] << endl;  
                }
           if(y==0.0 && x==0.0){
                dataFileC << z/ang2bohr  << '\t'  << delLogEps[3*ijk+2] << '\t' << rho_polar[ijk] << endl;  
                }
         }
       }
     }
     dataFileG.close();
     dataFileZ.close();
     dataFileA.close();
     dataFileB.close();
     dataFileC.close();
     dataFileD.close();
     dataFileE.close();
     dataFileF.close();
   }

   void PrintEpsilon(void) {
     int i,j,k,ijk,Nx=NPoints[0],Ny=NPoints[1],Nz=NPoints[2];
     ofstream dataFileXY("EpsilonXY.data");
     ofstream dataFileYZ("EpsilonYZ.data");
     ofstream dataFileZX("EpsilonZX.data");
     ofstream dataFileRXY("RhoElecXY.data");
     ofstream dataFileRYZ("RhoElecYZ.data");
     ofstream dataFileRZX("RhoElecZX.data");
     for (i=0;i<Nx;i++){
       for (k=0;k<Nz;k++){
         double x = valueOfX(i);
         double z = valueOfZ(k);
         for (j=0;j<Ny;j++){
           double y = valueOfY(j);
           ijk = get_multInd(i,j,k);
           if (y==0.0 && x == 0.0){ 
		dataFileZX << z << '\t'  << epsilon[ijk] << '\t' << rho_polar[ijk] << '\t' << phi_polar[ijk] << '\t' << phi_sccs[ijk] << endl;
	   }
	   if (x==0.0 && y==0.0){ 
                dataFileYZ << y << '\t'  << epsilon[ijk] << '\t' << rho_polar[ijk] << '\t' << phi_polar[ijk] << '\t' << phi_sccs[ijk] << endl;
           }
	   if (z==0.0 && y==0.0){ 
                dataFileXY << x << '\t'  << epsilon[ijk] << '\t' << rho_polar[ijk] << '\t' << phi_polar[ijk] << '\t' << phi_sccs[ijk] << endl;
           }
         }
       }
     }
     dataFileXY.close();
     dataFileYZ.close();
     dataFileZX.close();
     dataFileRXY.close();
     dataFileRYZ.close();
     dataFileRZX.close();
   }

// Solvent Aware Surface cavity =========================//
   void GetSCCS(); 
   void SolventAwareSurfaceCavity(); 
   double getProbe(double); 
// ======================================================//

   void Read_PEQS();
   void Read_PEQSgrid();
   void Read_EpsilonGrid();
   void DiscreetBulkCavity();
   void Allocate();  
   void Coordinates();
   void PrepareSolver();
   void BoundaryCondition();
   void InitRhoPhiEps(double*);
   void ElectronicRho(double*);
   void EvalDensity(double*);
   void NuclearRhoPhi();
   void ComputeDielectric(bool);
   void RigidVDWSurface();
   void RigidVDWInterface();
   void SphereCenter(double*);        
   void SphereRadius();
   void SphericalSurface(double*);
   void SphericalInterface(double*);
   void HybridAxes(double*);
   double VDWRadius(int);
   void GetCOM(double*);
   double PeriodicTable(int);
   void FinDif1stDeriv(double*,double*);
   void ComputeCharge(double*,double*,bool);
   bool MultiGridControl(int);
   double ConjGrad(int,int);
   void RestrictError(double*,double*,double*);
   void InterpolError(double*,double*);
   void ReferenceJob();
   void ReferenceRhoPhi();
   bool UpdateRhoIter();
   void MarcusRefJob();
   void MarcusNeqJob();
   void MarcusNeqRhoPhi();
   bool MarcusNeqRhoIter();
   void PekarRefJob();
   void PekarNeqJob();
   void PekarNeqRhoPhi();
   bool PekarNeqRhoIter();
   double SolvEnergy(double);
   void ReadRhoPhi();
   void WriteRhoPhi();

//CS: add some functions for Poisson-Boltzmann
   void PBJob();
   bool UpdateRhoIons(Electrolyte);
   void ComputeLambda();
   double SolvEnergyPB(Electrolyte el);
};

#endif /*_PEQS_CLASS_H */
