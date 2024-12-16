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
#include <Eigen/Core>
#include <unsupported/Eigen/Splines>

#define INVALID_INT  -99999
#define INVALID_REAL -99999.9

class PEQSclass {
private:

   enum PEQS_CavType {
      PEQS_CAVTYPE_NONE  =0,
      PEQS_CAVTYPE_VDW   =1,
      PEQS_CAVTYPE_SPHERE=2,
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
   double filter,sigma,solver_thresh,cg_thresh;
   double *XMin,*XMax,*delX,*XYZ,*Carts,*Center;
   double *vdWScale,*epsilon,*log_epsilon,*delLogEps,*delPhi,*delPhiFast,*delRhoElec;
   double *delPhiN, *delPhiE;
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


//===== Interpolation ============================//
   
   double *delCart; // ACG
   double *val_elec; // Do_IDW2
   double *mono, *dip, *quad, *octa, *hexa; //ACG
   double *fifth, *sixth, *seventh, *eight; //ACG
   double *gauss; // ACG
   double *Field, *rho_from_field; // electric field
   INTEGER *test; // Eigen Spline 
   double *RhoCart; // ACG
   double *valCart, *valACG, *phi_irr, *phiACG, *ACGField; // ACG
   double *val2, *val3, *delValCart, *valEField; //ACG
   double power, shift1, shift2, shift3;  // ACG
   double *xacg,*yacg,*zacg; // ACG
   double *val3f, *delVal, *valSf; // ACG
   INTEGER *occ;
   INTEGER *InSide;
   INTEGER *OutSide;
   INTEGER *MidZone;
   int cutoff, interp, count, MulOrder; //ACG
   int nGrid, potgrid, GrdPat; //ACG
   int inhouse, outhouse; // Do_IDW2
   double *efx, *efy, *efz, *efn, *valGf;
   int order,max_der,delta;
   int orderp1,halfordp1,orderp1t2,norder,nhalf,nth;
   double gd,gd2;

   double **table;
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

//=====================================================//


   int NZ()            const{return NPoints[2];}
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
   double hypertan(double x, double a, double b) const{return 1+tanh((x-a)/b);}

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
     double *origin=QAllocDouble(3);
     GetCOM(origin);
     int i,j,k,ijk,Nx=NPoints[0],Ny=NPoints[1],Nz=NPoints[2];
     ofstream dataFileZ("PEQS.data");
     ofstream dataFileA("Debug_x.data");
     ofstream dataFileB("Debug_y.data");
     ofstream dataFileC("Debug_z.data");
     ofstream dataFileG1("phiXY.data");
     ofstream dataFileG2("phiZX.data");
     ofstream dataFileG3("phiYZ.data");

     for (i=0;i<Nx;i++){
       for (j=0;j<Ny;j++){
         for (k=0;k<Nz;k++){
           ijk = get_multInd(i,j,k);
           double x = valueOfX(i);
           double y = valueOfY(j);
           double z = valueOfZ(k);
	   if (z==0.0) dataFileG1 << x/ang2bohr << '\t' << y/ang2bohr << '\t' << 
			phi_elec[ijk] << '\t' << gauss[ijk] << '\t' << valCart[ijk] << '\t' << val2[ijk] << '\t' << val3[ijk] << '\t' << 
			rho_elec[ijk] << '\t' << val_elec[ijk] << '\t' << abs(val2[ijk]-phi_elec[ijk]) << '\t' << abs(rho_elec[ijk]-val_elec[ijk]) << '\t' <<
			Field[3*ijk+2] << '\t' << delPhiE[3*ijk+2] << '\t' << rho_from_field[ijk] << valEField[3*ijk+2] << endl;
 
	   if (y==0.0) dataFileG2 << z/ang2bohr << '\t' << x/ang2bohr << '\t' <<
			phi_elec[ijk] << '\t' << gauss[ijk] << '\t' << valCart[ijk] << '\t' << val2[ijk] << '\t' << val3[ijk] << '\t' <<
                        rho_elec[ijk] << '\t' << val_elec[ijk] << '\t' << abs(val2[ijk]-phi_elec[ijk]) << '\t' << abs(rho_elec[ijk]-val_elec[ijk]) << endl;
 
	   if (x==0.0) dataFileG3 << y/ang2bohr << '\t' << z/ang2bohr << '\t' <<
			phi_elec[ijk] << '\t' << gauss[ijk] << '\t' << valCart[ijk] << '\t' << val2[ijk] << '\t' << val3[ijk] << '\t' <<
                        rho_elec[ijk] << '\t' << val_elec[ijk] << '\t' << abs(val2[ijk]-phi_elec[ijk]) << '\t' << abs(rho_elec[ijk]-val_elec[ijk]) << endl;

           if (z==0.0) dataFileZ << x/ang2bohr << '\t' << y/ang2bohr << '\t' 
               << rho_total[ijk] << '\t' << rho_elec[ijk]  << '\t' << rho_nuke[ijk]   << '\t' << '\t' 
               << rho_polar[ijk] << '\t' << phi_total[ijk] << '\t' << phi_elec[ijk]   << '\t' << abs(phi_elec[ijk]-valCart[ijk]) << '\t' 
               << phi_nuke[ijk]  << '\t' << phi_polar[ijk] << '\t' << delPhi[3*ijk+2] << '\t' 
               << epsilon[ijk]   << '\t' << delLogEps[3*ijk+2] << '\t' << rho_ions[ijk] << '\t' 
               << '\t' << rho_iter[ijk]  << '\t' << lambda[ijk] << '\t' << valCart[ijk] << endl; 
	   if(z==0.0 && y==0.0){
		dataFileA << x/ang2bohr  << '\t' << phi_elec[ijk]/val3[ijk] << '\t' << rho_iter[ijk] << '\t'
			<< valCart[ijk] << '\t' << phi_elec[ijk] << '\t' << val2[ijk] << '\t' << val3[ijk] << '\t' 
			<< rho_elec[ijk] << '\t' << val_elec[ijk] << '\t' << rho_polar[ijk] << '\t' << phi_nuke[ijk] << '\t' << rho_nuke[ijk] << '\t'
			<< delPhi[3*ijk] << '\t' << Field[3*ijk] << '\t' << delPhiE[3*ijk] << '\t' << delPhiN[3*ijk] << '\t'
			<< rho_from_field[ijk] << '\t' << valEField[3*ijk+0] << '\t' << valGf[3*ijk+0] << '\t' << valEField[3*ijk+0]/valGf[3*ijk+0] << '\t' << valSf[3*ijk] << endl;
                }
           if(x==0.0 && z==0.0){
		dataFileB << y/ang2bohr  << '\t' << phi_elec[ijk]/val3[ijk] << '\t' << rho_iter[ijk] << '\t'  
			<< valCart[ijk] << '\t' << phi_elec[ijk] << '\t' << val2[ijk] << '\t' << val3[ijk] <<  '\t' 
			<< rho_elec[ijk] << '\t' << val_elec[ijk] << '\t' << rho_polar[ijk] << '\t' << phi_nuke[ijk] << '\t' << rho_nuke[ijk] << '\t'
			<< delPhi[3*ijk+1] << '\t' << Field[3*ijk+1] << '\t' << delPhiE[3*ijk+1] << '\t' << delPhiN[3*ijk+1] << '\t'
		  	<< rho_from_field[ijk] << '\t' << valEField[3*ijk+1] << '\t' << valGf[3*ijk+1] << '\t' << valEField[3*ijk+1]/valGf[3*ijk+1] << '\t' << valSf[3*ijk+1] << endl;
                }
           if(y==0.0 && x==0.0){
		dataFileC << z/ang2bohr  << '\t' << phi_elec[ijk]/val3[ijk] << '\t' << rho_iter[ijk] << '\t'
			<< valCart[ijk] << '\t' << phi_elec[ijk] << '\t' << val2[ijk] << '\t' << val3[ijk] << '\t' 
			<< rho_elec[ijk] << '\t' << val_elec[ijk] << '\t' << rho_polar[ijk] << '\t' << phi_nuke[ijk] << '\t' << rho_nuke[ijk] << '\t'
			<< delPhi[3*ijk+2] << '\t' << Field[3*ijk+2] << '\t' << delPhiE[3*ijk+2] << '\t' << delPhiN[3*ijk+2] << '\t'
			<< rho_from_field[ijk] << '\t' << valEField[3*ijk+2] << '\t' << valGf[3*ijk+2] << '\t' << valEField[3*ijk+2]/valGf[3*ijk+2] << '\t' << valSf[3*ijk+2] << endl;
                }
         }
       }
     }
     dataFileG1.close();
     dataFileG2.close();
     dataFileG3.close();
     dataFileZ.close();
     dataFileA.close();
     dataFileB.close();
     dataFileC.close();
     QFree(origin);
   }

   void PrintPot(void){
     int i,j,k,ijk,Nx=NPoints[0],Ny=NPoints[1],Nz=NPoints[2];
     ofstream PotX("PotX.dat");
     ofstream PotY("PotY.dat");
     ofstream PotZ("PotZ.dat");

     for (i=0;i<Nx;i++){
       for (j=0;j<Ny;j++){
         for (k=0;k<Nz;k++){
           ijk = get_multInd(i,j,k);
           double x = valueOfX(i);
           double y = valueOfY(j);
           double z = valueOfZ(k);

                if(z==0.0 && y==0.0){
                PotX << x/ang2bohr  << '\t' << rho_elec[ijk] << '\t' << valCart[ijk] << '\t'
                       << rho_polar[ijk] << '\t' << phi_elec[ijk] << '\t' << val3[ijk] << endl;
                }
                if(x==0.0 && z==0.0){
                PotY << y/ang2bohr  << '\t' << rho_elec[ijk] << '\t' << valCart[ijk] << '\t'
                       << rho_polar[ijk] << '\t' << phi_elec[ijk] << '\t' << val3[ijk] << endl;
                }
                if(y==0.0 && x==0.0){
                PotZ << z/ang2bohr  << '\t' << rho_elec[ijk] << '\t' << valCart[ijk] << '\t'
                       << rho_polar[ijk] << '\t' << phi_elec[ijk] << '\t' << val3[ijk] << endl;
                }
          }
       }
     }
     PotX.close();
     PotY.close();
     PotZ.close();
   }

   void PrintPot2D(void){
     int i,j,k,ijk,Nx=NPoints[0],Ny=NPoints[1],Nz=NPoints[2];
     ofstream PotXY("PotXY.dat");
     ofstream PotYZ("PotYZ.dat");
     ofstream PotZX("PotZX.dat");

     for (i=0;i<Nx;i++){
       for (j=0;j<Ny;j++){
         for (k=0;k<Nz;k++){
           ijk = get_multInd(i,j,k);
           double x = valueOfX(i);
           double y = valueOfY(j);
           double z = valueOfZ(k);

                if(z==0.0){
                PotXY << x/ang2bohr  << '\t' << y/ang2bohr << '\t' << rho_elec[ijk] << '\t' << valCart[ijk] << '\t'
                       << rho_polar[ijk] << '\t' << phi_elec[ijk] << '\t' << val3[ijk] << endl;
                }
                if(x==0.0){
                PotYZ << y/ang2bohr  << '\t' << z/ang2bohr << '\t' << rho_elec[ijk] << '\t' << valCart[ijk] << '\t'
                       << rho_polar[ijk] << '\t' << phi_elec[ijk] << '\t' << val3[ijk] << endl;
                }
                if(y==0.0){
                PotZX << z/ang2bohr  << '\t' << x/ang2bohr << '\t' << rho_elec[ijk] << '\t' << valCart[ijk] << '\t'
                       << rho_polar[ijk] << '\t' << phi_elec[ijk] << '\t' << val3[ijk] << endl;
                }
          }
       }
     }
     PotXY.close();
     PotYZ.close();
     PotZX.close();
   }

   void Make_Atom_Centered_Grid(); // ACG
   void Get_Electronic_Density_on_ACG(double*); // ACG
   void Get_Electric_Field_on_ACG(double*); // ACG
   double get_coeff(double,double,double,double,double,double); //ACG
   void Inverse_Distance_Weighted_interpolation(); // ACG
   void Expanded_IDW(); // ACG
   void Do_Filter(double *, int); // ACG
   void Combine_IDW_and_MulPot(double*); //ACG
   double GetFact(int); //ACG
   double NPR(int, int); //ACG
   double GetMulPot(double *, double *); //ACG
   double make_gaussian(double, double, double); //ACG

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
   double SolvEnergy();
   void ReadRhoPhi();
   void WriteRhoPhi();

//CS: add some functions for Poisson-Boltzmann
   void PBJob();
   bool UpdateRhoIons(Electrolyte);
   void ComputeLambda();
   double SolvEnergyPB(Electrolyte el);
   int testMain();

};

// ================== Eigen test =============================================================//

class SplineFunction: public PEQSclass {
public:

  SplineFunction(Eigen::VectorXd const &x_vec,
                 Eigen::VectorXd const &y_vec)
    : x_min(x_vec.minCoeff()),
      x_max(x_vec.maxCoeff()),
  spline_(Eigen::SplineFitting<Eigen::Spline<double, 1>>::Interpolate(
                y_vec.transpose(),
  std::min<int>(x_vec.rows() - 1, 3),
                scaled_values(x_vec)))
  { }

  double operator()(double x) const {
    return spline_(scaled_value(x))(0);
  }

private:

   double scaled_value(double x) const {
    return (x - x_min) / (x_max - x_min);
  }

  Eigen::RowVectorXd scaled_values(Eigen::VectorXd const &x_vec) const {
    return x_vec.unaryExpr([this](double x) { return scaled_value(x); }).transpose();
  }

  double x_min;
  double x_max;

  Eigen::Spline<double, 1> spline_;

};
// ================== Eigen test =============================================================//


#endif /*_PEQS_CLASS_H */
