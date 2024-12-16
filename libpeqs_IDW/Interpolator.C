#include "peqs_class.h"
#include "qchem.h"
#include "rem_values.h"
#include "BasisSet.hh"
#include "OneEMtrx.hh"
#include "hirshfeld.h"
#include "$QCHEM_path/anlman/CM5.h"
#include "$QCHEM_path/anlman/makemesh.h"
#include "$QCHEM_path/libdftn/xcclass.h"
#include "$QCHEM_path/libgen/evlbasis.h"
#include "$QCHEM_path/liblas/liblas.h"
#include "ext_libs/fftw/include/fftw.h"
#include "ext_libs/fftw/include/rfftw.h"
#include <Eigen/Core>
#include <unsupported/Eigen/Splines>
#include <Eigen/Dense>
#include <unsupported/Eigen/Splines>
#include "$QCHEM_path/thirdparty/eigen-3.3.4/unsupported/Eigen/src/Splines/SplineFwd.h"
#include "$QCHEM_path/thirdparty/eigen-3.3.4/unsupported/Eigen/src/Splines/Spline.h"
#include "$QCHEM_path/thirdparty/eigen-3.3.4/unsupported/Eigen/src/Splines/SplineFitting.h"

extern "C" void cnvtmm(double*,int*);

void PEQSclass::Make_Atom_Centered_Grid()
{ 
  cout << "Creating atomic grid points " << endl;
  INTEGER NAtoms = rem_read(REM_NATOMS);
    BasisSet BasisDen(bSetMgr.crntCode());
    BasisSet BasisOrb(bSetMgr.crntCode());
    int bCodeDen = BasisDen.code();
    int grdTyp   = rem_read(REM_IGRDTY);
    int IPrint   = rem_read(REM_PRINT_DFT_GRID);
    int IBCode   = bCodeDen;
    SET_QC_FCOMMON(&IBCode);
    XCFunctional xcFunc = XCFunctional(XCFUNC_SCF);     // read from input
    
    XCAtoms    *xcatom;
    XCJobPara  *xcpara;
    XCBasisSet *basDen;
    MoleGrid   *mgrid;
    
    xcatom = new XCAtoms; 
    xcpara = new XCJobPara(xcFunc, *xcatom, 0);
    xcpara->nDen = 1;
    int nAtoms = xcatom->getNAtoms();
    double thresh = xcpara->thresh;
    basDen = new XCBasisSet(bCodeDen, nAtoms, thresh);
    xcatom->setSize(*basDen);
    mgrid = new MoleGrid(*xcatom, grdTyp, xcpara->nDrvNuc, thresh);
    
    
    int totalgrid(0),offset(0),offset2(0),offset3(0);
    double *jPts2, *jWts;
    int nBatch   = mgrid->getNBatch();
    int nDrvGrd = mgrid->getDrvGrd();
    double nThresh = mgrid->getThresh();
    int NPts = mgrid->getNPts();

    for (int ibat=0; ibat<nBatch; ibat++) {
      BatchGrid grid(*mgrid,ibat);
      nGrid  = grid.getNGrid();
      jPts2  = grid.getPts();
      jWts  = grid.getWts();

      
      offset = count;
      offset2 = totalgrid+offset;
      offset3 = offset2*3;
	FileMan(FM_WRITE,FILE_PEQS_GRID,FM_DP,nGrid*3,0,FM_CUR,jPts2);
        FileMan(FM_WRITE,FILE_ATM_BATCH_WTS,FM_DP,nGrid,0,FM_CUR,jWts);
	
      totalgrid += nGrid;
      if (ibat == (nBatch - 1)) count+= totalgrid;
    }

    delete xcatom, xcpara, basDen, mgrid;
}

void PEQSclass::Get_Electronic_Density_on_ACG(double *PAv)
{
  int i,j,k,ijk,l,m;
  double sum1=0.0,sum2=0.0,sum3=0.0,sumC1=0.0,sumC2=0.0;
  INTEGER NBasis = bSetMgr.crntShlsStats(STAT_NBASIS);
  INTEGER N2 = NBasis*NBasis;
  double *kPts2 = QAllocDouble(count*3);  
  double *kWts = QAllocDouble(count);   
  FileMan(FM_READ,FILE_PEQS_GRID,FM_DP,count*3,0,FM_BEG,kPts2);  
  FileMan(FM_READ,FILE_ATM_BATCH_WTS,FM_DP,count,0,FM_BEG,kWts); 

  double *jDMA = QAllocDouble(N2);
  double *jDMB = QAllocDouble(N2);
  FileMan(FM_READ,FILE_DENSITY_MATRIX,FM_DP,N2,0,FM_BEG,jDMA);
  FileMan(FM_READ,FILE_DENSITY_MATRIX,FM_DP,N2,N2,FM_BEG,jDMB);
  double *jChiFree, *aChi;
  aChi = QAllocDoubleWithInit(NBasis*count);
  EvlBasis(aChi,NBasis,kPts2,count);
   double *RhoA = QAllocDoubleWithInit(count);
   double *RhoB = QAllocDoubleWithInit(count);
   double *RhoACG = QAllocDoubleWithInit(count);
   double *unwtRhoACG = QAllocDoubleWithInit(count);
   VRload(RhoACG,count,0.0);
   VRload(unwtRhoACG,count,0.0);

   for(j = 0; j < count; j++){
      for (l = 0; l < NBasis; l++) {
        for (m = 0; m < NBasis; m++) {
        RhoA[j] += jDMA[l*NBasis+m]*aChi[l*count+j]*aChi[m*count+j];
        RhoB[j] += jDMB[l*NBasis+m]*aChi[l*count+j]*aChi[m*count+j];
        }
      }
   }
   for(j = 0; j < count; j++){
        sum2 += kWts[j]*RhoA[j];   
        sum2 += kWts[j]*RhoB[j];   
        RhoACG[j] = kWts[j]*(RhoA[j] + RhoB[j]);
	unwtRhoACG[j] = (RhoA[j] + RhoB[j]);
   }

//========= Electronic Potential calculations ==================================//
  qtime_t shackle08;
  double shackle09[3];
  shackle08 = QTimerOn();
  int cycles=MaxPPer;
  int acgPPer=cycles, offset=0;
  double *irregPts=QAllocDouble(3*cycles),*acg_phi=QAllocDouble(cycles);
  phi_irr=QAllocDouble(count);
  VRload(phi_irr,count,0.0);
  for(l=0;l<count;l+=cycles){
  if (cycles <= (count-l)) acgPPer = cycles;
  else acgPPer = (count-l);
    VRload(acg_phi,acgPPer,0.0);
    VRcopy(irregPts,&kPts2[(3*l)],3*acgPPer);   
    AOints(acg_phi,NULL,NULL,NULL,PAv,irregPts,NULL,NULL,&acgPPer,12);
    VRcopy(&phi_irr[l],acg_phi,acgPPer);
  }
  cout << "Potential have been calculated" << endl;

   QTimerOff(shackle09,shackle08);
   printf(" Potential calculation in ACG points completed in %2.2f (cpu) %2.2f (wall) seconds. \n",shackle09[0],shackle09[2]);

   QFree(jDMA), QFree(jDMB);
   QFree(RhoA), QFree(RhoB), QFree(RhoACG);
   QFree(aChi); 
   QFree(irregPts); QFree(acg_phi);
   QFree(kPts2); QFree(kWts);
   QFree(unwtRhoACG);
}

void PEQSclass::Inverse_Distance_Weighted_interpolation()
{
   double *xcg,*ycg,*zcg;
   double *hPts2=QAllocDouble(count*3);
   double *norm=QAllocDouble(NTotPts);
   double *cartval=QAllocDouble(NTotPts);

   occ=QAllocINTEGER(NTotPts);
   valACG = QAllocDouble(count); 
   phiACG = QAllocDouble(count); 
   int ijk,ibox,jbox,kbox,i,j,k,l,m,n=0,sum_occ=0;
   int Nx=NPoints[0],Ny=NPoints[1],Nz=NPoints[2];
   double gdx=(1.0/delX[0]);
   double gdy=(1.0/delX[1]);
   double gdz=(1.0/delX[2]);

   xacg=QAllocDouble(count);
   yacg=QAllocDouble(count);
   zacg=QAllocDouble(count);
   xcg=QAllocDouble(Nx); 
   ycg=QAllocDouble(Ny); 
   zcg=QAllocDouble(Nz);
   qtime_t shackle06;
   double shackle07[3];
   shackle06 = QTimerOn();

   INTEGER *idx_box=QAllocINTEGER(count*8);
   INTEGER *coeff_box=QAllocINTEGER(count*8);
   VRload(idx_box,count*8,0);
   VRload(coeff_box,count*8,0);
   VRload(valACG,count,0.0);
   VRload(phiACG,count,0.0);
   VRload(occ,NTotPts,0);
   VRload(valCart,NTotPts,0.0);
   VRload(norm,NTotPts,0.0);
   VRload(cartval,NTotPts,0.0);
   FileMan(FM_READ,FILE_PEQS_GRID,FM_DP,count*3,0,FM_BEG,hPts2);
   FileMan(FM_READ,FILE_PEQS_GRID,FM_DP,count,count*3+count+count*3+count,FM_BEG,phiACG);

   for (i=0,j=0,k=0; i<Nx,j<Ny,k<Nz; i++,j++,k++){
     xcg[i] = XMin[0] + i*delX[0];
     ycg[j] = XMin[1] + j*delX[1];
     zcg[k] = XMin[2] + k*delX[2];
   }
   for(l=0;l<count;l++){
     xacg[l]=hPts2[3*l];
     yacg[l]=hPts2[3*l+1];
     zacg[l]=hPts2[3*l+2];
   }

   int ind1,ind2,ind3,ind4,ind5,ind6,ind7,ind8;
   cout << "power = " << power << endl;
   for(l=0;l<count;l++){

     ibox= (int) floor((xacg[l]-XMin[0])*gdx);
     jbox= (int) floor((yacg[l]-XMin[1])*gdy);
     kbox= (int) floor((zacg[l]-XMin[2])*gdz);
     ijk = ibox*Ny*Nz+jbox*Nz+kbox;
     occ[ijk]++;

     ind1=ijk;
     ind2=ijk+1;
     ind3=ijk+Nz;
     ind4=ijk+Nz+1;
     ind5=ijk+Ny*Nz;
     ind6=ijk+Ny*Nz+1;
     ind7=ijk+Ny*Nz+Nz;
     ind8=ijk+Ny*Nz+Nz+1;

     coeff_box[l*8+0]=get_coeff(xcg[ibox],ycg[jbox],zcg[kbox],xacg[l],yacg[l],zacg[l]);
     coeff_box[l*8+1]=get_coeff(xcg[ibox],ycg[jbox],zcg[kbox+1],xacg[l],yacg[l],zacg[l]);
     coeff_box[l*8+2]=get_coeff(xcg[ibox],ycg[jbox+1],zcg[kbox],xacg[l],yacg[l],zacg[l]);
     coeff_box[l*8+3]=get_coeff(xcg[ibox],ycg[jbox+1],zcg[kbox+1],xacg[l],yacg[l],zacg[l]);
     coeff_box[l*8+4]=get_coeff(xcg[ibox+1],ycg[jbox],zcg[kbox],xacg[l],yacg[l],zacg[l]);
     coeff_box[l*8+5]=get_coeff(xcg[ibox+1],ycg[jbox],zcg[kbox+1],xacg[l],yacg[l],zacg[l]);
     coeff_box[l*8+6]=get_coeff(xcg[ibox+1],ycg[jbox+1],zcg[kbox],xacg[l],yacg[l],zacg[l]);
     coeff_box[l*8+7]=get_coeff(xcg[ibox+1],ycg[jbox+1],zcg[kbox+1],xacg[l],yacg[l],zacg[l]);

     idx_box[l*8+0]=ind1;
     idx_box[l*8+1]=ind2;
     idx_box[l*8+2]=ind3;
     idx_box[l*8+3]=ind4;
     idx_box[l*8+4]=ind5;
     idx_box[l*8+5]=ind6;
     idx_box[l*8+6]=ind7;
     idx_box[l*8+7]=ind8;
   }

   cout << "acg points have been printed" << endl;

   for(ijk=0;ijk<NTotPts;ijk++) sum_occ+=occ[ijk];

   for(l=0;l<count;l++){
     for(m=0;m<8;m++){
       norm[idx_box[l*8+m]]+=coeff_box[l*8+m];
       cartval[idx_box[l*8+m]]+=(coeff_box[l*8+m]*phiACG[l]);
     }
   }

   for(ijk=0;ijk<NTotPts;ijk++){
     if(norm[ijk]!=0.0) valCart[ijk]=cartval[ijk]/norm[ijk];
     else valCart[ijk]=0.0;
   }
   double intCart=0.0,sumCart=0.0;
   intCart=Integrator(valCart);
   for(ijk=0;ijk<NTotPts;ijk++) sumCart+=valCart[ijk];

   FinDif1stDeriv(valCart,delCart); 

   QTimerOff(shackle07,shackle06);
   printf(" Sorting and interpolation completed in %2.2f (cpu) %2.2f (wall) seconds. \n",shackle07[0],shackle07[2]);

   QFree(xcg), QFree(ycg), QFree(zcg);
   QFree(idx_box), QFree(coeff_box);
   QFree(hPts2), QFree(norm), QFree(cartval);
}

double PEQSclass::get_coeff(double x,double y,double z,double a,double b,double c)
{
   double coeff,distance=0.0,tol=0.0000001;
   double side=delX[0];
   double fact=0.0;
   distance = sqrt(((x-a)*(x-a))+((y-b)*(y-b))+((z-c)*(z-c)));
   if(distance>tol) coeff=(1.0/pow(distance,power));
   else coeff=(1.0/pow(tol,power));
   return coeff;
}

void PEQSclass::Expanded_IDW()
{
   double *xcg,*ycg,*zcg;
   double *hPts2=QAllocDouble(count*3);
   double *norm=QAllocDouble(NTotPts);
   double *cartval=QAllocDouble(NTotPts);  

   int incount,outcount,k1=0,k2=0;
   int counter;
   int neighDim = outhouse*outhouse*outhouse*8;

   occ=QAllocINTEGER(NTotPts);
   valACG = QAllocDouble(count); 
   phiACG = QAllocDouble(count); 
   int ijk,ibox,jbox,kbox,i,j,k,l,m,n=0,sum_occ=0;
   int Nx=NPoints[0],Ny=NPoints[1],Nz=NPoints[2];
   double gdx=(1.0/delX[0]);
   double gdy=(1.0/delX[1]);
   double gdz=(1.0/delX[2]);

   xacg=QAllocDouble(count);
   yacg=QAllocDouble(count);
   zacg=QAllocDouble(count);
   xcg=QAllocDouble(Nx);
   ycg=QAllocDouble(Ny);
   zcg=QAllocDouble(Nz);
   qtime_t shackle06;
   double shackle07[3];
   shackle06 = QTimerOn();

   INTEGER *idx_box=QAllocINTEGER(count*neighDim);
   INTEGER *coeff_box=QAllocINTEGER(count*neighDim);
   INTEGER *ind=QAllocINTEGER(neighDim);

   VRload(idx_box,count*neighDim,0);
   VRload(coeff_box,count*neighDim,0);
   VRload(ind,neighDim,0);
   VRload(xcg,Nx,0.0);
   VRload(ycg,Ny,0.0);
   VRload(zcg,Nz,0.0);
   VRload(xacg,count,0.0);
   VRload(yacg,count,0.0);
   VRload(zacg,count,0.0);
   VRload(valACG,count,0.0);
   VRload(phiACG,count,0.0);
   VRload(occ,NTotPts,0);
   VRload(valCart,NTotPts,0.0);
   VRload(norm,NTotPts,0.0);
   VRload(cartval,NTotPts,0.0);


   FileMan(FM_READ,FILE_PEQS_GRID,FM_DP,count*3,0,FM_BEG,hPts2);  // ACG points

   for (i=0,j=0,k=0; i<Nx,j<Ny,k<Nz; i++,j++,k++){
     xcg[i] = XMin[0] + i*delX[0];
     ycg[j] = XMin[1] + j*delX[1];
     zcg[k] = XMin[2] + k*delX[2];
   }
   cout << "acg pts to be printed." << endl;
   for(l=0;l<count;l++){
     xacg[l]=hPts2[3*l];
     yacg[l]=hPts2[3*l+1];
     zacg[l]=hPts2[3*l+2];
   cout << xacg[l] << '\t' << yacg[l] << '\t' << zacg[l] << '\t' <<  phi_irr[l] << endl;
   }


   double xmin=Carts[0],ymin=Carts[1],zmin=Carts[2];
   double xmax=xmin, ymax=ymin, zmax=zmin;
   double xcar,ycar,zcar;
   for(l=0;l<NAtoms;l++){
     xcar=Carts[3*l+0];
     ycar=Carts[3*l+1];
     zcar=Carts[3*l+2];
     xmin=min(xcar,xmin);
     ymin=min(ycar,ymin);
     zmin=min(zcar,zmin);
     xmax=max(xcar,xmax);
     ymax=max(ycar,ymax);
     zmax=max(zcar,zmax);     
   }
   double lowX=xmin-shift1,lowY=ymin-shift1,lowZ=zmin-shift1;
   double highX=xmax+shift1,highY=ymax+shift1,highZ=zmax+shift1;
   cout << "LowX = " << lowX/ang2bohr << '\t' << "HighX = " << highX/ang2bohr << endl;
   cout << "LowY = " << lowY/ang2bohr << '\t' << "HighY = " << highY/ang2bohr << endl;
   cout << "LowZ = " << lowZ/ang2bohr << '\t' << "HighZ = " << highZ/ang2bohr << endl;
//=================================================================//
   
   double axmin=xacg[0],aymin=yacg[0],azmin=zacg[0];
   double axmax=axmin, aymax=aymin, azmax=azmin;
   for(l=1;l<count;l++){
     axmin=min(xacg[l],axmin);
     aymin=min(yacg[l],aymin);
     azmin=min(zacg[l],azmin);
     axmax=max(xacg[l],axmax);
     aymax=max(yacg[l],aymax);
     azmax=max(zacg[l],azmax);
   }

   cout << "LowXACG = " << axmin/ang2bohr << '\t' << "HighXACG = " << axmax/ang2bohr << endl;
   cout << "LowYACG = " << aymin/ang2bohr << '\t' << "HighYACG = " << aymax/ang2bohr << endl;
   cout << "LowZACG = " << azmin/ang2bohr << '\t' << "HighZACG = " << azmax/ang2bohr << endl;

//====== assorting ACG pts in to inBatch and outBatch============//
   InSide=QAllocINTEGER(count);
   OutSide=QAllocINTEGER(count);
   VRload(InSide,count,0);
   VRload(OutSide,count,0);
   for(l=0;l<count;l++){
     ibox= (int) floor((xacg[l]-XMin[0])*gdx);
     jbox= (int) floor((yacg[l]-XMin[1])*gdy);
     kbox= (int) floor((zacg[l]-XMin[2])*gdz);
     ijk = ibox*Ny*Nz+jbox*Nz+kbox;
     
     occ[ijk]++;

     if((lowX <= XYZ[3*ijk] && XYZ[3*ijk] <= highX) &&
        (lowY <= XYZ[3*ijk+1] && XYZ[3*ijk+1] <= highY) &&
        (lowZ <= XYZ[3*ijk+2] && XYZ[3*ijk+2] <= highZ)){
	
	InSide[k1]=l;
        k1++;
     }
     else{
        OutSide[k2]=l;
        k2++;
     }
   incount=k1;
   outcount=k2;
   }


//=================================================================//
   int q1=0;
   for(l=0;l<incount;l++){
     ibox= (int) floor((xacg[InSide[l]]-XMin[0])*gdx);
     jbox= (int) floor((yacg[InSide[l]]-XMin[1])*gdy);
     kbox= (int) floor((zacg[InSide[l]]-XMin[2])*gdz);
     ijk = ibox*Ny*Nz+jbox*Nz+kbox;
     
     counter=0;
// loops in x,y,z should run from (1-n) to n, where n is no. of neigh boxes considered
     for(i=(1-inhouse);i<=inhouse;i++){
       for(j=(1-inhouse);j<=inhouse;j++){
	 for(k=(1-inhouse);k<=inhouse;k++){
	   ind[counter]=(ibox+i)*Ny*Nz+(jbox+j)*Nz+(kbox+k);

	   idx_box[InSide[l]*neighDim+counter]=ind[counter];
	   coeff_box[InSide[l]*neighDim+counter]=get_coeff(xcg[ibox+i],ycg[jbox+j],zcg[kbox+k],xacg[InSide[l]],yacg[InSide[l]],zacg[InSide[l]]);
	   counter++;
	 }
       }
     }
   }
   cout << "InBatch interpolation calculation done." << endl;
   int q2=0;
   for(l=0;l<outcount;l++){
     ibox= (int) floor((xacg[OutSide[l]]-XMin[0])*gdx);
     jbox= (int) floor((yacg[OutSide[l]]-XMin[1])*gdy);
     kbox= (int) floor((zacg[OutSide[l]]-XMin[2])*gdz);
     ijk = ibox*Ny*Nz+jbox*Nz+kbox;

     counter=0;
// loops in x,y,z should run from (1-n) to n, where n is no. of neigh boxes considered
     for(i=(1-outhouse);i<=outhouse;i++){
       for(j=(1-outhouse);j<=outhouse;j++){
         for(k=(1-outhouse);k<=outhouse;k++){
           ind[counter]=(ibox+i)*Ny*Nz+(jbox+j)*Nz+(kbox+k);

           idx_box[OutSide[l]*neighDim+counter]=ind[counter];
           coeff_box[OutSide[l]*neighDim+counter]=get_coeff(xcg[ibox+i],ycg[jbox+j],zcg[kbox+k],xacg[OutSide[l]],yacg[OutSide[l]],zacg[OutSide[l]]);
           counter++;
         }
       }
     }
   }
   cout << "OutBatch interpolation calculation done." << endl;
   
   for(l=0;l<count;l++){
     for(m=0;m<neighDim;m++){
       norm[idx_box[l*neighDim+m]]+=coeff_box[l*neighDim+m];
       cartval[idx_box[l*neighDim+m]]+=(coeff_box[l*neighDim+m]*phi_irr[l]);  // interpolatting epot
     }
   }
   cout << "Norm and cartval arrays done." << endl;

   for(ijk=0;ijk<NTotPts;ijk++){
     if(norm[ijk]!=0.0){ 
	valCart[ijk]=cartval[ijk]/norm[ijk];
     }
     else valCart[ijk]=0.0;	
   }
   cout << "valcart generation done." << endl;

   QTimerOff(shackle07,shackle06);
   printf(" Sorting and interpolation completed in %2.2f (cpu) %2.2f (wall) seconds. \n",shackle07[0],shackle07[2]);
 
   QFree(xcg); QFree(ycg); QFree(zcg);
   QFree(idx_box); QFree(coeff_box); QFree(ind);
   QFree(hPts2); QFree(norm); QFree(cartval);
}

//============ genrating a gaussian filter=========================//

void PEQSclass::Do_Filter(double *efn, int l)
{
   double *origin=QAllocDouble(3);
   GetCOM(origin);
   int ijk,i,j,k;
   double xvec,yvec,zvec;

   cout << "l = " << l << endl;
   double cnst=1/sqrt(2*Pi*filter*filter);
   VRload(gauss,NTotPts,0.0); 
   cout << "Filter Width = " << filter << endl;
   for(ijk=0; ijk<NTotPts; ijk++){
    get_3ind_from_multInd(i,j,k,ijk);
    xvec = valueOfX(i);
    yvec = valueOfY(j);
    zvec = valueOfZ(k);
    gauss[ijk] = cnst*cnst*cnst*exp(-(xvec*xvec+yvec*yvec+zvec*zvec)/(2*filter*filter));
   }
   QFree(origin);

//=================================================================//
//==========convolution========================================
   int nx=NPoints[0], ny=NPoints[1], nz=NPoints[2];
   fftw_real cr[NTotPts],icr[NTotPts];
   fftw_real zcr[NTotPts], yzcr[NTotPts],xyzcr[NTotPts];
   fftw_complex b1[NTotPts],b2[NTotPts],cc[NTotPts];

   int sum;
   double scale = 1.0/(NTotPts);

   rfftwnd_plan p, pinv;
   p = rfftw3d_create_plan(nx, ny, nz, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
   pinv = rfftw3d_create_plan(nx, ny, nz, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);

   rfftwnd_one_real_to_complex(p, efn, b1);  //target function
   rfftwnd_one_real_to_complex(p, gauss, b2);   // kernel function

   for(i=0; i<NTotPts; i++){
     cc[i].re = (b1[i].re * b2[i].re - b1[i].im * b2[i].im) * scale ;
     cc[i].im = (b1[i].re * b2[i].im + b1[i].im * b2[i].re) * scale ;
   }
  
  rfftwnd_one_complex_to_real(pinv, cc, cr);

  rfftwnd_destroy_plan(p);
  rfftwnd_destroy_plan(pinv);


  for (i=0;i<NPoints[0];i++){
  for (j=0;j<NPoints[1];j++){
  for (k=0;k<NPoints[2];k++){
  double tmp;
  int s=k;
  int d=(k+NPoints[2]/2)%NPoints[2];
  int src = i*NPoints[1]*NPoints[2] + j*NPoints[2] + s;
  int dst = i*NPoints[1]*NPoints[2] + j*NPoints[2] + d;
  tmp=cr[src];
  zcr[src]=cr[dst];
  zcr[dst]=tmp;
  }
  }
  }

  for (i=0;i<NPoints[0];i++){
  for (j=0;j<NPoints[1];j++){
  for (k=0;k<NPoints[2];k++){
  double tmp;
  int s=j;
  int d=(j+NPoints[1]/2)%NPoints[1];
  int src = i*NPoints[1]*NPoints[2] + s*NPoints[2] + k;
  int dst = i*NPoints[1]*NPoints[2] + d*NPoints[2] + k;
  tmp=zcr[src];
  yzcr[src]=zcr[dst];
  yzcr[dst]=tmp;
  }
  }
  }

  for (i=0;i<NPoints[0];i++){
  for (j=0;j<NPoints[1];j++){
  for (k=0;k<NPoints[2];k++){
  double tmp;
  int s=i;
  int d=(i+NPoints[0]/2)%NPoints[0];
  int src = s*NPoints[1]*NPoints[2] + j*NPoints[2] + k;
  int dst = d*NPoints[1]*NPoints[2] + j*NPoints[2] + k;
  tmp=yzcr[src];
  xyzcr[src]=yzcr[dst];
  xyzcr[dst]=tmp;
  }
  }
  }

  for(ijk=0;ijk<NTotPts;ijk++){
    valGf[3*ijk+l]=xyzcr[ijk]*delX[0]*delX[1]*delX[2]; 
  }
}


double PEQSclass::GetFact(int n)
{
   int fact=1;
   for(int i=1;i<=n;i++){
     fact *= i;
   }
   return fact;
}

double PEQSclass::NPR(int n, int l)
{
   double res=GetFact(l)/GetFact(l-n);
   return res;
}

double PEQSclass::GetMultipolePotential(double *coord, double *mom)
{
   double x=coord[0],y=coord[1],z=coord[2];
   double r2=x*x+y*y+z*z;
   double r=sqrt(r2);
   double denom=1.0/r;
   double pot=mom[0]*denom;
   int lx,ly,lz,l,k;

   for(l=1;l<=MulOrder;l++){
     k=l*(l+1)*(l+2)/6;
     denom /= r2;
     denom /= GetFact(l);
     for(lz=0;lz<=l;lz++){
       for(ly=0;ly<=(l-lz);ly++){
	 lx=l-ly-lz;
         if((lx!=0 && ly==0 && lz==0) || (lx==0 && ly!=0 && lz==0) || (lx==0 && ly==0 && lz!=0)){
           pot += mom[k]*pow(x,double(lx))*pow(y,double(ly))*pow(z,double(lz))*denom;	 
	   k++;
	 }
         else if((lx!=0 && ly!=0 && lz==0) || (lx==0 && ly!=0 && lz!=0) || (lx!=0 && ly==0 && lz!=0)){
           pot += NPR(2,l)*mom[k]*pow(x,double(lx))*pow(y,double(ly))*pow(z,double(lz))*denom;
           k++;
         }
	 else{
	   pot += mom[k]*pow(x,double(lx))*pow(y,double(ly))*pow(z,double(lz))*denom;
           k++;
         }
       }
     }
   }
   return pot;
}

void PEQSclass::Combine_IDW_and_MulPot(double *PAv)
{
   int i,j,k=0,l=0,ijk,ibcount,obcount;
   cout << "Inside CleanPhi() line 768" << endl;
   INTEGER *InnerBatch=QAllocINTEGER(NTotPts);
   INTEGER *OuterBatch=QAllocINTEGER(NTotPts);

   VRload(InnerBatch,NTotPts,0);
   VRload(OuterBatch,NTotPts,0);
   VRload(val3,NTotPts,0.0);
   cout << "Inside cleanup routine" << endl;

// Creating fake-envelop and  storing indices
   double xmin=Carts[0],ymin=Carts[1],zmin=Carts[2];
   double xmax=xmin, ymax=ymin, zmax=zmin;
   double xcar,ycar,zcar;
   for(l=0;l<NAtoms;l++){
     xcar=Carts[3*l+0];
     ycar=Carts[3*l+1];
     zcar=Carts[3*l+2];
     xmin=min(xcar,xmin);
     ymin=min(ycar,ymin);
     zmin=min(zcar,zmin);
     xmax=max(xcar,xmax);
     ymax=max(ycar,ymax);
     zmax=max(zcar,zmax);
   }
   double lowX=xmin-shift3,lowY=ymin-shift3,lowZ=zmin-shift3;
   double highX=xmax+shift3,highY=ymax+shift3,highZ=zmax+shift3;
   cout << "LowX = " << lowX/ang2bohr << '\t' << "HighX = " << highX/ang2bohr << endl;
   cout << "LowY = " << lowY/ang2bohr << '\t' << "HighY = " << highY/ang2bohr << endl;
   cout << "LowZ = " << lowZ/ang2bohr << '\t' << "HighZ = " << highZ/ang2bohr << endl;
   
   for(ijk=0;ijk<NTotPts;ijk++){
     if((lowX <= XYZ[3*ijk] && XYZ[3*ijk] <= highX) && 
        (lowY <= XYZ[3*ijk+1] && XYZ[3*ijk+1] <= highY) && 
        (lowZ <= XYZ[3*ijk+2] && XYZ[3*ijk+2] <= highZ)){
        InnerBatch[k]=ijk;
	k++;
     }
     else{
	OuterBatch[l]=ijk;
	l++;
     }
   }
   ibcount=k;
   obcount=l;

   cout << "Sorting into IB and OB done" << endl;
   cout << "Inside envelop = " << ibcount << endl;
   cout << "Outside envelop = " << obcount << endl; 

// Copying interpolated-phi-values inside the envelop
   
   for(i=0;i<ibcount;i++){
      val3[InnerBatch[i]]=valCart[InnerBatch[i]];
   }
   
//======== MOMENT CALCULATION ========================================//
   cout << "Will start moments calculation" << endl;
   int Unrestricted = rem_read(REM_JUSTAL);
   int NB2    = rem_read(REM_NB2);
   int NB2car = rem_read(REM_NB2CAR);
   int NBas   = bSetMgr.crntShlsStats(STAT_NBASIS);
   int NBas2  = NBas*NBas;
   int NDen = (Unrestricted) ? 2 : 1;
   int NMoments=LFuncC(0,MulOrder);
   double *moments=QAllocDouble(NMoments);
   VRload(moments,NMoments,0.0);
   double *Pv=QAllocDouble(NB2car*NDen);
   double *denMat=QAllocDouble(NBas2*NDen);

   cout << "NB2 = " << NB2 << endl;
   cout << "NB2car = " << NB2car << endl;
   cout << "NBas = " << NBas << endl;
   cout << "NDen = " << NDen << endl;
   cout << "MulOrder = " << MulOrder << endl;
   cout << "Unrestricted = " << Unrestricted << endl;
   cout << "NMoments = " << NMoments << endl;

   FileMan(FM_READ,FILE_DENSITY_MATRIX,FM_DP,NBas2*NDen,0,FM_BEG,denMat);
   VRload(Pv,NB2car*NDen,0.0);
   ScaV2M(denMat,Pv,1,0);
   if(Unrestricted){
      ScaV2M(denMat+NBas2,Pv+NB2car,1,0);
      VRadd2(Pv,Pv+NB2car,NB2);
   }
   else VRscale(Pv,NB2,2.0);
   MultipoleMoments(moments,Pv,1,MulOrder,0); // 0 means not including nuclear moments
   cout << "Finished multipole moments calculation" << endl;

//======== Automated Multipole Potential ========================//
   double *rX=QAllocDouble(3);
   double *com=QAllocDouble(3);
   double *dist=QAllocDouble(3);
   GetCOM(com);
   cout << "COM is = " << com[0] << '\t' << com[1] << '\t' << com[2] << endl;
   for(int jj=0;jj<obcount;jj++){
     VRcopy(rX,&XYZ[3*(OuterBatch[jj])],3);
     VRsub(dist,rX,com,3);     
     val3[OuterBatch[jj]] = GetMulPot(dist,moments);      // this is the actual line
   }

   QFree(rX);
   QFree(com);
//===============================================================//

   QFree(Pv), QFree(denMat), QFree(moments);
   cout << "outerBatch potential generation done" << endl; 
   QFree(InnerBatch), QFree(OuterBatch);   
}

double PEQSclass::make_gaussian(double x,double y,double z)
{
   double value1,value2,value3;

   value1=exp(-(x*x+y*y+z*z));
   value2=(1.0/(1+ (x*x+ y*y + z*z)));
   value3=(x+y+z)*(exp(-(x*x+y*y+z*z)));

   return value1;
}

