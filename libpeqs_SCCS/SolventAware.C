#include "peqs_class.h"
#include <armadillo>
#include <algorithm>
#include "ext_libs/fftw/include/fftw.h"
#include "ext_libs/fftw/include/rfftw.h"


void PEQSclass::SolventAwareSurfaceCavity()
{
  int i,j,k,l,ijk;
  int c1=0,c2=0,c3=0;
  double logRmax = log(rho_max);
  double logRmin = log(rho_min);
  double logEsolv = log(eps_solvent);
  double p1,p2,ratio,logRho;
  double *xc=QAllocDouble(NPoints[0]);
  double *yc=QAllocDouble(NPoints[1]);
  double *zc=QAllocDouble(NPoints[2]);
  
  FinDif1stDeriv(rho_elec,delRhoElec);

  cout << "Entered Gygi Cavity" << endl;
  cout << "NTotPts = " << NTotPts << endl;
  cout << "eps_vacuum = " << eps_vacuum << endl;
  cout << "eps_solvent = " << eps_solvent << endl;
  VRload(eps,NTotPts,0.0);


  for(ijk=0;ijk<NTotPts;ijk++){
  if(rho_elec[ijk]>=rho_max) {
	ifunc[ijk] = 1.0;
        c1 = c1+1;
        }
  else if(rho_elec[ijk]<=rho_min) {
	ifunc[ijk] = 0.0;
        c2 = c2+1;
        }

  else if(rho_elec[ijk]>rho_min && rho_elec[ijk]<rho_max){
    c3 = c3+1;
    logRho = log(rho_elec[ijk]);
    ratio = (logRmax-logRho)/(logRmax-logRmin);
    p1 = logEsolv*ratio;
    p2 = (logEsolv/(2*Pi))*sin(2*Pi*ratio);
    ifunc[ijk] = 1.0-((p1-p2)/logEsolv);
    }
  }
//==========Gygi+Rigid==========================================//
   cout << "VDW_shift = " << vdwShift << endl;
   for(ijk=0;ijk<NTotPts;ijk++){
     get_3ind_from_multInd(i,j,k,ijk);
     double x=valueOfX(i),y=valueOfY(j),z=valueOfZ(k);
     double vdw=1.0,h=1.0;
     for(l=0;l<NAtoms;l++){
       double xcar=Carts[3*l+0],ycar=Carts[3*l+1],zcar=Carts[3*l+2];
       double xcar2=xcar*xcar,ycar2=ycar*ycar,zcar2=zcar*zcar;
       int Atom=AtNo[l];
       vdw=vdwScale*ang2bohr*(VDWRadius(Atom)+vdwShift);
       double rX=x-xcar,rY=y-ycar,rZ=z-zcar;
       double dist2=rX*rX+rY*rY+rZ*rZ,dist = sqrt(dist2);
       h *= 0.5*(1.0 + erf((dist-vdw)/erf_delta));
     }
     eps[ijk]=(ifunc[ijk]-1.0)*h + 1.0;
   }
//================================================================//


//==========probe function======================================== 
  double *rX=QAllocDouble(3);
  double *com=QAllocDouble(3);
  double *dist=QAllocDouble(3);
  double *pro=QAllocDouble(NTotPts);
  VRload(pro,NTotPts,0.0);
  VRload(npro,NTotPts,0.0);
  double d,pro_sum=0.0;
  GetCOM(com);

  for(ijk=0;ijk<NTotPts;ijk++){
    VRload(rX,3,0.0);
    VRcopy(rX,&XYZ[3*ijk],3);
    VRsub(dist,rX,com,3);
    d=sqrt(rX[0]*rX[0]+rX[1]*rX[1]+rX[2]*rX[2]);
    pro[ijk]=getProbe(d);
    pro_sum += pro[ijk];
  }

  for(ijk=0;ijk<NTotPts;ijk++){
    npro[ijk] = pro[ijk]/pro_sum;
  }

//==========convolution========================================
  int nx=NPoints[0], ny=NPoints[1], nz=NPoints[2];

  fftw_real cr[NTotPts],icr[NTotPts];
  fftw_real zcr[NTotPts], yzcr[NTotPts],xyzcr[NTotPts];
  fftw_complex b1[NTotPts],b2[NTotPts],cc[NTotPts];
  double scale = 1.0/(NTotPts);

  rfftwnd_plan p, pinv;
  p = rfftw3d_create_plan(nx, ny, nz, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
  pinv = rfftw3d_create_plan(nx, ny, nz, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);

  rfftwnd_one_real_to_complex(p, eps, b1);
  rfftwnd_one_real_to_complex(p, npro, b2); 

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

  VRload(FillFrac,NTotPts,0.0);
  VRload(t_eta,NTotPts,0.0);
  VRload(mod_sr,NTotPts,0.0);
  VRload(tst_sr,NTotPts,0.0);

  for(i=0;i<nx;i++){
   for(j=0;j<ny;j++){
    for(k=0;k<nz;k++){
    ijk = i*ny*nz + j*nz + k;  
     FillFrac[ijk] = xyzcr[ijk];
     t_eta[ijk] = 0.5*(1+erf((FillFrac[ijk]-FillZero)/eta_delta));
     mod_sr[ijk] = eps[ijk]+(1-eps[ijk])*t_eta[ijk];
    }
   }
  }

  for(ijk=0;ijk<NTotPts;ijk++){
    epsilon[ijk]=exp(log(eps_solvent)*(1-mod_sr[ijk]));
  }

  double QiVol=0.0, QfVol=0.0, FillVol=0.0, cell;

  FillVol = Integrator(t_eta);
  QiVol = Integrator(ifunc);
  QfVol = Integrator(mod_sr);
  cout << "Fill volm = " << FillVol << endl;
  cout << "Qi-vol = " << QiVol << endl;
  cout << "Qf-vol = " << QfVol << endl;

  QFree(rX);
  QFree(com);
  QFree(dist);
  QFree(pro);

  PrintEpsilon();
  Print();
}

double PEQSclass::getProbe(double dist)
{
  double probe=0.0;
  probe = 0.5*erfc((dist-(alfa*Rsol))/zhi_delta);
  return probe; 
}

/*
  Evaluate basis function on atom centered grid and
  compute SCCS polarization response for augmenting to
  Gas phase Fock matrix
*/

void PEQSclass::GetSCCS()
{
  int i,j,k,l,ijk;
  double logRmax = log(rho_max);
  double logRmin = log(rho_min);
  double logEsolv = log(eps_solvent);
  double p1,p2,ratio,logRho;
  double dot,OneByRho;
  double xm,ym,zm;
  VRload(phi_sccs,NTotPts,0.0);

  for(ijk=0;ijk<NTotPts;ijk++){
    dot=0.0;
    for (l=0;l<3;l++) dot += delPhi[3*ijk+l];
    if(rho_elec[ijk]>rho_min && rho_elec[ijk]<rho_max){
      logRho = log(rho_elec[ijk]);
      ratio = (logRmax-logRho)/(logRmax-logRmin);
      OneByRho = (1.0/rho_elec[ijk]);
      p1 = (logEsolv/(logRmax-logRmin));
      p2 = cos(2*Pi*ratio)-1.0;
      phi_sccs[ijk] = (-1.0/(8*Pi))*epsilon[ijk]*p1*p2*OneByRho*dot*dot;
    }
  }

  int NGrid;
  Veps=QAllocDouble(NBasis*NBasis);
  VRload(Veps,NBasis*NBasis,0.0);

  qtime_t shackle11,shackle13;
  double shackle12[3];
  double shackle14[3];
  shackle11 = QTimerOn();

  for(ijk=0;ijk<NTotPts;ijk+=MaxPPer){
    if(MaxPPer <= (NTotPts-ijk)) NGrid=MaxPPer;
    else NGrid=(NTotPts-ijk);
    double *Chi=QAllocDouble(NBasis*NGrid);
    double *PDT2=QAllocDouble(NBasis*NBasis); // for liblas
    double *cPts=QAllocDouble(3*NGrid);
    VRload(Chi,NBasis*NGrid,0.0);
    VRload(PDT2,NBasis*NBasis,0.0); // for liblas
    VRcopy(cPts,&XYZ[(3*ijk)],3*NGrid);

    EvlBasis(Chi,NBasis,cPts,NGrid);
    double *Res=QAllocDouble(NBasis*NGrid); // VRlib
    VRload(Res,NBasis*NGrid,0.0); // VRlib

    double *C=QAllocDouble(NGrid);
    double *F=QAllocDouble(NGrid);
    double *P=QAllocDouble(NGrid);

    for (i=0;i<NBasis;i++){
      VRload(P,0.0,NGrid);
      VRmultelmt(P,phi_sccs+ijk,Chi+i*NGrid,NGrid);
      VRcopy(&Res[i*NGrid],P,NGrid);
    }

    AtimsB(PDT2,Res,Chi,NBasis,NBasis,NGrid,NBasis,NBasis,NBasis,3); // liblas
    VRadd(Veps,Veps,PDT2,NBasis*NBasis); // liblas

    QFree(Chi);
    QFree(cPts);
    QFree(Res);
    QFree(C),QFree(F),QFree(P);
    QFree(PDT2);
  }

  VRscale(Veps,NBasis*NBasis,volumeElmnt());
  QTimerOff(shackle12,shackle11);
     printf("\tPEqS:  Computed the Veps-Quad in %2.2f (cpu) %2.2f (wall) sec. \n",
           shackle12[0],shackle12[2]);
  double sum1=0.0,sum2=0.0;
  for(i=0;i<NBasis*NBasis;i++){
    sum1 += Veps[i];
  }
}
