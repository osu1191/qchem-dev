#include "peqs_class.h"
#include <armadillo>
#include <algorithm>
#include "ext_libs/fftw/include/fftw.h"
#include "ext_libs/fftw/include/rfftw.h"

void PEQSclass::ComputeDielectric(bool recompute)
{
  VRload(epsilon,NTotPts,eps_vacuum);
  if (peqs_print > 0)
     printf("\tPEqS:  Solvent diel = %6.3f\n",eps_solvent);
  dielectric = eps_solvent;
  if (noneq_state == 1 || (noneq_state == 0 && recompute)){
    if (peqs_print > 0)
       printf("\tPEqS:  Optical diel = %6.3f\n",eps_infi);
    dielectric = eps_infi;
  }


  if (CavityType == PEQS_CAVTYPE_NONE)
     printf("\tPEqS:  Performing an electrostatic calculation in vacuum.\n");
  else if (CavityType == PEQS_CAVTYPE_VDW)
  {
     if (peqs_print > 1)
     {
        printf("\tPEqS:  Constructing diel function; switching region = %.4f bohr.\n",erf_delta);
        if (vdwType == PEQS_VDW_UNSCALED)
           printf("\tPEqS:  Switching functions use unscaled vdW radii.\n");
        else if (vdwType == PEQS_VDW_SCALED)
           printf("\tPEqS:  Switching functions use vdW radii scaled by %.3f.\n",vdwScale);
        else if (vdwType == PEQS_VDW_SHIFTED || vdwType == PEQS_VDW_HYBRID)
           printf("\tPEqS:  Switching functions use vdW radii scaled by %.3f, shifted by %.3f bohr.\n",
             vdwScale,vdwShift);
        if (vdwType == PEQS_VDW_HYBRID)
           printf("\tPEqS:  Using the hybrid cavity algorithm.\n");
     }
     if (!Interface)
        RigidVDWSurface();
     else
        RigidVDWInterface();
  }
  else if (CavityType == PEQS_CAVTYPE_SPHERE)
  {
    if (FindCenter == 2) SphereCenter(Center);
    else if (FindCenter == 3) GetCOM(Center);
    if (AutoRadius) SphereRadius();
    if (peqs_print > 1)
    {
       printf("\tPEqS:  Using a spherical cavity, R = %.3f bohr.\n",r_cavity);
       printf("\tPEqS:  Cavity center = (%.3f, %.3f, %.3f) bohr.\n",Center[0],Center[1],Center[2]);
       if (tanh_cut1 > 0.0)
          printf("\tPEqS:  Cutoff of %.3f bohr added to cavity radius.\n",tanh_cut1);
       printf("\tPEqS:  Dielectric interpolated over %.3f bohr\n",tanh_cut2);
       printf("\tPEqS:  Smoothing function scale factor = %.3f bohr^(-1)\n",tanh_scale);
    }
    if (!Interface)
       SphericalSurface(Center);
    else
       SphericalInterface(Center);
  }
  else if(CavityType == PEQS_CAVTYPE_GYGI)
  {
    cout << "Using Gygi Cavity based on electronic density" << endl;
    if(!Interface) GygiBulkCavity();
    GetSCCS();
  }
  else if(CavityType == PEQS_CAVTYPE_ARBI){  // changed today
    cout << "User input dielectric is being used for arbitrary cavity" << endl; // changed today
    if (!Interface)
    DiscreetBulkCavity();  // changed today
  }

  for (int ijk=0;ijk<NTotPts;ijk++) 
     log_epsilon[ijk] = log(epsilon[ijk]);
  if (CavityType != PEQS_CAVTYPE_SPHERE || Interface)
     FinDif1stDeriv(log_epsilon,delLogEps);
}

void PEQSclass::DiscreetBulkCavity()
{
  double *pts = QAllocDouble(4*NTotPts);
  FileMan(FM_READ,FILE_PEQS_GRID,FM_DP,4*NTotPts,0,FM_BEG,pts);
  for (int i=3,ijk=0;i<4*NTotPts,ijk<NTotPts;i+=4,ijk++){
//        cout << "pore thik kore nebo" << endl;
        epsilon[ijk]=pts[i];
        }
}

void PEQSclass::RigidVDWSurface()
{
  int i,j,k,ijk,l;
  for (ijk=0;ijk<NTotPts;ijk++){
    get_3ind_from_multInd(i,j,k,ijk);
    double x=valueOfX(i),x2=x*x,y=valueOfY(j),y2=y*y,z=valueOfZ(k),z2=z*z,grid_dist2=(x2+y2+z2);
    double d_0=1.0,h = 1.0;
    for (l=0;l<NAtoms;l++){
      double xcar=Carts[3*l+0],ycar=Carts[3*l+1],zcar=Carts[3*l+2];
      double xcar2=xcar*xcar,ycar2=ycar*ycar,zcar2=zcar*zcar,cart_dist2=xcar2+ycar2+zcar2;
      int AtomicNo = AtNo[l];
      if (vdwType == 1) d_0 = ang2bohr*VDWRadius(AtomicNo);
      else if (vdwType == 2) d_0 = vdWScale[l]*ang2bohr*VDWRadius(AtomicNo);
      else d_0 = vdWScale[l]*ang2bohr*(VDWRadius(AtomicNo)+vdwShift);
      //add probe radius to effectively use solvent accessible surface (SAS)
      if (probe_radius > 0 && vdwType == 1)
         d_0 += probe_radius;
      else if(probe_radius > 0.0 && vdwType != 1)
         std::cout << "Attempt to construct solvent accessible surface not implemented for scaled vdW radii." << std::endl;
      double rX=x-xcar,rY=y-ycar,rZ= z-zcar;
      double dist2=rX*rX+rY*rY+rZ*rZ,dist = sqrt(dist2);
      h *= 0.5*(1.0 + erf((dist-d_0)/erf_delta));
    }
    if (vdwType < 4) epsilon[ijk] = (dielectric-eps_vacuum)*h + eps_vacuum;
    else{
      double *Axes=QAllocDouble(3);
      HybridAxes(Axes);
      double a=Axes[0],a2=a*a,b=Axes[1],b2=b*b,c=Axes[2],c2=c*c;
      double ellip=(x2/a2)+(y2/b2)+(z2/c2);
      double r_hybrid2=r_hybrid*r_hybrid;
      if (grid_dist2 < r_hybrid2 && ellip < 1.0) epsilon[ijk] = eps_vacuum;
      else epsilon[ijk] = (dielectric-eps_vacuum)*h + eps_vacuum;
      QFree(Axes);
    }
  }
  VRload(eps,NTotPts,0.0);
  PrintEpsilon();
}

void PEQSclass::RigidVDWInterface()
{
  int i,j,k,ijk,l;
  if (z_direction > 0){
    for (ijk=0;ijk<NTotPts;ijk++){
      get_3ind_from_multInd(i,j,k,ijk);
      double x=valueOfX(i),x2=x*x,y=valueOfY(j),y2=y*y,z=valueOfZ(k),z2=z*z,grid_dist2=(x2+y2+z2);
      double d_0=1.0,rmid_2=r_interface*0.5,h=1.0;
      for (l=0;l<NAtoms;l++){
        double xcar=Carts[3*l+0],ycar=Carts[3*l+1],zcar=Carts[3*l+2];
        double xcar2=xcar*xcar,ycar2=ycar*ycar,zcar2=zcar*zcar,cart_dist2=xcar2+ycar2+zcar2;
        int AtomicNo = AtNo[l];
        if (vdwType == 1) d_0 = ang2bohr*VDWRadius(AtomicNo);
        else if (vdwType == 2) d_0 = vdWScale[l]*ang2bohr*VDWRadius(AtomicNo);
        else d_0 = vdWScale[l]*ang2bohr*(VDWRadius(AtomicNo)+vdwShift);
        double rX=x-xcar,rY=y-ycar,rZ=z-zcar;
        double dist2=rX*rX+rY*rY+rZ*rZ,dist=sqrt(dist2);
        h *= 0.5*(1.0 + erf((dist-d_0)/erf_delta));
      }
      if (vdwType < 4) epsilon[ijk] = (dielectric-eps_vacuum)*h + eps_vacuum;
      else{
        double *Axes=QAllocDouble(3);
        HybridAxes(Axes);
        double a=Axes[0],a2=a*a,b=Axes[1],b2=b*b,c=Axes[2],c2=c*c;
        double ellip=(x2/a2)+(y2/b2)+(z2/c2);
        double r_hybrid2=r_hybrid*r_hybrid;
        if (grid_dist2 < r_hybrid2 && ellip < 1.0) epsilon[ijk] = eps_vacuum;
        else epsilon[ijk] = (dielectric-eps_vacuum)*h + eps_vacuum;
        QFree(Axes);
      }
      if (z > gibbsDS+rmid_2) epsilon[ijk] = eps_vacuum;
      else if (z >= gibbsDS-rmid_2 && z <= gibbsDS+rmid_2 && epsilon[ijk] > 0.75*dielectric)
        epsilon[ijk] = 0.5*(eps_vacuum+dielectric)
                     + 0.5*(eps_vacuum-dielectric)*tanh(interface_scale*(z-gibbsDS));
    }
  }
  else if (z_direction < 0){
    for (ijk=0;ijk<NTotPts;ijk++){
      get_3ind_from_multInd(i,j,k,ijk);
      double x=valueOfX(i),x2=x*x,y=valueOfY(j),y2=y*y,z=valueOfZ(k),z2=z*z,grid_dist2=(x2+y2+z2);
      double d_0=1.0,rmid_2=r_interface*0.5,h=1.0;
      for (l=0;l<NAtoms;l++){
        double xcar=Carts[3*l+0],ycar=Carts[3*l+1],zcar=Carts[3*l+2];
        double xcar2=xcar*xcar,ycar2=ycar*ycar,zcar2=zcar*zcar,cart_dist2=xcar2+ycar2+zcar2;
        int AtomicNo = AtNo[l];
        if (vdwType == 1) d_0 = ang2bohr*VDWRadius(AtomicNo);
        else if (vdwType == 2) d_0 = vdWScale[l]*ang2bohr*VDWRadius(AtomicNo);
        else d_0 = vdWScale[l]*ang2bohr*(VDWRadius(AtomicNo)+vdwShift);
        double rX=x-xcar,rY=y-ycar,rZ=z-zcar;
        double dist2=rX*rX+rY*rY+rZ*rZ,dist=sqrt(dist2);
        h *= 0.5*(1.0 + erf((dist-d_0)/erf_delta));
      }
      if (vdwType < 4) epsilon[ijk] = (dielectric-eps_vacuum)*h + eps_vacuum;
      else{
        double *Axes=QAllocDouble(3);
        HybridAxes(Axes);
        double a=Axes[0],a2=a*a,b=Axes[1],b2=b*b,c=Axes[2],c2=c*c;
        double ellip=(x2/a2)+(y2/b2)+(z2/c2);
        double r_hybrid2=r_hybrid*r_hybrid;
        if (grid_dist2 < r_hybrid2 && ellip < 1.0) epsilon[ijk] = eps_vacuum;
        else epsilon[ijk] = (dielectric-eps_vacuum)*h + eps_vacuum;
        QFree(Axes);
      }
      if (z < gibbsDS-rmid_2) epsilon[ijk] = eps_vacuum;
      else if (z <= gibbsDS+rmid_2 && z >= gibbsDS-rmid_2 && epsilon[ijk] > 0.75*dielectric)
        epsilon[ijk] = 0.5*(eps_vacuum+dielectric)
                     + 0.5*(eps_vacuum-dielectric)*tanh(interface_scale*(gibbsDS-z));
    }
  }
}

void PEQSclass::GygiBulkCavity()
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

void PEQSclass::SphericalSurface(double *Center)
{
  double dist,dist2,rX,rY,rZ;
  double rmid=r_cavity+tanh_cut1+0.5*tanh_cut2,rtotal=r_cavity+tanh_cut1,rtotal2=rtotal*rtotal;
  int i,j,k,ijk,l;

  for (ijk=0;ijk<NTotPts;ijk++){
    get_3ind_from_multInd(i,j,k,ijk);
    rX = valueOfX(i)-Center[0];
    rY = valueOfY(j)-Center[1];
    rZ = valueOfZ(k)-Center[2];
    double z = valueOfZ(k);
    dist2 = rX*rX + rY*rY + rZ*rZ;
    dist = sqrt(dist2);
    epsilon[ijk] = 0.5*(dielectric+eps_vacuum)
                 + 0.5*(dielectric-eps_vacuum)*tanh(tanh_scale*(dist-rmid));
    delLogEps[3*ijk+0] = ((0.5*(dielectric-eps_vacuum))/epsilon[ijk])
                       * (1.0-tanh(tanh_scale*(dist-rmid))*tanh(tanh_scale*(dist-rmid)))
                       *(tanh_scale*rX/dist);
    delLogEps[3*ijk+1] = ((0.5*(dielectric-eps_vacuum))/epsilon[ijk])
                       * (1.0-tanh(tanh_scale*(dist-rmid))*tanh(tanh_scale*(dist-rmid)))
                       * (tanh_scale*rY/dist);
    delLogEps[3*ijk+2] = ((0.5*(dielectric-eps_vacuum))/epsilon[ijk])
                       * (1.0-tanh(tanh_scale*(dist-rmid))*tanh(tanh_scale*(dist-rmid)))
                       * (tanh_scale*rZ/dist);
  }
}

void PEQSclass::SphericalInterface(double *Center)
{
  double dist,dist2,rX,rY,rZ;
  double rmid=r_cavity+tanh_cut1+0.5*tanh_cut2,rtotal=r_cavity+tanh_cut1,rtotal2=rtotal*rtotal;
  double rmid_2=r_interface*0.5;
  int i,j,k,ijk,l;

  for (ijk=0;ijk<NTotPts;ijk++){
    get_3ind_from_multInd(i,j,k,ijk);
    rX = valueOfX(i)-Center[0];
    rY = valueOfY(j)-Center[1];
    rZ = valueOfZ(k)-Center[2];
    double z = valueOfZ(k);
    dist2 = rX*rX + rY*rY + rZ*rZ;
    dist = sqrt(dist2);
    if (z_direction > 0){
      if (z > gibbsDS+rmid_2)
        epsilon[ijk] = eps_vacuum;
      else if (z >= gibbsDS-rmid_2 && z <= gibbsDS+rmid_2){
        if (dist2 > rtotal2){
          epsilon[ijk] = 0.5*(eps_vacuum+dielectric)
                       + 0.5*(eps_vacuum-dielectric)*tanh(interface_scale*(z-gibbsDS));
        }
        else epsilon[ijk] = eps_vacuum;
      }
      else if (z < gibbsDS-rmid_2){
        epsilon[ijk] = 0.5*(dielectric+eps_vacuum)
                     + 0.5*(dielectric-eps_vacuum)*tanh(tanh_scale*(dist-rmid));
      }
    }
    else if (z_direction < 0){
      if (z < gibbsDS-rmid_2)
        epsilon[ijk] = eps_vacuum;
      else if (z <= gibbsDS+rmid_2 && z >= gibbsDS-rmid_2){
        if (dist2 > rtotal2){
          epsilon[ijk] = 0.5*(eps_vacuum+dielectric)
                       + 0.5*(eps_vacuum-dielectric)*tanh(interface_scale*(gibbsDS-z));
        }
        else epsilon[ijk] = eps_vacuum;
      }
      else if (z > gibbsDS+rmid_2){
        epsilon[ijk] = 0.5*(dielectric+eps_vacuum)
                     + 0.5*(dielectric-eps_vacuum)*tanh(tanh_scale*(dist-rmid));
      }
    }
  }
}

void PEQSclass::SphereCenter(double* Center)
{
  int l,m;
  double *xmin=QAllocDouble(3),*xmax=QAllocDouble(3);
  VRload(xmin,3,-INVALID_REAL); // note that INVALID_REAL is negative
  VRload(xmax,3, INVALID_REAL);
  for (l=0;l<NAtoms;l++){
    for(m=0;m<3;m++){
      if (xmax[m] < Carts[3*l+m]) xmax[m] = Carts[3*l+m];
      if (xmin[m] > Carts[3*l+m]) xmin[m] = Carts[3*l+m];
    }
  }
  for (m=0;m<3;m++) Center[m] = 0.5*(xmax[m]+xmin[m]);
  QFree(xmin);QFree(xmax);
}

void PEQSclass::SphereRadius()
{
  int l,m;
  double *xmin=QAllocDouble(3),*xmax=QAllocDouble(3);
  r_cavity = INVALID_REAL;
  VRload(xmin,3,-INVALID_REAL); // note that INVALID_REAL is negative
  VRload(xmax,3, INVALID_REAL);
  for (l=0;l<NAtoms;l++){
    for(m=0;m<3;m++){
      if (xmax[m] < Carts[3*l+m]) xmax[m] = Carts[3*l+m];
      if (xmin[m] > Carts[3*l+m]) xmin[m] = Carts[3*l+m];      
    }
  }
  for (m=0;m<3;m++){
    double length = 0.5*(xmax[m]-xmin[m]);
    if (r_cavity < length) r_cavity = length; 
  }
  QFree(xmin);QFree(xmax);
}

void PEQSclass::HybridAxes(double *Axes)
{
  int l,m;
  double *xmin=QAllocDouble(3),*xmax=QAllocDouble(3),d_0=1.0;
  VRload(xmin,3,-INVALID_REAL); // note that INVALID_REAL is negative
  VRload(xmax,3, INVALID_REAL);
  for (l=0;l<NAtoms;l++){
    int AtomicNo = AtNo[l];
    d_0 = vdWScale[l]*ang2bohr*(VDWRadius(AtomicNo)+vdwShift);
    double pos_ext=d_0-2.0*erf_delta,neg_ext=-d_0+2.0*erf_delta;
    for(m=0;m<3;m++){
      if (Carts[3*l+m] >= 0.0 && (xmax[m] < (Carts[3*l+m]+pos_ext)))
        xmax[m] = Carts[3*l+m]+pos_ext;
      if (Carts[3*l+m] < 0.0  && (xmin[m] > (Carts[3*l+m]+neg_ext)))
        xmin[m] = Carts[3*l+m]+neg_ext;
    }
  }
  for (m=0;m<3;m++){
    if (xmax[m] >= -1.0*xmin[m]) Axes[m] = -1.0*xmin[m];
    else Axes[m] = xmax[m];
  }
  QFree(xmin);QFree(xmax);
}
