#include "peqs_class.h"

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
//  else if(CavityType == PEQS_CAVTYPE_ARBI){
//    printf("\tPEqS:  User input dielectric is being used for arbitrary cavity\n");
//    ArbitraryCavity();
//  }

  for (int ijk=0;ijk<NTotPts;ijk++) 
     log_epsilon[ijk] = log(epsilon[ijk]);
  if (CavityType != PEQS_CAVTYPE_SPHERE || Interface)
     FinDif1stDeriv(log_epsilon,delLogEps);
}

/*
void PEQSclass::ArbitraryCavity()
{
  double *pts = QAllocDouble(4*NTotPts);
  FileMan(FM_READ,FILE_PEQS_GRID,FM_DP,4*NTotPts,0,FM_BEG,pts);
  for(int i=3,ijk=0;i<4*NTotPts,ijk<NTotPts;i+=4,ijk++){
    epsilon[ijk]=pts[i];
  }
}
*/

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
