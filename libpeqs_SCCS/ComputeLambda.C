#include "peqs_class.h"

void PEQSclass::ComputeLambda(){
//simplest case for spherical cavity following Fisicaro et al. doi: 10.1063/1.4939125
//lambda = (epsilon(r) -1)/(epsilon_sol -1)
if(CavityType == PEQS_CAVTYPE_SPHERE){
  int ijk;
  #pragma omp for
  for(ijk=0;ijk<NTotPts;ijk++)
  {
    lambda[ijk] = (epsilon[ijk] -1.0)/(eps_solvent - 1.0);
  }
}
else if(CavityType == PEQS_CAVTYPE_VDW){
//for vdw cavity we define lambda directly in terms of an error function just as epsilon
//this allows for an inclusion of the Stern layer
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
      h *= 0.5*(1.0 + erf((dist-d_0-stern_thickness)/erf_delta));
    }
    lambda[ijk] = h;
  }
}


}
