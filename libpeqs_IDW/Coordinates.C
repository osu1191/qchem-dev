#include "peqs_class.h"

void PEQSclass::Coordinates()
{
  int i,j;
  get_carts(NULL,&Carts,&AtNo,&NAtoms);
  if (CavityType == PEQS_CAVTYPE_VDW)
  {
    for (i=0;i<NAtoms;i++){
      if (AtNo[i] == 1) vdWScale[i] = 1.0;
      else vdWScale[i] = vdwScale;
    }
  }
  if (multigrid && ((NPoints[0]-1) % pow(2,(mg_order-1)) != 0
                ||  (NPoints[1]-1) % pow(2,(mg_order-1)) != 0
                ||  (NPoints[2]-1) % pow(2,(mg_order-1)) != 0))
  {
    printf(" PEqS grid problem\n");
    printf(" Multi-grid order (MG) = %d\n",mg_order);
    printf(" (Nx, Ny, Nz) = (%d, %d, %d)\n",NPoints[0],NPoints[1],NPoints[2]);
    QCrash(" For PEqS multi-grid calculations, (N-1) must be divisible by 2^(MG-1)");
  }
  if (mg_order < 0 || mg_order > 4)
    QCrash(" The multi-grid order must be = 1, 2, 3, or 4");

  for (i=1; i<mg_order; i++){
    NPoints[3*i+0] = (NPoints[3*(i-1)+0]-1)/2 + 1;
    NPoints[3*i+1] = (NPoints[3*(i-1)+1]-1)/2 + 1;
    NPoints[3*i+2] = (NPoints[3*(i-1)+2]-1)/2 + 1;
  }
  for (i=0; i<mg_order; i++){
    NTotPts_mg[i] = 1;
    for (j=0; j<NumDim; j++){
       NTotPts_mg[i] *= NPoints[3*i+j];
       delX[3*i+j] = (XMax[j]-XMin[j])/(double((NPoints[3*i+j]-1)));
    }
  }

  if (XYZ == NULL)
     XYZ = QAllocDouble(NTotPts*NumDim);
  else
     QCrash(" Attempted to create an existing XYZ list");
  for (int ijk=0; ijk<NTotPts; ijk++){
    int i,j,k;
    get_3ind_from_multInd(i,j,k,ijk);
    XYZ[NumDim*ijk+0] = valueOfX(i);
    XYZ[NumDim*ijk+1] = valueOfY(j);
    XYZ[NumDim*ijk+2] = valueOfZ(k);
//  cout << valueOfX(i) << '\t' << valueOfY(j) << '\t' << valueOfZ(k) << endl;
  }
}
