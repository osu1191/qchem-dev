#include "peqs_class.h"

void PEQSclass::GetCOM(double *COM)
{
  int ijk,AtomicNo;
  double AtMass=0.0,TotMass=0.0;
  VRload(COM,3,0.0);
  for (ijk=0;ijk<NAtoms;ijk++){
    AtomicNo = AtNo[ijk];
    AtMass = PeriodicTable(AtomicNo);
    COM[0] += AtMass*Carts[3*ijk+0];
    COM[1] += AtMass*Carts[3*ijk+1];
    COM[2] += AtMass*Carts[3*ijk+2];
    TotMass += AtMass;
  }
  VRscale(COM,3,(1.0/TotMass));
}

