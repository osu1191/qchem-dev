#include "peqs_class.h"
#include "rem_values.h"
#include "../libscrf/libscrf/pcm/pcm.h"
//
// Atomic vdW radii for cavity construction.
// Updated to allow options other than Bondi.
// JMH (2024)
//

double PEQSclass::VDWRadius(int Z)
{ // value is returned in Angstroms
  double r;
  int itype = rem_read(REM_PEQS_VDW_RADII);
  if (itype == BONDI)
     r = BondiRadius(Z);
  else if (itype == UFF) 
     r = UFFRadius(Z);
  else if (itype == READ)
     r = UserDefdRadius2(Z);
  else
  {
     printf(" PEqS vdW radii type %d are selected\n",itype);
     QCrash("Unrecognized vdW radius choice in PEQSclass::VDWRadius");
  }
  
  //printf("\tPEqS:  Use vdW radius r = %.3f A for Z = %d\n",r,Z);
  return r;
}

