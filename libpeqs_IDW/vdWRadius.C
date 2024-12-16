#include "peqs_class.h"

double PEQSclass::VDWRadius(int AtomicNo)
{
  double VDWRadius=1.0;
  if (AtomicNo == 1) VDWRadius = 1.10;
  else if (AtomicNo == 2) VDWRadius = 1.40;
  else if (AtomicNo == 3) VDWRadius = 1.81;
  else if (AtomicNo == 4) VDWRadius = 1.53;
  else if (AtomicNo == 5) VDWRadius = 1.92;
  else if (AtomicNo == 6) VDWRadius = 1.70;
  else if (AtomicNo == 7) VDWRadius = 1.55;
  else if (AtomicNo == 8) VDWRadius = 1.52;
  else if (AtomicNo == 9) VDWRadius = 1.47;
  else if (AtomicNo == 10) VDWRadius = 1.54;
  else if (AtomicNo == 11) VDWRadius = 2.27;
  else if (AtomicNo == 12) VDWRadius = 1.73;
  else if (AtomicNo == 13) VDWRadius = 1.84;
  else if (AtomicNo == 14) VDWRadius = 2.10;
  else if (AtomicNo == 15) VDWRadius = 1.80;
  else if (AtomicNo == 16) VDWRadius = 1.80;
  else if (AtomicNo == 17) VDWRadius = 1.75;
  else if (AtomicNo == 18) VDWRadius = 1.88;
  else if (AtomicNo == 19) VDWRadius = 2.75;
  else if (AtomicNo == 20) VDWRadius = 2.31;
  else if (AtomicNo == 22) VDWRadius = 1.588;
  else if (AtomicNo == 29) VDWRadius = 1.40;
  else if (AtomicNo == 30) VDWRadius = 1.39;
  else if (AtomicNo == 31) VDWRadius = 1.87;
  else if (AtomicNo == 32) VDWRadius = 2.11;
  else if (AtomicNo == 33) VDWRadius = 1.85;
  else if (AtomicNo == 34) VDWRadius = 1.90;
  else if (AtomicNo == 35) VDWRadius = 1.85;
  else if (AtomicNo == 36) VDWRadius = 2.02;
  else if (AtomicNo == 37) VDWRadius = 3.03;
  else if (AtomicNo == 38) VDWRadius = 2.49;
  else if (AtomicNo == 47) VDWRadius = 1.72;
  else if (AtomicNo == 49) VDWRadius = 1.93;
  else if (AtomicNo == 50) VDWRadius = 2.17;
  else if (AtomicNo == 51) VDWRadius = 2.06;
  else if (AtomicNo == 52) VDWRadius = 2.06;
  else if (AtomicNo == 53) VDWRadius = 1.98;
  else if (AtomicNo == 54) VDWRadius = 2.16;
  else if (AtomicNo == 55) VDWRadius = 3.43;
  else if (AtomicNo == 56) VDWRadius = 2.68;
  else if (AtomicNo == 78) VDWRadius = 1.75;
  else if (AtomicNo == 79) VDWRadius = 1.66;
  else if (AtomicNo == 81) VDWRadius = 1.96;
  else if (AtomicNo == 82) VDWRadius = 2.02;
  else if (AtomicNo == 83) VDWRadius = 2.07;
  else if (AtomicNo == 84) VDWRadius = 1.97;
  else if (AtomicNo == 85) VDWRadius = 2.02;
  else if (AtomicNo == 86) VDWRadius = 2.30;
  else if (AtomicNo == 87) VDWRadius = 4.38;
  else if (AtomicNo == 88) VDWRadius = 2.83;
  else QCrash(" Requested van der Waals radius not available!");
  return VDWRadius;
}
