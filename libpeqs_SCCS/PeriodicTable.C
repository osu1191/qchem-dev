#include "peqs_class.h"

double PEQSclass::PeriodicTable(int AtomicNo)
{
  double AtomicMass=0.0;
  if (AtomicNo == 1) AtomicMass = 1.00784;
  else if (AtomicNo == 2) AtomicMass = 4.002602;
  else if (AtomicNo == 3) AtomicMass = 6.938;
  else if (AtomicNo == 4) AtomicMass = 9.012183;
  else if (AtomicNo == 5) AtomicMass = 10.806;
  else if (AtomicNo == 6) AtomicMass = 12.0096;
  else if (AtomicNo == 7) AtomicMass = 14.00643;
  else if (AtomicNo == 8) AtomicMass = 15.99903;
  else if (AtomicNo == 9) AtomicMass = 18.998403;
  else if (AtomicNo == 10) AtomicMass = 20.1797;
  else if (AtomicNo == 11) AtomicMass = 22.989769;
  else if (AtomicNo == 12) AtomicMass = 24.304;
  else if (AtomicNo == 13) AtomicMass = 26.9815385;
  else if (AtomicNo == 14) AtomicMass = 28.084;
  else if (AtomicNo == 15) AtomicMass = 30.973761998;
  else if (AtomicNo == 16) AtomicMass = 32.059;
  else if (AtomicNo == 17) AtomicMass = 35.446;
  else if (AtomicNo == 18) AtomicMass = 39.948;
  else if (AtomicNo == 19) AtomicMass = 39.0983;
  else if (AtomicNo == 20) AtomicMass = 40.078;
  else if (AtomicNo == 21) AtomicMass = 44.955908;
  else if (AtomicNo == 22) AtomicMass = 47.867;
  else if (AtomicNo == 23) AtomicMass = 50.9415;
  else if (AtomicNo == 24) AtomicMass = 51.9961;
  else if (AtomicNo == 25) AtomicMass = 54.938044;
  else if (AtomicNo == 26) AtomicMass = 55.845;
  else if (AtomicNo == 27) AtomicMass = 58.933194;
  else if (AtomicNo == 28) AtomicMass = 58.6934;
  else if (AtomicNo == 29) AtomicMass = 63.546;
  else if (AtomicNo == 30) AtomicMass = 65.38;
  else if (AtomicNo == 31) AtomicMass = 69.723;
  else if (AtomicNo == 32) AtomicMass = 72.630;
  else if (AtomicNo == 33) AtomicMass = 74.921595;
  else if (AtomicNo == 34) AtomicMass = 78.971;
  else if (AtomicNo == 35) AtomicMass = 79.901;
  else if (AtomicNo == 36) AtomicMass = 83.798;
  else if (AtomicNo == 37) AtomicMass = 85.4678;
  else if (AtomicNo == 38) AtomicMass = 87.62;
  else if (AtomicNo == 39) AtomicMass = 88.90584;
  else if (AtomicNo == 40) AtomicMass = 91.224;
  else if (AtomicNo == 41) AtomicMass = 92.90637;
  else if (AtomicNo == 42) AtomicMass = 95.95;
  else if (AtomicNo == 43) {AtomicMass = 98.0; cout << " Warning, need to add and specify which isotope to use!" << endl;}
  else if (AtomicNo == 44) AtomicMass = 101.07;
  else if (AtomicNo == 45) AtomicMass = 102.90550;
  else if (AtomicNo == 46) AtomicMass = 106.42;
  else if (AtomicNo == 47) AtomicMass = 107.8682;
  else if (AtomicNo == 48) AtomicMass = 112.414;
  else if (AtomicNo == 49) AtomicMass = 114.818;
  else if (AtomicNo == 50) AtomicMass = 118.710;
  else if (AtomicNo == 51) AtomicMass = 121.760;
  else if (AtomicNo == 52) AtomicMass = 127.60;
  else if (AtomicNo == 53) AtomicMass = 126.90447;
  else if (AtomicNo == 54) AtomicMass = 131.293;
  else QCrash(" Requested atomic mass not available!");
  return AtomicMass;
}
