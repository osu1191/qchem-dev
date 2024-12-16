#include "Electrolyte.h"
#include <cmath>

Ion::Ion(int charge):
  charge_(charge) {

}

Ion::Ion(int charge, double effective_radius):
  charge_(charge),effective_radius_(effective_radius) {
  calculate_c_max();
}

Ion::~Ion(){
}

void Ion::calculate_c_max(){
  c_max_ = 0.74*3.0/(4.0*M_PI*std::pow(effective_radius_,3.0));
}


Electrolyte::Electrolyte(std::vector<Ion> Ions, double concentration):
  Ions_(Ions), concentration_(concentration) {
  concentration_au_ = concentration*conv_conc_au_;
  calculate_c_comb();
//  std::cout << "Volume fraction occupied by ions "<< 2.0*concentration_au_/c_comb_ << " ." << std::endl;
}

Electrolyte::Electrolyte(){
}

Electrolyte::~Electrolyte(){
}

void Electrolyte::calculate_c_comb(){
//right now implemented only for electrolyte ions of equal size
  c_comb_ = Ions_[0].get_c_max();
}

void Electrolyte::set_concentration(double conc){
  concentration_ = conc;
   concentration_au_ = conc*conv_conc_au_;
}

void Electrolyte::set_ions(std::vector<Ion> Ions){
  Ions_ = Ions;
}
