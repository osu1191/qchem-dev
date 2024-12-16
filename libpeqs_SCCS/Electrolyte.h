#include <vector>
#include <iostream>

class Ion{

private:
 int charge_;
 double effective_radius_ = 0.0; 
 double c_max_ = 0.0;

public:
 //Constructor
 Ion(int charge);

 //alternative Constructor
 Ion(int charge, double effective_radius);

 //Destructor
 ~Ion();

 int get_charge(){return charge_;}
 double get_effective_radius(){return effective_radius_;}
 double get_c_max(){return c_max_;}
 void calculate_c_max();
};


class Electrolyte{

private:
 double conv_conc_au_ = 8.92389387e-5;
 double concentration_;
 double concentration_au_;
 double c_comb_;

public:

 //Constructors
 Electrolyte(std::vector<Ion> Ions, double concentration);
 Electrolyte();

 //Destructor
 ~Electrolyte();

 void calculate_c_comb();

 double get_concentration(){return concentration_;}
 double get_concentration_au(){return concentration_au_;}
 double get_number_ions(){return Ions_.size();}
 double get_c_comb(){return c_comb_;}

 void set_concentration(double conc);
 void set_ions(std::vector<Ion>);

 std::vector<Ion> Ions_;
};

