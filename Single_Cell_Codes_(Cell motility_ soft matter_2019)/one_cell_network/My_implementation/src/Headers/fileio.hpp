// Output functions
// Store length of the lattice
#include "globals.hpp"

void OutputLx(std::string frac, std::string job_index) {

    std::ofstream info;
    info.exceptions(std::ofstream::failbit | std::ofstream::badbit);

    info.open(out_dir + "f" + frac + "id" + job_index + "Lx.txt");

    info.close();
}
