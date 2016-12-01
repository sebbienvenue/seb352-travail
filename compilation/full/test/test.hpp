# include <cstdlib>
# include <cmath>
# include <ctime>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <sstream>
# include <string>
//# include "test.hpp"
using namespace std;

int seb ( int spatialDIM, int Npoints, int avoid );
int sobol ( int argc, char *argv[] );
int i8_bit_hi1 ( long long int n );
int i8_bit_lo0 ( long long int n );
void i8_sobol ( int dim_num, long long int *seed, double quasi[ ] );
double *i8_sobol_generate ( int m, int n, int skip );
void r8mat_write ( string output_filename, int m, int n, double table[] );
void timestamp ( );


