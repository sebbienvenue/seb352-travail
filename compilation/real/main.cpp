# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "sobol.hpp"

int main () {

int m, n, seed =1;
int a=0, b=1;
int *cseed;
double t[12];
float tab[6*2*2];

cseed=&seed;
n=4;
m=6;

timestamp ( );
r8mat_write("/home/bienvenue/Documents/seb352-travail/compilation/real/ok.txt", m, n, t);
cout << "\n2\n";
cout << i8_uniform(a, b, cseed) << "\n";

tab=i4_sobol_generate(m, n, 1);

return 0;
}
