//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
//#include <iostream.h>
#include <time.h>
#include <string>
//#include <stdarg.h>
#include "Distributions.h"
#include "SimSeqs.h"

long NR;
double para1, para2;
LDSqBase*numbers;

void DoDist(Distribution*dist, ofstream &out) {
  double nr=0, n=0;
  double m1=0, m2=0, m3=0, m4=0;
  clock_t time;

  time=clock();
  for (int i=1; i<=NR; i++) {
    numbers->NextElement(&n,1);
    nr=dist->Transform(n, para1, para2);
//    nr=n;

    m1 += nr;
    m2 += nr*nr;
    m3 += nr*nr*nr;
    m4 += nr*nr*nr*nr;
  }
  time=clock()-time;
  // Mn denote the mean of x^n
  double M1=m1/(double)NR;
  double M2=m2/(double)NR;
  double M3=m3/(double)NR;
  double M4=m4/(double)NR;

//  m2[j]/NR-(m[j]*m[j])/(1.*NR*NR));
  out<<"{"<<(double)time/CLOCKS_PER_SEC<<", \""<<dist->DistributionName()<<"\", param1="<<para1<<", param2="<<para2<<", "<<
      M1<<", "<<
      (double)(NR)/(double)(NR-1)*(M2-M1*M1)<<", "<<
      M3-3*M1*M2 + 2*M1*M1*M1<<", "<<
      M4-4*M3*M1 + 6*M2*M1*M1- 3*M1*M1*M1*M1<<"}}"<<flush;
  cout<<"   {"<<(double)time/CLOCKS_PER_SEC<<", \""<<dist->DistributionName()<<"\", param1="<<para1<<", param2="<<para2<<", "<<
      M1<<", "<<
      (double)(NR)/(double)(NR-1)*(M2-M1*M1)<<", "<<
      M3-3*M1*M2 + 2*M1*M1*M1<<", "<<
      M4-4*M3*M1 + 6*M2*M1*M1- 3*M1*M1*M1*M1<<"}}"<<endl<<flush;
//      (double)(m4/(double)NR - (4.*m3*m1)/(double)(NR*NR) +
//           (6.*m2*m1*m1)/(double)(NR*NR*NR) - (3.*m1*m1*m1*m1)/(double)(NR*NR*NR*NR))<<"}"<<flush;
}


int main(int argc, char*argv[]) {
  NR=100000;
  para1=0;
  para2=1;
  if (argc>1) {
    sscanf(argv[1], " %li ", &NR);
  }
  if (argc>2) {
    sscanf(argv[2], " %lg ", &para1);
  }
  if (argc>3) {
    sscanf(argv[3], " %lg ", &para2);
  }
  cout<<"Testing with NR="<<NR<<" for param1="<<para1<<" and para2="<<para2<<". \n";
  cout<<"If you don't like these values, give them as parameters to ";
  cout<<"the programm, e.g. DistTest.x 1000000 0 1\n";


  numbers=(LDSqBase*) new LDSqMonteCarlo();
  ofstream out("DistTest.m");
  out<<"{\n  ";
  Distribution*dist;
  for (int i=1; i<=6; i++) {
    dist=CreateDistribution(i, para1, para2);
    cout<<"Testing "<<dist->DistributionName()<<"\n";
    DoDist(dist, out);
    delete dist;
    out<<",\n  "<<flush;
  }
  out.close();

  return 0;
}
