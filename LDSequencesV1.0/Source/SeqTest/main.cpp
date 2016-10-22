#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream.h>
#include <time.h>
#include <string.h>
#include <stdarg.h>
#include "SimSeqs.h"

long NR;
long DIM;
clock_t timebeforeinit;

void DoSequence(LDSqBase*Numbers, FILE*out, char*name) {
  double buffer[DIM];
  double m[DIM], m2[DIM], m3[DIM], m4[DIM];
  double mean[DIM], vari[DIM], skew[DIM], kurt[DIM];
  int j, i;
  clock_t time;

  // initialize with 0
  for (i=0; i<DIM; i++) {
    buffer[i]=0;
    m[i]=0;
    m2[i]=0;
    m3[i]=0;
    m4[i]=0;
  }

  time=clock();
  for (i=1; i<=NR; i++) {
    Numbers->NextElement(&buffer[0], DIM);

    for (j=0; j<DIM; j++) {
      m[j] +=buffer[j];
      m2[j]+=buffer[j]*buffer[j];
      m3[j]+=buffer[j]*buffer[j]*buffer[j];
      m4[j]+=buffer[j]*buffer[j]*buffer[j]*buffer[j];
    }
  }
  time=clock()-time;
  timebeforeinit=clock()-timebeforeinit;

  fprintf(out, "{%s, %f, %f", name, (double)((double)time/(double)CLOCKS_PER_SEC), (double)((double)timebeforeinit/(double)CLOCKS_PER_SEC));

  for (j=0; j<DIM; j++) {
    mean[j]=m[j]/NR;
    vari[j]=(m2[j]-(mean[j]*mean[j]*NR))/(NR-1);
    skew[j]=(m3[j]/NR - 3.*m2[j]*m[j]/(1.*NR*NR) + 2.*m[j]*m[j]*m[j]/(1.*NR*NR*NR)) * (1.*NR*NR)/((NR-1.)*(NR-2.));
    kurt[j]=(m4[j]/NR - 4.*m3[j]*m[j]/(1.*NR*NR) +
                 6.*m2[j]*m[j]*m[j]/(1.*NR*NR*NR) - 3.*m[j]*m[j]*m[j]*m[j]/(1.*NR*NR*NR*NR));
    kurt[j]=( NR*(1.*NR*NR-2.*NR+3.)*kurt[j] - 3.*NR*(2.*NR-3.)*vari[j]*vari[j]*(NR-1)*(NR-1)/(NR*NR) ) / ((NR-1.)*(NR-2.)*(NR-3.));
    kurt[j]=kurt[j]/(vari[j]*vari[j]) -3;
  }

  fprintf(out, ",\n   {\"mean\"");
  for (j=0; j<DIM; j++) fprintf(out, ", %f", mean[j]-0.5);
  fprintf(out, "}");
  fflush(out);

  fprintf(out, ",\n   {\"variance\"");
  for (j=0; j<DIM; j++) fprintf(out, ", %f", vari[j]-1./12.);
  fprintf(out, "}");
  fflush(out);

  fprintf(out, ",\n   {\"skewness\"");
  for (j=0; j<DIM; j++) fprintf(out, ", %f", skew[j]-0.);
  fprintf(out, "}");
  fflush(out);

  fprintf(out, ",\n   {\"kurtosis\"");
  for (j=0; j<DIM; j++) fprintf(out, ", %f", kurt[j]-(-1.2));
  fprintf(out, "}");
  fflush(out);

  fprintf(out, "\n}");
  fflush(out);

}


int main(int argc, char*argv[]) {
  NR=10000;
  DIM=10;
  if (argc>2) {
    sscanf(argv[2], " %li ", &NR);
  }
  if (argc>1) {
    sscanf(argv[1], " %li ", &DIM);
  }
  cout<<"Testing for dim="<<DIM<<" and nr="<<NR<<". \n";
  cout<<"If you don't like these values, give them as parameters to ";
  cout<<"the programm, e.g. SeqTest 15 1000000\n";

  FILE*out=fopen("NumberTest.m", "w+");
  LDSqBase*Numbers;
  fprintf(out, "{\n  ");
  long bases[DIM];
  ReadPrimes(&bases[0], DIM);

  cout<<"Testing ran1\n";
  timebeforeinit=clock();
  Numbers=(LDSqBase*)new LDSqMonteCarlo(MC1, -12345663);
  DoSequence(Numbers, out, "ran1");
  free(Numbers);
  fprintf(out, ",\n  ");
  fflush(out);

  cout<<"Testing ran2\n";
  timebeforeinit=clock();
  Numbers=(LDSqBase*)new LDSqMonteCarlo(MC2, -12345663);
  DoSequence(Numbers, out, "ran2");
  free(Numbers);
  fprintf(out, ",\n  ");
  fflush(out);

  cout<<"Testing Halton\n";
  timebeforeinit=clock();
  Numbers=(LDSqBase*)new LDSqHalton(&bases[0], DIM, NR, BGENAU, GENAU);
  DoSequence(Numbers, out, "Halton");
  free(Numbers);
  fprintf(out, ",\n  ");
  fflush(out);

  cout<<"Testing (0,s) Nets\n";
  timebeforeinit=clock();
  Numbers=(LDSqBase*)new LDSqNetz(&bases[0], NextPrime(DIM), DIM, NR, 1);
  DoSequence(Numbers, out, "Netz");
  free(Numbers);
  fprintf(out, ",\n  ");
  fflush(out);

  cout<<"Testing Niederreiter\n";
  timebeforeinit=clock();
  Numbers=(LDSqBase*)new LDSqNiederreiter(NULL, NextPrime(DIM), DIM, NR);
  DoSequence(Numbers, out, "Niederreiter");
  free(Numbers);
  fprintf(out, ",\n  ");
  fflush(out);

  cout<<"Testing Sobol\n";
  timebeforeinit=clock();
  Numbers=(LDSqBase*)new LDSqSobol(&bases[0], DIM, NR);
  DoSequence(Numbers, out, "Sobol");
  free(Numbers);
  fprintf(out, ",\n  ");
  fflush(out);

  cout<<"Testing Faure\n";
  timebeforeinit=clock();
  Numbers=(LDSqBase*)new LDSqFaure(&bases[0], DIM, NR);
  DoSequence(Numbers, out, "Faure");
  free(Numbers);
  fprintf(out, ",\n  ");
  fflush(out);

  cout<<"Testing Atanassov\n";
  timebeforeinit=clock();
  Numbers=(LDSqBase*)new LDSqAtanassov(&bases[0], DIM, NR);
  DoSequence(Numbers, out, "Atanassov");
  free(Numbers);
  fprintf(out, ",\n  ");
  fflush(out);

  cout<<"Testing Hammersley\n";
  timebeforeinit=clock();
  Numbers=(LDSqBase*)new LDSqHammersley(&bases[0], DIM, NR);
  DoSequence(Numbers, out, "Hammersley");
  free(Numbers);


  fprintf(out, "\n}");
  fflush(out);
  fclose(out);
  return 0;
}
