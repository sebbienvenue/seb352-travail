int seb ( int spatialDIM, int Npoints, int avoid )

{
  int i;
  int m;
  ostringstream m_ostring;
  int n;
  ostringstream n_ostring;
  string output_filename;
  double *r;
  int skip;


  m = spatialDIM;
  n = Npoints;
  skip = avoid;

//  Compute the data.
//
  r = i8_sobol_generate ( m, n, skip );
//
//  Write it to a file.
//
  m_ostring << m;
  n_ostring << n;

  output_filename = "sobol_" + m_ostring.str ( ) + "_" 
    + n_ostring.str ( ) + ".txt";

  r8mat_write ( output_filename, m, n, r );

  cout << "\n";
  cout << "  The data was written to the file \"" 
    << output_filename << "\".\n";
//
//  Terminate.
//
  delete [] r;

  cout << "\n";
  cout << "SOBOL_DATASET:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}


swig -c++ -python test.i
g++ -O2 -fPIC -c changed.cpp
g++ -O2 -fPIC -c test_wrap.cxx `pkg-config --cflags python`
g++ -shared changed.o test_wrap.o -o _test.so

