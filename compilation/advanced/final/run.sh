#!bin/sh
clear

swig -c++ -python sobol.i
g++ -O2 -fPIC -c changed.cpp
g++ -O2 -fPIC -c sobol_wrap.cxx `pkg-config --cflags python`
g++ -shared changed.o sobol_wrap.o -o _sobol.so
