#!bin/sh
clear

swig -c++ -python test.i
g++ -O2 -fPIC -c changed.cpp
g++ -O2 -fPIC -c test_wrap.cxx `pkg-config --cflags python`
g++ -shared changed.o test_wrap.o -o _test.so
