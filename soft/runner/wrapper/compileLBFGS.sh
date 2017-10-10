#!/bin/bash
source settings.sh

f2py3 --f90exec=${compiler} \
      --fcompiler=${compiler_type} \
      --f90flags="${fflags}" \
      --opt="${opt}" \
      -lopenblas \
      -c lbfgs.pyf \
      lbfgs.f
