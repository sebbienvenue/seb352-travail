#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import collections
import copy
import sys

from mpi4py import MPI
import numpy as np

import lbfgs
from RuNNer import wrapper as RuNNer

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# input.data files formatting
LAT_FMT = 'lattice {lat[0]:16.8f} {lat[1]:16.8f} {lat[2]:16.8f}'
ATOM_FMT = 'atom    {xyz[0]:16.8f} {xyz[1]:16.8f} {xyz[2]:16.8f} {symbol:2s} 0 {nrg:16.8f} {forces[0]:16.8f} {forces[1]:16.8f} {forces[2]:16.8f}'
NRG_FMT = 'energy {nrg:17.8f}'

# mapping of atomic symbol to number
symbolToNumber = {
    'O'  :  8,
    'Cu' : 29,
    'Zn' : 30
}
# reverse mapping
numberToAtom = {}
for k, v in symbolToNumber.items(): numberToAtom[v] = k

def init(fileName):
    '''Initializes the calculation by distributing the initial structure and initializing the RuNNer wrapper.'''
    if rank == 0:
        elems, kinds, lat_count, natoms, xyz, zelem = set(), [], 0, 0, [], []
        periodic, lattice = np.array(0), np.zeros((3, 3))
        with open(fileName) as fh:
            for line in fh:
                if 'atom' in line:
                    fields = line.split()
                    elems.add(fields[4])
                    kinds.append(int(fields[-1]))
                    xyz.append(list(map(float, fields[1:4])))
                    zelem.append(symbolToNumber[fields[4]])
                    natoms += 1
                elif 'end' in line: break
                elif 'lattice' in line:
                    fields = line.split()
                    lattice[:,lat_count] = list(map(float, fields[1:4]))
                    periodic = np.array(1)
                    lat_count += 1

        natoms = np.array(natoms)
        nelems = np.array(len(elems))
        xyz = np.array(xyz, order='F').transpose()
        zelem = np.array(zelem, dtype=np.int32)

        for var in (natoms, nelems, periodic, lattice, xyz, zelem): comm.Bcast(var)

        RuNNer(natoms, nelems)
        return xyz, zelem, lattice, periodic, kinds

    else:
        natoms, nelems = np.array(0), np.array(0)
        periodic, lattice = np.array(0), np.zeros((3, 3))

        for var in (natoms, nelems, periodic, lattice): comm.Bcast(var)

        xyz = np.zeros((3, natoms))
        zelem = np.zeros(natoms, dtype=np.int32)
        comm.Bcast(xyz)
        comm.Bcast(zelem)

        RuNNer(natoms, nelems)
        return xyz, zelem, lattice, periodic, None

def main():
    xyz, zelem, lattice, periodic, kinds = init(sys.argv[1])

    do_continue = np.array(1)
    if rank == 0:
        # some constants for lbfgs; TODO: make them adjustable
        max_num_iter = 200
        m = 10
        eps = 1.0e-6
        xtol = 1.0e-16

        # some constants needed
        w = np.zeros(xyz.size * (2 * m + 1) + 2 * m, dtype=np.float64)
        diagco = False
        diag = np.zeros(xyz.size, dtype=np.float64)
        iprint = np.zeros(2, dtype=np.int32); iprint[0] = -1
        
        # set up constraints
        mask = np.zeros(xyz.shape, dtype=np.float64, order='F')
        for i, kind in enumerate(kinds):
            if kind == 1: mask[:,i] = [1, 1, 1]

        # flag to control the optimization
        iflag = np.array(0)

        with open('optimization.data', 'w') as fh:
            for i in range(max_num_iter):
                comm.Bcast(do_continue)
                comm.Bcast(xyz)
                nrg, forces, atomicNrgs, ex_atoms, bond_underrun = RuNNer.getenergyandforces(xyz, zelem, lattice, periodic)
                writeInputData(fh, lattice, xyz, zelem, atomicNrgs, forces, nrg)
                print('step {}\nenergy: {}'.format(i + 1, nrg))
                if np.any(ex_atoms): print('extrapolated atoms: {}'.format(str(ex_atoms.nonzero()[0] + 1)))
                if bond_underrun: print('bond underrun occured')
                x, g = xyz.ravel(), (-forces * mask).ravel()
                lbfgs.lbfgs(m, x, nrg, g, diagco, diag, iprint, eps, xtol, w, iflag)
                if iflag == 1: xyz = x.reshape(xyz.shape)
                else: break
        do_continue = np.array(0)
        comm.Bcast(do_continue)
    else:
        while True:
            comm.Bcast(do_continue)
            if do_continue == 0: break
            comm.Bcast(xyz)
            RuNNer.getenergyandforces(xyz, zelem, lattice, periodic)

def writeInputData(fh, lattice, xyz, zelem, atomicNrgs, forces, nrg):
    '''Creates an input.data formatted file.'''
    print('begin', file=fh)
    for i in range(3): print(LAT_FMT.format(lat=lattice[i]), file=fh)
    for i in range(xyz.shape[1]):
        print(ATOM_FMT.format(xyz=xyz[:,i], symbol=numberToAtom[zelem[i]], nrg=atomicNrgs[i], forces=forces[:,i]), file=fh)
    print(NRG_FMT.format(nrg=nrg), file=fh)
    print('end', file=fh)
    fh.flush()

if __name__ == '__main__': main()

# vim: set et ft=python sw=4 ts=4 tw=132:
