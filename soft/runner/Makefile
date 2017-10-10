# Makefile for RuNNer, JB 2008

#VERSION = 0.65 since 2011 versions not used anymore

# general settings
FC     = gfortran
#FC     = gfortran
MPIFC  = mpif90

###############################################
# MACHINE SPECIFIC SETTINGS
###############################################

# basel:
FFLAGS        = -O3 -ffree-line-length-none $(LIB)
#FFLAGS        = -O3 -132 -free -warn unused -check all -traceback -CB -g -p $(LIB)
#FFLAGS        = -O3 -132 -Bstatic $(LIB)
#FFLAGS        = -O3 -132 -free -warn unused -check all -traceback -CB -g -p $(LIB)
FLAGS_MPI    = -O3 -132 $(LIB_MPI) -I/opt/openmpi/1.3.3/include
#LIB           = -L/opt/intel/mkl/10.0.5.025/lib/em64t -lguide -lmkl_lapack -lmkl_core -lmkl_intel_thread -lmkl_intel_lp64 -lpthread
#LIB_MPI       = -L/opt/intel/mkl/10.0.5.025/lib/em64t -lguide -lmkl_lapack -lmkl_core -lmkl_intel_thread -lmkl_intel_lp64

#LIB     = -L/opt/intel/mkl/10.2.5.035/lib/em64t -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack -liomp5
#LIB           = -L/opt/intel/mkl/10.0.5.025/lib/em64t -lguide -lmkl_em64t -lpthread
#LIB_MPI       = -L//opt/intel/mkl/10.0.5.025/lib/em64t -lguide -lmkl_em64t

# breda:
#FFLAGS        = -O3 -132 $(LIB)
#FFLAGS        = -O3 -132 -free -check all -traceback -CB -g -p -inline_debug_info $(LIB)
#FFLAGS_MPI    = -O3 -132 $(LIB_MPI) -I/opt/openmpi/include
#LIB           = -L/opt/intel/mkl721/lib/em64t -lguide -lmkl_em64t -lpthread
#LIB_MPI       = -L/opt/intel/mkl721/lib/em64t -lguide -lmkl_em64t
#LIB_SCALAPACK = -L/opt/intel/mkl721/lib/em64t -lguide -lmkl_em64t
#FFLAGS        = -O3 -132 $(LIB)
#FFLAGS        = -O3 -132 -free -check all -traceback -CB -g -p -inline_debug_info $(LIB)
#FFLAGS_MPI    = -O3 -132 $(LIB_MPI) -I/opt/openmpi/include
#LIB           = -L/opt/intel/mkl721/lib/em64t -lguide -lmkl_em64t -lpthread
#LIB_MPI       = -L/opt/intel/mkl721/lib/em64t -lguide -lmkl_em64t
#LIB_SCALAPACK = -L/opt/intel/mkl721/lib/em64t -lguide -lmkl_em64t

# lugano: (old laptop)
#FFLAGS        = -O3 -132 $(LIB)
#LIB           = -L/opt/intel/mkl/10.0.1.014/lib/32 -lmkl_ia32 -lmkl -lmkl_lapack -liomp5 -lpthread

# tamaro:
#FFLAGS        = -O3 -132 $(LIB)
#FFLAGS        = -O3 -132 -free -check all -traceback -CB -g -p -inline_debug_info $(LIB)
#FFLAGS        = -O3 -132 -free -check all -traceback -CB -g -p -inline_debug_info $(LIB) -I/opt/openmpi/1.2.6/include
#FFLAGS_MPI    = -O3 -132 $(LIB_MPI) -I/opt/openmpi/1.2.6/include
#FFLAGS_MPI    = -O3 -132 -free -check all -traceback -CB -g -p -inline_debug_info $(LIB_MPI) -I/opt/openmpi/1.2.6/include
#FFLAGS_SCALAPACK = -O3 -132 $(LIB_SCALAPACK) -I/opt/openmpi/1.2.6/include
#LIB              = -L/opt/intel/mkl/10.0.5.025/lib/em64t -lguide -lmkl -lpthread
#LIB_MPI          = -L/opt/intel/mkl/10.0.5.025/lib/em64t -lguide -lmkl
#LIB_SCALAPACK    = -L/opt/intel/mkl/10.0.5.025/lib/em64t -lmkl_solver -lmkl_cdft -lmkl_core -lmkl_blacs_openmpi_lp64 -lmkl_scalapack -lmkl_intel_thread -liomp5 -lguide -lmkl

# generoso:
#FFLAGS     = -O3 -132 $(LIB)
#FFLAGS     = -O3 -132 -free -check all -traceback -CB -g -p -inline_debug_info  $(LIB)
#FFLAGS     = -O3 -132 -warn unused  -free -check all -traceback -CB -g -p  $(LIB)
#FFLAGS_MPI = -O3 -132 -free -check all -traceback -CB -g -p -inline_debug_info $(LIB_MPI) -I/cvos/shared/apps/openmpi/intel/64/1.4/include
#FFLAGS_MPI = -O3 -132 $(LIB_MPI) -I/cvos/shared/apps/openmpi/intel/64/1.4/include
#FFLAGS_MPI = -O3 -132 -free -check all -traceback -CB -g -p -inline_debug_info $(LIB_MPI) -I/cvos/shared/apps/openmpi/intel/64/1.4/include
#FFLAGS_MPI = -O3 -132 $(LIB_MPI) -I/cvos/shared/apps/openmpi/intel/64/1.4/include
#FFLAGS_SCALAPACK = -O3 -132 $(LIB_SCALAPACK)
#LIB     = -L/cvos/shared/apps/intel/mkl/10.2.2.025/lib/em64t -lguide -lmkl_lapack -lmkl_core -lmkl_intel_thread -lmkl_intel_lp64
#LIB_MPI = -L/cvos/shared/apps/intel/mkl/10.2.2.025/lib/em64t -lguide -lmkl_lapack -lmkl_core -lmkl_intel_thread -lmkl_intel_lp64
#LIB_SCALAPACK = -L/cvos/shared/apps/intel/mkl/10.2.2.025/lib/em64t

# bre:
#FFLAGS     = -O3 -132 $(LIB) -I/opt/openmpi/intel/64/1.4/include
#FFLAGS_MPI = -O3 -132 -free -check all -traceback -CB -g -p -inline_debug_info $(LIB_MPI) -I/opt/openmpi/intel/64/1.4/include/
#LIB     = -L/opt/intel/mkl/10.2.2.025/lib/em64t -lguide -lmkl_lapack -lmkl_core -lmkl_intel_thread -lmkl_intel_lp64
#LIB_MPI     = -L/opt/intel/mkl/10.2.2.025/lib/em64t -lguide -lmkl_lapack -lmkl_core -lmkl_intel_thread -lmkl_intel_lp64
#

# runner:
#FFLAGS     = -O3 -132 $(LIB)
#FFLAGS     = -O3 -132 -free -warn unused -check all -traceback -CB -g -p  $(LIB)
#FFLAGS_MPI = -O3 -132 -free $(LIB_MPI) -I/cvos/shared/apps/openmpi/intel/64/1.4/include
#FFLAGS_MPI = -O3 -132 -free -check all -traceback -CB -g -p -inline_debug_info $(LIB_MPI) -I/cvos/shared/apps/openmpi/intel/64/1.4/include
#FFLAGS_SCALAPACK = -O3 -132 $(LIB_SCALAPACK)
#LIB     = -L/opt/intel/Compiler/11.1/072/mkl/lib/em64t -lmkl_lapack -lmkl_core -lmkl_intel_lp64 -lmkl_sequential
#########LIB     = -L/opt/intel/mkl/10.2.5.035/lib/em64t -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack -liomp5
#LIB_MPI = -L/opt/intel/Compiler/11.1/072/mkl/lib/em64t -lmkl_lapack -lmkl_core -lmkl_intel_lp64
#LIB_SCALAPACK = -L/opt/intel/Compiler/11.1/072/mkl/lib/em64t

# gfortran and Macbook
#FFLAGS     = -O3 -framework vecLib $(LIB)
#FFLAGS     = -O3 -Wall -Wextra -pedantic -fcheck=all -Wsurprising -fbacktrace -fimplicit-none -framework vecLib $(LIB)
#FFLAGS     = -O3 -fcheck=all -Wsurprising -fdump-fortran-original -fbacktrace -fimplicit-none -pedantic -framework vecLib $(LIB)
#FLAGS_MPI  = -O3 $(LIB_MPI) -I/usr/local/include
#LIB     = -framework vecLib

###############################################
# NOTES ON COMPILER OPTIONS
###############################################
# -warn alignments
# -warn interfaces
# -warn unused
# -warn usage

SUBROUTINES = main.o \
	abstime.o \
	addatoms.o \
	addvdw.o \
	adjustkalman.o \
	adjustkalman_short.o \
	adjustkalman_elec.o \
	allocatestructures.o \
	analyzecorrmat.o \
	analyzeinput.o \
	analyzeweights.o \
	atomsymfunction1.o \
	atomsymfunction2.o \
	atomsymfunction3.o \
	atomsymfunction3Andi.o \
	atomsymfunction4.o \
	atomsymfunction5.o \
	atomsymfunction6.o \
	atomsymfunction8.o \
	atomsymfunction9.o \
	calcfunctions.o \
	calcrmse_charge.o \
	calcrmse_hextoff.o \
	calcrmse_energy.o \
	calcrmse_forces.o \
	calconecharge.o \
	calconefunction_atomic.o \
	calconenn.o \
	calconefunction_pair.o \
        calconefunction_hextoff_mode1.o \
        calconefunction_hextoff_mode3.o \
	calconeshort.o \
	calconeshortpair.o \
	calconeshort_para.o \
	calconeshort_parapair.o \
	checkelement.o \
	checkextrapolation.o \
	checkextrapolationpair.o \
	checkfiles.o \
	checkfunction.o \
	checkinputnn.o \
	checkonestructure.o \
	checkstructures.o \
	cleanup.o \
	comperrfct.o \
	compileinfo.o \
	coulomb_para.o \
	countpoints.o \
	dataclustering.o \
	dcomperrfct.o \
	dcopy.o \
	ddot.o \
	deallocatestructures.o \
	debugweights.o \
	dgemm.o \
	distribute_fittingoptions.o \
	distribute_globaloptions.o \
	distribute_predictionoptions.o \
	dspmv.o \
	dspr.o \
	dsymm.o \
	electrostatic.o \
	enforcecharge.o \
	environmentanalysis.o \
	erroranalysis.o \
	ewaldenergies_para.o \
	ewaldreal.o \
	ewaldreal_mode3.o \
	ewaldrecip.o \
	ewaldrecip_mode3.o \
	ewaldself.o \
	ewaldself_mode3.o \
	findcontradictions.o \
	fit_density.o \
	fit_hexton.o \
        fit_hextoff.o \
	fit_overlap.o \
	fitting_atomic.o \
	fitting_electrostatic.o \
	fitting_nntb.o \
	fitting_short_atomic.o \
	fitting_short_pair.o \
	generalinfo.o \
	getallelectrostatic.o \
	getallforces.o \
	getallshortforces.o \
	getallshortforcespair.o \
	getatom.o \
	getatom_para.o \
	getatomicoutput_para.o \
	gethextoffoutput.o \
	getatomsymfunctions.o \
	getavcharge.o \
	getcharges.o \
	getchargesatomic.o \
	getcoulombforces.o \
	getcoulombforcesone.o \
	getcoulombdchargedxyz.o \
	getdchargedsfunc.o \
	getdchargedsfunc_para.o \
	getdchargedxyz.o \
	getddnodes_values.o \
	getdeshortdsfunc.o \
	getdeshortdw.o \
	getdepairdsfunc.o \
	getdepairdw.o \
	getdeshortdw_para.o \
	getdfshortdw.o \
	getdfpairdw.o \
	getdfshortdw_para.o \
	getdfpairdw_para.o \
	getdimensions.o \
	getdnnoutdw.o \
	getdnodes_values.o \
	get_dnodes_values_dg.o \
	getenergystatistics.o \
	gethamiltonianstatistics.o \
	geterror_atomic.o \
	geterror_elec.o \
	geterror_hextoff.o \
	geterror_short_atomic.o \
	geterror_short_pair.o \
	geteshort.o \
	geteshortpair.o \
	getewald.o \
	getewaldcorr.o \
	getewaldenergy.o \
	getewaldstress.o \
	getfilenames.o \
	getfilenames_elec.o \
	getfilenamespair.o \
	getfixedweights.o \
        getforces.o \
		getforcescaling.o \
        getforcestatistics.o \
        getglobalscaling.o \
	gethextoff.o \
	getidx.o \
	getkaldims.o \
	getkaldims_elec.o \
	getkaldims_short.o \
        getkaldims_hextoff.o \
	getkalmanmatrices_short.o \
	getkalmanmatrices_elec.o \
        getkalmanmatrices_hextoff.o \
	getlocalscaling.o \
	getlistdim.o \
	getmaxnumneighbors.o \
	getneighboridxatomic_para.o \
	getneighboridxatomic.o \
	getnumatomsall.o \
	getneighborsatomic.o \
	getneighborsatomic_para.o \
	getnumpairs.o \
        getnumtriplets.o \
	getonededb.o \
	getonededw.o \
	getoneshortforce_para.o \
	getoneshortforcepair_para.o \
	getonestructure.o \
	getpointindex.o \
	getreclat.o \
	getrmse.o \
	getrmse_hextoff.o \
	getrmse_elec.o \
	getrmse_short.o \
	getscale.o \
	getscale_hextoff.o \
	getscreenfunctionforelectrostatics.o \
	getsensitivity_atomic.o \
	getsensitivity_short.o \
	getsensitivity_elec.o \
	getshortatomic.o \
	getshortenergies_para.o \
	getshortenergies_parapair.o \
	getshortforces.o \
	getshortforcespair.o \
        getshortforces_parapair.o\
	getstructure_mode3.o \
	getstructures.o \
	getstructures_mixed.o \
	getsymmetryfunctions.o \
	getvdwparams.o \
	getvolume.o \
	getwconstraintidx.o \
	getweightfilenames.o \
	getweightfilenames_elec.o \
	getwshift.o \
	initialcorrmatrix.o \
	initialization.o \
	initializecounters.o \
	initialweights.o \
	initialweights_elec.o \
	initmode3.o \
	initnn.o \
	inputnndefaults.o \
	invert3x3.o \
	lsame.o \
	mode1.o \
	mode2.o \
	mpifitdistribution.o \
	neighbor.o \
	neighbor_para.o \
	neighbor_para_short.o \
	neighborpair_para.o \
        neighbortriplet_para.o \
	nguyenwidrowweights.o \
	nuccharge.o \
	optimize_atomic.o \
	optimize_ewald.o \
	optimize_hextoff.o \
	optimize_short_combined.o \
	optimize_short_combinedpair.o \
	paircount.o \
	pairsymfunction1.o \
	pairsymfunction2.o \
	pairsymfunction3.o \
	pairsymfunction4.o \
	pairsymfunction5.o \
	pairsymfunction6.o \
        hextoffsymfunction1.o \
	precondition_atomic.o \
	precondition_electrostatic.o \
	precondition_short_atomic.o \
	precondition_short_pair.o \
	predict.o \
	sockets.o\
	predict_ipi.o \
	predictelec.o \
	predictionshortatomic.o \
	predictionshortpair.o \
	predictionelectrostatic.o \
	preparemd.o \
	printdateandtime.o \
	printinputnn.o \
	printtimings.o \
	random.o \
	readatomenergies.o \
	readbasis.o \
	readcharges.o \
	readelementlayersatomic.o \
	readelementlayerspair.o \
	readforces.o \
	readforces_mixed.o \
	readfunctions.o \
	readfunctions_hextoff.o \
	readfunctions_mixed.o \
	readfunctions_hextoff_mixed.o \
	readinput.o \
	readkeywords.o \
        readonestructure_ham.o \
	readonestructure.o \
	readotf.o \
	readscale.o \
	readscale_nntb.o \
	readstructures.o \
	readsymfunctionatomic.o \
	readsymfunctionelementatomic.o \
	readsymfunctionglobalatomic.o \
	readsymfunctionpair.o \
	readsymfunctionelementpair.o \
	readsymfunctionglobalpair.o \
        readsymfunctionhextoff.o \
	readweights.o \
	removeatoms.o \
	scaledsfunc.o \
	scaledsfuncpair.o \
	scaledsfunc_para.o \
	scaledsfuncpair_para.o \
	scalesym.o \
	scalesymone.o \
	scalesympair.o \
	scalesym_hextoff.o \
	scalesym_para.o \
        scalesympair_para.o\
	scalesymfit_para.o \
	scalesymfit_parapair.o \
	select.o \
	setglobalactivation.o \
	shuffleweightsshortatomic.o \
	sortint.o \
	sort2realint.o \
	sortchargeerror.o \
	sortelements.o \
	sorteshorterror.o \
	sortforceerror.o \
	sorthextoffsymfunctions.o \
	sortpairsymfunctions.o \
	sortsymfunctions.o \
	splitcoulombforces.o \
	structurecount.o \
	systematicweights.o \
	translate.o \
	updateforces_ed.o \
	updateforces_edv2.o \
	updatekalman.o \
	updatekalman_ed.o \
	updatekalman_para.o \
	updateshort.o \
	updatesteepest.o \
	writeepochtime.o \
	writefitstat.o \
	writefitstat_short.o \
	writefitstat_elec.o \
	writeheader.o \
	writeinputdata.o \
	writekalman.o \
	writekalman_elec.o \
	writekalman_short.o \
	writeoptfit_atomic.o \
	writeoptfit_elec.o \
	writeoptfit_short.o \
	writeoptweights.o \
	writescale.o \
	writescale_nntb.o \
	writeatomicsymfunctions.o \
	writepairsymfunctions.o \
        writehextoffsymfunctions.o \
	writetmpweights.o \
	writeweights.o \
	writepdb.o \
	writepov.o \
	writepwscf.o \
	writexyz.o \
	xerbla.o \
	zerotime.o

MPI_DUMMY_SUBROUTINES = mpi_dummy_routines.o

MODULES      = fileunits.o \
               predictionoptions.o \
               nnflags.o \
               nnconstants.o \
               globaloptions.o \
               basismod.o \
               inputnncounters.o \
               mode1options.o \
               nnshort_atomic.o \
               nnewald.o \
               nnshort_pair.o \
               nnham.o \
               runneripi.o \
               saturation.o \
               structures.o \
               symfunctions.o \
               timings.o \
               fittingoptions.o

MODULES_MPI = mpi_module.o

MODULES_MPI_DUMMY = mpi_dummy.o

# begin make

all:
	cat compile.msg

# compile modules
fileunits.o: fileunits.f90
	rm fileunits.mod ; \
	$(FC) $(FFLAGS_MPI) -c fileunits.f90

fittingoptions.o: fittingoptions.f90
	rm fittingoptions.mod ; \
	$(FC) $(FFLAGS_MPI) -c fittingoptions.f90

nnflags.o: nnflags.f90
	rm nnflags.mod ; \
	$(FC) $(FFLAGS_MPI) -c nnflags.f90

nnconstants.o: nnconstants.f90
	rm nnconstants.mod ; \
	$(FC) $(FFLAGS_MPI) -c nnconstants.f90

globaloptions.o: globaloptions.f90
	rm globaloptions.mod ; \
	$(FC) $(FFLAGS_MPI) -c globaloptions.f90

basismod.o: basismod.f90
	rm basismod.mod ; \
	$(FC) $(FFLAGS_MPI) -c basismod.f90

inputnncounters.o: inputnncounters.f90
	rm inputnncounters.mod ; \
	$(FC) $(FFLAGS_MPI) -c inputnncounters.f90

mode1options.o: mode1options.f90
	rm mode1options.mod ; \
	$(FC) $(FFLAGS_MPI) -c mode1options.f90

nnshort_atomic.o: nnshort_atomic.f90
	rm nnshort_atomic.mod ; \
	$(FC) $(FFLAGS_MPI) -c nnshort_atomic.f90

nnewald.o: nnewald.f90
	rm nnewald.mod ; \
	$(FC) $(FFLAGS_MPI) -c nnewald.f90

nnshort_pair.o: nnshort_pair.f90
	rm nnshort_pair.mod ; \
	$(FC) $(FFLAGS_MPI) -c nnshort_pair.f90

nnham.o: nnham.f90
	rm nnham.mod ; \
	$(FC) $(FFLAGS_MPI) -c nnham.f90

runneripi.o: runneripi.f90
	rm runneripi.mod ; \
	$(FC) $(FFLAGS_MPI) -c runneripi.f90

saturation.o: saturation.f90
	rm saturation.mod ; \
	$(FC) $(FFLAGS_MPI) -c saturation.f90

structures.o: structures.f90
	rm structures.mod ; \
	$(FC) $(FFLAGS_MPI) -c structures.f90

mpi_module.o: mpi_module.f90
	rm mpi_module.mod ; \
	$(MPIFC) $(FFLAGS_MPI) -c mpi_module.f90

symfunctions.o: symfunctions.f90
	rm symfunctions.mod ; \
	$(FC) $(FFLAGS_MPI) -c symfunctions.f90

timings.o: timings.f90
	rm timings.mod ; \
	$(FC) $(FFLAGS_MPI) -c timings.f90

predictionoptions.o: predictionoptions.f90
	rm predicitonoptions.mod ; \
	$(FC) $(FFLAGS_MPI) -c predictionoptions.f90

# compile MPI modules

# compile MPI dummy modules

mpi_dummy.o: mpi_dummy.f90
	rm mpi_dummy.mod ; \
	$(FC) $(FFLAGS) -c mpi_dummy.f90

#        $(MPIFC) $(F90FLAGS) $(CBFLAGS) $(MPIF) -c mpi_module.f90

# compile subroutines in alphabetic order

main.o: main.f90
	$(FC) $(FFLAGS) -c main.f90

abstime.o: abstime.f90
	$(FC) $(FFLAGS) -c abstime.f90

addatoms.o: addatoms.f90
	$(FC) $(FFLAGS) -c addatoms.f90

addvdw.o: addvdw.f90
	$(FC) $(FFLAGS) -c addvdw.f90

adjustkalman.o: adjustkalman.f90
	$(FC) $(FFLAGS) -c adjustkalman.f90

adjustkalman_short.o: adjustkalman_short.f90
	$(FC) $(FFLAGS) -c adjustkalman_short.f90

adjustkalman_elec.o: adjustkalman_elec.f90
	$(FC) $(FFLAGS) -c adjustkalman_elec.f90

allocatestructures.o: allocatestructures.f90
	$(FC) $(FFLAGS) -c allocatestructures.f90

analyzecorrmat.o: analyzecorrmat.f90
	$(FC) $(FFLAGS) -c analyzecorrmat.f90

analyzeinput.o: analyzeinput.f90
	$(FC) $(FFLAGS) -c analyzeinput.f90

analyzeweights.o: analyzeweights.f90
	$(FC) $(FFLAGS) -c analyzeweights.f90

atomsymfunction1.o: atomsymfunction1.f90
	$(FC) $(FFLAGS) -c atomsymfunction1.f90

atomsymfunction2.o: atomsymfunction2.f90
	$(FC) $(FFLAGS) -c atomsymfunction2.f90

atomsymfunction3.o: atomsymfunction3.f90
	$(FC) $(FFLAGS) -c atomsymfunction3.f90

atomsymfunction3Andi.o: atomsymfunction3Andi.f90
	$(FC) $(FFLAGS) -c atomsymfunction3Andi.f90

atomsymfunction4.o: atomsymfunction4.f90
	$(FC) $(FFLAGS) -c atomsymfunction4.f90

atomsymfunction5.o: atomsymfunction5.f90
	$(FC) $(FFLAGS) -c atomsymfunction5.f90

atomsymfunction6.o: atomsymfunction6.f90
	$(FC) $(FFLAGS) -c atomsymfunction6.f90

atomsymfunction8.o: atomsymfunction8.f90
	$(FC) $(FFLAGS) -c atomsymfunction8.f90

atomsymfunction9.o: atomsymfunction9.f90
	$(FC) $(FFLAGS) -c atomsymfunction9.f90

calcfunctions.o: calcfunctions.f90
	$(FC) $(FFLAGS) -c calcfunctions.f90

calcrmse_charge.o: calcrmse_charge.f90
	$(FC) $(FFLAGS) -c calcrmse_charge.f90

calcrmse_hextoff.o: calcrmse_hextoff.f90
	$(FC) $(FFLAGS) -c calcrmse_hextoff.f90

calcrmse_energy.o: calcrmse_energy.f90
	$(FC) $(FFLAGS) -c calcrmse_energy.f90

calcrmse_forces.o: calcrmse_forces.f90
	$(FC) $(FFLAGS) -c calcrmse_forces.f90

calconecharge.o: calconecharge.f90
	$(FC) $(FFLAGS) -c calconecharge.f90

calconefunction_atomic.o: calconefunction_atomic.f90
	$(FC) $(FFLAGS) -c calconefunction_atomic.f90

calconenn.o: calconenn.f90
	$(FC) $(FFLAGS) -c calconenn.f90

calconefunction_pair.o: calconefunction_pair.f90
	$(FC) $(FFLAGS) -c calconefunction_pair.f90

calconefunction_hextoff_mode1.o: calconefunction_hextoff_mode1.f90
	$(FC) $(FFLAGS) -c calconefunction_hextoff_mode1.f90

calconefunction_hextoff_mode3.o: calconefunction_hextoff_mode3.f90
	$(FC) $(FFLAGS) -c calconefunction_hextoff_mode3.f90

calconeshort.o: calconeshort.f90
	$(FC) $(FFLAGS) -c calconeshort.f90

calconeshortpair.o: calconeshortpair.f90
	$(FC) $(FFLAGS) -c calconeshortpair.f90

calconeshort_para.o: calconeshort_para.f90
	$(FC) $(FFLAGS) -c calconeshort_para.f90

calconeshort_parapair.o: calconeshort_parapair.f90
	$(FC) $(FFLAGS) -c calconeshort_parapair.f90

checkelement.o: checkelement.f90
	$(FC) $(FFLAGS) -c checkelement.f90

checkextrapolation.o: checkextrapolation.f90
	$(FC) $(FFLAGS) -c checkextrapolation.f90

checkextrapolationpair.o: checkextrapolationpair.f90
	$(FC) $(FFLAGS) -c checkextrapolationpair.f90

checkfiles.o: checkfiles.f90
	$(FC) $(FFLAGS) -c checkfiles.f90

checkfunction.o: checkfunction.f90
	$(FC) $(FFLAGS) -c checkfunction.f90

checkinputnn.o: checkinputnn.f90
	$(FC) $(FFLAGS) -c checkinputnn.f90

checkonestructure.o: checkonestructure.f90
	$(FC) $(FFLAGS) -c checkonestructure.f90

checkstructures.o: checkstructures.f90
	$(FC) $(FFLAGS) -c checkstructures.f90

cleanup.o: cleanup.f90
	$(FC) $(FFLAGS) -c cleanup.f90

comperrfct.o: comperrfct.f90
	$(FC) $(FFLAGS) -c comperrfct.f90

	rm compileinfo.o;
compileinfo.o: compileinfo.f90
	./compileinfo.csh ; cat compileinfoheader compileinfomiddle compileinfotail > compileinfo.f90 ; rm -f compileinfomiddle ; $(FC) $(FFLAGS) -c compileinfo.f90

coulomb_para.o: coulomb_para.f90
	$(FC) $(FFLAGS) -c coulomb_para.f90

countpoints.o: countpoints.f90
	$(FC) $(FFLAGS) -c countpoints.f90

dataclustering.o: dataclustering.f90
	$(FC) $(FFLAGS) -c dataclustering.f90

dcomperrfct.o: dcomperrfct.f90
	$(FC) $(FFLAGS) -c dcomperrfct.f90

dcopy.o: dcopy.f90
	$(FC) $(FFLAGS) -c dcopy.f90

ddot.o: ddot.f90
	$(FC) $(FFLAGS) -c ddot.f90

deallocatestructures.o: deallocatestructures.f90
	$(FC) $(FFLAGS) -c deallocatestructures.f90

debugweights.o: debugweights.f90
	$(FC) $(FFLAGS) -c debugweights.f90

dgemm.o: dgemm.f90
	$(FC) $(FFLAGS) -c dgemm.f90

distribute_fittingoptions.o: distribute_fittingoptions.f90
	$(FC) $(FFLAGS) -c distribute_fittingoptions.f90

distribute_globaloptions.o: distribute_globaloptions.f90
	$(FC) $(FFLAGS) -c distribute_globaloptions.f90

distribute_predictionoptions.o: distribute_predictionoptions.f90
	$(FC) $(FFLAGS) -c distribute_predictionoptions.f90

dspmv.o: dspmv.f90
	$(FC) $(FFLAGS) -c dspmv.f90

dspr.o: dspr.f90
	$(FC) $(FFLAGS) -c dspr.f90

dsymm.o: dsymm.f90
	$(FC) $(FFLAGS) -c dsymm.f90

electrostatic.o: electrostatic.f90
	$(FC) $(FFLAGS) -c electrostatic.f90

enforcecharge.o: enforcecharge.f90
	$(FC) $(FFLAGS) -c enforcecharge.f90

environmentanalysis.o: environmentanalysis.f90
	$(FC) $(FFLAGS) -c environmentanalysis.f90

erroranalysis.o: erroranalysis.f90
	$(FC) $(FFLAGS) -c erroranalysis.f90

ewaldenergies_para.o: ewaldenergies_para.f90
	$(FC) $(FFLAGS) -c ewaldenergies_para.f90

ewaldreal.o: ewaldreal.f90
	$(FC) $(FFLAGS) -c ewaldreal.f90

ewaldreal_mode3.o: ewaldreal_mode3.f90
	$(FC) $(FFLAGS) -c ewaldreal_mode3.f90

ewaldrecip.o: ewaldrecip.f90
	$(FC) $(FFLAGS) -c ewaldrecip.f90

ewaldrecip_mode3.o: ewaldrecip_mode3.f90
	$(FC) $(FFLAGS) -c ewaldrecip_mode3.f90

ewaldself.o: ewaldself.f90
	$(FC) $(FFLAGS) -c ewaldself.f90

ewaldself_mode3.o: ewaldself_mode3.f90
	$(FC) $(FFLAGS) -c ewaldself_mode3.f90

findcontradictions.o: findcontradictions.f90
	$(FC) $(FFLAGS) -c findcontradictions.f90

fit_density.o: fit_density.f90
	$(FC) $(FFLAGS) -c fit_density.f90

fit_hexton.o: fit_hexton.f90
	$(FC) $(FFLAGS) -c fit_hexton.f90

fit_hextoff.o: fit_hextoff.f90
	$(FC) $(FFLAGS) -c fit_hextoff.f90

fit_overlap.o: fit_overlap.f90
	$(FC) $(FFLAGS) -c fit_overlap.f90

fitting_atomic.o: fitting_atomic.f90
	$(FC) $(FFLAGS) -c fitting_atomic.f90

fitting_electrostatic.o: fitting_electrostatic.f90
	$(FC) $(FFLAGS) -c fitting_electrostatic.f90

fitting_nntb.o: fitting_nntb.f90
	$(FC) $(FFLAGS) -c fitting_nntb.f90

fitting_short_atomic.o: fitting_short_atomic.f90
	$(FC) $(FFLAGS) -c fitting_short_atomic.f90

fitting_short_pair.o: fitting_short_pair.f90
	$(FC) $(FFLAGS) -c fitting_short_pair.f90

generalinfo.o: generalinfo.f90
	$(FC) $(FFLAGS) -cpp -c generalinfo.f90

getallelectrostatic.o: getallelectrostatic.f90
	$(FC) $(FFLAGS) -c getallelectrostatic.f90

getallforces.o: getallforces.f90
	$(FC) $(FFLAGS) -c getallforces.f90

getallshortforces.o: getallshortforces.f90
	$(FC) $(FFLAGS) -c getallshortforces.f90

getallshortforcespair.o: getallshortforcespair.f90
	$(FC) $(FFLAGS) -c getallshortforcespair.f90

getatom.o: getatom.f90
	$(FC) $(FFLAGS) -c getatom.f90

getatom_para.o: getatom_para.f90
	$(FC) $(FFLAGS) -c getatom_para.f90

getatomicoutput_para.o: getatomicoutput_para.f90
	$(FC) $(FFLAGS) -c getatomicoutput_para.f90

gethextoffoutput.o: gethextoffoutput.f90
	$(FC) $(FFLAGS) -c gethextoffoutput.f90

getatomsymfunctions.o: getatomsymfunctions.f90
	$(FC) $(FFLAGS) -c getatomsymfunctions.f90

getavcharge.o: getavcharge.f90
	$(FC) $(FFLAGS) -c getavcharge.f90

getcharges.o: getcharges.f90
	$(FC) $(FFLAGS) -c getcharges.f90

getchargesatomic.o: getchargesatomic.f90
	$(FC) $(FFLAGS) -c getchargesatomic.f90

getcoulombforces.o: getcoulombforces.f90
	$(FC) $(FFLAGS) -c getcoulombforces.f90

getcoulombforcesone.o: getcoulombforcesone.f90
	$(FC) $(FFLAGS) -c getcoulombforcesone.f90

getcoulombdchargedxyz.o: getcoulombdchargedxyz.f90
	$(FC) $(FFLAGS) -c getcoulombdchargedxyz.f90

getdchargedsfunc.o: getdchargedsfunc.f90
	$(FC) $(FFLAGS) -c getdchargedsfunc.f90

getdchargedsfunc_para.o: getdchargedsfunc_para.f90
	$(FC) $(FFLAGS) -c getdchargedsfunc_para.f90

getdchargedxyz.o: getdchargedxyz.f90
	$(FC) $(FFLAGS) -c getdchargedxyz.f90

getddnodes_values.o: getddnodes_values.f90
	$(FC) $(FFLAGS) -c getddnodes_values.f90

getdeshortdsfunc.o: getdeshortdsfunc.f90
	$(FC) $(FFLAGS) -c getdeshortdsfunc.f90

getdeshortdw.o: getdeshortdw.f90
	$(FC) $(FFLAGS) -c getdeshortdw.f90

getdepairdsfunc.o: getdepairdsfunc.f90
	$(FC) $(FFLAGS) -c getdepairdsfunc.f90

getdepairdw.o: getdepairdw.f90
	$(FC) $(FFLAGS) -c getdepairdw.f90

getdeshortdw_para.o: getdeshortdw_para.f90
	$(FC) $(FFLAGS) -c getdeshortdw_para.f90

getdfshortdw.o: getdfshortdw.f90
	$(FC) $(FFLAGS) -c getdfshortdw.f90

getdfpairdw.o: getdfpairdw.f90
	$(FC) $(FFLAGS) -c getdfpairdw.f90

getdfshortdw_para.o: getdfshortdw_para.f90
	$(FC) $(FFLAGS) -c getdfshortdw_para.f90

getdfpairdw_para.o: getdfpairdw_para.f90
	$(FC) $(FFLAGS) -c getdfpairdw_para.f90

getdimensions.o: getdimensions.f90
	$(FC) $(FFLAGS) -c getdimensions.f90

getdnnoutdw.o: getdnnoutdw.f90
	$(FC) $(FFLAGS) -c getdnnoutdw.f90

getdnodes_values.o: getdnodes_values.f90
	$(FC) $(FFLAGS) -c getdnodes_values.f90

get_dnodes_values_dg.o: get_dnodes_values_dg.f90
	$(FC) $(FFLAGS) -c get_dnodes_values_dg.f90

getenergystatistics.o: getenergystatistics.f90
	$(FC) $(FFLAGS) -c getenergystatistics.f90

gethamiltonianstatistics.o: gethamiltonianstatistics.f90
	$(FC) $(FFLAGS) -c gethamiltonianstatistics.f90

geterror_atomic.o: geterror_atomic.f90
	$(FC) $(FFLAGS) -c geterror_atomic.f90

geterror_elec.o: geterror_elec.f90
	$(FC) $(FFLAGS) -c geterror_elec.f90

geterror_hextoff.o: geterror_hextoff.f90
	$(FC) $(FFLAGS) -c geterror_hextoff.f90

geterror_short_atomic.o: geterror_short_atomic.f90
	$(FC) $(FFLAGS) -c geterror_short_atomic.f90

geterror_short_pair.o: geterror_short_pair.f90
	$(FC) $(FFLAGS) -c geterror_short_pair.f90

geteshort.o: geteshort.f90
	$(FC) $(FFLAGS) -c geteshort.f90

geteshortpair.o: geteshortpair.f90
	$(FC) $(FFLAGS) -c geteshortpair.f90

getewald.o: getewald.f90
	$(FC) $(FFLAGS) -c getewald.f90

getewaldcorr.o: getewaldcorr.f90
	$(FC) $(FFLAGS) -c getewaldcorr.f90

getewaldenergy.o: getewaldenergy.f90
	$(FC) $(FFLAGS) -c getewaldenergy.f90

getewaldstress.o: getewaldstress.f90
	$(FC) $(FFLAGS) -c getewaldstress.f90

getfilenames.o: getfilenames.f90
	$(FC) $(FFLAGS) -c getfilenames.f90

getfilenames_elec.o: getfilenames_elec.f90
	$(FC) $(FFLAGS) -c getfilenames_elec.f90

getfilenamespair.o: getfilenamespair.f90
	$(FC) $(FFLAGS) -c getfilenamespair.f90

getfixedweights.o: getfixedweights.f90
	$(FC) $(FFLAGS) -c getfixedweights.f90

getforces.o: getforces.f90
	$(FC) $(FFLAGS) -c getforces.f90

getforcescaling.o: getforcescaling.f90
	$(FC) $(FFLAGS) -c getforcescaling.f90

getforcestatistics.o: getforcestatistics.f90
	$(FC) $(FFLAGS) -c getforcestatistics.f90

getglobalscaling.o: getglobalscaling.f90
	$(FC) $(FFLAGS) -c getglobalscaling.f90

gethextoff.o: gethextoff.f90
	$(FC) $(FFLAGS) -c gethextoff.f90

getidx.o: getidx.f90
	$(FC) $(FFLAGS) -c getidx.f90

getkaldims.o: getkaldims.f90
	$(FC) $(FFLAGS) -c getkaldims.f90

getkaldims_elec.o: getkaldims_elec.f90
	$(FC) $(FFLAGS) -c getkaldims_elec.f90

getkaldims_short.o: getkaldims_short.f90
	$(FC) $(FFLAGS) -c getkaldims_short.f90

getkaldims_hextoff.o: getkaldims_hextoff.f90
	$(FC) $(FFLAGS) -c getkaldims_hextoff.f90

getkalmanmatrices_short.o: getkalmanmatrices_short.f90
	$(FC) $(FFLAGS) -c getkalmanmatrices_short.f90

getkalmanmatrices_elec.o: getkalmanmatrices_elec.f90
	$(FC) $(FFLAGS) -c getkalmanmatrices_elec.f90

getkalmanmatrices_hextoff.o: getkalmanmatrices_hextoff.f90
	$(FC) $(FFLAGS) -c getkalmanmatrices_hextoff.f90

getlistdim.o: getlistdim.f90
	$(FC) $(FFLAGS) -c getlistdim.f90

getlocalscaling.o: getlocalscaling.f90
	$(FC) $(FFLAGS) -c getlocalscaling.f90

getreclat.o: getreclat.f90
	$(FC) $(FFLAGS) -c getreclat.f90

getrmse.o: getrmse.f90
	$(FC) $(FFLAGS) -c getrmse.f90

getrmse_hextoff.o: getrmse_hextoff.f90
	$(FC) $(FFLAGS) -c getrmse_hextoff.f90

getrmse_elec.o: getrmse_elec.f90
	$(FC) $(FFLAGS) -c getrmse_elec.f90

getrmse_short.o: getrmse_short.f90
	$(FC) $(FFLAGS) -c getrmse_short.f90

getmaxnumneighbors.o: getmaxnumneighbors.f90
	$(FC) $(FFLAGS) -c getmaxnumneighbors.f90

getneighboridxatomic_para.o: getneighboridxatomic_para.f90
	$(FC) $(FFLAGS) -c getneighboridxatomic_para.f90

getneighboridxatomic.o: getneighboridxatomic.f90
	$(FC) $(FFLAGS) -c getneighboridxatomic.f90

getnumatomsall.o: getnumatomsall.f90
	$(FC) $(FFLAGS) -c getnumatomsall.f90

getneighborsatomic.o: getneighborsatomic.f90
	$(FC) $(FFLAGS) -c getneighborsatomic.f90

getneighborsatomic_para.o: getneighborsatomic_para.f90
	$(FC) $(FFLAGS) -c getneighborsatomic_para.f90

getnumpairs.o: getnumpairs.f90
	$(FC) $(FFLAGS) -c getnumpairs.f90

getnumtriplets.o: getnumtriplets.f90
	$(FC) $(FFLAGS) -c getnumtriplets.f90

getonededb.o: getonededb.f90
	$(FC) $(FFLAGS) -c getonededb.f90

getonededw.o: getonededw.f90
	$(FC) $(FFLAGS) -c getonededw.f90

getoneshortforce_para.o: getoneshortforce_para.f90
	$(FC) $(FFLAGS) -c getoneshortforce_para.f90

getoneshortforcepair_para.o: getoneshortforcepair_para.f90
	$(FC) $(FFLAGS) -c getoneshortforcepair_para.f90

getonestructure.o: getonestructure.f90
	$(FC) $(FFLAGS) -c getonestructure.f90

getpointindex.o: getpointindex.f90
	$(FC) $(FFLAGS) -c getpointindex.f90

getscale.o: getscale.f90
	$(FC) $(FFLAGS) -c getscale.f90

getscale_hextoff.o: getscale_hextoff.f90
	$(FC) $(FFLAGS) -c getscale_hextoff.f90

getscreenfunctionforelectrostatics.o: getscreenfunctionforelectrostatics.f90
	$(FC) $(FFLAGS) -c getscreenfunctionforelectrostatics.f90

getsensitivity_atomic.o: getsensitivity_atomic.f90
	$(FC) $(FFLAGS) -c getsensitivity_atomic.f90

getsensitivity_short.o: getsensitivity_short.f90
	$(FC) $(FFLAGS) -c getsensitivity_short.f90

getsensitivity_elec.o: getsensitivity_elec.f90
	$(FC) $(FFLAGS) -c getsensitivity_elec.f90

getshortatomic.o: getshortatomic.f90
	$(FC) $(FFLAGS) -c getshortatomic.f90

getshortenergies_para.o: getshortenergies_para.f90
	$(FC) $(FFLAGS) -c getshortenergies_para.f90

getshortenergies_parapair.o: getshortenergies_parapair.f90
	$(FC) $(FFLAGS) -c getshortenergies_parapair.f90

getshortforces.o: getshortforces.f90
	$(FC) $(FFLAGS) -c getshortforces.f90

getshortforcespair.o: getshortforcespair.f90
	$(FC) $(FFLAGS) -c getshortforcespair.f90

getshortforces_parapair.o: getshortforces_parapair.f90
	$(FC) $(FFLAGS) -c getshortforces_parapair.f90

getstructure_mode3.o: getstructure_mode3.f90
	$(FC) $(FFLAGS) -c getstructure_mode3.f90

getstructures.o: getstructures.f90
	$(FC) $(FFLAGS) -c getstructures.f90

getstructures_mixed.o: getstructures_mixed.f90
	$(FC) $(FFLAGS) -c getstructures_mixed.f90

getsymmetryfunctions.o: getsymmetryfunctions.f90
	$(FC) $(FFLAGS) -c getsymmetryfunctions.f90

getvdwparams.o: getvdwparams.f90
	$(FC) $(FFLAGS) -c getvdwparams.f90

getvolume.o: getvolume.f90
	$(FC) $(FFLAGS) -c getvolume.f90

getwconstraintidx.o: getwconstraintidx.f90
	$(FC) $(FFLAGS) -c getwconstraintidx.f90

getweightfilenames.o: getweightfilenames.f90
	$(FC) $(FFLAGS) -c getweightfilenames.f90

getweightfilenames_elec.o: getweightfilenames_elec.f90
	$(FC) $(FFLAGS) -c getweightfilenames_elec.f90

getwshift.o: getwshift.f90
	$(FC) $(FFLAGS) -c getwshift.f90

initialcorrmatrix.o: initialcorrmatrix.f90
	$(FC) $(FFLAGS) -c initialcorrmatrix.f90

initialization.o: initialization.f90
	$(FC) $(FFLAGS) -c initialization.f90

initialweights.o: initialweights.f90
	$(FC) $(FFLAGS) -c initialweights.f90

initialweights_elec.o: initialweights_elec.f90
	$(FC) $(FFLAGS) -c initialweights_elec.f90

initializecounters.o: initializecounters.f90
	$(FC) $(FFLAGS) -c initializecounters.f90

initmode3.o: initmode3.f90
	$(FC) $(FFLAGS) -c initmode3.f90

initnn.o: initnn.f90
	$(FC) $(FFLAGS) -c initnn.f90

inputnndefaults.o: inputnndefaults.f90
	$(FC) $(FFLAGS) -c inputnndefaults.f90

invert3x3.o: invert3x3.f90
	$(FC) $(FFLAGS) -c invert3x3.f90

mpi_dummy_routines.o: mpi_dummy_routines.f90
	$(FC) $(FFLAGS) -c mpi_dummy_routines.f90

lsame.o: lsame.f90
	$(FC) $(FFLAGS) -c lsame.f90

mode1.o: mode1.f90
	$(FC) $(FFLAGS) -c mode1.f90

mode2.o: mode2.f90
	$(FC) $(FFLAGS) -c mode2.f90

mpifitdistribution.o: mpifitdistribution.f90
	$(FC) $(FFLAGS) -c mpifitdistribution.f90

neighbor.o: neighbor.f90
	$(FC) $(FFLAGS) -c neighbor.f90

neighbor_para.o: neighbor_para.f90
	$(FC) $(FFLAGS) -c neighbor_para.f90

neighbor_para_short.o: neighbor_para_short.f90
	$(FC) $(FFLAGS) -c neighbor_para_short.f90

neighborpair_para.o: neighborpair_para.f90
	$(FC) $(FFLAGS) -c neighborpair_para.f90

neighbortriplet_para.o: neighbortriplet_para.f90
	$(FC) $(FFLAGS) -c neighbortriplet_para.f90

nguyenwidrowweights.o: nguyenwidrowweights.f90
	$(FC) $(FFLAGS) -c nguyenwidrowweights.f90

nuccharge.o: nuccharge.f90
	$(FC) $(FFLAGS) -c nuccharge.f90

structurecount.o: structurecount.f90
	$(FC) $(FFLAGS) -c structurecount.f90

optimize_atomic.o: optimize_atomic.f90
	$(FC) $(FFLAGS) -c optimize_atomic.f90

optimize_ewald.o: optimize_ewald.f90
	$(FC) $(FFLAGS) -c optimize_ewald.f90

optimize_hextoff.o: optimize_hextoff.f90
	$(FC) $(FFLAGS) -c optimize_hextoff.f90

optimize_short_combined.o: optimize_short_combined.f90
	$(FC) $(FFLAGS) -c optimize_short_combined.f90

optimize_short_combinedpair.o: optimize_short_combinedpair.f90
	$(FC) $(FFLAGS) -c optimize_short_combinedpair.f90

paircount.o: paircount.f90
	$(FC) $(FFLAGS) -c paircount.f90

pairsymfunction1.o: pairsymfunction1.f90
	$(FC) $(FFLAGS) -c pairsymfunction1.f90

pairsymfunction2.o: pairsymfunction2.f90
	$(FC) $(FFLAGS) -c pairsymfunction2.f90

pairsymfunction3.o: pairsymfunction3.f90
	$(FC) $(FFLAGS) -c pairsymfunction3.f90

pairsymfunction4.o: pairsymfunction4.f90
	$(FC) $(FFLAGS) -c pairsymfunction4.f90

pairsymfunction5.o: pairsymfunction5.f90
	$(FC) $(FFLAGS) -c pairsymfunction5.f90

pairsymfunction6.o: pairsymfunction6.f90
	$(FC) $(FFLAGS) -c pairsymfunction6.f90

hextoffsymfunction1.o: hextoffsymfunction1.f90
	$(FC) $(FFLAGS) -c hextoffsymfunction1.f90

precondition_atomic.o: precondition_atomic.f90
	$(FC) $(FFLAGS) -c precondition_atomic.f90

precondition_electrostatic.o: precondition_electrostatic.f90
	$(FC) $(FFLAGS) -c precondition_electrostatic.f90

precondition_short_atomic.o: precondition_short_atomic.f90
	$(FC) $(FFLAGS) -c precondition_short_atomic.f90

precondition_short_pair.o: precondition_short_pair.f90
	$(FC) $(FFLAGS) -c precondition_short_pair.f90

predict.o: predict.f90
	$(FC) $(FFLAGS) -c predict.f90

sockets.o: sockets.c
	$(CC) $(CFLAGS) -c -o sockets.o sockets.c 

predict_ipi.o: predict_ipi.f90
	$(FC) $(FFLAGS) -c predict_ipi.f90

predictelec.o: predictelec.f90
	$(FC) $(FFLAGS) -c predictelec.f90

predictionshortatomic.o: predictionshortatomic.f90
	$(FC) $(FFLAGS) -c predictionshortatomic.f90

predictionshortpair.o: predictionshortpair.f90
	$(FC) $(FFLAGS) -c predictionshortpair.f90

predictionelectrostatic.o: predictionelectrostatic.f90
	$(FC) $(FFLAGS) -c predictionelectrostatic.f90

preparemd.o: preparemd.f90
	$(FC) $(FFLAGS) -c preparemd.f90

printdateandtime.o: printdateandtime.f90
	$(FC) $(FFLAGS) -c printdateandtime.f90

printinputnn.o: printinputnn.f90
	$(FC) $(FFLAGS) -c printinputnn.f90

printtimings.o: printtimings.f90
	$(FC) $(FFLAGS) -c printtimings.f90

random.o: random.f90
	$(FC) $(FFLAGS) -c random.f90

readatomenergies.o: readatomenergies.f90
	$(FC) $(FFLAGS) -c readatomenergies.f90

readbasis.o: readbasis.f90
	$(FC) $(FFLAGS) -c readbasis.f90

readcharges.o: readcharges.f90
	$(FC) $(FFLAGS) -c readcharges.f90

readelementlayersatomic.o: readelementlayersatomic.f90
	$(FC) $(FFLAGS) -c readelementlayersatomic.f90

readelementlayerspair.o: readelementlayerspair.f90
	$(FC) $(FFLAGS) -c readelementlayerspair.f90

readforces.o: readforces.f90
	$(FC) $(FFLAGS) -c readforces.f90

readforces_mixed.o: readforces_mixed.f90
	$(FC) $(FFLAGS) -c readforces_mixed.f90

readfunctions.o: readfunctions.f90
	$(FC) $(FFLAGS) -c readfunctions.f90

readfunctions_hextoff.o: readfunctions_hextoff.f90
	$(FC) $(FFLAGS) -c readfunctions_hextoff.f90

readfunctions_mixed.o: readfunctions_mixed.f90
	$(FC) $(FFLAGS) -c readfunctions_mixed.f90

readfunctions_hextoff_mixed.o: readfunctions_hextoff_mixed.f90
	$(FC) $(FFLAGS) -c readfunctions_hextoff_mixed.f90

readinput.o: readinput.f90
	$(FC) $(FFLAGS) -c readinput.f90

readkeywords.o: readkeywords.f90
	$(FC) $(FFLAGS) -c readkeywords.f90

readonestructure_ham.o: readonestructure_ham.f90
	$(FC) $(FFLAGS) -c readonestructure_ham.f90

readonestructure.o: readonestructure.f90
	$(FC) $(FFLAGS) -c readonestructure.f90

readotf.o: readotf.f90
	$(FC) $(FFLAGS) -c readotf.f90

readscale.o: readscale.f90
	$(FC) $(FFLAGS) -c readscale.f90

readscale_nntb.o: readscale_nntb.f90
	$(FC) $(FFLAGS) -c readscale_nntb.f90

readstructures.o: readstructures.f90
	$(FC) $(FFLAGS) -c readstructures.f90

readsymfunctionatomic.o: readsymfunctionatomic.f90
	$(FC) $(FFLAGS) -c readsymfunctionatomic.f90

readsymfunctionelementatomic.o: readsymfunctionelementatomic.f90
	$(FC) $(FFLAGS) -c readsymfunctionelementatomic.f90

readsymfunctionglobalatomic.o: readsymfunctionglobalatomic.f90
	$(FC) $(FFLAGS) -c readsymfunctionglobalatomic.f90

readsymfunctionpair.o: readsymfunctionpair.f90
	$(FC) $(FFLAGS) -c readsymfunctionpair.f90

readsymfunctionelementpair.o: readsymfunctionelementpair.f90
	$(FC) $(FFLAGS) -c readsymfunctionelementpair.f90

readsymfunctionglobalpair.o: readsymfunctionglobalpair.f90
	$(FC) $(FFLAGS) -c readsymfunctionglobalpair.f90

readsymfunctionhextoff.o: readsymfunctionhextoff.f90
	$(FC) $(FFLAGS) -c readsymfunctionhextoff.f90

readweights.o: readweights.f90
	$(FC) $(FFLAGS) -c readweights.f90

removeatoms.o: removeatoms.f90
	$(FC) $(FFLAGS) -c removeatoms.f90

scaledsfunc.o: scaledsfunc.f90
	$(FC) $(FFLAGS) -c scaledsfunc.f90

scaledsfuncpair.o: scaledsfuncpair.f90
	$(FC) $(FFLAGS) -c scaledsfuncpair.f90

scaledsfunc_para.o: scaledsfunc_para.f90
	$(FC) $(FFLAGS) -c scaledsfunc_para.f90

scaledsfuncpair_para.o: scaledsfuncpair_para.f90
	$(FC) $(FFLAGS) -c scaledsfuncpair_para.f90

scalesym.o: scalesym.f90
	$(FC) $(FFLAGS) -c scalesym.f90

scalesymone.o: scalesymone.f90
	$(FC) $(FFLAGS) -c scalesymone.f90

scalesympair.o: scalesympair.f90
	$(FC) $(FFLAGS) -c scalesympair.f90

scalesym_hextoff.o: scalesym_hextoff.f90
	$(FC) $(FFLAGS) -c scalesym_hextoff.f90

scalesym_para.o: scalesym_para.f90
	$(FC) $(FFLAGS) -c scalesym_para.f90

scalesympair_para.o: scalesympair_para.f90
	$(FC) $(FFLAGS) -c scalesympair_para.f90


scalesymfit_para.o: scalesymfit_para.f90
	$(FC) $(FFLAGS) -c scalesymfit_para.f90

scalesymfit_parapair.o: scalesymfit_parapair.f90
	$(FC) $(FFLAGS) -c scalesymfit_parapair.f90

select.o: select.f90
	$(FC) $(FFLAGS) -c select.f90

setglobalactivation.o: setglobalactivation.f90
	$(FC) $(FFLAGS) -c setglobalactivation.f90

shuffleweightsshortatomic.o: shuffleweightsshortatomic.f90
	$(FC) $(FFLAGS) -c shuffleweightsshortatomic.f90

sortint.o: sortint.f90
	$(FC) $(FFLAGS) -c sortint.f90

sort2realint.o: sort2realint.f90
	$(FC) $(FFLAGS) -c sort2realint.f90

sortchargeerror.o: sortchargeerror.f90
	$(FC) $(FFLAGS) -c sortchargeerror.f90

sortelements.o: sortelements.f90
	$(FC) $(FFLAGS) -c sortelements.f90

sorteshorterror.o: sorteshorterror.f90
	$(FC) $(FFLAGS) -c sorteshorterror.f90

sortforceerror.o: sortforceerror.f90
	$(FC) $(FFLAGS) -c sortforceerror.f90

sorthextoffsymfunctions.o: sorthextoffsymfunctions.f90
	$(FC) $(FFLAGS) -c sorthextoffsymfunctions.f90

sortpairsymfunctions.o: sortpairsymfunctions.f90
	$(FC) $(FFLAGS) -c sortpairsymfunctions.f90

sortsymfunctions.o: sortsymfunctions.f90
	$(FC) $(FFLAGS) -c sortsymfunctions.f90

splitcoulombforces.o: splitcoulombforces.f90
	$(FC) $(FFLAGS) -c splitcoulombforces.f90

systematicweights.o: systematicweights.f90
	$(FC) $(FFLAGS) -c systematicweights.f90

translate.o: translate.f90
	$(FC) $(FFLAGS) -c translate.f90

updateforces_ed.o: updateforces_ed.f90
	$(FC) $(FFLAGS) -c updateforces_ed.f90

updateforces_edv2.o: updateforces_edv2.f90
	$(FC) $(FFLAGS) -c updateforces_edv2.f90

updatekalman.o: updatekalman.f90
	$(FC) $(FFLAGS) -c updatekalman.f90

updatekalman_ed.o: updatekalman_ed.f90
	$(FC) $(FFLAGS) -c updatekalman_ed.f90

updatekalman_para.o: updatekalman_para.f90
	$(FC) $(FFLAGS) -c updatekalman_para.f90

updateshort.o: updateshort.f90
	$(FC) $(FFLAGS) -c updateshort.f90

updatesteepest.o: updatesteepest.f90
	$(FC) $(FFLAGS) -c updatesteepest.f90

writeepochtime.o: writeepochtime.f90
	$(FC) $(FFLAGS) -c writeepochtime.f90

writefitstat.o: writefitstat.f90
	$(FC) $(FFLAGS) -c writefitstat.f90

writefitstat_short.o: writefitstat_short.f90
	$(FC) $(FFLAGS) -c writefitstat_short.f90

writefitstat_elec.o: writefitstat_elec.f90
	$(FC) $(FFLAGS) -c writefitstat_elec.f90

writeheader.o: writeheader.f90
	$(FC) $(FFLAGS) -c writeheader.f90

writeinputdata.o: writeinputdata.f90
	$(FC) $(FFLAGS) -c writeinputdata.f90

writekalman.o: writekalman.f90
	$(FC) $(FFLAGS) -c writekalman.f90

writekalman_elec.o: writekalman_elec.f90
	$(FC) $(FFLAGS) -c writekalman_elec.f90

writekalman_short.o: writekalman_short.f90
	$(FC) $(FFLAGS) -c writekalman_short.f90

writeoptfit_atomic.o: writeoptfit_atomic.f90
	$(FC) $(FFLAGS) -c writeoptfit_atomic.f90

writeoptfit_elec.o: writeoptfit_elec.f90
	$(FC) $(FFLAGS) -c writeoptfit_elec.f90

writeoptfit_short.o: writeoptfit_short.f90
	$(FC) $(FFLAGS) -c writeoptfit_short.f90

writeoptweights.o: writeoptweights.f90
	$(FC) $(FFLAGS) -c writeoptweights.f90

writescale.o: writescale.f90
	$(FC) $(FFLAGS) -c writescale.f90

writescale_nntb.o: writescale_nntb.f90
	$(FC) $(FFLAGS) -c writescale_nntb.f90

writeatomicsymfunctions.o: writeatomicsymfunctions.f90
	$(FC) $(FFLAGS) -c writeatomicsymfunctions.f90

writepairsymfunctions.o: writepairsymfunctions.f90
	$(FC) $(FFLAGS) -c writepairsymfunctions.f90

writehextoffsymfunctions.o: writehextoffsymfunctions.f90
	$(FC) $(FFLAGS) -c writehextoffsymfunctions.f90

writetmpweights.o: writetmpweights.f90
	$(FC) $(FFLAGS) -c writetmpweights.f90

writeweights.o: writeweights.f90
	$(FC) $(FFLAGS) -c writeweights.f90

writepdb.o: writepdb.f90
	$(FC) $(FFLAGS) -c writepdb.f90

writepov.o: writepov.f90
	$(FC) $(FFLAGS) -c writepov.f90

writepwscf.o: writepwscf.f90
	$(FC) $(FFLAGS) -c writepwscf.f90

writexyz.o: writexyz.f90
	$(FC) $(FFLAGS) -c writexyz.f90

xerbla.o: xerbla.f90
	$(FC) $(FFLAGS) -c xerbla.f90

zerotime.o: zerotime.f90
	$(FC) $(FFLAGS) -c zerotime.f90

serial: ./RuNNer.serial.x
./RuNNer.serial.x: $(MODULES_MPI_DUMMY) $(MODULES) $(SUBROUTINES) $(MPI_DUMMY_SUBROUTINES)
	$(FC) $(FFLAGS) -o $@ \
               $(SUBROUTINES) $(MODULES_MPI_DUMMY) $(MODULES) $(MPI_DUMMY_SUBROUTINES) $(LIB) ; rm -f compileinfo.o

mpi: ./RuNNer.mpi.x
./RuNNer.mpi.x: $(MODULES_MPI) $(MODULES) $(SUBROUTINES)
	$(MPIFC) $(FFLAGS_MPI) -o $@ \
               $(SUBROUTINES) $(MODULES_MPI) $(MODULES) $(LIB) ; rm -f compileinfo.o

scalapack: ./RuNNer.scalapack.x
./RuNNer.scalapack.x: $(MODULES) $(MODULES_MPI) $(SUBROUTINES)
	$(MPIFC) $(FFLAGS_SCALAPACK) -o $@ \
               $(SUBROUTINES) $(MODULES) $(MODULES_MPI) ; rm -f compileinfo.o

clean:
	rm -f *.o *.mod
