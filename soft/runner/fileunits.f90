      module fileunits 

      implicit none

!! input files
      integer, parameter :: nnunit          = 20 ! input.nn
      integer, parameter :: dataunit        = 21 ! input.data

!! input and output files
      integer, parameter :: scaleunit       = 22 ! scaling.data, scalinge.data 
      integer, parameter :: runit           = 23 ! trainforces.data, testforces.data
      integer, parameter :: trainfunit      = 24 ! trainforces.data
      integer, parameter :: testfunit       = 25 ! testforces.data
      integer, parameter :: trainfeunit     = 26 ! trainforcese.data
      integer, parameter :: testfeunit      = 27 ! testforcese.data

      integer, parameter :: symunit         = 28 ! function.data, testing.data 
      integer, parameter :: tymunit         = 29 ! testing.data 
      integer, parameter :: symeunit        = 30 ! functione.data, testinge.data 
      integer, parameter :: tymeunit        = 31 ! testinge.data 
      integer, parameter :: symdensunit     = 60 ! functiond.data, testingd.data 
      integer, parameter :: tymdensunit     = 61 ! testingd.data 
      integer, parameter :: symhextoffunit  = 62 ! functionx.data, testingx.data 
      integer, parameter :: tymhextoffunit  = 63 ! testingx.data 
      integer, parameter :: symoverunit     = 64 ! functiono.data, testingo.data 
      integer, parameter :: tymoverunit     = 65 ! testingo.data 
      integer, parameter :: symhextonunit   = 66 ! functionx.data, testingx.data 
      integer, parameter :: tymhextonunit   = 67 ! testingx.data 
      integer, parameter :: symhamunit      = 68 ! functionx.data, testingx.data
      integer, parameter :: tymhamunit      = 69 ! testingx.data 
      integer, parameter :: otfunit         = 70 ! input.otf 

      integer, parameter :: teststructunit  = 32 ! teststruct.data 
      integer, parameter :: trainstructunit = 33 ! trainstruct.data 

!! output files
      integer, parameter :: ounit           = 6 ! standard out (or file runner.out, if ounit.ne.6)
      integer, parameter :: pdbunit         = 34 ! structure.pdb
      integer, parameter :: xyzunit         = 35 ! structure.xyz
      integer, parameter :: pwunit          = 36 ! structure.pw
      integer, parameter :: povunit         = 37 ! povray.ini, input.pov 
      integer, parameter :: nnaunit         = 38 ! nnatoms.out 
      integer, parameter :: nnpunit         = 39 ! nnpairs.out 
      integer, parameter :: outunit         = 40 ! output.data 
      integer, parameter :: nnsunit         = 41 ! nnstress.out 
      integer, parameter :: nnfunit         = 42 ! nnforces.out 
      integer, parameter :: nneunit         = 43 ! energy.out 
      integer, parameter :: debugunit       = 44 ! debug.out 
!! element-specific files
      integer, parameter :: wunit           = 45 ! weights.XXX.data, weightse.XXX.data, weights.XXX.XXX.data
      integer, parameter :: woutunit        = 46 ! XXXXXX.short.XXX.out, XXXXXX.ewald.XXX.out, XXXXXX.weights.XXX.XXX.out
      integer, parameter :: wtmpunit        = 47 ! tmpweights.000.out, tmpweightse.XXX.out
      integer, parameter :: kalunit         = 48 ! kalman.short.XXX.data or kalman.pair.XXX.XXX.data 
      integer, parameter :: kaleunit        = 49 ! kalman.elec.XXX.data 
      integer, parameter :: kalcunit        = 50 ! kalmanc.data 
!! epoch-specific output files
      integer, parameter :: trainpxunit     = 51 ! trainpoints.XXXXXX.out
      integer, parameter :: testpxunit      = 52 ! testpoints.XXXXXX.out
      integer, parameter :: trainqxunit     = 53 ! traincharges.XXXXXX.out 
      integer, parameter :: testqxunit      = 54 ! testcharges.XXXXXX.out 
      integer, parameter :: trainfxunit     = 55 ! trainforces.XXXXXX.out
      integer, parameter :: testfxunit      = 56 ! testforces.XXXXXX.out

      integer, parameter :: nnmdunit        = 57 ! nnmd.in 
      integer, parameter :: chargeunit      = 58 ! charges.in 


      end module fileunits 

