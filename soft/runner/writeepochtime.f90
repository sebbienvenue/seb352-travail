!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fitting.f90
!!
      subroutine writeepochtime(countepoch)
!!
      use fileunits
      use nnflags
      use globaloptions
      use fittingoptions
      use timings
!!
      implicit none
!!
      integer countepoch
!!
      write(ounit,'(a42,i5,x,f8.2)')' EPOCHTIME full epoch                ',countepoch,timeepoch/60.d0
      write(ounit,'(a42,i5,x,f8.2)')'   EPOCHTIME geterror                ',countepoch,timegeterror/60.d0
      write(ounit,'(a42,i5,x,f8.2)')'   EPOCHTIME mixing points           ',countepoch,timemix/60.d0
      timemix       = 0.0d0
      write(ounit,'(a42,i5,x,f8.2)')'   EPOCHTIME IO                      ',countepoch,timeio/60.d0
      timeio        = 0.0d0
      if(lshort)then
        write(ounit,'(a42,i5,x,f8.2)')'   EPOCHTIME short range fit         ',countepoch,timeshortfit/60.d0
        timeshortfit  = 0.0d0
        write(ounit,'(a42,i5,x,f8.2)')'     EPOCHTIME energy fit            ',countepoch,timeefitting/60.d0
        timeefitting  = 0.0d0
        write(ounit,'(a42,i5,x,f8.2)')'       EPOCHTIME energy error        ',countepoch,timeeerror/60.d0
        timeeerror = 0.0d0
        write(ounit,'(a42,i5,x,f8.2)')'       EPOCHTIME deshortdw           ',countepoch,timedeshortdw/60.d0
        timedeshortdw = 0.0d0
        write(ounit,'(a42,i5,x,f8.2)')'       EPOCHTIME energyupdate        ',countepoch,timeeupdate/60.d0
        timeeupdate   = 0.0d0
        if(luseforces)then
          write(ounit,'(a42,i5,x,f8.2)')'     EPOCHTIME force fit             ',countepoch,timeffitting/60.d0
          timeffitting  = 0.0d0
          write(ounit,'(a42,i5,x,f8.2)')'       EPOCHTIME force symfunc       ',countepoch,timefsym/60.d0
          timefsym      = 0.0d0
          write(ounit,'(a42,i5,x,f8.2)')'       EPOCHTIME forceerror          ',countepoch,timeferror/60.d0
          timeferror    = 0.0d0
          write(ounit,'(a42,i5,x,f8.2)')'       EPOCHTIME dfshortdw           ',countepoch,timedfshortdw/60.d0
          timedfshortdw = 0.0d0
          write(ounit,'(a42,i5,x,f8.2)')'       EPOCHTIME forceupdate         ',countepoch,timefupdate/60.d0
          timefupdate   = 0.0d0
          if(lrepeate)then
            write(ounit,'(a42,i5,x,f8.2)')'     EPOCHTIME repeated energy fit   ',countepoch,timeefittingrepeat/60.d0
            timeefittingrepeat = 0.0d0
            write(ounit,'(a42,i5,x,f8.2)')'       EPOCHTIME energy error repeat ',countepoch,timeeerrorrepeat/60.d0
            timeeerrorrepeat = 0.0d0
            write(ounit,'(a42,i5,x,f8.2)')'       EPOCHTIME deshortdw repeat    ',countepoch,timedeshortdwrepeat/60.d0
            timedeshortdwrepeat = 0.0d0
            write(ounit,'(a42,i5,x,f8.2)')'       EPOCHTIME energyupdate repeat ',countepoch,timeeupdaterepeat/60.d0
            timeeupdaterepeat = 0.0d0
          endif ! lrepeate
        endif ! luseforces
      endif ! lshort
      if(lelec.and.(nn_type_elec.eq.1))then
        write(ounit,'(a42,i5,x,f8.2)')'   EPOCHTIME electrostatic fit       ',countepoch,timeelecfit/60.d0
        timeelecfit  = 0.0d0
        write(ounit,'(a42,i5,x,f8.2)')'     EPOCHTIME charge error          ',countepoch,timeqerror/60.d0
        timeqerror = 0.0d0
        write(ounit,'(a42,i5,x,f8.2)')'     EPOCHTIME dqdw                  ',countepoch,timedqdw/60.d0
        timedqdw = 0.0d0
        write(ounit,'(a42,i5,x,f8.2)')'     EPOCHTIME charge update         ',countepoch,timequpdate/60.d0
        timequpdate   = 0.0d0
      endif ! lelec
!! 
      return
!!
      end
