!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: - electrostatic.f90  
!!            - electrostatic_para.f90   
!!            - getcoulombforces.f90
!!            - getcoulombforcesone.f90
!!            - getcoulombforces_para.f90
!!            - splitcoulombforces.f90
!!
      subroutine getscreenfunctionforelectrostatics(rij,&
                    fscreen,fscreenderiv,drijdxyz)
!!
      use globaloptions
!!
      implicit none
!!
      real*8, intent(in)        :: rij,drijdxyz
      real*8, intent(out)       :: fscreen, fscreenderiv
      real*8, parameter         :: pi = 3.141592654d0
!!
      if(rij.gt.rscreen_cut) then      
        fscreen=1.0d0
        fscreenderiv=0.0d0
      elseif(rij.lt.rscreen_onset) then      
        fscreen=0.0d0
        fscreenderiv=0.0d0
      else
        fscreen=0.5d0*(1.0d0-dcos(pi*(rij-rscreen_onset)&
          /(rscreen_cut-rscreen_onset)))
        fscreenderiv=0.5d0*dsin(pi*(rij-rscreen_onset)/(rscreen_cut-rscreen_onset))* &
                     pi/(rscreen_cut-rscreen_onset)*drijdxyz
      endif
!!
      return
      end
!!
