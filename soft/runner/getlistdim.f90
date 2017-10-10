!######################################################################                                                                                                                                              
! This routine is part of                                                                                                                                                                                            
! RuNNer - Ruhr-University Neural Network Energy Representation                                                                                                                                                      
! (c) Dr. Joerg Behler 2008                                                                                                                                                                                          
!######################################################################                                                                                                                                              

!! called by: 
!! - main.f90
!!
      subroutine getlistdim()
!!
      use nnflags
      use globaloptions
      use fileunits
!!
      implicit none
!!
      integer listdimtemp1
      integer listdimtemp2
      integer listdimtemp3
      integer listdimtemp4
!!
!! set a more clever value for listdim here
      listdim=0
      if(lshort.and.(nn_type_short.eq.1))then
        listdimtemp1=800*max_num_atoms
!! make sure here that in parallel cases neighbor arrays are not larger than necessary
        if(mode.eq.3)then
          listdimtemp1=200*nblock
        endif
        listdim=max(listdim,listdimtemp1)
      endif
      if(lshort.and.(nn_type_short.eq.2))then
        listdimtemp2=400*max_num_pairs
!! make sure here that in parallel cases neighbor arrays are not larger than necessary
        if(mode.eq.3)then
          listdimtemp2=200*nblock
        endif
        listdim=max(listdim,listdimtemp2)
      endif
!! separate electrostatic NN:
!!      if(lelec.and.(nn_type_elec.eq.1))then ! this does not work because also for fixed charges the Ewald sum needs neighbor lists
      if(lelec)then
        listdimtemp3=400*max_num_atoms
!! make sure here that in parallel cases neighbor arrays are not larger than necessary
        if(mode.eq.3)then
          listdimtemp3=200*nblock
        endif
        listdim=max(listdim,listdimtemp3)
      endif
!! Hamiltonian NN:
      if(lnntb)then
        listdimtemp4=200*max_num_triplets
        if(mode.eq.3)then
          listdimtemp4=200*nblock
        endif
        listdim=max(listdim,listdimtemp4)
      endif
!! 
!! TODO: lsta can be reduced in size for prediction mode (at least)
!!    if(mode.eq.3)then
!!      lstadim=ceiling(dble(nblock) / dble(mpisize))
!!    else
!!      lstadim=max_num_atoms
!!    endif

      return
      end
