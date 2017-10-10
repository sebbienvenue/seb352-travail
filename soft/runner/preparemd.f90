!######################################################################                                                                                                                                
! This routine is part of                                                                                                                                                                              
! RuNNer - Ruhr-University Neural Network Energy Representation                                                                                                                                        
! (c) Dr. Joerg Behler 2008                                                                                                                                                                            
!######################################################################                                                                                                                                
!!                                                                    
!! called by:                                                                                                                                                                                          
!! - predict.f90
!!
      subroutine preparemd()
!!
      use nnflags
      use globaloptions
      use fileunits
      use predictionoptions
!!
      implicit none
!!
      write(ounit,*)'ERROR: Writing nnmd.in is not yet implemented'
      stop
!!
      open(nnmdunit,file='nnmd.in',form='formatted',status='replace')
      rewind(nnmdunit)
!!
      write(nnmdunit)nn_type_short
      write(nnmdunit)lshort
      write(nnmdunit)lelec
!!
!! write short range NN parameters, atomic case
      if(lshort.and.(nn_type_short.eq.1))then

      endif
!!
!! write short range NN parameters, pair case
      if(lshort.and.(nn_type_short.eq.2))then

      endif
!!
      if(lelec.and.(nn_type_elec.eq.1))then

      endif
!!
!!
      close(nnmdunit)
!!
      return
      end 
