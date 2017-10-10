!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: 
!! - fitting.f90
!! - fittingpair.f90
!!
      subroutine readforces_mixed(npoints,&
        ntrain,block_counter,pointindex,num_atoms_all)
!!
      use fileunits
      use globaloptions
      use structures
!!
      implicit none
!!
      integer ntrain                                             ! in
      integer block_counter                                      ! in
      integer pointindex(ntrain)                                 ! in
      integer num_atoms_all(ntrain)                              ! in
      integer npoints                                            ! in
      integer idummy                                             ! internal
      integer i1,i2,i3                                           ! internal
      integer jcount                                             ! internal
!!
      real*8 forces_local(3,max_num_atoms)                       ! internal
      real*8 rdummy                                              ! internal
!!    
!!      write(ounit,*)'readforces_mixed npoints ',npoints,block_counter
!!      do i1=block_counter,block_counter+npoints-1
!!        write(ounit,'(a,2i6)')'JB pointindex ',i1,pointindex(i1)
!!      enddo
!!   
!      totalforce_list(:,:,:)=0.0d0
      shortforce_list(:,:,:)=0.0d0
!!
      open(trainfunit,file='trainforces.data',form='formatted',status='old')
      rewind(trainfunit)
!!'
!! here we read only a subset of strutures
      do i1=1,ntrain
!! read forces of one structure into temporary array
        read(trainfunit,*)idummy
        do i2=1,num_atoms_all(i1)
          read(trainfunit,*)(forces_local(i3,i2),i3=1,3)
        enddo ! i2
        jcount=block_counter
        do i2=1,npoints
!!          write(ounit,*)'JBlooking for ',i2,pointindex(jcount)
          if(i1.eq.pointindex(jcount))then
!!            write(ounit,'(a,i5,a,i5,a,i5)')'JBreadforces structure ',i1,' is number ',jcount,' and goes to field ',i2
!! copy the forces of this structure to the shortforce_list array
            do i3=1,num_atoms_all(i1)
              shortforce_list(:,i3,i2)=forces_local(:,i3)
            enddo ! i3
          endif
          jcount=jcount+1
        enddo ! i2
      enddo ! i1
!!
      close(trainfunit)
!!
      return
      end


