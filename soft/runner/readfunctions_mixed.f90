!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine

!! Purpose:
!! get the information for the arrays
!! - symfunction_list
!! - zelem_list
!! - zelemp_list
!! - num_atoms_list
!! - num_pairs_list
!! - totalcharge_list
!! - totalenergy_list
!! - shortenergy_list
!! - elecenergy_list
!! for a block on npoints structures, the data is passed to RuNNer via module structures

!! called by:
!!
      subroutine readfunctions_mixed(npoints,ndim1,ndim2,&
         ntrain,block_counter,pointindex,iswitch,&
         maxnum_funcvalues_local,num_funcvalues_local,&
         symfunction_list_local)
!!
      use fileunits
      use globaloptions
      use structures
!!
      implicit none
!!     
      integer ndim1                                                 ! in (nelem or npairs)
      integer ndim2                                                 ! in (max_num_atoms or max_num_pairs)
      integer npoints                                               ! in
      integer maxnum_funcvalues_local                               ! in
      integer num_funcvalues_local(ndim1)                           ! in
      integer i1,i2,i3                                              ! internal
      integer ntrain                                                ! in
      integer iswitch                                               ! in
      integer block_counter                                         ! in
      integer pointindex(ntrain)                                    ! in
      integer jcount                                                ! internal
      integer idummy                                                ! internal
      integer num_atoms_local                                       ! internal
      integer num_pairs_local                                       ! internal
      integer zelem_local(max_num_atoms)                            ! internal
      integer zelemp_local(2,max_num_pairs)                           ! internal
!!
      real*8 symfunction_list_local(maxnum_funcvalues_local,ndim2,nblock)  ! out
      real*8 symfunction_local(maxnum_funcvalues_local,ndim2)       ! internal 
      real*8 totalcharge_local                                      ! internal
      real*8 totalenergy_local                                      ! internal
      real*8 shortenergy_local                                      ! internal
      real*8 elecenergy_local                                       ! internal
!!
!!
!!      write(ounit,*)'readfunctions_mixed npoints ',npoints,block_counter
!!
      if(iswitch.eq.0)then ! atom case
        open(symunit,file='function.data',form='formatted',status='old')
      elseif(iswitch.eq.1)then
        open(symunit,file='functione.data',form='formatted',status='old')
      elseif(iswitch.eq.2)then ! pair case
        open(symunit,file='function.data',form='formatted',status='old')
      else
        write(ounit,*)'ERROR: unknown iswitch in readfunctions_mixed'
        stop
      endif
      rewind(symunit)
!!
!! here we read only a subset of strutures
      do i1=1,ntrain
!! read symfunctions and zelem of one structure into temporary array
        if((iswitch.eq.0).or.(iswitch.eq.1))then ! atomic case
          read(symunit,*)num_atoms_local
          do i2=1,num_atoms_local
            read(symunit,*)zelem_local(i2)
            backspace(symunit)
            read(symunit,*)idummy,&
              (symfunction_local(i3,i2),i3=1,num_funcvalues_local(elementindex(zelem_local(i2))))
          enddo ! i2
        elseif(iswitch.eq.2)then ! pair case
          read(symunit,*)num_atoms_local,num_pairs_local
          do i2=1,num_pairs_local
            read(symunit,*)zelemp_local(1,i2),zelemp_local(2,i2)
            backspace(symunit)
            read(symunit,*)idummy,idummy,&
              (symfunction_local(i3,i2),i3=1,&
                num_funcvalues_local(pairindex(zelemp_local(1,i2),zelemp_local(2,i2))))
          enddo ! i2
        endif ! iswitch
        read(symunit,*) totalcharge_local,totalenergy_local,&
                      shortenergy_local,elecenergy_local 
!! check if we need this structure now
!! loop over all pointindex values of this block of data
        jcount=block_counter
        do i2=1,npoints
          if(i1.eq.pointindex(jcount))then
!!            write(ounit,'(a,i5,a,i5,a,i5)')'JBreadfunctions structure ',i1,' is number ',jcount,' and goes to field ',i2
            if((iswitch.eq.0).or.(iswitch.eq.1))then
              num_atoms_list(i2)=num_atoms_local
              zelem_list(i2,:)=zelem_local(:)
              symfunction_list_local(:,:,i2)=symfunction_local(:,:)
            elseif(iswitch.eq.2)then ! pair case
              num_atoms_list(i2)=num_atoms_local
              num_pairs_list(i2)=num_pairs_local
              zelemp_list(1,i2,:)=zelemp_local(1,:)
              zelemp_list(2,i2,:)=zelemp_local(2,:)
              symfunction_list_local(:,:,i2)=symfunction_local(:,:)
            endif ! iswitch
            totalenergy_list(i2)=totalenergy_local
            totalcharge_list(i2)=totalcharge_local
            shortenergy_list(i2)=shortenergy_local
            elecenergy_list(i2) =elecenergy_local
          endif
          jcount=jcount+1
        enddo ! i2
      enddo ! i1
!!
      close(symunit)
!!
      return
      end
