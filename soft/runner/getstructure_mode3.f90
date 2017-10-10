!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose: Read structure and derive all structure-related data for mode 3
!!
      subroutine getstructure_mode3(num_atoms,zelem,&
        num_atoms_element,lattice,xyzstruct,&
        totalenergy,totalcharge,totalforce,atomenergy,atomcharge,&
        elementsymbol,lperiodic&
        )
!!
      use fileunits
      use mpi_mod
      use globaloptions
!!
      implicit none
!!
      integer zelem(max_num_atoms)                           ! out 
      integer num_atoms                                      ! out 
      integer num_atoms_element(nelem)                       ! out 
      integer i1,i2                                          ! internal
!!
!! DFT data (not necessarily provided in predicition mode)
      real*8 totalcharge                                     ! out 
      real*8 totalenergy                                     ! out 
      real*8 totalforce(3,max_num_atoms)                     ! out 
      real*8 atomcharge(max_num_atoms)                       ! out 
      real*8 atomenergy(max_num_atoms)                       ! out 
      real*8 lattice(3,3)                                    ! out 
      real*8 xyzstruct(3,max_num_atoms)                      ! out 
!!
      character*2 elementsymbol(max_num_atoms)               ! out 
!!
      logical lperiodic                                      ! out 
!!
!!=====================================================================
!! initializations for structure
!!=====================================================================
      lattice(:,:)        =0.0d0
      xyzstruct(:,:)      =0.0d0
      num_atoms           =0
      num_atoms_element(:)=0
      zelem(:)            =0
      elementsymbol(:)    ='  '
      lperiodic           =.false.
!!=====================================================================
!! initializations for DFT data (usually are not given)
!!=====================================================================
      totalcharge    =0.0d0
      totalenergy    =0.0d0
      atomcharge(:)  =0.0d0
      atomenergy(:)  =0.0d0
      totalforce(:,:)=0.0d0
!!
!!=====================================================================
!! read and distribute structure
!!=====================================================================
      if(mpirank.eq.0)then
        open(dataunit,file='input.data',form='formatted',status='old')
          rewind(dataunit)
          call readonestructure(num_atoms,&
            zelem,num_atoms_element,lattice,&
            totalcharge,totalenergy,atomcharge,atomenergy,xyzstruct,&
            totalforce,elementsymbol,lperiodic)
        close(dataunit)
      endif
!!
!!=====================================================================
!! write structure to output file 
!!=====================================================================
      if(mpirank.eq.0)then
        write(ounit,'(a,i5,a)')' Structure with ',num_atoms,' atoms in input.data:'
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'begin'
        if(lperiodic)then
          write(ounit,'(a,3f14.8)')' lattice ',lattice(1,1),lattice(1,2),lattice(1,3)
          write(ounit,'(a,3f14.8)')' lattice ',lattice(2,1),lattice(2,2),lattice(2,3)
          write(ounit,'(a,3f14.8)')' lattice ',lattice(3,1),lattice(3,2),lattice(3,3)
        endif
        do i1=1,num_atoms
          write(ounit,'(a,3f16.9,x,a2,x,5f16.9)')' atom ',&
            xyzstruct(1,i1),xyzstruct(2,i1),xyzstruct(3,i1),elementsymbol(i1),&
            atomcharge(i1),atomenergy(i1),totalforce(1,i1),totalforce(2,i1),totalforce(3,i1)
        enddo
        write(ounit,'(a,f18.8)')' energy ',totalenergy
        write(ounit,'(a,f18.8)')' charge ',totalcharge
        write(ounit,*)'end'
        write(ounit,*)'-------------------------------------------------------------'
      endif ! 'mpirank.eq.0
!!
!!=====================================================================
!! distribute full structure arrays independent of system size to all processes
!!=====================================================================
      call mpi_bcast(num_atoms,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(lattice,9,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(num_atoms_element,nelem,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(totalcharge,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(totalenergy,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(lperiodic,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(zelem,max_num_atoms,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(atomcharge,max_num_atoms,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(atomenergy,max_num_atoms,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(xyzstruct,3*max_num_atoms,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(totalforce,3*max_num_atoms,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(elementsymbol,max_num_atoms,mpi_character,0,mpi_comm_world,mpierror)
!!
      return
      end
