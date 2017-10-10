!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose: determine num_pairs for a given structure

!! called by: 
!! - paircount.f90
!! - predict.f90
!!
      subroutine getnumpairs(num_atoms,num_pairs,zelem,&
        cutoff,lattice,xyzstruct,lperiodic)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer i1,i2,i3                      ! internal
      integer flag_store                    ! internal
      integer i,j                           ! internal
      integer indc                          ! internal
      integer na,nb,nc                      ! internal
      integer n1,n2,n3                      ! internal
      integer num_atoms                     ! in 
      integer n_count                       ! internal
      integer num_pairs                     ! out

      integer zelem(max_num_atoms)          ! in
      integer lsta(2,max_num_atoms)         ! out
      integer lstc(listdim)                 ! out
      integer lste(listdim)                 ! out

      real*8  lstb(listdim,4)               ! out
      real*8  absaxb,absaxc,absbxc          ! internal
      real*8  axb(3),axc(3),bxc(3)          ! internal
      real*8  cutoff,rr                     ! in
      real*8  lattice(3,3)                  ! in
      real*8  proja,projb,projc             ! internal
      real*8  xrel1,xrel2,xrel3             ! internal
      real*8  xtemp,ytemp,ztemp             ! internal
      real*8  xyzstruct(3,max_num_atoms)    ! in

      real*8  mult
      logical lperiodic
!!
      mult =10000000.d0
      if(lperiodic) then
!! determine how often we have to multiply the cell in which direction
      na = 0
      nb = 0
      nc = 0
!! calculate the norm vectors of the 3 planes
      axb(1)=lattice(1,2)*lattice(2,3)-lattice(1,3)*lattice(2,2)
      axb(2)=lattice(1,3)*lattice(2,1)-lattice(1,1)*lattice(2,3)
      axb(3)=lattice(1,1)*lattice(2,2)-lattice(1,2)*lattice(2,1)
      absaxb=sqrt(axb(1)**2+axb(2)**2+axb(3)**2)
      axb(1)=axb(1)/absaxb
      axb(2)=axb(2)/absaxb
      axb(3)=axb(3)/absaxb
      axc(1)=lattice(1,2)*lattice(3,3)-lattice(1,3)*lattice(3,2)
      axc(2)=lattice(1,3)*lattice(3,1)-lattice(1,1)*lattice(3,3)
      axc(3)=lattice(1,1)*lattice(3,2)-lattice(1,2)*lattice(3,1)
      absaxc=sqrt(axc(1)**2+axc(2)**2+axc(3)**2)
      axc(1)=axc(1)/absaxc
      axc(2)=axc(2)/absaxc
      axc(3)=axc(3)/absaxc
      bxc(1)=lattice(2,2)*lattice(3,3)-lattice(2,3)*lattice(3,2)
      bxc(2)=lattice(2,3)*lattice(3,1)-lattice(2,1)*lattice(3,3)
      bxc(3)=lattice(2,1)*lattice(3,2)-lattice(2,2)*lattice(3,1)
      absbxc=sqrt(bxc(1)**2+bxc(2)**2+bxc(3)**2)
      bxc(1)=bxc(1)/absbxc
      bxc(2)=bxc(2)/absbxc
      bxc(3)=bxc(3)/absbxc
!! calculate the projections
      proja=lattice(1,1)*bxc(1)+lattice(1,2)*bxc(2)+lattice(1,3)*bxc(3)
      projb=lattice(2,1)*axc(1)+lattice(2,2)*axc(2)+lattice(2,3)*axc(3)
      projc=lattice(3,1)*axb(1)+lattice(3,2)*axb(2)+lattice(3,3)*axb(3)
      proja=abs(proja)
      projb=abs(projb)
      projc=abs(projc)
!! determine how often we have to multiply the cell in which direction
      na = 0
      nb = 0
      nc = 0
      do while(dble(na)*proja.le.cutoff)
        na=na+1
      enddo
      do while(dble(nb)*projb.le.cutoff)
        nb=nb+1
      enddo
      do while(dble(nc)*projc.le.cutoff)
        nc=nc+1
      enddo
      endif

!! check if we have pairs

!!
!! CHANGE ANDI: GFORTRAN: gfortran forbids .eq. for logicals, use .eqv. instead.
!!                        ifort allows both.
!!
     !if((num_atoms.eq.1).and.(lperiodic.eq..false.))then
      if((num_atoms.eq.1).and.(lperiodic.eqv..false.))then
!! END CHANGE

         write(*,*)'ERROR: nonperiodic structure with 1 atom does not have pairs!'
         stop !'
       endif

!! determine lsta,lstb,lstc,lste
      if(lperiodic)then
        indc=0
        do i=1,num_atoms  !! loop over atoms for which we want the neighbors
          indc=indc+1
          lsta(1,i)=indc !+1
          do j=1,num_atoms !! loop over all potentially neighboring atoms
            do n1=-na,na
              do n2=-nb,nb
                do n3=-nc,nc
                  if((n1.eq.0).and.(n2.eq.0).and.(n3.eq.0).and.(i.eq.j)) then !!avoid self interaction
                  else
                    xtemp=xyzstruct(1,j)+ dble(n1)*lattice(1,1) + dble(n2)*lattice(2,1) + dble(n3)*lattice(3,1)
                    ytemp=xyzstruct(2,j)+ dble(n1)*lattice(1,2) + dble(n2)*lattice(2,2) + dble(n3)*lattice(3,2)
                    ztemp=xyzstruct(3,j)+ dble(n1)*lattice(1,3) + dble(n2)*lattice(2,3) + dble(n3)*lattice(3,3)
                    xrel1=xtemp-xyzstruct(1,i)
                    xrel2=ytemp-xyzstruct(2,i)
                    xrel3=ztemp-xyzstruct(3,i)
                    rr=xrel1**2+xrel2**2+xrel3**2
                    if(rr.le.rmin)then
                      write(ounit,*)'Distance B/W Atoms in the System < rmin'
                      stop
                    endif
                    if(rr.le.cutoff**2) then
                      indc=indc+1
                      if(indc.gt.listdim)then
                        write(ounit,*)'Redimension lstb'
                        stop
                      endif
                      lstb(indc,1)=xtemp
                      lstb(indc,2)=ytemp
                      lstb(indc,3)=ztemp
                      lstb(indc,4)=sqrt(rr)
                      lstc(indc)  =j          ! identification of atom
                      lste(indc)  =zelem(j)   ! nuclear charge of atom
                    endif
                  endif
                enddo ! n3
              enddo ! n2
            enddo ! n1
          enddo ! j
          lsta(2,i)=indc
        enddo ! i
!!
      else ! .not. lperiodic
        indc=0
        do i=1,num_atoms
          lsta(1,i)=indc+1
          do j=1,num_atoms
            if(i.eq.j) then
            else
              xtemp=xyzstruct(1,j)
              ytemp=xyzstruct(2,j)
              ztemp=xyzstruct(3,j)
              xrel1=xtemp-xyzstruct(1,i)
              xrel2=ytemp-xyzstruct(2,i)
              xrel3=ztemp-xyzstruct(3,i)
              rr=xrel1**2+xrel2**2+xrel3**2
              if(rr.le.cutoff**2) then
                if(rr.le.rmin)then
                  write(ounit,*)'Distance B/W Atoms in the System < rmin'
                  stop
                endif
                indc=indc+1
                if(indc.gt.listdim)then
                  write(ounit,*)'Redimension lstb'
                  stop
                endif
                lstb(indc,1)=xtemp
                lstb(indc,2)=ytemp
                lstb(indc,3)=ztemp
                lstb(indc,4)=sqrt(rr)
                lstc(indc)=j        ! identification of atom
                lste(indc)=zelem(j) ! nuclear charge of atom
              endif
            endif
          enddo   ! j
          lsta(2,i)=indc
        enddo     ! i
        n_count = num_atoms
      endif ! lperiodic
!!      write(*,*) 'Get num pairs'
!!      write(*,*) 'lsta', lsta(1,1), lsta(2,1)
!!      write(*,*) 'lsta', lsta(1,2), lsta(2,2)
!!      write(*,*) 'lsta', lsta(1,3), lsta(2,3)
!!
!! find pairs 
      num_pairs=1 ! use as counter and index here
      do i2=1,num_atoms ! loop over all atoms
        do i3=lsta(1,i2),lsta(2,i2) ! loop over all neighbors of atom i2
          flag_store=3
          if(zelem(i2).lt.lste(i3))then     !! heteronuclear pair
            flag_store=2 ! keep
          elseif(zelem(i2).eq.lste(i3))then !! homonuclear pair
            if(idnint(dble(mult*xyzstruct(1,i2))).lt.idnint(dble(mult*lstb(i3,1))))then
              flag_store=2 ! keep
            elseif(idnint(dble(mult*xyzstruct(1,i2))).eq.idnint(dble(mult*lstb(i3,1))))then
              if(idnint(dble(mult*xyzstruct(2,i2))).lt.idnint(dble(mult*lstb(i3,2))))then 
                flag_store=2 ! keep
              elseif(idnint(dble(mult*xyzstruct(2,i2))).eq.idnint(dble(mult*lstb(i3,2))))then
                if(idnint(dble(mult*xyzstruct(3,i2))).lt.idnint(dble(mult*lstb(i3,3))))then
                  flag_store=2 ! keep
                else
                  flag_store=1 ! don't keep
                endif 
              else
                flag_store=1 ! don't keep
              endif 
            endif  
          else ! zelem(i2).gt.lste(i3) => don't keep 
            flag_store=1
          endif
          if(flag_store.eq.2)then ! keep information of this pair
            num_pairs                   = num_pairs + 1
          endif
        enddo
      enddo    
      num_pairs=num_pairs-1 ! correct counter
!!
      return
      end
