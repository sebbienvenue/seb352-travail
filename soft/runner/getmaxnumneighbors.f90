!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine

!! Goal: determine max_num_neighbors as dimension for large array dsfuncdxyz
!! for atomic case

!! called by:
!!
      subroutine getmaxnumneighbors(n_start,n_end,&
        num_atoms,max_num_atoms,max_num_neighbors_local,&
        maxcutoff_local,lattice,xyzstruct,lperiodic)
!!
      use fileunits
!!
      implicit none
!!
      integer i,j                                            ! internal
      integer icount                                         ! internal
      integer max_num_neighbors_local                        ! out
      integer max_num_atoms                                  ! in        
      integer n_start                                        ! in
      integer n_end                                          ! in
      integer na,nb,nc                                       ! internal
      integer n1,n2,n3                                       ! internal
      integer num_atoms                                      ! in
      integer lsta_local(2)                                  ! internal
      integer indc                                           ! internal
!!
      real*8 maxcutoff_local                                 ! in
      real*8 lattice(3,3)                                    ! in
      real*8 xyzstruct(3,max_num_atoms)                      ! in
      real*8 axb(3),axc(3),bxc(3)                            ! internal
      real*8 absaxb,absaxc,absbxc                            ! internal
      real*8 proja,projb,projc                               ! internal
      real*8 rr                                              ! internal
      real*8 xtemp,ytemp,ztemp                               ! internal
      real*8 xrel1,xrel2,xrel3                               ! internal
!!
      logical lperiodic                                      ! in
!!
      max_num_neighbors_local=0
!!
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
        do while(dble(na)*proja.le.maxcutoff_local)
          na=na+1
        enddo
        do while(dble(nb)*projb.le.maxcutoff_local)
          nb=nb+1
        enddo
        do while(dble(nc)*projc.le.maxcutoff_local)
          nc=nc+1
        enddo
      endif
!!
      icount=0
      if(lperiodic)then
        indc=0
!! loop over atoms for which we want the neighbors
        do i=n_start,n_end
          lsta_local(1)=indc+1
!! loop over all potentially neighboring atoms
          do j=1,num_atoms
            do n1=-na,na
              do n2=-nb,nb
                do n3=-nc,nc
!! avoid interaction of atom with itself
                  if((n1.eq.0).and.(n2.eq.0).and.(n3.eq.0).and.(i.eq.j)) then
                  else
                    xtemp=xyzstruct(1,j)&
      +dble(n1)*lattice(1,1)+dble(n2)*lattice(2,1)+dble(n3)*lattice(3,1)
                    ytemp=xyzstruct(2,j)&
      +dble(n1)*lattice(1,2)+dble(n2)*lattice(2,2)+dble(n3)*lattice(3,2)
                    ztemp=xyzstruct(3,j)&
      +dble(n1)*lattice(1,3)+dble(n2)*lattice(2,3)+dble(n3)*lattice(3,3)
                    xrel1=xtemp-xyzstruct(1,i)
                    xrel2=ytemp-xyzstruct(2,i)
                    xrel3=ztemp-xyzstruct(3,i)
                    rr=xrel1**2+xrel2**2+xrel3**2
                    if(rr.le.maxcutoff_local**2) then
                      indc=indc+1
                    endif
                  endif
                enddo ! n3
              enddo ! n2
            enddo ! n1
          enddo ! j
          lsta_local(2)=indc
          max_num_neighbors_local=max(max_num_neighbors_local,lsta_local(2)-lsta_local(1)+1)
        enddo ! i
!!
      else ! .not. lperiodic
        indc=0
        do i=n_start,n_end
          lsta_local(1)=indc+1
!! loop over all potentially neighboring atoms
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
              if(rr.le.maxcutoff_local**2) then
                indc=indc+1
              endif
            endif
          enddo ! j
          lsta_local(2)=indc
          max_num_neighbors_local=max(max_num_neighbors_local,lsta_local(2)-lsta_local(1)+1)
        enddo ! i
      endif
!!
      return
      end
