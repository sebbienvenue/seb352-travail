!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose: Find the neighbors of the pairs in a structure

!! FIXME: This subroutine is not yet parallel

!! called by: 
!! - calconepairfunction_para.f90
!!
      subroutine neighbortriplet_para(num_atoms,num_pairs,zelem,&
            xyzstruct,cutoff,lperiodic,tripletstore)
!!          lsta,lstb,cutoff,lattice,xyzstruct,&
!!          pairs_rr,pairs_atom,pairs_xyz,pairs_charge,&
!!          pair_lsta,pair_lstb,pair_lstc,pair_lste,lperiodic)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer i1,i2,i3,i4                   ! internal
      integer index_store,flag_store        ! internal
      integer i,j,k,l,m                     ! internal
      integer indc,indd                     ! internal
      integer na,nb,nc                      ! internal
      integer n1,n2,n3,n4,n5,n6             ! internal
      integer num_atoms                     ! in 
      integer n_count                       ! internal
      integer num_pairs                     ! out
      integer num_triplets                  ! out
      integer zelem(max_num_atoms)          ! in
      integer lsta(2,max_num_atoms)         ! out
      integer lsta2(2,max_num_atoms*(max_num_atoms-1))        ! out
      integer lstc(listdim)                 ! out
      integer lstc2(listdim)                ! out
      integer lste(listdim)                 ! internal
      integer lste2(listdim)                ! internal
      real*8  lstb(listdim,4)               ! out
      real*8  lstb2(listdim,5)              ! out
      integer tripletstore((max_num_atoms*(max_num_atoms-1)),3)
      real*8 triplet_rr((max_num_atoms*(max_num_atoms-1)),3)
      integer triplet_charge((max_num_atoms*(max_num_atoms-1)),3)
      real*8 triplet_xyz((max_num_atoms*(max_num_atoms-1)),3,3)
      integer pair_lsta(listdim,2)          ! out
      integer pair_lstc(listdim)            ! out
      integer pair_lste(listdim)            ! out
      integer pairs_atom(listdim,2)         ! out
      integer pairs_charge(2,listdim)       ! out
      integer temp_lstc(listdim)            ! internal
      integer temp_c(listdim)               ! internal 

      real*8 absaxb,absaxc,absbxc          ! internal
      real*8 axb(3),axc(3),bxc(3)          ! internal
      real*8 cutoff
      real*8 rr,rr1,rr2                        ! internal
      real*8  lattice(3,3)                  ! in
      real*8  proja,projb,projc             ! internal
      real*8  rr_r(2)                       ! internal
      real*8  xrel1,xrel2,xrel3             ! internal
      real*8  xrel12,xrel22,xrel32          ! internal
      real*8  xtemp,ytemp,ztemp             ! internal
      real*8  xtemp2,ytemp2,ztemp2          ! internal
      real*8  xx_r,yy_r,zz_r                ! internal
      real*8  xyzstruct(3,max_num_atoms)    ! in

      real*8  pairs_rr(listdim)             ! out
      real*8  pairs_xyz(listdim,2,3)        ! out
      real*8  pair_lstb(listdim,5)          ! out

      real*8  temp_x(listdim)               ! internal
      real*8  temp_y(listdim)               ! internal
      real*8  temp_z(listdim)               ! internal
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

!! expand system to some neighboring boxes in case of periodic system
      if(lperiodic)then
        n_count = 0
        do i=1,num_atoms
          do n1=-na,na
            do n2=-nb,nb
              do n3=-nc,nc
                n_count = n_count + 1      
                temp_x(n_count)=xyzstruct(1,i)+ dble(n1)*lattice(1,1) + dble(n2)*lattice(2,1) + dble(n3)*lattice(3,1)
                temp_y(n_count)=xyzstruct(2,i)+ dble(n1)*lattice(1,2) + dble(n2)*lattice(2,2) + dble(n3)*lattice(3,2)
                temp_z(n_count)=xyzstruct(3,i)+ dble(n1)*lattice(1,3) + dble(n2)*lattice(2,3) + dble(n3)*lattice(3,3)
                temp_c(n_count)=zelem(i)
                temp_lstc(n_count)=i
              enddo
            enddo
          enddo
        enddo
      endif ! lperiodic

!! determine lsta,lstb,lstc,lste
      if(lperiodic)then
        indc=0
        indd=0
        do i=1,num_atoms  !! loop over atoms for which we want the neighbors
          lsta(1,i)=indc+1
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
                    if(rr.le.rmin**2)then
                      write(ounit,*)'ERROR: Distance B/W Atoms in the System < rmin'
                      write(*,*) 'rmin', rmin, 'rr', rr
                      stop
                    endif
                    if(rr.le.cutoff**2) then
                      indc=indc+1
                      if(indc.gt.listdim)then
                         write(ounit,*)'neighborpair_para: redimension lstb'
                         stop
                      endif
                      lstb(indc,1)=xtemp
                      lstb(indc,2)=ytemp
                      lstb(indc,3)=ztemp
                      lstb(indc,4)=dsqrt(rr)
                      lstc(indc)  =j          ! identification of atom
                      lste(indc)  =zelem(j)   ! nuclear charge of atom
                      !!
                      !! here we need to nest the j loop, and do it again but this time we need to check i,j,k are not the same atom
                      !! it also means the storage arrays need modification
                      !!
                      lsta2(1,indc) = indd+1
                      do k=1,num_atoms
                        do n4=-na,na
                          do n5=-nb,nb
                            do n6=-nc,nc
                              if((n4.eq.0).and.(n5.eq.0).and.(n6.eq.0).and.((k.eq.j).or.(k.eq.i))) then !!avoid self interaction
                              else
                                xtemp2=xyzstruct(1,k) + dble(n1)*lattice(1,1) + dble(n2)*lattice(2,1) + dble(n3)*lattice(3,1)
                                ytemp2=xyzstruct(2,k) + dble(n1)*lattice(1,2) + dble(n2)*lattice(2,2) + dble(n3)*lattice(3,2)
                                ztemp2=xyzstruct(3,k) + dble(n1)*lattice(1,3) + dble(n2)*lattice(2,3) + dble(n3)*lattice(3,3)
                                xrel1=xtemp2-xyzstruct(1,i)
                                xrel2=ytemp2-xyzstruct(2,i)
                                xrel3=ztemp2-xyzstruct(3,i)
                                rr1=xrel1**2+xrel2**2+xrel3**2
                                xrel12=xtemp2-xyzstruct(1,j)
                                xrel22=ytemp2-xyzstruct(2,j)
                                xrel32=ztemp2-xyzstruct(3,j)
                                rr2=xrel12**2+xrel22**2+xrel32**2
                                if((rr.le.cutoff**2).or.(rr2.le.cutoff**2)) then
                                  if((rr.le.rmin).or.(rr2.le.rmin)) then
                                    write(ounit,*)'Distance B/W Atoms in the System < rmin 2'
                                    stop
                                  endif
                                  indd=indd+1
                                  if(indd.gt.listdim)then
                                    write(ounit,*)'Redimension lstb2'
                                    stop
                                  endif
                                  lstb2(indd,1)=xtemp2
                                  lstb2(indd,2)=ytemp2
                                  lstb2(indd,3)=ztemp2
                                  lstb2(indd,4)=dsqrt(rr1)
                                  lstb2(indd,5)=dsqrt(rr2)
                                  lstc2(indd)=k        ! identification of atom
                                  lste2(indd)=zelem(k) ! nuclear charge of atom
                                endif
                              endif
                            enddo
                          enddo
                        enddo
                      enddo !k
                      lsta2(2,indc)=indd
                    endif
                  endif
                enddo ! n3
              enddo ! n2
            enddo ! n1
          enddo ! j
          lsta(2,i)=indc
!!          write(ounit,*)'atom ',i
!!          write(ounit,*)'lsta ',lsta(1,i),lsta(2,i)
!!          do i1=lsta(1,i),lsta(2,i)
!!            write(ounit,*)'lstb ',i1,lstb(i1,4)
!!          enddo
        enddo ! i
!!
      else ! .not. lperiodic
!!
        write(*,*) 'Neighbours search begins here'
        indc=0
        indd=0
        do i=1,num_atoms
          lsta(1,i)=indc+1
          do j=1,num_atoms
            if(i.eq.j) then
            else
              xrel1=xyzstruct(1,j)-xyzstruct(1,i)
              xrel2=xyzstruct(2,j)-xyzstruct(2,i)
              xrel3=xyzstruct(3,j)-xyzstruct(3,i)
              rr=xrel1**2+xrel2**2+xrel3**2
!! check if neighbor is close
              if(rr.le.cutoff**2) then
                if(rr.le.rmin**2)then
                  write(ounit,*)'Distance B/W Atoms in the System < rmin'
                  write(*,*) 'rmin', rmin, 'rr', rr, 'cutoff', cutoff
                  write(*,*) i, j
                  write(*,*) xyzstruct(:,j)
                  write(*,*) xyzstruct(:,i)
                  stop
                endif
                indc=indc+1
                if(indc.gt.listdim)then
                  write(ounit,*)'Redimension lstb'
                  stop
                endif
!! store information of neighbor
                lstb(indc,1)=xyzstruct(1,j)
                lstb(indc,2)=xyzstruct(2,j)
                lstb(indc,3)=xyzstruct(3,j)
                lstb(indc,4)=dsqrt(rr)
                lstc(indc)=j        ! identification of atom
                lste(indc)=zelem(j) ! nuclear charge of atom
                !!
                !! Nest the loop in here.
                !!
                lsta2(1,indc) = indd+1
                do k=1,num_atoms
                  if((k.eq.j).or.(k.eq.i)) then
                  else
                    xtemp2=xyzstruct(1,k)
                    ytemp2=xyzstruct(2,k)
                    ztemp2=xyzstruct(3,k)
                    xrel1=xtemp2-xyzstruct(1,i)
                    xrel2=ytemp2-xyzstruct(2,i)
                    xrel3=ztemp2-xyzstruct(3,i)
                    rr1=xrel1**2+xrel2**2+xrel3**2
                    xrel12=xtemp2-xyzstruct(1,j)
                    xrel22=ytemp2-xyzstruct(2,j)
                    xrel32=ztemp2-xyzstruct(3,j)
                    rr2=xrel12**2+xrel22**2+xrel32**2
                    if((rr1.le.cutoff**2).or.(rr2.le.cutoff**2)) then
                      if((rr1.le.rmin).or.(rr2.le.rmin)) then
                        write(ounit,*)'Distance B/W Atoms in the System < rmin 2'
                        stop
                      endif
                      indd=indd+1
                      if(indd.gt.listdim)then
                        write(ounit,*)'Redimension lstb2'
                        stop
                      endif
                      lstb2(indd,1)=xtemp2
                      lstb2(indd,2)=ytemp2
                      lstb2(indd,3)=ztemp2
                      lstb2(indd,4)=dsqrt(rr1)
                      lstb2(indd,5)=dsqrt(rr2)
                      lstc2(indd)=k        ! identification of atom
                      lste2(indd)=zelem(k) ! nuclear charge of atom
                    endif
                  endif
                enddo !k
                lsta2(2,indc)=indd
              endif
            endif
          enddo ! j
          temp_c(i)=zelem(i)
          lsta(2,i)=indc
          temp_lstc(i)=i
        enddo ! i
!!
        do i=1,num_atoms
          temp_x(i) = xyzstruct(1,i)
          temp_y(i) = xyzstruct(2,i)
          temp_z(i) = xyzstruct(3,i)
        enddo
        n_count = num_atoms
      endif ! lperiodic
!!
!! find triplets
!!

      num_triplets=1 ! use as counter and index here
      do i2=1,num_atoms ! loop over all atoms
        do i3=lsta(1,i2),lsta(2,i2) ! loop over all neighbors of atom i2
          do i4=lsta2(1,i3),lsta2(2,i3) ! loop over all the neighbours of the possible pairs of i2,i3
            flag_store=3
           !! write(*,*) 'Triplets'
           !! write(*,*) i2,lstc(i3),lstc2(i4)
           !! write(*,*) zelem(i2), lste(i3), lste2(i4)
            if(lste(i3).lt.lste2(i4)) then
              if(zelem(i2).lt.lste(i3)) then !! heteronuclear trimer
                flag_store=2 ! keep
                tripletstore(num_triplets,1) = i2
                tripletstore(num_triplets,2) = lstc(i3)
                tripletstore(num_triplets,3) = lstc2(i4)
                triplet_rr(num_triplets,1) = lstb(i3,4)
                triplet_rr(num_triplets,2) = lstb2(i4,4)
                triplet_rr(num_triplets,3) = lstb2(i4,5)
                triplet_charge(num_triplets,1) = zelem(i2)
                triplet_charge(num_triplets,2) = lste(i3)
                triplet_charge(num_triplets,3) = lste2(i4)
                triplet_xyz(num_triplets,1,1) = xyzstruct(1,i2)
                triplet_xyz(num_triplets,1,2) = xyzstruct(1,i2)
                triplet_xyz(num_triplets,1,3) = xyzstruct(1,i2)
                triplet_xyz(num_triplets,2,1) = lstb(i3,1)
                triplet_xyz(num_triplets,2,2) = lstb(i3,2)
                triplet_xyz(num_triplets,2,3) = lstb(i3,3)
                triplet_xyz(num_triplets,3,1) = lstb2(i4,1)
                triplet_xyz(num_triplets,3,2) = lstb2(i4,2)
                triplet_xyz(num_triplets,3,3) = lstb2(i4,3)
              elseif(zelem(i2).eq.lste(i3)) then !! homo nuclear pair with a hetero triplet atom of the AAB type
                if(idnint(dble(mult*xyzstruct(1,i2))).lt.idnint(dble(mult*lstb(i3,1))))then !! x < x < _ case -2
                  flag_store=2 ! keep
                  tripletstore(num_triplets,1) = i2
                  tripletstore(num_triplets,2) = lstc(i3)
                  tripletstore(num_triplets,3) = lstc2(i4)
                  triplet_rr(num_triplets,1) = lstb(i3,4)
                  triplet_rr(num_triplets,2) = lstb2(i4,4)
                  triplet_rr(num_triplets,3) = lstb2(i4,5)
                  triplet_charge(num_triplets,1) = zelem(i2)
                  triplet_charge(num_triplets,2) = lste(i3)
                  triplet_charge(num_triplets,3) = lste2(i4)
                  triplet_xyz(num_triplets,1,1) = xyzstruct(1,i2)
                  triplet_xyz(num_triplets,1,2) = xyzstruct(1,i2)
                  triplet_xyz(num_triplets,1,3) = xyzstruct(1,i2)
                  triplet_xyz(num_triplets,2,1) = lstb(i3,1)
                  triplet_xyz(num_triplets,2,2) = lstb(i3,2)
                  triplet_xyz(num_triplets,2,3) = lstb(i3,3)
                  triplet_xyz(num_triplets,3,1) = lstb2(i4,1)
                  triplet_xyz(num_triplets,3,2) = lstb2(i4,2)
                  triplet_xyz(num_triplets,3,3) = lstb2(i4,3)
                elseif(idnint(dble(mult*xyzstruct(1,i2))).eq.idnint(dble(mult*lstb(i3,1))))then !! x = x < _
                  if(idnint(dble(mult*xyzstruct(2,i2))).lt.idnint(dble(mult*lstb(i3,2))))then !! y < y < _ case -1
                   flag_store=2 ! keep
                   tripletstore(num_triplets,1) = i2
                   tripletstore(num_triplets,2) = lstc(i3)
                   tripletstore(num_triplets,3) = lstc2(i4)
                   triplet_rr(num_triplets,1) = lstb(i3,4)
                   triplet_rr(num_triplets,2) = lstb2(i4,4)
                   triplet_rr(num_triplets,3) = lstb2(i4,5)
                   triplet_charge(num_triplets,1) = zelem(i2)
                   triplet_charge(num_triplets,2) = lste(i3)
                   triplet_charge(num_triplets,3) = lste2(i4)
                   triplet_xyz(num_triplets,1,1) = xyzstruct(1,i2)
                   triplet_xyz(num_triplets,1,2) = xyzstruct(1,i2)
                   triplet_xyz(num_triplets,1,3) = xyzstruct(1,i2)
                   triplet_xyz(num_triplets,2,1) = lstb(i3,1)
                   triplet_xyz(num_triplets,2,2) = lstb(i3,2)
                   triplet_xyz(num_triplets,2,3) = lstb(i3,3)
                   triplet_xyz(num_triplets,3,1) = lstb2(i4,1)
                   triplet_xyz(num_triplets,3,2) = lstb2(i4,2)
                   triplet_xyz(num_triplets,3,3) = lstb2(i4,3)
                  elseif(idnint(dble(mult*xyzstruct(2,i2))).eq.idnint(dble(mult*lstb(i3,2))))then !! y < y < _
                    if(idnint(dble(mult*xyzstruct(3,i2))).lt.idnint(dble(mult*lstb(i3,3))))then !! z < z < _ case 0
                      flag_store=2 ! keep
                      tripletstore(num_triplets,1) = i2
                      tripletstore(num_triplets,2) = lstc(i3)
                      tripletstore(num_triplets,3) = lstc2(i4)
                      triplet_rr(num_triplets,1) = lstb(i3,4)
                      triplet_rr(num_triplets,2) = lstb2(i4,4)
                      triplet_rr(num_triplets,3) = lstb2(i4,5)
                      triplet_charge(num_triplets,1) = zelem(i2)
                      triplet_charge(num_triplets,2) = lste(i3)
                      triplet_charge(num_triplets,3) = lste2(i4)
                      triplet_xyz(num_triplets,1,1) = xyzstruct(1,i2)
                      triplet_xyz(num_triplets,1,2) = xyzstruct(1,i2)
                      triplet_xyz(num_triplets,1,3) = xyzstruct(1,i2)
                      triplet_xyz(num_triplets,2,1) = lstb(i3,1)
                      triplet_xyz(num_triplets,2,2) = lstb(i3,2)
                      triplet_xyz(num_triplets,2,3) = lstb(i3,3)
                      triplet_xyz(num_triplets,3,1) = lstb2(i4,1)
                      triplet_xyz(num_triplets,3,2) = lstb2(i4,2)
                      triplet_xyz(num_triplets,3,3) = lstb2(i4,3)
                    else
                      flag_store=1
                    endif
                  else
                    flag_store=1
                  endif
                else
                  flag_store=1
                endif
              else
                flag_store=1
              endif
            elseif(lste(i3).eq.lste2(i4)) then
              if(zelem(i2).lt.lste(i3)) then !! homo nuclear pair with a hetero triplet atom of the type ABB
                if(idnint(dble(mult*lstb(i3,1))).lt.idnint(dble(mult*lstb2(i4,1))))then
                  flag_store=2 ! keep
                  tripletstore(num_triplets,1) = i2
                  tripletstore(num_triplets,2) = lstc(i3)
                  tripletstore(num_triplets,3) = lstc2(i4)
                  triplet_rr(num_triplets,1) = lstb(i3,4)
                  triplet_rr(num_triplets,2) = lstb2(i4,4)
                  triplet_rr(num_triplets,3) = lstb2(i4,5)
                  triplet_charge(num_triplets,1) = zelem(i2)
                  triplet_charge(num_triplets,2) = lste(i3)
                  triplet_charge(num_triplets,3) = lste2(i4)
                  triplet_xyz(num_triplets,1,1) = xyzstruct(1,i2)
                  triplet_xyz(num_triplets,1,2) = xyzstruct(2,i2)
                  triplet_xyz(num_triplets,1,3) = xyzstruct(3,i2)
                  triplet_xyz(num_triplets,2,1) = lstb(i3,1)
                  triplet_xyz(num_triplets,2,2) = lstb(i3,2)
                  triplet_xyz(num_triplets,2,3) = lstb(i3,3)
                  triplet_xyz(num_triplets,3,1) = lstb2(i4,1)
                  triplet_xyz(num_triplets,3,2) = lstb2(i4,2)
                  triplet_xyz(num_triplets,3,3) = lstb2(i4,3)
                elseif(idnint(dble(mult*lstb(i3,1))).eq.idnint(dble(mult*lstb2(i4,1))))then
                  if(idnint(dble(mult*lstb(i3,2))).lt.idnint(dble(mult*lstb2(i4,2))))then
                    flag_store=2 ! keep
                    tripletstore(num_triplets,1) = i2
                    tripletstore(num_triplets,2) = lstc(i3)
                    tripletstore(num_triplets,3) = lstc2(i4)
                    triplet_rr(num_triplets,1) = lstb(i3,4)
                    triplet_rr(num_triplets,2) = lstb2(i4,4)
                    triplet_rr(num_triplets,3) = lstb2(i4,5)
                    triplet_charge(num_triplets,1) = zelem(i2)
                    triplet_charge(num_triplets,2) = lste(i3)
                    triplet_charge(num_triplets,3) = lste2(i4)
                    triplet_xyz(num_triplets,1,1) = xyzstruct(1,i2)
                    triplet_xyz(num_triplets,1,2) = xyzstruct(2,i2)
                    triplet_xyz(num_triplets,1,3) = xyzstruct(3,i2)
                    triplet_xyz(num_triplets,2,1) = lstb(i3,1)
                    triplet_xyz(num_triplets,2,2) = lstb(i3,2)
                    triplet_xyz(num_triplets,2,3) = lstb(i3,3)
                    triplet_xyz(num_triplets,3,1) = lstb2(i4,1)
                    triplet_xyz(num_triplets,3,2) = lstb2(i4,2)
                    triplet_xyz(num_triplets,3,3) = lstb2(i4,3)
                  elseif(idnint(dble(mult*lstb(i3,2))).eq.idnint(dble(mult*lstb2(i4,2))))then
                    if(idnint(dble(mult*lstb(i3,3))).lt.idnint(dble(mult*lstb2(i4,3))))then
                      flag_store=2 ! keep
                      tripletstore(num_triplets,1) = i2
                      tripletstore(num_triplets,2) = lstc(i3)
                      tripletstore(num_triplets,3) = lstc2(i4)
                      triplet_rr(num_triplets,1) = lstb(i3,4)
                      triplet_rr(num_triplets,2) = lstb2(i4,4)
                      triplet_rr(num_triplets,3) = lstb2(i4,5)
                      triplet_charge(num_triplets,1) = zelem(i2)
                      triplet_charge(num_triplets,2) = lste(i3)
                      triplet_charge(num_triplets,3) = lste2(i4)
                      triplet_xyz(num_triplets,1,1) = xyzstruct(1,i2)
                      triplet_xyz(num_triplets,1,2) = xyzstruct(2,i2)
                      triplet_xyz(num_triplets,1,3) = xyzstruct(3,i2)
                      triplet_xyz(num_triplets,2,1) = lstb(i3,1)
                      triplet_xyz(num_triplets,2,2) = lstb(i3,2)
                      triplet_xyz(num_triplets,2,3) = lstb(i3,3)
                      triplet_xyz(num_triplets,3,1) = lstb2(i4,1)
                      triplet_xyz(num_triplets,3,2) = lstb2(i4,2)
                      triplet_xyz(num_triplets,3,3) = lstb2(i4,3)
                    else
                      flag_store=1 ! don't keep
                    endif
                  else
                    flag_store=1
                  endif
                else
                  flag_store=1
                endif
              elseif(zelem(i2).eq.lste(i3)) then !! homo nuclear triplet
                if(idnint(dble(mult*lstb(i3,1))).lt.idnint(dble(mult*lstb2(i4,1))))then    !!  _ < x < x
                  if(idnint(dble(mult*xyzstruct(1,i2))).lt.idnint(dble(mult*lstb(i3,1))))then !! x < x < x case 1
                    flag_store=2 ! keep
                    tripletstore(num_triplets,1) = i2
                    tripletstore(num_triplets,2) = lstc(i3)
                    tripletstore(num_triplets,3) = lstc2(i4)
                    triplet_rr(num_triplets,1) = lstb(i3,4)
                    triplet_rr(num_triplets,2) = lstb2(i4,4)
                    triplet_rr(num_triplets,3) = lstb2(i4,5)
                    triplet_charge(num_triplets,1) = zelem(i2)
                    triplet_charge(num_triplets,2) = lste(i3)
                    triplet_charge(num_triplets,3) = lste2(i4)
                    triplet_xyz(num_triplets,1,1) = xyzstruct(1,i2)
                    triplet_xyz(num_triplets,1,2) = xyzstruct(2,i2)
                    triplet_xyz(num_triplets,1,3) = xyzstruct(3,i2)
                    triplet_xyz(num_triplets,2,1) = lstb(i3,1)
                    triplet_xyz(num_triplets,2,2) = lstb(i3,2)
                    triplet_xyz(num_triplets,2,3) = lstb(i3,3)
                    triplet_xyz(num_triplets,3,1) = lstb2(i4,1)
                    triplet_xyz(num_triplets,3,2) = lstb2(i4,2)
                    triplet_xyz(num_triplets,3,3) = lstb2(i4,3)
                  elseif(idnint(dble(mult*xyzstruct(1,i2))).eq.idnint(dble(mult*lstb(i3,1))))then !! x = x < x
                    if(idnint(dble(mult*xyzstruct(2,i2))).lt.idnint(dble(mult*lstb(i3,2))))then !! y < yx < x case 2
                      flag_store=2 ! keep
                      tripletstore(num_triplets,1) = i2
                      tripletstore(num_triplets,2) = lstc(i3)
                      tripletstore(num_triplets,3) = lstc2(i4)
                      triplet_rr(num_triplets,1) = lstb(i3,4)
                      triplet_rr(num_triplets,2) = lstb2(i4,4)
                      triplet_rr(num_triplets,3) = lstb2(i4,5)
                      triplet_charge(num_triplets,1) = zelem(i2)
                      triplet_charge(num_triplets,2) = lste(i3)
                      triplet_charge(num_triplets,3) = lste2(i4)
                      triplet_xyz(num_triplets,1,1) = xyzstruct(1,i2)
                      triplet_xyz(num_triplets,1,2) = xyzstruct(2,i2)
                      triplet_xyz(num_triplets,1,3) = xyzstruct(3,i2)
                      triplet_xyz(num_triplets,2,1) = lstb(i3,1)
                      triplet_xyz(num_triplets,2,2) = lstb(i3,2)
                      triplet_xyz(num_triplets,2,3) = lstb(i3,3)
                      triplet_xyz(num_triplets,3,1) = lstb2(i4,1)
                      triplet_xyz(num_triplets,3,2) = lstb2(i4,2)
                      triplet_xyz(num_triplets,3,3) = lstb2(i4,3)
                    elseif(idnint(dble(mult*xyzstruct(2,i2))).eq.idnint(dble(mult*lstb(i3,2))))then !! yx = yx < x
                      if(idnint(dble(mult*xyzstruct(3,i2))).lt.idnint(dble(mult*lstb(i3,3))))then !! z < zyx < x case 3
                        flag_store=2 ! keep
                        tripletstore(num_triplets,1) = i2
                        tripletstore(num_triplets,2) = lstc(i3)
                        tripletstore(num_triplets,3) = lstc2(i4)
                        triplet_rr(num_triplets,1) = lstb(i3,4)
                        triplet_rr(num_triplets,2) = lstb2(i4,4)
                        triplet_rr(num_triplets,3) = lstb2(i4,5)
                        triplet_charge(num_triplets,1) = zelem(i2)
                        triplet_charge(num_triplets,2) = lste(i3)
                        triplet_charge(num_triplets,3) = lste2(i4)
                        triplet_xyz(num_triplets,1,1) = xyzstruct(1,i2)
                        triplet_xyz(num_triplets,1,2) = xyzstruct(2,i2)
                        triplet_xyz(num_triplets,1,3) = xyzstruct(3,i2)
                        triplet_xyz(num_triplets,2,1) = lstb(i3,1)
                        triplet_xyz(num_triplets,2,2) = lstb(i3,2)
                        triplet_xyz(num_triplets,2,3) = lstb(i3,3)
                        triplet_xyz(num_triplets,3,1) = lstb2(i4,1)
                        triplet_xyz(num_triplets,3,2) = lstb2(i4,2)
                        triplet_xyz(num_triplets,3,3) = lstb2(i4,3)
                      else
                        flag_store=1 ! keep
                      endif
                    else
                      flag_store=1
                    endif
                  else
                    flag_store=1
                  endif
                elseif(idnint(dble(mult*lstb(i3,1))).eq.idnint(dble(mult*lstb2(i4,1))))then    !!  _ < x = x
                  if(idnint(dble(mult*xyzstruct(1,i2))).lt.idnint(dble(mult*lstb(i3,1))))then    !!  x < x = x
                    if(idnint(dble(mult*lstb(i3,2))).lt.idnint(dble(mult*lstb2(i4,2))))then    !!  _ < y < y case 4
                      flag_store=2 ! keep
                      tripletstore(num_triplets,1) = i2
                      tripletstore(num_triplets,2) = lstc(i3)
                      tripletstore(num_triplets,3) = lstc2(i4)
                      triplet_rr(num_triplets,1) = lstb(i3,4)
                      triplet_rr(num_triplets,2) = lstb2(i4,4)
                      triplet_rr(num_triplets,3) = lstb2(i4,5)
                      triplet_charge(num_triplets,1) = zelem(i2)
                      triplet_charge(num_triplets,2) = lste(i3)
                      triplet_charge(num_triplets,3) = lste2(i4)
                      triplet_xyz(num_triplets,1,1) = xyzstruct(1,i2)
                      triplet_xyz(num_triplets,1,2) = xyzstruct(2,i2)
                      triplet_xyz(num_triplets,1,3) = xyzstruct(3,i2)
                      triplet_xyz(num_triplets,2,1) = lstb(i3,1)
                      triplet_xyz(num_triplets,2,2) = lstb(i3,2)
                      triplet_xyz(num_triplets,2,3) = lstb(i3,3)
                      triplet_xyz(num_triplets,3,1) = lstb2(i4,1)
                      triplet_xyz(num_triplets,3,2) = lstb2(i4,2)
                      triplet_xyz(num_triplets,3,3) = lstb2(i4,3)
                    elseif(idnint(dble(mult*lstb(i3,2))).eq.idnint(dble(mult*lstb2(i4,2))))then    !!  _ < y = y
                      if(idnint(dble(mult*lstb(i3,3))).lt.idnint(dble(mult*lstb2(i4,3))))then    !!  _ < z < z case 5
                        flag_store=2 ! keep
                        tripletstore(num_triplets,1) = i2
                        tripletstore(num_triplets,2) = lstc(i3)
                        tripletstore(num_triplets,3) = lstc2(i4)
                        triplet_rr(num_triplets,1) = lstb(i3,4)
                        triplet_rr(num_triplets,2) = lstb2(i4,4)
                        triplet_rr(num_triplets,3) = lstb2(i4,5)
                        triplet_charge(num_triplets,1) = zelem(i2)
                        triplet_charge(num_triplets,2) = lste(i3)
                        triplet_charge(num_triplets,3) = lste2(i4)
                        triplet_xyz(num_triplets,1,1) = xyzstruct(1,i2)
                        triplet_xyz(num_triplets,1,2) = xyzstruct(2,i2)
                        triplet_xyz(num_triplets,1,3) = xyzstruct(3,i2)
                        triplet_xyz(num_triplets,2,1) = lstb(i3,1)
                        triplet_xyz(num_triplets,2,2) = lstb(i3,2)
                        triplet_xyz(num_triplets,2,3) = lstb(i3,3)
                        triplet_xyz(num_triplets,3,1) = lstb2(i4,1)
                        triplet_xyz(num_triplets,3,2) = lstb2(i4,2)
                        triplet_xyz(num_triplets,3,3) = lstb2(i4,3)
                      else
                        flag_store=1
                      endif
                    else
                      flag_store=1
                    endif
                  elseif(idnint(dble(mult*xyzstruct(1,i2))).eq.idnint(dble(mult*lstb(i3,1))))then    !!  x = x = x
                    if(idnint(dble(mult*lstb(i3,2))).lt.idnint(dble(mult*lstb2(i4,2))))then    !!  _ < y < y
                      if(idnint(dble(mult*xyzstruct(2,i2))).lt.idnint(dble(mult*lstb(i3,2))))then !! y < y < y case 6
                        flag_store=2 ! keep
                        tripletstore(num_triplets,1) = i2
                        tripletstore(num_triplets,2) = lstc(i3)
                        tripletstore(num_triplets,3) = lstc2(i4)
                        triplet_rr(num_triplets,1) = lstb(i3,4)
                        triplet_rr(num_triplets,2) = lstb2(i4,4)
                        triplet_rr(num_triplets,3) = lstb2(i4,5)
                        triplet_charge(num_triplets,1) = zelem(i2)
                        triplet_charge(num_triplets,2) = lste(i3)
                        triplet_charge(num_triplets,3) = lste2(i4)
                        triplet_xyz(num_triplets,1,1) = xyzstruct(1,i2)
                        triplet_xyz(num_triplets,1,2) = xyzstruct(2,i2)
                        triplet_xyz(num_triplets,1,3) = xyzstruct(3,i2)
                        triplet_xyz(num_triplets,2,1) = lstb(i3,1)
                        triplet_xyz(num_triplets,2,2) = lstb(i3,2)
                        triplet_xyz(num_triplets,2,3) = lstb(i3,3)
                        triplet_xyz(num_triplets,3,1) = lstb2(i4,1)
                        triplet_xyz(num_triplets,3,2) = lstb2(i4,2)
                        triplet_xyz(num_triplets,3,3) = lstb2(i4,3)
                      elseif(idnint(dble(mult*xyzstruct(2,i2))).eq.idnint(dble(mult*lstb(i3,2))))then !! y = y < y
                        if(idnint(dble(mult*xyzstruct(3,i2))).lt.idnint(dble(mult*lstb(i3,3))))then !! z < z < _ case 7
                          flag_store=2 ! keep
                          tripletstore(num_triplets,1) = i2
                          tripletstore(num_triplets,2) = lstc(i3)
                          tripletstore(num_triplets,3) = lstc2(i4)
                          triplet_rr(num_triplets,1) = lstb(i3,4)
                          triplet_rr(num_triplets,2) = lstb2(i4,4)
                          triplet_rr(num_triplets,3) = lstb2(i4,5)
                          triplet_charge(num_triplets,1) = zelem(i2)
                          triplet_charge(num_triplets,2) = lste(i3)
                          triplet_charge(num_triplets,3) = lste2(i4)
                          triplet_xyz(num_triplets,1,1) = xyzstruct(1,i2)
                          triplet_xyz(num_triplets,1,2) = xyzstruct(2,i2)
                          triplet_xyz(num_triplets,1,3) = xyzstruct(3,i2)
                          triplet_xyz(num_triplets,2,1) = lstb(i3,1)
                          triplet_xyz(num_triplets,2,2) = lstb(i3,2)
                          triplet_xyz(num_triplets,2,3) = lstb(i3,3)
                          triplet_xyz(num_triplets,3,1) = lstb2(i4,1)
                          triplet_xyz(num_triplets,3,2) = lstb2(i4,2)
                          triplet_xyz(num_triplets,3,3) = lstb2(i4,3)
                        else
                          flag_store=1
                        endif
                      else
                        flag_store=1
                      endif
                    elseif(idnint(dble(mult*lstb(i3,2))).lt.idnint(dble(mult*lstb2(i4,2))))then    !!  _ < y = y
                      if(idnint(dble(mult*xyzstruct(2,i2))).lt.idnint(dble(mult*lstb(i3,2))))then !! y < y = y
                        if(idnint(dble(mult*lstb(i3,3))).lt.idnint(dble(mult*lstb2(i4,3))))then !! _ < z < z case 8
                          flag_store=2 ! keep
                          tripletstore(num_triplets,1) = i2
                          tripletstore(num_triplets,2) = lstc(i3)
                          tripletstore(num_triplets,3) = lstc2(i4)
                          triplet_rr(num_triplets,1) = lstb(i3,4)
                          triplet_rr(num_triplets,2) = lstb2(i4,4)
                          triplet_rr(num_triplets,3) = lstb2(i4,5)
                          triplet_charge(num_triplets,1) = zelem(i2)
                          triplet_charge(num_triplets,2) = lste(i3)
                          triplet_charge(num_triplets,3) = lste2(i4)
                          triplet_xyz(num_triplets,1,1) = xyzstruct(1,i2)
                          triplet_xyz(num_triplets,1,2) = xyzstruct(2,i2)
                          triplet_xyz(num_triplets,1,3) = xyzstruct(3,i2)
                          triplet_xyz(num_triplets,2,1) = lstb(i3,1)
                          triplet_xyz(num_triplets,2,2) = lstb(i3,2)
                          triplet_xyz(num_triplets,2,3) = lstb(i3,3)
                          triplet_xyz(num_triplets,3,1) = lstb2(i4,1)
                          triplet_xyz(num_triplets,3,2) = lstb2(i4,2)
                          triplet_xyz(num_triplets,3,3) = lstb2(i4,3)
                        else
                          flag_store=1
                        endif
                      elseif(idnint(dble(mult*xyzstruct(2,i2))).eq.idnint(dble(mult*lstb(i3,2))))then !! y = y = y
                        if(idnint(dble(mult*lstb(i3,3))).lt.idnint(dble(mult*lstb2(i4,3))))then !! _ < z < z
                          if(idnint(dble(mult*xyzstruct(3,i2))).lt.idnint(dble(mult*lstb(i3,3))))then !! z < z < z case 9
                            flag_store=2 ! keep
                            tripletstore(num_triplets,1) = i2
                            tripletstore(num_triplets,2) = lstc(i3)
                            tripletstore(num_triplets,3) = lstc2(i4)
                            triplet_rr(num_triplets,1) = lstb(i3,4)
                            triplet_rr(num_triplets,2) = lstb2(i4,4)
                            triplet_rr(num_triplets,3) = lstb2(i4,5)
                            triplet_charge(num_triplets,1) = zelem(i2)
                            triplet_charge(num_triplets,2) = lste(i3)
                            triplet_charge(num_triplets,3) = lste2(i4)
                            triplet_xyz(num_triplets,1,1) = xyzstruct(1,i2)
                            triplet_xyz(num_triplets,1,2) = xyzstruct(2,i2)
                            triplet_xyz(num_triplets,1,3) = xyzstruct(3,i2)
                            triplet_xyz(num_triplets,2,1) = lstb(i3,1)
                            triplet_xyz(num_triplets,2,2) = lstb(i3,2)
                            triplet_xyz(num_triplets,2,3) = lstb(i3,3)
                            triplet_xyz(num_triplets,3,1) = lstb2(i4,1)
                            triplet_xyz(num_triplets,3,2) = lstb2(i4,2)
                            triplet_xyz(num_triplets,3,3) = lstb2(i4,3)
                          else
                            flag_store=1
                          endif
                        else
                          flag_store=1
                        endif
                      else
                        flag_store=1
                      endif
                    else
                      flag_store=1
                    endif
                  else
                    flag_store=1
                  endif
                else
                  flag_store=1
                endif
              else
                flag_store=1
              endif
            else
              flag_store=1
            endif
            if(flag_store.eq.2) then
              num_triplets                   = num_triplets + 1
            endif
          enddo
        enddo
      enddo
      num_triplets=num_triplets-1 ! correct counter
!!      write(*,*) max_num_triplets

!!      do 12=1,num_triplets
!!        triplet_atom(i2,1) = tripletstore(i2,1)
!!        triplet_atom(i2,2) = tripletstore(i2,2)
!!        triplet_atom(i2,3) = tripletstore(i2,3)
!!
!!        triplet_rr(i2,1) = lstb
!!      enddo
!!      write(*,*) 'Triplet Store'
!!      write(*,*) 'num triplets', num_triplets
      
!!      write(*,*) triplet_charge(1,1), triplet_charge(1,2)
!!      STOP
!!      do i2=1,num_triplets
!!        write(*,*) tripletstore(i2,1), tripletstore(i2,2), tripletstore(i2,3)
!!        write(*,*) triplet_rr(i2,1), triplet_rr(i2,2), triplet_rr(i2,3)
!!        write(*,*) triplet_charge(i2,:)
!!        write(*,*) triplet_xyz(i2,:,:) 
!!      enddo
     


!!      STOP









!!
      return
      end
