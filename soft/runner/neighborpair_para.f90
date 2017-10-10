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
      subroutine neighborpair_para(num_atoms,num_pairs,zelem,&
          lsta,lstb,cutoff,lattice,xyzstruct,&
          pairs_rr,pairs_atom,pairs_xyz,pairs_charge,&
          pair_lsta,pair_lstb,pair_lstc,pair_lste,lperiodic)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer i1,i2,i3                      ! internal
      integer index_store,flag_store        ! internal
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
      integer lste(listdim)                 ! internal 

      integer pair_lsta(listdim,2)          ! out
      integer pair_lstc(listdim)            ! out
      integer pair_lste(listdim)            ! out
      integer pairs_atom(listdim,2)         ! out
      integer pairs_charge(2,listdim)       ! out
      integer temp_lstc(listdim)            ! internal
      integer temp_c(listdim)               ! internal 

      real*8 lstb(listdim,4)               ! out
      real*8 absaxb,absaxc,absbxc          ! internal
      real*8 axb(3),axc(3),bxc(3)          ! internal
      real*8 cutoff
      real*8 rr                            ! internal
      real*8  lattice(3,3)                  ! in
      real*8  proja,projb,projc             ! internal
      real*8  rr_r(2)                       ! internal
      real*8  xrel1,xrel2,xrel3             ! internal
      real*8  xtemp,ytemp,ztemp             ! internal
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
!! find pairs 
!! This is a bit complicated because we don't want to count any pair twice
      num_pairs=1
      do i2=1,num_atoms ! loop over all atoms
        do i3=lsta(1,i2),lsta(2,i2) ! loop over all neighbors of this atom
          flag_store=3
!!          write(ounit,*)i2,lstc(i3),'charges ',zelem(i2),lste(i3)
          if(zelem(i2).lt.lste(i3))then  !! heteronuclear pair 
            flag_store=2
          elseif(zelem(i2).eq.lste(i3))then !! homonuclear pair
            if(idnint(dble(mult*xyzstruct(1,i2))).lt.idnint(dble(mult*lstb(i3,1)))) then
              flag_store=2
            elseif(idnint(dble(mult*xyzstruct(1,i2))).eq.idnint(dble(mult*lstb(i3,1))))then
              if(idnint(dble(mult*xyzstruct(2,i2))).lt.idnint(dble(mult*lstb(i3,2))))then
                flag_store=2
              elseif(idnint(dble(mult*xyzstruct(2,i2))).eq.idnint(dble(mult*lstb(i3,2))))then
                if(idnint(dble(mult*xyzstruct(3,i2))).lt.idnint(dble(mult*lstb(i3,3))))then
                  flag_store=2
                else
                  flag_store=1
                endif ! inner most
              else
                flag_store=1
              endif ! next most
            else
              flag_store=1
            endif ! next most
          else ! zelem(i2).gt.lste(i3) => rejected pair
            flag_store=1
          endif
          if(flag_store.eq.2)then ! keep this pair
!!            write(ounit,*)'Keeping pair'
!! store number of both atoms in pair
            pairs_atom(num_pairs,1)     = i2          ! pair index A
            pairs_atom(num_pairs,2)     = lstc(i3)    ! pair index B
!! store interatomic distance in pair
            pairs_rr(num_pairs)         = lstb(i3,4)  ! distance b/w A&B
!! store nuclear charges in pair
            pairs_charge(1,num_pairs)   = zelem(i2)   ! charge A
            pairs_charge(2,num_pairs)   = lste(i3)    ! charge B
!!            write(ounit,*)'pairs_charge ',num_pairs,pairs_charge(1,num_pairs),pairs_charge(2,num_pairs)
!! store positions of both atoms in pair
            pairs_xyz(num_pairs,1,1)    = xyzstruct(1,i2)   ! x of A
            pairs_xyz(num_pairs,1,2)    = xyzstruct(2,i2)   ! y of A
            pairs_xyz(num_pairs,1,3)    = xyzstruct(3,i2)   ! z of A
            pairs_xyz(num_pairs,2,1)    = lstb(i3,1)  ! x of B
            pairs_xyz(num_pairs,2,2)    = lstb(i3,2)  ! y of B
            pairs_xyz(num_pairs,2,3)    = lstb(i3,3)  ! z of B
            num_pairs                   = num_pairs + 1
          elseif(flag_store.eq.3)then
            write(ounit,*)'ERROR in neighborpair_para '
            stop
          endif
        enddo ! i3
      enddo ! i2   
      num_pairs = num_pairs - 1 ! correct counter

!! find neighbors of these pairs 
      index_store = 1
!!      do i1 = n_start,n_end ! loop over all pairs of this process
      do i1 = 1,num_pairs ! loop over all pairs 
        pair_lsta(i1,1) = index_store ! index of first neighbor of this pair
        do i3 =1,n_count ! loop over all atoms that could be neighbors 
          do i2 = 1,2 ! get distance of atom i3 to both atoms in pair 
            xx_r     =  pairs_xyz(i1,i2,1) - temp_x(i3)
            yy_r     =  pairs_xyz(i1,i2,2) - temp_y(i3) 
            zz_r     =  pairs_xyz(i1,i2,3) - temp_z(i3)
            rr_r(i2) = sqrt(xx_r**2 + yy_r**2 + zz_r**2)
          enddo ! i2
!! check if i3 is close to at least one of both atoms in the pair
          if((rr_r(1).lt.cutoff).or.(rr_r(2).lt.cutoff))then
!! check if i3 is not too close
            if((rr_r(1).gt.rmin).and.(rr_r(2).gt.rmin))then
!! store position of neighbor
              pair_lstb(index_store,1) =  temp_x(i3)    ! nbr x
              pair_lstb(index_store,2) =  temp_y(i3)    ! nbr y
              pair_lstb(index_store,3) =  temp_z(i3)    ! nbr z
!! store both distances of neighbor
              pair_lstb(index_store,4) =  rr_r(1)       ! dist b/w nbr & A
              pair_lstb(index_store,5) =  rr_r(2)       ! dist b/w nbr & B
!! store nuclear charge of neighbor
              pair_lste(index_store)   =  temp_c(i3)    ! nbr charge 
!! store atom number of neighbor
              pair_lstc(index_store)   =  temp_lstc(i3) ! nbr index
!! store index of last neighbor of this pair 
              pair_lsta(i1,2)          =  index_store
              index_store              =  index_store + 1
            endif
          endif
        enddo ! i3
      enddo ! i1
!!
      return
      end
