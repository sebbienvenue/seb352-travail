!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - calconefunction_para.f90
!!
      subroutine atomsymfunction1(i1,i2,iindex,natoms,natomsdim,nelem,&
        max_num_atoms,max_num_neighbors,invneighboridx,&
        jcount,listdim,lsta,lstc,lste,symelement_local,maxnum_funcvalues_local,&
        lstb,funccutoff_local,pi,xyzstruct,symfunction,dsfuncdxyz,strs,&
        ldoforces,ldostress)
!!
      use fileunits
!!
      implicit none
!!
      integer i1,i2,i3
      integer listdim
      integer nelem 
      integer natomsdim 
      integer natoms 
      integer iindex
      integer max_num_atoms 
      integer max_num_neighbors 
      integer maxnum_funcvalues_local 
      integer symelement_local(maxnum_funcvalues_local,2,nelem)          ! in
      integer lsta(2,max_num_atoms)
      integer lstc(listdim)
      integer lste(listdim)
      integer n
      integer jcount 
      integer invneighboridx(natoms,max_num_atoms)
!!
      real*8 strs(3,3,maxnum_funcvalues_local,natomsdim)
      real*8 lstb(listdim,4)
      real*8 funccutoff_local(maxnum_funcvalues_local,nelem)
      real*8 rij
      real*8 symfunction(maxnum_funcvalues_local,natomsdim)
      real*8 xyzstruct(3,max_num_atoms)
      real*8 deltaxj,deltayj,deltazj
      real*8 drijdxi, drijdyi, drijdzi
      real*8 drijdxj, drijdyj, drijdzj
      real*8 temp1
      real*8 pi
      real*8 dsfuncdxyz(maxnum_funcvalues_local,natomsdim,0:max_num_neighbors,3) 
!!
      logical ldoforces
      logical ldostress
!!
      do i3=lsta(1,i1),lsta(2,i1) ! loop over neighbors
        if(symelement_local(i2,1,iindex).eq.lste(i3))then
!!
          rij=lstb(i3,4)
          if(rij.le.funccutoff_local(i2,iindex))then
            n=lstc(i3)
!!
            symfunction(i2,i1)=symfunction(i2,i1)&
              +0.5d0*(dcos(pi*rij/funccutoff_local(i2,iindex))+1.d0)
!!
            deltaxj=-1.d0*(xyzstruct(1,jcount)-lstb(i3,1))
            deltayj=-1.d0*(xyzstruct(2,jcount)-lstb(i3,2))
            deltazj=-1.d0*(xyzstruct(3,jcount)-lstb(i3,3))
            drijdxi=-deltaxj/rij
            drijdyi=-deltayj/rij
            drijdzi=-deltazj/rij
            drijdxj=-1.d0*drijdxi
            drijdyj=-1.d0*drijdyi
            drijdzj=-1.d0*drijdzi
!!
            if(ldoforces)then
              temp1=-0.5d0*dsin(pi*rij/funccutoff_local(i2,iindex))&
                *pi/funccutoff_local(i2,iindex)
!! Calculation of derivatives for forces 
!! dsfunc/dx
              dsfuncdxyz(i2,i1,invneighboridx(i1,jcount),1)=dsfuncdxyz(i2,i1,invneighboridx(i1,jcount),1)&
                +(temp1*drijdxi)
              dsfuncdxyz(i2,i1,invneighboridx(i1,n),1) =dsfuncdxyz(i2,i1,invneighboridx(i1,n),1)&
                +(temp1*drijdxj)
!! dsfunc/dy
              dsfuncdxyz(i2,i1,invneighboridx(i1,jcount),2)=dsfuncdxyz(i2,i1,invneighboridx(i1,jcount),2)&
                +(temp1*drijdyi)
              dsfuncdxyz(i2,i1,invneighboridx(i1,n),2) =dsfuncdxyz(i2,i1,invneighboridx(i1,n),2)&
                +(temp1*drijdyj)
!! dsfunc/dz
              dsfuncdxyz(i2,i1,invneighboridx(i1,jcount),3)=dsfuncdxyz(i2,i1,invneighboridx(i1,jcount),3)&
                +(temp1*drijdzi)
              dsfuncdxyz(i2,i1,invneighboridx(i1,n),3) =dsfuncdxyz(i2,i1,invneighboridx(i1,n),3)&
                +(temp1*drijdzj)
!!
              if(ldostress)then
!! Calculation of derivatives for stress
!! dsfunc/dx
                strs(1,1,i2,i1)=strs(1,1,i2,i1)&
                  +deltaxj*(drijdxj*temp1)
                strs(2,1,i2,i1)=strs(2,1,i2,i1)&
                  +deltayj*(drijdxj*temp1)
                strs(3,1,i2,i1)=strs(3,1,i2,i1)&
                  +deltazj*(drijdxj*temp1)
!! dsfunc/dy
                strs(1,2,i2,i1)=strs(1,2,i2,i1)&
                  +deltaxj*(drijdyj*temp1)
                strs(2,2,i2,i1)=strs(2,2,i2,i1)&
                  +deltayj*(drijdyj*temp1)
                strs(3,2,i2,i1)=strs(3,2,i2,i1)&
                  +deltazj*(drijdyj*temp1)
!! dsfunc/dz
                strs(1,3,i2,i1)=strs(1,3,i2,i1)&
                  +deltaxj*(drijdzj*temp1)
                strs(2,3,i2,i1)=strs(2,3,i2,i1)&
                  +deltayj*(drijdzj*temp1)
                strs(3,3,i2,i1)=strs(3,3,i2,i1)&
                  +deltazj*(drijdzj*temp1)
!!
              endif ! ldostress 
            endif ! ldoforces
          endif ! rij.le.funccutoff_local(i2,iindex)
        endif ! (symelement_local(i2,1,iindex).eq.lste(i3))then
      enddo ! i3

      return
      end
