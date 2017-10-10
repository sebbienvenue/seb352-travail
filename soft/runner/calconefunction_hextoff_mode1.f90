!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose: calculate the pair symmetry functions for one structure

!! FIXME: parallel case probably not yet working

!! called by:
!! - optimize_short_combinedpair.f90
!! - prediction_hextoff.f90 
!! - calcpairfunctions.f90
!! - getallshortforcespair.f90
!!
!! CMH
!! Unlike the other versions, this section requires the preprocessing of the structures for there
!! to be a prediction performed/analysis of the structure for training.
!! For training we have preprocessed already by having all training examples in the correct orientation
!! but for mode 3 this will be different. 
!!
!!
      subroutine calconefunction_hextoff_mode1(n_start,&
         istruct,num_atoms,zelem,xyzstruct,symfunction_hextoff_local)
!!
      use fileunits
      use globaloptions
      use symfunctions
      use nnham
      use nnflags
      use nnconstants
!!
      implicit none
!!
      integer n_start                                                          ! in
      integer i1,i2,i3                                                         ! internal
      integer num_atoms                                                        ! in
      integer istruct                                                          ! in, number of structure  
      integer zelem(max_num_atoms)                                             ! in

!!      integer triplet_charge(2,listdim)                                        ! in
      integer nuc1,nuc2,nuc3                                                   ! internal
      integer itriplet                                                         ! internal

      real*8 dsfuncdxyz_hextoff(maxnum_funcvalues_hextoff,max_num_atoms,ntriplets,3) ! out
      real*8 symfunction_hextoff_local(maxnum_funcvalues_hextoff)              ! in   check new dimension
      real*8 xyzstruct(3,max_num_atoms)                                        ! in

      integer oatom,zatom,tatom                                                 ! internal

      real*8 deltaxij,deltayij,deltazij                                        ! internal
      real*8 deltaxik,deltayik,deltazik                                        ! internal
      real*8 deltaxjk,deltayjk,deltazjk                                        ! internal
      real*8 rij                                                               ! internal
      real*8 rjk                                                               ! internal
      real*8 rik                                                               ! internal
      real*8 dexpijdxi,dexpijdyi,dexpijdzi
      real*8 dexpijdxj,dexpijdyj,dexpijdzj
      real*8 dexpijdxk,dexpijdyk,dexpijdzk
      real*8 dexpjkdxi,dexpjkdyi,dexpjkdzi
      real*8 dexpjkdxj,dexpjkdyj,dexpjkdzj
      real*8 dexpjkdxk,dexpjkdyk,dexpjkdzk
      real*8 dexpikdxi,dexpikdyi,dexpikdzi
      real*8 dexpikdxj,dexpikdyj,dexpikdzj
      real*8 dexpikdxk,dexpikdyk,dexpikdzk
      real*8 expij,expjk,expik,expijk
      real*8 drijdxi, drijdyi, drijdzi                                         ! internal
      real*8 drijdxj, drijdyj, drijdzj                                         ! internal
      real*8 drijdxk, drijdyk, drijdzk                                         ! internal
      real*8 drikdxi, drikdyi, drikdzi
      real*8 drikdxj, drikdyj, drikdzj
      real*8 drikdxk, drikdyk, drikdzk
      real*8 drjkdxi, drjkdyi, drjkdzi
      real*8 drjkdxj, drjkdyj, drjkdzj
      real*8 drjkdxk, drjkdyk, drjkdzk
      real*8 dfcutijdxi,dfcutijdyi,dfcutijdzi
      real*8 dfcutijdxj,dfcutijdyj,dfcutijdzj
      real*8 dfcutijdxk,dfcutijdyk,dfcutijdzk
      real*8 dfcutikdxi,dfcutikdyi,dfcutikdzi
      real*8 dfcutikdxj,dfcutikdyj,dfcutikdzj
      real*8 dfcutikdxk,dfcutikdyk,dfcutikdzk
      real*8 dfcutjkdxi,dfcutjkdyi,dfcutjkdzi
      real*8 dfcutjkdxj,dfcutjkdyj,dfcutjkdzj
      real*8 dfcutjkdxk,dfcutjkdyk,dfcutjkdzk
      real*8 fcutij
      real*8 fcutik
      real*8 fcutjk
      real*8 temp1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!! initialization
      symfunction_hextoff_local(:)       = 0.0d0
      dsfuncdxyz_hextoff(:,:,:,:)= 0.0d0
!!      strs_hextoff(:,:,:,:)      = 0.0d0
!!      lrmin                   =.true.
      oatom = 0
      zatom = 0
      tatom = 0

      
!!
!! Prepare values used for calculations of the symmetry functions and abbreviations
!!
!!      write(*,*) hextoff_training_triplet
!!      write(*,*) tripletindex(hextoff_training_triplet(1),hextoff_training_triplet(2),hextoff_training_triplet(3))
!!      itriplet = tripletindex(hextoff_training_triplet(1),hextoff_training_triplet(2),hextoff_training_triplet(3))
      itriplet = 1 ! always in modes 1 and 2 
      nuc1=hextoff_training_triplet(1)
      nuc2=hextoff_training_triplet(2)
      nuc3=hextoff_training_triplet(3)
!!      write(*,*) element(elementindex(elemtriplet(itriplet,1))),&
!!                 element(elementindex(elemtriplet(itriplet,2))),&
!!                 element(elementindex(elemtriplet(itriplet,3))),&
!!                 elementindex(elemtriplet(itriplet,1)),&
!!                 elementindex(elemtriplet(itriplet,2)),&
!!                 elementindex(elemtriplet(itriplet,3)),&
!!                 elemtriplet(itriplet,1),&
!!                 elemtriplet(itriplet,2),&
!!                 elemtriplet(itriplet,3)
!!
!!  Determine origin atom, z axis atom and triplet atom (structure may have been read in in a irregular manner).
!!
      do i1=1,num_atoms
        if(((abs(xyzstruct(1,i1))).le.(0.0E-10))&
          .and.(abs((xyzstruct(2,i1))).le.(0.0E-10))&
          .and.(abs((xyzstruct(3,i1))).le.(0.0E-10))) then
          oatom = i1
        elseif((abs(xyzstruct(1,i1)).le.0.0E-14)&
          .and.(abs(xyzstruct(2,i1)).le.0.0E-14)&
          .and.(xyzstruct(3,i1).gt.0.5)) then
          if(zatom.ne.0) then
            if(xyzstruct(3,i1).lt.xyzstruct(3,zatom)) then
              tatom = zatom
              zatom = i1        
            endif
          else
            zatom = i1
          endif
        else
          tatom = i1
        endif
      enddo
!!
!!   Calculate the distance in the pair atoms
!!
!!      write(*,*) 'atom order in structure',oatom,zatom,tatom                      
      deltaxij = (xyzstruct(1,oatom) - xyzstruct(1,zatom))
      deltayij = (xyzstruct(2,oatom) - xyzstruct(2,zatom))
      deltazij = (xyzstruct(3,oatom) - xyzstruct(3,zatom))
      rij = deltaxij**2 + deltayij**2 + deltazij**2
      rij = sqrt(rij)
      deltaxik = (xyzstruct(1,oatom) - xyzstruct(1,tatom))
      deltayik = (xyzstruct(2,oatom) - xyzstruct(2,tatom))
      deltazik = (xyzstruct(3,oatom) - xyzstruct(3,tatom))
      rik = deltaxik**2 + deltayik**2 + deltazik**2
      rik = sqrt(rik)
      deltaxjk = (xyzstruct(1,zatom) - xyzstruct(1,tatom))
      deltayjk = (xyzstruct(2,zatom) - xyzstruct(2,tatom))
      deltazjk = (xyzstruct(3,zatom) - xyzstruct(3,tatom))
      rjk = deltaxjk**2 + deltayjk**2 + deltazjk**2
      rjk = sqrt(rjk)
!!
!!   Distances have been precalculated so now calculate the symmetry functions
!!

!! loop over all symmetry functions
      do i1=1,num_funcvalues_hextoff(itriplet)
!!        write(ounit,*)'Calculating symfunction ',i1
        if(function_type_hextoff(i1,itriplet).eq.1)then 
!!          write(*,*) 'function 1 starts'
          if(rij.le.funccutoff_hextoff(i1,itriplet))then
            drijdxi=(xyzstruct(1,oatom)-xyzstruct(1,zatom))/rjk
            drijdyi=(xyzstruct(2,oatom)-xyzstruct(2,zatom))/rjk
            drijdzi=(xyzstruct(3,oatom)-xyzstruct(3,zatom))/rjk
            drijdxj=-1.d0*drijdxi
            drikdyj=-1.d0*drijdyi
            drijdzj=-1.d0*drijdzi
            drijdxk=0.0d0
            drijdyk=0.0d0
            drijdzk=0.0d0
!!
            fcutij=0.5d0*(dcos(pi*rij/funccutoff_hextoff(i1,itriplet))+1.d0)
            temp1=0.5d0*(-dsin(pi*rij/funccutoff_hextoff(i1,itriplet)))*&
                        (pi/funccutoff_hextoff(i1,itriplet))
            dfcutijdxi=temp1*drijdxi
            dfcutijdyi=temp1*drijdyi
            dfcutijdzi=temp1*drijdzi
            dfcutijdxj=-1.d0*dfcutijdxi
            dfcutijdyj=-1.d0*dfcutijdyi
            dfcutijdzj=-1.d0*dfcutijdzi
            dfcutijdxk=0.0d0
            dfcutijdyk=0.0d0
            dfcutijdzk=0.0d0
!!
            expij=(dexp(-1.d0*eta_hextoff(i1,itriplet)* (rij)**2))
!!
            dexpijdxi=(-2.0d0*expij*(xyzstruct(1,oatom)-xyzstruct(1,zatom))*eta_hextoff(i1,itriplet))
            dexpijdyi=(-2.0d0*expij*(xyzstruct(2,oatom)-xyzstruct(2,zatom))*eta_hextoff(i1,itriplet))
            dexpijdzi=(-2.0d0*expij*(xyzstruct(3,oatom)-xyzstruct(3,zatom))*eta_hextoff(i1,itriplet))
            dexpijdxj= -1.0d0*dexpijdxi
            dexpijdyj= -1.0d0*dexpijdyi
            dexpijdzj= -1.0d0*dexpijdzi
            dexpijdxk=  0.0d0
            dexpijdyk=  0.0d0
            dexpijdzk=  0.0d0
          endif ! rjk
!!
          symfunction_hextoff_local(i1)=symfunction_hextoff_local(i1)+ expij * fcutij
!!          write(*,*) symfunction_hextoff_local(i1), expij, rij
!!          write(*,*) symfunction_hextoff_local(:)
!!          write(*,*) num_funcvalues_hextoff(itriplet), maxnum_funcvalues_hextoff

        elseif(function_type_hextoff(i1,itriplet).eq.2)then
!!          write(*,*) 'function 2 starts'
          if(nuc1.eq.nuc2)then
            if((rij.le.funccutoff_hextoff(i1,itriplet)).and.&
               (rik.le.funccutoff_hextoff(i1,itriplet)).and.&
               (rjk.le.funccutoff_hextoff(i1,itriplet)))then
!!
               fcutij=0.5d0*(dcos(pi*rij/funccutoff_hextoff(i1,itriplet))+1.d0)
               fcutik=0.5d0*(dcos(pi*rik/funccutoff_hextoff(i1,itriplet))+1.d0)
               fcutjk=0.5d0*(dcos(pi*rjk/funccutoff_hextoff(i1,itriplet))+1.d0)
               expijk=(dexp(-1.d0*eta_hextoff(i1,itriplet)* (rik+rjk)**2))
               symfunction_hextoff_local(i1)=symfunction_hextoff_local(i1)+ expijk * fcutij *(fcutik+fcutjk)
!!          write(*,*) symfunction_hextoff_local(i1), expij, rij
!!          write(*,*) symfunction_hextoff_local(:)
!!          write(*,*) num_funcvalues_hextoff(itriplet), maxnum_funcvalues_hextoff

            endif
          else
            if((rij.le.funccutoff_hextoff(i1,itriplet)).and.&
               (rik.le.funccutoff_hextoff(i1,itriplet)).and.&
               (rjk.le.funccutoff_hextoff(i1,itriplet)))then
!!
               fcutij=0.5d0*(dcos(pi*rij/funccutoff_hextoff(i1,itriplet))+1.d0)
               fcutik=0.5d0*(dcos(pi*rik/funccutoff_hextoff(i1,itriplet))+1.d0)
               fcutjk=0.5d0*(dcos(pi*rjk/funccutoff_hextoff(i1,itriplet))+1.d0)
               expijk=(dexp(-1.d0*eta_hextoff(i1,itriplet)* rik**2))
               symfunction_hextoff_local(i1)=symfunction_hextoff_local(i1)+ expijk * fcutij * fcutik
            endif
          endif       
!!
        elseif(function_type_hextoff(i1,itriplet).eq.3)then
!!          write(*,*) 'function 3 starts'
          if(nuc1.eq.nuc2)then
            if((rij.le.funccutoff_hextoff(i1,itriplet)).and.&
               (rik.le.funccutoff_hextoff(i1,itriplet)).and.&
               (rjk.le.funccutoff_hextoff(i1,itriplet)))then
!!
               fcutij=0.5d0*(dcos(pi*rij/funccutoff_hextoff(i1,itriplet))+1.d0)
               fcutik=0.5d0*(dcos(pi*rik/funccutoff_hextoff(i1,itriplet))+1.d0)
               fcutjk=0.5d0*(dcos(pi*rjk/funccutoff_hextoff(i1,itriplet))+1.d0)
               expijk=(dexp(-1.d0*eta_hextoff(i1,itriplet)* (rij-rjk)**2))
               symfunction_hextoff_local(i1)=symfunction_hextoff_local(i1)+ expijk * fcutij* (fcutik+fcutjk)
!!          write(*,*) symfunction_hextoff_local(i1), expij, rij
!!          write(*,*) symfunction_hextoff_local(:)
!!          write(*,*) num_funcvalues_hextoff(itriplet), maxnum_funcvalues_hextoff
            endif
          else
            if((rij.le.funccutoff_hextoff(i1,itriplet)).and.&
               (rik.le.funccutoff_hextoff(i1,itriplet)).and.&
               (rjk.le.funccutoff_hextoff(i1,itriplet)))then
!!
               fcutij=0.5d0*(dcos(pi*rij/funccutoff_hextoff(i1,itriplet))+1.d0)
               fcutik=0.5d0*(dcos(pi*rik/funccutoff_hextoff(i1,itriplet))+1.d0)
               fcutjk=0.5d0*(dcos(pi*rjk/funccutoff_hextoff(i1,itriplet))+1.d0)
               expijk=(dexp(-1.d0*eta_hextoff(i1,itriplet)* rjk**2))
               symfunction_hextoff_local(i1)=symfunction_hextoff_local(i1)+ expijk * fcutij * fcutjk
!!          write(*,*) symfunction_hextoff_local(i1), expij, rij
!!          write(*,*) symfunction_hextoff_local(:)
!!          write(*,*) num_funcvalues_hextoff(itriplet), maxnum_funcvalues_hextoff
            endif
          endif
          
!!
!!             call hextoffsymfunction3(i1,i2,jcount,num_atoms,listdim,iindex,&
!!                jindex,&
!!                pair_lstc,maxnum_funcvalues_hextoff,npairs,ntriplets,max_num_pairs,&
!!                max_num_atoms,pairs_atom,pair_lste,symelement_hextoff,&
!!                pair_lsta,pair_lstb,pairs_xyz,eta_hextoff,funccutoff_hextoff,strs_hextoff,&
!!                dsfuncdxyz_hextoff,symfunction_hextoff,pi,ldoforces_local,ldostress)
!!
        else
          write(ounit,*)'ERROR mode1: Unknown hextoff symmetry function type ',&
            function_type_hextoff(i1,itriplet)
          stop
        endif              
      enddo
!!
!!      write(*,*) symfunction_hextoff_local(:)
!!      write(*,*) 'End calconefunction_hextoff'
      return
      end
