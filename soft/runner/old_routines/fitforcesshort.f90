!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fitting.f90
!!
      subroutine fitforcesshort(ounit,punit,nblock,npoints,idx,&
        nelem,num_weightsshort,max_num_atoms,num_funcvalues,&
        num_functions,nucelem,listdim,nforcegroup,&
        kseed,paramode,nnindex,function_type,&
        num_layersshort,windex,zelem_list,maxnodes_short,&
        num_atoms_list,nodes_short,elementindex,optmode,&
        kalmanthresholdf,kalmanlambda,kalmannue,&
        symfunction_list,weights_short,lattice_list,&
        rmse_force_s_ref,xyzstruct_list,xyzforce_list,&
        corrmatrix_list,kalgainmat_list,&
        funccutoff,maxcutoff,eta,zeta,rshift,lambda,&
        minvalue,maxvalue,forcernd,&
        actfunc_short,&
        lperiodic_list,lscalesym,ldebug)
!!
      use kalmandims
!!
      implicit none
!!
      integer i0,i1,i2,i3,i4
!!
      integer ounit                                    ! in
      integer punit                                    ! in
      integer nblock                                   ! in
      integer npoints                                  ! in
      integer idx(nblock)                              ! in
      integer nelem                                    ! in
      integer num_weightsshort                         ! in
      integer max_num_atoms                            ! in
      integer num_funcvalues                           ! in
      integer num_functions                            ! in
      integer num_layersshort                          ! in
      integer windex(2*num_layersshort)                ! in
      integer zelem_list(nblock,max_num_atoms)         ! in 
      integer zelem(max_num_atoms)                     ! internal 
      integer maxnodes_short                           ! in
      integer num_atoms_list(nblock)                   ! in
      integer num_atoms                                ! internal
      integer nodes_short(0:num_layersshort)           !
      integer elementindex(102)                        ! in
      integer optmode                                  ! in
      integer nnindex(num_functions,nelem,nelem)       ! in
      integer function_type(num_functions)             ! in
      integer nucelem(nelem)                           ! in
      integer listdim                                  ! in
      integer num_atoms_element(nelem)                 ! internal 
      integer nforcegroup                              ! in
      integer nforce                                   ! internal
      integer kseed                                    !
      integer paramode                                 ! in
!!
      real*8 kalmanthresholdf                                                !
      real*8 kalmanlambda(nelem)                                             !
      real*8 kalmannue                                                      !
      real*8 symfunction_list(num_funcvalues,max_num_atoms,nblock)          ! in
      real*8 symfunctiondummy(num_funcvalues,max_num_atoms)                 ! internal
      real*8 weights_short(num_weightsshort,nelem)                          ! in/out
      real*8 weights(num_weightsshort)                                      ! internal
      real*8 rmse_force_s_ref                                               ! in
      real*8 xyzstruct_list(3,max_num_atoms,nblock)                         ! in 
      real*8 xyzstruct(3,max_num_atoms)                                     ! internal 
      real*8 xyzforce_list(nblock,max_num_atoms,3)                          ! in
      real*8 xyzforce(max_num_atoms,3)                                      ! internal
      real*8 nnxyzforce(max_num_atoms,3)                                    ! internal
      real*8 corrmatrix_list(corrdim,nelem)                                 !
      real*8 kalgainmat_list(kaldim,nelem)                                  !
      real*8 lattice_list(3,3,nblock)                                       ! in
      real*8 minvalue(nelem,num_funcvalues)                                 ! in
      real*8 maxvalue(nelem,num_funcvalues)                                 ! in
      real*8 dsfuncdxyz(num_funcvalues,max_num_atoms,max_num_atoms,3)       ! internal
      real*8 strs(3,3,num_funcvalues,max_num_atoms)                         ! internal dummy
      real*8 funccutoff(num_functions)                                      ! in
      real*8 eta(num_functions)                                             ! in
      real*8 rshift(num_functions)                                          ! in
      real*8 lambda(num_functions)                                          ! in
      real*8 zeta(num_functions)                                            ! in
      real*8 maxcutoff                                                      ! in
      real*8 deshortdsfunc(max_num_atoms,num_funcvalues)                    ! dummy here 
      real*8 abserror                                                       ! internal
      real*8 error                                                          ! internal
      real*8 dfshortdw(num_weightsshort)                                    ! internal
      real*8 dfshortdwsum(num_weightsshort,nelem)                           ! internal
      real*8 errorsum(nelem)                                                ! internal
      real*8 ran0,z                                                         ! internal
      real*8 forcernd                                                       ! in
!!
      character*1 actfunc_short(num_layersshort)                            ! in
!!
      logical lperiodic_list(nblock)                                        ! in
      logical lperiodic                                                     ! internal
      logical lscalesym                                                     ! in
      logical ldebug                                                        ! in
!!
!! initializations
      xyzforce(:,:)        = 0.0d0
      nnxyzforce(:,:)      = 0.0d0 
      dfshortdw(:)         = 0.0d0
      dfshortdwsum(:,:)    = 0.0d0
      errorsum(:)          = 0.0d0
      nforce               = 0
      num_atoms_element(:) = 0
!!
!!
!! loop over all training structures 
      do i1=1,npoints
!!
!!        write(ounit,*)'fitforcesshort point ',i1
!!
!! get the arrays for training point i1 
        num_atoms        = num_atoms_list(idx(i1))
        zelem(:)         = zelem_list(idx(i1),:)
        lperiodic        = lperiodic_list(idx(i1))
        xyzforce(:,:)    = xyzforce_list(idx(i1),:,:)
!!
!! get dsfuncdxyz:
!! this is independent of the weights and needs to be done only once for each point i1
!! in principle dsfuncdxyz could be precalculated, but needs too much storage on disk
!! Caution: symmetry functions from this subroutine cannot be used because they are not scaled
        call calconefunction(ounit,num_functions,num_funcvalues,&
          max_num_atoms,nelem,num_atoms,zelem,listdim,&
          nnindex,elementindex,&
          function_type,lattice_list(1,1,idx(i1)),&
          xyzstruct_list(1,1,idx(i1)),xyzforce,symfunctiondummy,maxcutoff,&
          funccutoff,eta,rshift,lambda,zeta,dsfuncdxyz,strs,&
          lperiodic,.true.,.false.,ldebug)
!!
!! scale dsfuncdxyz
        if(lscalesym)then
          call scaledsfunc(ounit,max_num_atoms,num_funcvalues,&
          nelem,num_atoms,elementindex,minvalue,maxvalue,&
          zelem,dsfuncdxyz,strs,&
          ldebug)
        endif ! lscalesym
!!
!! calculate the number of atoms per element, because we don't want to 
!! update weights for elements not being present in current example
        num_atoms_element(:) = 0
        do i2=1,num_atoms
          do i3=1,3 ! loop over x,y and z
            num_atoms_element(elementindex(zelem(i2)))=&
              num_atoms_element(elementindex(zelem(i2)))+1
          enddo ! i3
        enddo ! i2
!!
!! calculate the error of the forces individually for each atom and component
        do i2=1,num_atoms
          do i3=1,3 ! loop over x,y and z
!!
!! decide randomly if this force should be used for update
            z=ran0(kseed)
            if(z.le.forcernd)then ! use force for weight update
!!
!! calculate the NN forces nnxyzforce for current weight values
              nnxyzforce(:,:)=0.0d0
!!
!! this subroutine calculates only the one force that is needed (specified by i2,i3), not the full array nnxyzforce
              call getoneshortforce(ounit,i2,i3,num_funcvalues,&
                max_num_atoms,num_atoms,num_weightsshort,&
                nelem,num_layersshort,nodes_short,maxnodes_short,&
                zelem,elementindex,windex,&
                weights_short,symfunction_list(1,1,idx(i1)),&
                dsfuncdxyz,nnxyzforce,actfunc_short,ldebug)
!!
!! calculate the error for the current force component
              abserror=abs(xyzforce(i2,i3)-nnxyzforce(i2,i3))
              error=xyzforce(i2,i3)-nnxyzforce(i2,i3)
              error=-1.d0*error          ! this must be used!
!!              write(ounit,'(a,f20.10)')'error forceupdate ',error
!!
!!
              if(abserror.gt.kalmanthresholdf*rmse_force_s_ref)then
!! calculate the derivative of the short range forces with respect to the weights
!! dfshortdw is for one atom only, but depends on all atoms
                nforce=nforce+1
!!
!! calculate dfshortdw
                dfshortdw(:)=0.0d0
                call getdfshortdw(ounit,nelem,i2,i3,num_weightsshort,&
                  zelem,elementindex,max_num_atoms,num_atoms,&
                  maxnodes_short,&
                  windex,num_layersshort,num_funcvalues,nodes_short,&
                  symfunction_list(1,1,idx(i1)),dsfuncdxyz,&
                  weights_short,dfshortdw,&
                  actfunc_short,&
                  ldebug)
!!
                dfshortdwsum(:,elementindex(zelem(i2)))=&
                  dfshortdwsum(:,elementindex(zelem(i2)))+dfshortdw(:)
                errorsum(elementindex(zelem(i2)))=&
                  errorsum(elementindex(zelem(i2)))+error
!!
                if((mod(nforce,nforcegroup).eq.0).or.&
                  ((i2.eq.num_atoms).and.(i3.eq.3).and.(i1.eq.npoints)))then
!!
                  do i4=1,nelem
!!
                    if(num_atoms_element(i4).gt.0)then
!!
!! JB: general note on normalization of dfshortdw and error: 
!! It seems to work best if both are divided by num_atoms_element
!!
!! normalize by the number of atoms in this element to reduce the importance of the forces with respect to the energies, should we do this?
                      dfshortdw(:)=dfshortdwsum(:,i4)/dble(num_atoms_element(i4)) ! CHECK!
!!                      dfshortdw(:)=dfshortdwsum(:,i4) ! CHECK!
!! if we average over several force derivatives we should normalize here 
                      dfshortdw(:)=dfshortdw(:)/dble(nforce) ! CHECK!
!!
                      error=errorsum(i4)/dble(num_atoms_element(i4)) ! CHECK!
!!                      error=errorsum(i4) ! CHECK!
!! if we average over several force derivatives we should normalize here 
                      error=error/dble(nforce) ! CHECK!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! update the weights:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      if(optmode.eq.1)then ! Kalman filter
!!
                        call updatekalman(ounit,kaldim,corrdim,&
                          num_weightsshort,kalmanlambda(i4),kalmannue,&
                          weights_short(1,i4),dfshortdw,&
                          kalgainmat_list(1,i4),corrmatrix_list(1,i4),error,ldebug)
!!
                      elseif(optmode.eq.2)then ! conjugate gradient'
                        write(ounit,*)'CG optimization not yet implemented in fitforcesshort' 
                        stop !'
!!
                      elseif(optmode.eq.3)then ! steepest descent'
                        weights(:)   =weights_short(:,elementindex(zelem(i2)))
!!
                        call updatesteepest(ounit,&
                          num_weightsshort,weights,dfshortdw,error,ldebug)
!!
                        weights_short(:,elementindex(zelem(i2)))=weights(:)
                      endif ! optmode
!!
                    endif !(num_atoms_element(i4).gt.0)then
!!
                  enddo ! i4
!!
!! reinitializations for next update
                  nforce           = 0
                  dfshortdwsum(:,:)= 0.0d0
                  errorsum(:)      = 0.0d0
!!
                endif !(mod(nforce,nforcegroup).eq.0)then
!!
              else ! point is not used for short update 
!!
              endif !(abserror.gt.kalmanthresholdf*rmse_force_s_ref)
!!            
            else ! don't use force for weight update
!!
            endif ! z.le.forcernd      ! use force for weight update
!!
          enddo ! i3 loop over x,y and z
!!
        enddo ! i2 loop over atoms
!!
      enddo ! i1 loop over training structures
!!
      return
      end
