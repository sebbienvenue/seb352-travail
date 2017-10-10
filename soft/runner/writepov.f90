!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: 
!! - prediction.f90
!!
!! execution of povray with
!! 'povray povray.ini +FP'
!! 'povray povray.ini +FP +L/usr/share/povray-3.5/include'
!!
      subroutine writepov(num_atoms,zelem,&
           lattice,xyzstruct,lperiodic)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer num_atoms                ! in
      integer zelem(max_num_atoms)     ! in
      integer i1,i2                    ! internal
!! 
      real*8 lattice(3,3)               ! in
      real*8 xyzstruct(3,max_num_atoms) ! in
      real*8 xyztemp(3,max_num_atoms)   ! internal
      real*8 cmsx                       ! internal
      real*8 cmsy                       ! internal
      real*8 cmsz                       ! internal
      real*8 a1                         ! internal
      real*8 a2                         ! internal
      real*8 a3                         ! internal
      real*8 b1                         ! internal
      real*8 b2                         ! internal
      real*8 b3                         ! internal
      real*8 c1                         ! internal
      real*8 c2                         ! internal
      real*8 c3                         ! internal
      real*8 radcyl                     ! internal
      real*8 radbond                    ! internal
      real*8 r                          ! internal
      real*8 rminimal                       ! internal
      real*8 xtemp                      ! internal
      real*8 ytemp                      ! internal
      real*8 ztemp                      ! internal
      real*8 rcut                       ! internal
      real*8 radatom(102)               ! internal
!!
      character*7 atomcolor(102)               ! internal
!!
      logical lperiodic                ! in
!!
!! initializations
      cmsx=0.0d0
      cmsy=0.0d0
      cmsz=0.0d0
      radcyl=0.2d0  ! thickness of box
      radbond=0.3d0 ! thickness of bonds 
      rcut=4.0d0    ! bond criterion in Bohr 
      radatom(:)=1.5d0 ! same radius for all elements
      atomcolor(:)='colorSi'
      xyztemp(:,:)=xyzstruct(:,:) ! internal copy for shifting
!! abbreviations
      a1=lattice(1,1)
      a2=lattice(1,2)
      a3=lattice(1,3)
      b1=lattice(2,1)
      b2=lattice(2,2)
      b3=lattice(2,3)
      c1=lattice(3,1)
      c2=lattice(3,2)
      c3=lattice(3,3)
!!
!! calculate the center of mass
      if(lperiodic)then
        do i1=1,3
          cmsx=cmsx+lattice(i1,1)/2.0d0
          cmsy=cmsy+lattice(i1,2)/2.0d0
          cmsz=cmsz+lattice(i1,3)/2.0d0
        enddo
      else ! not periodic
        do i1=1,num_atoms
          cmsx=cmsx+xyztemp(1,i1)
          cmsy=cmsy+xyztemp(2,i1)
          cmsz=cmsz+xyztemp(3,i1)
        enddo  ! i1
        cmsx=cmsx/dble(num_atoms)
        cmsy=cmsy/dble(num_atoms)
        cmsz=cmsz/dble(num_atoms)
      endif ! lperiodic
!!
!! remove the center of mass
      do i1=1,num_atoms
        xyztemp(1,i1)=xyztemp(1,i1)-cmsx
        xyztemp(2,i1)=xyztemp(2,i1)-cmsy
        xyztemp(3,i1)=xyztemp(3,i1)-cmsz
      enddo
!!
      open(povunit,file='povray.ini',form='formatted',status='replace')
      write(povunit,*)'Input_File_Name=runner.pov'
      write(povunit,*)' +A0.1'
      write(povunit,*)' +UL'
      write(povunit,*)' +UV'
      write(povunit,*)' +Q9'
      write(povunit,*)' +W600 +H400'
      close(povunit)
!!
      open(povunit,file='runner.pov',form='formatted',status='replace')
!!
!! write header for povray
      write(povunit,*)'#include "colors.inc"'
!!      write(povunit,*)'#include "myheader.inc"'
!!      write(povunit,*)'#include "mycolors.inc"'
!!
!! write general header
      write(povunit,*)'//Camera definition. Default uses RIGHT hand system (negative right vector).'
      write(povunit,*) 'camera' !'
      write(povunit,*) '{ orthographic'
      write(povunit,*) 'location  <250,30,30>'
      write(povunit,*) 'direction <-1,0,0>'
      write(povunit,*) 'look_at   <0,0,0>'
      write(povunit,*) 'sky       <0,1,0>'
      write(povunit,*) 'up        <0,1,0>'
      write(povunit,*) '// right     <1.33333333333333, 0.0, 0.0>'
      write(povunit,*) 'angle 10'
      write(povunit,*) '}'
!!
      write(povunit,*) '//  Object: PointLightObject - Light 1'
      write(povunit,*) 'light_source  {  <100,50,50> color White'
      write(povunit,*) '}'
      write(povunit,*) '//  Object: PointLightObject - Light 2'
      write(povunit,*) 'light_source { <50,-10,-10> color White'
      write(povunit,*) '}'
!!
!!      write(povunit,*) '//texture applied to all atoms and bond cylinders'
!!      write(povunit,*) '#default {'
!!      write(povunit,*) 'texture {'
!!      write(povunit,*) '  finish'
!!      write(povunit,*) ' {  diffuse 0.6 brilliance 1 ambient 0.2 '
!!      write(povunit,*) 'reflection 0 phong 2 phong_size 90 '
!!      write(povunit,*) ' crand 0' 
!!      write(povunit,*) ' specular 0  roughness 0.05'
!!      write(povunit,*) ' }'
!!      write(povunit,*) '}}'
!!      write(povunit,*) ' '
!!
!! define the atomic colors and atomic radii
      radatom(1)  =1.0d0
      atomcolor(1)='colorH '
      write(povunit,*) '#declare colorH ='
      write(povunit,*) 'texture {'
      write(povunit,*) 'pigment { color rgb<0.75,0.75,0.75> }'
      write(povunit,*) 'finish {  diffuse 0.6 brilliance 1 ambient 0.2 '
      write(povunit,*) 'reflection 0 phong 2 phong_size 90 crand 0 '
      write(povunit,*) 'specular 0  roughness 0.05 }}'
      write(povunit,*) ' '
!!
      radatom(8)  =1.0d0
      atomcolor(8)='colorO' 
      write(povunit,*) '#declare colorO ='
      write(povunit,*) 'texture {'
      write(povunit,*) 'pigment { color rgb<0.9,0.2,0.2> }'
      write(povunit,*) 'finish {  diffuse 0.6 brilliance 1 ambient 0.2 '
      write(povunit,*) 'reflection 0 phong 2 phong_size 90 crand 0 '
      write(povunit,*) 'specular 0  roughness 0.05 }}'
      write(povunit,*) ' '
!!
      radatom(14)  =1.0d0
      atomcolor(14)='colorSi'
      write(povunit,*) '#declare colorSi ='
      write(povunit,*) 'texture {'
      write(povunit,*) 'pigment { color rgb<0.75,0.75,0.25> }'
      write(povunit,*) 'finish {  diffuse 0.6 brilliance 1 ambient 0.2 '
      write(povunit,*) 'reflection 0 phong 2 phong_size 90 crand 0 '
      write(povunit,*) 'specular 0  roughness 0.05 }}'
      write(povunit,*) ' '
!!
      radatom(29)  =1.0d0
      atomcolor(29)='colorCu'
      write(povunit,*) '#declare colorCu ='
      write(povunit,*) 'texture {'
      write(povunit,*) 'pigment { color rgb<0.5,0.0,0.0> }'
      write(povunit,*) 'finish {  diffuse 0.6 brilliance 1 ambient 0.2 '
      write(povunit,*) 'reflection 0 phong 2 phong_size 90 crand 0 '
      write(povunit,*) 'specular 0  roughness 0.05 }}'
      write(povunit,*) ' '
!!
      radatom(30)  =1.0d0
      atomcolor(30)='colorZn'
      write(povunit,*) '#declare colorZn ='
      write(povunit,*) 'texture {'
      write(povunit,*) 'pigment { color rgb<0.75,0.75,0.75> }'
      write(povunit,*) 'finish {  diffuse 0.6 brilliance 1 ambient 0.2 '
      write(povunit,*) 'reflection 0 phong 2 phong_size 90 crand 0 '
      write(povunit,*) 'specular 0  roughness 0.05 }}'
      write(povunit,*) ' '
!!
!! write union
!!      write(povunit,*)'#declare unitcell = object{union{'

      if(lperiodic)then
!! write the box boundary in povray file
!! 0-a
      write(povunit,*)' '
      write(povunit,*)'// simulation box edge '
      write(povunit,*)'object {cylinder{'
      write(povunit,*)'<',-cmsx,',',-cmsy,',',-cmsz,'>'
      write(povunit,*)'<',a1-cmsx,',',a2-cmsy,',',a3-cmsz,'>'
      write(povunit,*)radcyl,'}'
      write(povunit,*)'texture{pigment{color White}}'
      write(povunit,*)'scale< 1.0,1.0,1.0>translate<0.0,0.0,0.0>}'
!! 0-b
      write(povunit,*)' '
      write(povunit,*)'// simulation box edge '
      write(povunit,*)'object {cylinder{'
      write(povunit,*)'<',-cmsx,',',-cmsy,',',-cmsz,'>'
      write(povunit,*)'<',b1-cmsx,',',b2-cmsy,',',b3-cmsz,'>'
      write(povunit,*)radcyl,'}'
      write(povunit,*)'texture{pigment{color White}}'
      write(povunit,*)'scale< 1.0,1.0,1.0>translate<0.0,0.0,0.0>}'
!! 0-c
      write(povunit,*)' '
      write(povunit,*)'// simulation box edge '
      write(povunit,*)'object {cylinder{'
      write(povunit,*)'<',-cmsx,',',-cmsy,',',-cmsz,'>'
      write(povunit,*)'<',c1-cmsx,',',c2-cmsy,',',c3-cmsz,'>'
      write(povunit,*)radcyl,'}'
      write(povunit,*)'texture{pigment{color White}}'
      write(povunit,*)'scale< 1.0,1.0,1.0>translate<0.0,0.0,0.0>}'
!! a - a+b
      write(povunit,*)' '
      write(povunit,*)'// simulation box edge '
      write(povunit,*)'object {cylinder{'
      write(povunit,*)'<',a1-cmsx,',',a2-cmsy,',',a3-cmsz,'>'
      write(povunit,*)'<',a1+b1-cmsx,',',a2+b2-cmsy,',',a3+b3-cmsz,'>'
      write(povunit,*)radcyl,'}'
      write(povunit,*)'texture{pigment{color White}}'
      write(povunit,*)'scale< 1.0,1.0,1.0>translate<0.0,0.0,0.0>}'
!! a - a+c
      write(povunit,*)' '
      write(povunit,*)'// simulation box edge '
      write(povunit,*)'object {cylinder{'
      write(povunit,*)'<',a1-cmsx,',',a2-cmsy,',',a3-cmsz,'>'
      write(povunit,*)'<',a1+c1-cmsx,',',a2+c2-cmsy,',',a3+c3-cmsz,'>'
      write(povunit,*)radcyl,'}'
      write(povunit,*)'texture{pigment{color White}}'
      write(povunit,*)'scale< 1.0,1.0,1.0>translate<0.0,0.0,0.0>}'
!! b - a+b
      write(povunit,*)' '
      write(povunit,*)'// simulation box edge '
      write(povunit,*)'object {cylinder{'
      write(povunit,*)'<',b1-cmsx,',',b2-cmsy,',',b3-cmsz,'>'
      write(povunit,*)'<',a1+b1-cmsx,',',a2+b2-cmsy,',',a3+b3-cmsz,'>'
      write(povunit,*)radcyl,'}'
      write(povunit,*)'texture{pigment{color White}}'
      write(povunit,*)'scale< 1.0,1.0,1.0>translate<0.0,0.0,0.0>}'
!! c - a+c
      write(povunit,*)' '
      write(povunit,*)'// simulation box edge '
      write(povunit,*)'object {cylinder{'
      write(povunit,*)'<',c1-cmsx,',',c2-cmsy,',',c3-cmsz,'>'
      write(povunit,*)'<',a1+c1-cmsx,',',a2+c2-cmsy,',',a3+c3-cmsz,'>'
      write(povunit,*)radcyl,'}'
      write(povunit,*)'texture{pigment{color White}}'
      write(povunit,*)'scale< 1.0,1.0,1.0>translate<0.0,0.0,0.0>}'
!! c - b+c
      write(povunit,*)' '
      write(povunit,*)'// simulation box edge '
      write(povunit,*)'object {cylinder{'
      write(povunit,*)'<',c1-cmsx,',',c2-cmsy,',',c3-cmsz,'>'
      write(povunit,*)'<',b1+c1-cmsx,',',b2+c2-cmsy,',',b3+c3-cmsz,'>'
      write(povunit,*)radcyl,'}'
      write(povunit,*)'texture{pigment{color White}}'
      write(povunit,*)'scale< 1.0,1.0,1.0>translate<0.0,0.0,0.0>}'
!! b - c+b
      write(povunit,*)' '
      write(povunit,*)'// simulation box edge '
      write(povunit,*)'object {cylinder{'
      write(povunit,*)'<',b1-cmsx,',',b2-cmsy,',',b3-cmsz,'>'
      write(povunit,*)'<',c1+b1-cmsx,',',c2+b2-cmsy,',',c3+b3-cmsz,'>'
      write(povunit,*)radcyl,'}'
      write(povunit,*)'texture{pigment{color White}}'
      write(povunit,*)'scale< 1.0,1.0,1.0>translate<0.0,0.0,0.0>}'
!! a+b - a+b+c
      write(povunit,*)' '
      write(povunit,*)'// simulation box edge '
      write(povunit,*)'object {cylinder{'
      write(povunit,*)'<',a1+b1-cmsx,',',a2+b2-cmsy,',',a3+b3-cmsz,'>'
      write(povunit,*)'<',a1+c1+b1-cmsx,',',a2+c2+b2-cmsy,',&
      ',a3+c3+b3-cmsz,'>'
      write(povunit,*)radcyl,'}'
      write(povunit,*)'texture{pigment{color White}}'
      write(povunit,*)'scale< 1.0,1.0,1.0>translate<0.0,0.0,0.0>}'
!! a+c - a+b+c
      write(povunit,*)' '
      write(povunit,*)'// simulation box edge '
      write(povunit,*)'object {cylinder{'
      write(povunit,*)'<',a1+c1-cmsx,',',a2+c2-cmsy,',',a3+c3-cmsz,'>'
      write(povunit,*)'<',a1+c1+b1-cmsx,',',a2+c2+b2-cmsy,',&
      ',a3+c3+b3-cmsz,'>'
      write(povunit,*)radcyl,'}'
      write(povunit,*)'texture{pigment{color White}}'
      write(povunit,*)'scale< 1.0,1.0,1.0>translate<0.0,0.0,0.0>}'
!! b+c - a+b+c
      write(povunit,*)' '
      write(povunit,*)'// simulation box edge '
      write(povunit,*)'object {cylinder{'
      write(povunit,*)'<',b1+c1-cmsx,',',b2+c2-cmsy,',',b3+c3-cmsz,'>'
      write(povunit,*)'<',a1+c1+b1-cmsx,',',a2+c2+b2-cmsy,',&
      ',a3+c3+b3-cmsz,'>'
      write(povunit,*)radcyl,'}'
      write(povunit,*)'texture{pigment{color White}}'
      write(povunit,*)'scale< 1.0,1.0,1.0>translate<0.0,0.0,0.0>}'
      write(povunit,*)' '
      endif ! lperiodic

!! write atoms in povray file
      do i1=1,num_atoms
        write(povunit,*)'// Atom '
        write(povunit,'(a8,x,f8.3,x,a1,x,f8.3,x,a1,x,f8.3,x,a2)')&
        'sphere{<',xyztemp(1,i1),',',&
                  xyztemp(2,i1),',',xyztemp(3,i1),'>,'
        write(povunit,'(f8.3,x,a9,x,a8,x,a2)')&
                     radatom(zelem(i1)),'texture {',&
                     atomcolor(zelem(i1)),'}}'
      enddo
!!
            write(povunit,*)' '
!!
!! calculated interatomic distances
      rminimal=1000.0d0
      do i1=1,num_atoms
        xtemp=xyztemp(1,i1)
        ytemp=xyztemp(2,i1)
        ztemp=xyztemp(3,i1)
        do i2=1,num_atoms
        if (i1.ne.i2) then
          rminimal=1000.0d0
!! original box
          r=(xtemp-xyztemp(1,i2))**2 &
          + (ytemp-xyztemp(2,i2))**2 &
          + (ztemp-xyztemp(3,i2))**2
          r=dsqrt(r)
          rminimal=min(r,rminimal)
          if(r.le.rcut) then
            write(povunit,*)'// bond '
            write(povunit,*)'cylinder{'
            write(povunit,'(a1,x,f8.3,x,a1,x,f8.3,x,a1,x,f8.3,x,a2)')&
                          '<',xtemp,',',ytemp,',',ztemp,'>,'
            write(povunit,'(a1,x,f8.3,x,a1,x,f8.3,x,a1,x,f8.3,x,a2)')&
                          '<',(xtemp+xyztemp(1,i2))/2.0d0,',',&
                              (ytemp+xyztemp(2,i2))/2.0d0,',',&
                              (ztemp+xyztemp(3,i2))/2.0d0,'>,'
            write(povunit,'(f8.3,x,a10,x,a8,x,a2)')&
              radbond,' texture {',atomcolor(zelem(i1)),'}}'

            write(povunit,*)'cylinder{'
            write(povunit,'(a1,x,f8.3,x,a1,x,f8.3,x,a1,x,f8.3,x,a2)')&
                          '<',(xtemp+xyztemp(1,i2))/2.0d0,&
                          ',',(ytemp+xyztemp(2,i2))/2.0d0,&
                          ',',(ztemp+xyztemp(3,i2))/2.0d0,'>,'
            write(povunit,'(a1,x,f8.3,x,a1,x,f8.3,x,a1,x,f8.3,x,a2)')&
                          '<',xyztemp(1,i2),',',&
                              xyztemp(2,i2),',',&
                              xyztemp(3,i2),'>,'
            write(povunit,'(f8.3,x,a10,x,a8,x,a2)')&
              radbond,' texture {',atomcolor(zelem(i2)),'}}'
          endif
        endif ! i1.ne.i2
        enddo ! i2
      enddo ! i1
!!
!! now the bonds to atoms in neighbored boxes should be drawn
!!
!!      write(povunit,*)'}}' ! for union
      close(povunit)
!!
      return
      end
