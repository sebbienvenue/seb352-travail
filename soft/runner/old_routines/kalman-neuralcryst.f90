subroutine kalman-neuralcryst( aw,nn,maxn,lay,lambda,nue,corr,&
     gaink,h,t,exam,nodes,error,x,mode)

  use definitions

  ! performs weights changes based on the kalman filter
  ! algorithm. this routine is only valid for a single
  ! output node.

  ! Arguments
  character*7, intent( in ) :: mode

  integer, intent ( in ) :: lay,nn(0:lay),maxn
  integer, intent ( in ) :: exam,nodes

  real ( quadro ), intent ( in )    :: h(nodes,nn(lay))
  real ( quadro ), intent ( in )    :: t(nn(lay),exam)
  real ( quadro ), intent ( in )    :: nue
  real ( quadro ), intent ( inout ) :: lambda
  real ( quadro ), intent ( inout ) :: aw(lay,0:maxn,maxn) 
  real ( quadro ), intent ( inout ) :: corr(nodes*(nodes+1)/2)
  real ( quadro ), intent ( inout ) :: gaink(nodes)
  real ( quadro ), intent ( in )    :: x(nn(0),exam)
  real ( quadro ), intent ( in )    :: error 

  ! Local
  integer :: i,j,k,l,m,n,nrun,i1,i2

  real ( quadro ) :: hp,inv,diff,mul
  real ( quadro ) :: coh(nodes)

  ! --------------------------------------------------------------------

!  write(*,*)'ekfvec starts'
  ! scalar
  diff = error  

  ! forgetting schedule
  hp    = 1.d0 / lambda

  ! Calculation of lambda^-1*P(n-1)*H(n)
!  write(*,*)'ekfvec calls dspmv'
  call dspmv('l',nodes,1.d0,corr,h(:,1),1,0.d0,coh,1)

  ! Calculation of K(n)=lambda^(-1)*P(n-1)*H(n)*
  ! [ I + lambda^(-1)*H^T(n)*P(n-1)*H(n) ]^(-1)

  ! 1. Matrix Inversion including weight matrix
  ! ddot=Calculation of H^T(n) * P(n-1)*H(n)
  inv = 1.d0 / ( lambda + ddot(nodes,h(:,1),1,coh,1) )

  gaink(:) = inv * coh(:) 

  ! Calculation of P(n) = lambda^(-1)*[P(n-1)-K(n)*H^T(n)*P(n-1)]
  ! with (P*H)^T=H^T*P^T=H^T*P, P symmetric
!  write(*,*)'ekfvec calls dsprsl'
  call dsprsl('l',nodes,-hp*inv,coh,1,corr,hp)

  ! Weight update : Weights(n)=Weights(n-1) + K(n)*[d(n)-h(i(n))]
!  write(*,*)'ekfvec calls update'
  call update(lay,nn,maxn,nodes,aw,diff,gaink,mode)

  ! -- damping factor for kalman filter
  lambda = nue * lambda + 1.d0 - nue

end subroutine ekfvec
! -*- f90 -*-

subroutine update(lay,nn,maxn,nodes,aw,diff,vect,mode)

  use definitions

  character*7, intent ( in ) :: mode

  integer, intent ( in ) :: lay, nn(0:lay), maxn, nodes

  real ( quadro ), intent ( inout ) :: aw(lay,0:maxn,maxn)
  real ( quadro ), intent ( in )    :: diff
  real ( quadro ), intent ( in )    :: vect(nodes)

  ! local variables 
  integer nrun, n, l, k, m, i1, i2, iout

  ! --------------------------------------------------------------------

  nrun = 0  

  do m=1,lay
     do i1=1,nn(m)
        ! aw(m,0,i1)
        nrun = nrun+1
        aw(m,0,i1) = aw(m,0,i1) + diff * vect(nrun) 
        do i2=1,nn(m-1)
           ! aw(m,i2,i1)
           nrun = nrun+1
           aw(m,i2,i1) = aw(m,i2,i1) + diff * vect(nrun) 
        enddo
     enddo
  enddo

end subroutine update
