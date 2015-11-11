SUBROUTINE KernelFtn(NP, N, xx, W, dW, NW, h)
!
!==============================================================================
IMPLICIT NONE
!== External variables ========================================================
! 
INTEGER :: NP, N, NW(NP,NP)
REAL*8 :: xx(NP,2), W(NP,NP), dW(NP,NP,2), h
! NW - Flag and corresponding particle data
! NW( ,1) = number of particles for a particle
! NW( ,2) = 1st particle for a particle
!== Internal variables ========================================================
! 
INTEGER :: i, j, k, D
REAL*8 :: pi, z(NP), C, grad(2)
!
! preliminary setting
pi = 3.141592654
k = 0
do i=1, NP
   z(i) = sqrt( (xx(i,1) - xx(N,1))**2 + (xx(i,2) - xx(N,2))**2 )/h
end do
z(N) = 1.e6
!
! 2-dimensional case
D = 2
C = 10./7./pi
!
! loop start
do i=1, NP
   if (z(i) < 1.0) then
      k = k + 1
      NW(N,k+1) = i
      grad(1) = (xx(N,1) - xx(i,1))/z(i)/h**2
      grad(2) = (xx(N,2) - xx(i,2))/z(i)/h**2
      W(N,i) = C*(1. - 1.5*z(i)**2 + 0.75*z(i)**3)/h**D
      dW(N,i,1) = C*(-3.*z(i) + 2.25*z(i)**2)*grad(1)/h**D
      dW(N,i,2) = C*(-3.*z(i) + 2.25*z(i)**2)*grad(2)/h**D
   else if (z(i) < 2.0) then
      k = k + 1
      NW(N,k+1) = i
      grad(1) = (xx(N,1) - xx(i,1))/z(i)/h**2
      grad(2) = (xx(N,2) - xx(i,2))/z(i)/h**2
      W(N,i) = 0.25*C*(2. - z(i))**3/h**D
      dW(N,i,1) = -0.75*C*(2. - z(i))**2*grad(1)/h**D
      dW(N,i,2) = -0.75*C*(2. - z(i))**2*grad(2)/h**D
   else
      W(N,i) = 0.0
      dW(N,i,1) = 0.0
      dW(N,i,2) = 0.0
   end if
end do
NW(N,1) = k
RETURN
END SUBROUTINE KernelFtn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subroutine for estimating smoothing length
SUBROUTINE LengthH(NP, MatDat, xx, h)
!
!==============================================================================
IMPLICIT NONE
!== External variables ========================================================
! 
INTEGER :: NP
REAL*8 :: MatDat(20,1), xx(NP,2), h(NP)
!== Internal variables ========================================================
! 
INTEGER :: i, j
REAL*8 :: a, z(NP)
! 
! h = min |x - x'| * a
a = MatDat(9,1)
do i=1, NP
   do j=1, NP
      z(j) = sqrt( (xx(i,1) - xx(j,1))**2 + (xx(i,2) - xx(j,2))**2)
   end do
   z(i) = 1.e6
h(i) = minval(z) * a
end do
RETURN
END SUBROUTINE LengthH



