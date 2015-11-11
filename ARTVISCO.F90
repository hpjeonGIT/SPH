SUBROUTINE ArtVisco(NP, N, NW, MatDat, xx, xv, h, rho, pi)
!
!==============================================================================
IMPLICIT NONE
!== External variables ========================================================
! 
INTEGER :: NP, N, NW(NP,NP)
REAL*8 :: MatDat(20,1), xx(NP,2), xv(NP,2), h(NP), rho(NP), pi(NP,NP)
!== Internal variables ========================================================
! 
INTEGER :: i, j
REAL*8 :: epsilon, cbar, alpha, beta, vx, hbar, x2, mu, rhobar
!==============================================================================
! common value is used
epsilon = MatDat(5,1)
!
! material property
cbar = MatDat(1,1)
alpha = MatDat(6,1)
beta = MatDat(7,1)
!
! 
do i=1,NW(N,1)
   j = NW(N,i+1)
   vx = (xv(N,1) - xv(j,1))*(xx(N,1) - xx(j,1)) + &
        (xv(N,2) - xv(j,2))*(xx(N,2) - xx(j,2))
   if (vx < 0.0) then
      ! 
      ! two particles are adjacent 
      hbar = (h(N) + h(j))/2.0
      x2 = (xx(N,1) - xx(j,1))**2 + (xx(N,2) - xx(j,2))**2
      mu = hbar*vx/(x2 + epsilon*hbar**2)
      rhobar = (rho(N) + rho(j))/2.0
      !
      ! edited - originally divided by rhobar
      pi(N,j) = (beta*mu**2 - alpha*cbar*mu)/rhobar
   else
      !
      ! two particles are away
      pi(N,j) = 0.0
   end if
end do
RETURN
END SUBROUTINE ArtVisco
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subroutine for estimating time increment
SUBROUTINE Courant(NP, MatDat, xv, c, h, dt)
!
!==============================================================================
IMPLICIT NONE
!== External variables ========================================================
! 
INTEGER :: NP
REAL*8 :: MatDat(20,1), xv(NP,2), h(NP), c, dt
!== Internal variables ========================================================
! 
INTEGER :: i
REAL*8 :: w, s(NP), criteria(NP)
!
! Initialization
w = MatDat(8,1)
! 
do i=1,NP
   s(i) = sqrt(xv(i,1)**2 + xv(i,2)**2)
end do
!
criteria = h*w/(c+s)
dt = 0.3*minval(criteria)
RETURN
END SUBROUTINE Courant











