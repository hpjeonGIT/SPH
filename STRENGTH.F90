SUBROUTINE ElastPlast(NP, N, MatDat, xx, xv, xm, rho, W, dW, NW, Stress, dt)
!
!==============================================================================
IMPLICIT NONE
!== External variables ========================================================
! 
INTEGER :: NP, N, NW(NP,NP)
REAL*8 :: MatDat(20,1), xx(NP,2), xv(NP,2), xm(NP), rho(NP), W(NP,NP), &
     dW(NP,NP,2), Stress(NP,4), dt
!== Internal variables ========================================================
! 
INTEGER :: i, j, k
REAL*8 :: E, nu, Y, rate(2,2), rateth, R(2), erate(4,1), C(4,4), Jaumann(4,1),&
     press, S(4), dSdt(4), f, J2, temp, onepi
!
! Initialization
E = MatDat(2,1) ! elastic modulus
nu = MatDat(3,1)  ! poisson ratio
Y = MatDat(4,1) ! Yield limit
onepi = 3.141592654
rate = 0.0
rateth = 0.0
C = 0.0
!
! Summation of velocity/kernel terms
do j=1, NW(N,1)
   k = NW(N,j+1)
   rate(1,1) = rate(1,1) + &
        xm(k)*(xv(k,1) - xv(N,1))*dW(N,k,1)/rho(k)/onepi/abs(xx(k,1))
   rate(2,2) = rate(2,2) + &
        xm(k)*(xv(k,2) - xv(N,2))*dW(N,k,2)/rho(k)/onepi/abs(xx(k,1))
   rate(1,2) = rate(1,2) + &
        xm(k)*(xv(k,1) - xv(N,1))*dW(N,k,2)/rho(k)/onepi/abs(xx(k,1))
   rate(2,1) = rate(2,1) + &
        xm(k)*(xv(k,2) - xv(N,2))*dW(N,k,1)/rho(k)/onepi/abs(xx(k,1))
   rateth = rateth + xm(k)*xv(k,1)*W(N,K)/rho(k)/onepi/abs(xx(k,1))/xx(k,1)
end do
!
! Rotation tensor
R(1) = 0.5*(rate(1,2) - rate(2,1))
R(2) = -R(1)
!
! strain rate (1-r 2-theta 3-z 4-shear)
erate(1,1) = rate(1,1)
erate(3,1) = rate(2,2)
erate(4,1) = rate(1,2) + rate(2,1)
erate(2,1) = rateth
!
! Elastic compliance
temp = E/(1.+nu)/(1.-2.*nu)
C(1,1) = temp*(1.-nu)
C(1,2) = temp*nu
C(2,1) = C(1,2)
C(2,2) = C(1,1)
C(2,3) = C(1,2)
C(3,2) = C(1,2)
C(3,3) = C(1,1)
C(4,4) = temp*(1.-2.*nu)/2.
!
! Jaumann rate
Jaumann = matmul(C, erate)
!
! Stress rate (1-r 2-z 3-theta 4-shear)
dSdt(1) = Jaumann(1,1) + 2.*R(1)*Stress(N,4)
dSdt(2) = Jaumann(3,1) + 2.*R(2)*Stress(N,4)
dSdt(4) = Jaumann(4,1) + R(1)*Stress(N,2) + R(2)*Stress(N,1)
dSdt(3) = Jaumann(2,1)
!
! Stress update
do i=1,4
   Stress(N,i) = Stress(N,i) + dSdt(i)*dt
end do
!
! press
press = -(Stress(N,1) + Stress(N,2) + Stress(N,3))/3.0
!
! dev. stress
do i=1,3
   S(i) = Stress(N,i) + press
end do
S(4) = Stress(N,4)
!
! J2 invariant
J2 = S(1)**2 + S(2)**2 + S(3)**2 + 2.*S(4)**2
!
! Yield check and modify - perfect plastic
f = 1. 
f = min(sqrt(2./3./J2)*Y,f)
!
! dev. stress update
S = S*f
!
! Stress update
do i=1,3
   Stress(N,i) = S(i) - press
end do
Stress(N,4) = S(4)
!
RETURN
END SUBROUTINE ElastPlast
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Elast(NP, N, MatDat, xx, xv, xm, rho, W, dW, NW, Stress, dt)
!
!==============================================================================
IMPLICIT NONE
!== External variables ========================================================
! 
INTEGER :: NP, N, NW(NP,NP)
REAL*8 :: MatDat(20,1), xx(NP,2), xv(NP,2), xm(NP), rho(NP), W(NP,NP), dW(NP,NP,2), &
     Stress(NP,4), dt
!== Internal variables ========================================================
! 
INTEGER :: i, j, k
REAL*8 :: E, nu, Y, rate(2,2), rateth, R(2), erate(4,1), C(4,4), Jaumann(4,1),&
     press, S(4), dSdt(4), f, J2, temp, onepi
!
! Initialization
E = MatDat(2,1) ! elastic modulus
nu = MatDat(3,1)  ! poisson ratio
Y = MatDat(4,1) ! Yield limit
onepi = 3.141592654*2.0
rate = 0.0
rateth = 0.0
C = 0.0
!
! Summation of velocity/kernel terms
do j=1, NW(N,1)
   k = NW(N,j+1)
   rate(1,1) = rate(1,1) + &
        xm(k)*(xv(k,1) - xv(N,1))*dW(N,k,1)/rho(k)/onepi/xx(k,1)
   rate(2,2) = rate(2,2) + &
        xm(k)*(xv(k,2) - xv(N,2))*dW(N,k,2)/rho(k)/onepi/xx(k,1)
   rate(1,2) = rate(1,2) + &
        xm(k)*(xv(k,1) - xv(N,1))*dW(N,k,2)/rho(k)/onepi/xx(k,1)
   rate(2,1) = rate(2,1) + &
        xm(k)*(xv(k,2) - xv(N,2))*dW(N,k,1)/rho(k)/onepi/xx(k,1)
   rateth = rateth + xm(k)*xv(k,1)*W(N,K)/rho(k)/onepi/xx(k,1)**2
end do
!
! Rotation tensor
R(1) = 0.5*(rate(1,2) - rate(2,1))
R(2) = -R(1)
!
! strain rate (1-r 2-theta 3-z 4-shear)
erate(1,1) = rate(1,1)
erate(3,1) = rate(2,2)
erate(4,1) = rate(1,2) + rate(2,1)
erate(2,1) = rateth
!
! Elastic compliance
temp = E/(1.+nu)/(1.-2.*nu)
C(1,1) = temp*(1.-nu)
C(1,2) = temp*nu
C(2,1) = C(1,2)
C(2,2) = C(1,1)
C(2,3) = C(1,2)
C(3,2) = C(1,2)
C(3,3) = C(1,1)
C(4,4) = temp*(1.-2*nu)/2.
!
! Jaumann rate
Jaumann = matmul(C, erate)
!
! Stress rate (1-r 2-z 3-theta 4-shear)
dSdt(1) = Jaumann(1,1) + 2.*R(1)*Stress(N,4)
dSdt(2) = Jaumann(3,1) + 2.*R(2)*Stress(N,4)
dSdt(4) = Jaumann(4,1) + R(1)*Stress(N,2) + R(2)*Stress(N,1)
dSdt(3) = Jaumann(2,1)
!
! Stress update
do i=1,4
   Stress(N,i) = Stress(N,i) + dSdt(i)*dt
end do
!
RETURN
END SUBROUTINE Elast






