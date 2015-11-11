SUBROUTINE EOS(NP, N, MatDat, xx, xv, xm, rho, rho0, W, dW, NW, Stress, Td, E)
!
!==============================================================================
IMPLICIT NONE
!== External variables ========================================================
! 
INTEGER :: NP, N, NW(NP,NP)
REAL*8 :: MatDat(20,1), xx(NP,2), xv(NP,2), xm(NP), rho(NP), rho0(NP), &
     W(NP,NP), dW(NP,NP,2), Stress(NP,4), Td, E(NP)
!== Internal variables ========================================================
! 
INTEGER :: i, j, k
REAL*8 :: press, twopi, A, B, R1, R2, omega, Em0, Cd, pH, eta, a0, b0, c0, &
     Gamma, C, S
!
! Initialization
twopi = 3.141592654*2.0
C = MatDat(1,1)
Gamma = MatDat(10,1)
S = MatDat(11,1)
A = MatDat(12,1)
B = MatDat(13,1)
R1 = MatDat(14,1)
R2 = MatDat(15,1)
omega = MatDat(16,1)
Em0 = MatDat(17,1)
Cd = MatDat(18,1) 
!
!
eta = rho(N)/rho0(N) - 1.
a0 = rho0(N)*C**2
b0 = a0*(1. + 2.*(S - 1.))
c0 = a0*(2.*(S - 1.) + 3.*(S - 1.)**2 )
!
!
if (eta > 0.) then
   pH = a0*eta + b0*eta**2 + c0*eta**3
else
   pH = a0*eta
end if
pH = pH*1.e9
!
! press
if (sqrt(xx(N,1)**2 + xx(N,2)**2).lt. (Td*Cd)) then
   press = A*(1. - omega*rho(N)/R1/rho0(N))*exp(-R1*rho0(N)/rho(N)) + &
        B*(1. - omega*rho(N)/R2/rho0(N))*exp(-R2*rho0(N)/rho(N)) + &
        omega*rho(N)/rho0(N)
else
   press = 0.0 !(1. - .5*Gamma*eta)*pH + Gamma*rho(N)*E(N)
end if
!
! Stress update
do i=1,3
   Stress(N,i) = -press*1.e9
end do
!
RETURN
END SUBROUTINE EOS

