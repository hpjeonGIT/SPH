SUBROUTINE Solver(NP, MatDat, xm, xx, xv, xa, TimeToStop, TimeInter)
!
!==============================================================================
IMPLICIT NONE
!== External variables ========================================================
! 
INTEGER :: NP
REAL*8 :: MatDat(20,1), xm(NP), xx(NP,2), xv(NP,2), xa(NP,2), &
     TimeToStop, TimeInter
!== Internal variables ========================================================
! 
INTEGER :: i, j, k, Interval, NW(NP,NP), Nstep
REAL*8 :: W(NP,NP), dW(NP,NP,2), h(NP), rho(NP), drhodt(NP), rho0(NP), &
     pi(NP,NP), D(NP), Stress(NP,4), dSdt(4), E(NP), dE(NP)
!
REAL*8 :: dt, TimeElapsed, c, mu, ee, kk, Elastic, nu, J2, Y, f, temp(3), &
     onepi
!== Forehead of loop ==========================================================
! 
c = MatDat(1,1) ! speed of sound
Elastic = MatDat(2,1) ! elastic modulus
nu = MatDat(3,1)  ! poisson ratio
Y = MatDat(4,1) ! Yield limit
mu = Elastic/2./(1. + nu) ! shear elastic modulus
kk = Elastic/3./(1. - 2.*nu) ! bulk modulus
Stress = 0.0   ! stress initialization
E = 0.0   ! Internal energy initialization
TimeElapsed = 0.0 ! total elapsed time
Interval = 0 ! parameter for output processing
Nstep = 0 ! Step number
onepi = 3.1415692
dt = 0.0
rho0 = 1632
rho = 1632
!
! Initial state postprocess
call Postprocess(NP, xv, xx, Stress, Interval, NW)
open(23,file="time.dat")
write(*,2000)
!
!== Estimate smoothed length of each particle =================================
!
10  drhodt = 0.0
!== Predict position and velocity =============================================
xx = xx + xv*dt + 0.5*xa*dt**2
xv = xv + 0.5*dt*xa
Call LengthH(NP, MatDat, xx, h)
!
!==============================================================================
! 1st loop
do i=1, NP
   !== Kernel function estimation =============================================
   ! 
   Call KernelFtn(NP, i, xx, W, dW, NW, h(i))
   !== Density estimation =====================================================
   ! 
   do j=1, NW(i,1)
      k = NW(i,1+j)
      drhodt(i) = drhodt(i) + xm(k)*abs( (xv(i,1) - xv(k,1))*dW(i,k,1) + &
           (xv(i,2) - xv(k,2))*dW(i,k,2) )/rho(k)/onepi/abs(xx(k,1))
   end do
   drhodt(i) = drhodt(i)*rho(i)
!   do j=1, NW(i,1)
!      k = NW(i,1+j)
!      rho(i) = rho(i) + xm(k)*W(i,k)/xx(k,1)/onepi
!   end do
end do
!
! Estimate time increment
Call Courant(NP, MatDat, xv, c, h, dt)
if (dt < 1.0E-10 )  stop "Time increment is too small"
TimeElapsed = TimeElapsed + dt
rho = rho + drhodt*dt
Nstep = Nstep + 1
do i=1,NP
   !== Artificial viscosity estimation ========================================
   ! 
   Call ArtVisco(NP, i, NW, MatDat, xx, xv, h, rho, pi)
   !== Constitutive relation ==================================================
   ! Elastic-perfectly plastic constitutive law
   Call EOS(NP, i, MatDat, xx, xv, xm, rho, rho0, W, dW, NW, Stress, TimeElapsed, E)
end do
!== Acceleration and internal energy estimation ===============================
! 
xa = 0.0  ! acceleration initialization
dE = 0.0
do i=1, NP
   do j=1, NW(i,1)
      k = NW(i,j+1)
      xa(i,1) = xa(i,1) + ( &
           (Stress(i,1)/rho(i)**2 + Stress(k,1)/rho(k)**2-pi(i,k))*dW(i,k,1) +&
           (Stress(i,4)/rho(i)**2 + Stress(k,4)/rho(k)**2)*dW(i,k,2) + &
           (Stress(k,1) - Stress(k,3))*W(i,k)/xx(k,1)/rho(k)**2 &
           )*xm(k)/onepi/abs(xx(k,1))
      xa(i,2) = xa(i,2) + ( &
           (Stress(i,2)/rho(i)**2 + Stress(k,2)/rho(k)**2-pi(i,k))*dW(i,k,2) +&
           (Stress(i,4)/rho(i)**2 + Stress(k,4)/rho(k)**2)*dW(i,k,1) + &
           Stress(k,4)*W(i,k)/xx(k,1)/rho(k)**2 &
           )*xm(k)/onepi/abs(xx(k,1))
      dE(i) = dE(i) + xm(k) * ( &
           (xv(i,1)-xv(k,1))*(Stress(i,1)/rho(i)**2+0.5*pi(i,k))*dW(i,k,1) + &
           (xv(i,1)-xv(k,1))*(Stress(i,4)/rho(i)**2+0.5*pi(i,k))*dW(i,k,2) + &
           (xv(i,2)-xv(k,2))*(Stress(i,4)/rho(i)**2+0.5*pi(i,k))*dW(i,k,1) + &
           (xv(i,2)-xv(k,2))*(Stress(i,2)/rho(i)**2+0.5*pi(i,k))*dW(i,k,2) )
   end do
end do
!== kinematic condition =======================================================
!
Call Kinematic(NP, xa, xv, xx, xm)
!== Correct - particle velocity ===============================================
! 
xv = xv + 0.5 * xa * dt
E = E + dE * dt
!== Print output ==============================================================
!
write(23,*) TimeElapsed, Stress(91,1), Stress(91,2), Stress(91,3), Stress(91,4)
if(TimeElapsed > TimeInter*Interval) then
   call Postprocess(NP, xv, xx, Stress, Interval, NW)
   write(*, 2001) Nstep, TimeElapsed, dt
2000 FORMAT("No. Step     Interim output (sec)      Increment of time (sec)")
2001 FORMAT(I6,14X, F8.6, 21X, F8.6)
end if
!
! check time to stop
if (TimeElapsed < TimeToStop) then
   go to 10
end if
close(23)
!
Return
END SUBROUTINE Solver








