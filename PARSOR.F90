SUBROUTINE Parsor(NP, MatDat, xm, xx, xv, TimeToStop, TimeInter)
!
!==============================================================================
IMPLICIT NONE
!==============================================================================
! External variables
INTEGER :: NP
REAL*8 :: xm(NP), xx(NP,2), xv(NP,2), TimeToStop, TimeInter, MatDat(20,1)
!==============================================================================
! Internal variables
INTEGER :: i, j
!==============================================================================
! 
! Initial coordinate of particle
do j=1, 20
   do i=1, 40
      xx((j-1)*40+i, 1) = real(i)*0.05 -1.025
      xx((j-1)*40+i, 2) = real(j)*0.05
      xm((j-1)*40+i) = 1632*3.1415*abs(xx((j-1)*40+i, 1))*0.05*0.05
   end do
end do
xv = 0.0
!
TimeToStop = 0.0004
TimeInter =  0.00005
!
! Material Data - TNT
MatDat(1,1) = 5.4e3  ! speed of sound
MatDat(2,1) = 200.e9 ! elastic modulus
MatDat(3,1) = 0.3    ! poisson ratio
MatDat(4,1) = 200.0E6! Yield lgimit
MatDat(5,1) = 0.01   ! epsilon - for artificial viscosity - 0.01
MatDat(6,1) = 0.5    ! alpha - for artificial viscosity - 0.5
MatDat(7,1) = 0.5    ! beta - for artificial viscosity - 0.5 
MatDat(8,1) = 0.3    ! w - for courant time criterion
MatDat(9,1) = 2.0    ! a - for smoothing length
MatDat(10,1) = 1.81  ! Gamma - hugoniot data
MatDat(11,1) = 1.8   ! S - hugoniot data
MatDat(12,1) = 524.4089 ! A - JWL
MatDat(13,1) = 4.900052 ! B - JWL
MatDat(14,1) = 4.579    ! R1 - JWL
MatDat(15,1) = 0.85     ! R2 - JWL
MatDat(16,1) = 0.23     ! w - JWL
MatDat(17,1) = 7.1      ! Em0 - JWL
MatDat(18,1) = 7.07E3   ! Cd - detonation velocity
RETURN
END SUBROUTINE Parsor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subroutine for estimating smoothing length
SUBROUTINE Kinematic(NP, xa, xv, xx, xm)
!
!==============================================================================
IMPLICIT NONE
!== External variables ========================================================
! 
INTEGER :: NP
REAL*8 :: xa(NP,2), xv(NP,2), xx(NP,2), xm(NP)
!== Internal variables ========================================================
! 
INTEGER :: i 
!
RETURN
END SUBROUTINE Kinematic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subroutine for estimating smoothing length
SUBROUTINE Postprocess(NP, xv, xx, S, Interval, NW)
!
!==============================================================================
IMPLICIT NONE
!== External variables ========================================================
! 
INTEGER :: NP, Interval, NW(NP,NP)
REAL*8 :: xv(NP,2), xx(NP,2), S(NP,4)
!== Internal variables ========================================================
! 
INTEGER :: i 
CHARACTER(LEN=15) :: FileName, File2
CHARACTER(LEN=10) :: NUMBER
!
!
FileName="output0000.dat"
File2   ="stress0000.dat"
DATA NUMBER/'0123456789'/
!
! File Name set
if(int(Interval/10) < 1) then
   FileName(10:10) = Number(Interval+1:Interval+1)
   File2(10:10) = Number(Interval+1:Interval+1)
else
   FileName(9:9) = Number(int(Interval/10)+1:int(Interval/10)+1)
   FileName(10:10) = Number(Interval - int(Interval/10)*10 +1 : &
        Interval - int(Interval/10)*10 + 1)
   File2(9:9) = Number(int(Interval/10)+1:int(Interval/10)+1)
   File2(10:10) = Number(Interval - int(Interval/10)*10 +1 : &
        Interval - int(Interval/10)*10 + 1)
end if
!
! File open
open(UNIT=16, FILE=FileName)
open(UNIT=17, FILE=File2)
!
! Write file
do i=1,NP
   write(16,2001) xx(i,1), xx(i,2), i
end do
do i=1,20
   write(17,2002) xx(i,1), NW(i,1)
end do
2001 FORMAT(E15.6, 2X, E15.6, 2X, I4)
2002 FORMAT(E15.6, 2X, I4)
!
! File close
close(16)
close(17)
Interval = Interval + 1
RETURN
END SUBROUTINE Postprocess

