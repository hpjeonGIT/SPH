PROGRAM SPH2D
!
! Test program of SPH 2D 
! 2001.08.15. HPJeon.
!
!==============================================================================
!
! Debug report
!
!==============================================================================
!
  IMPLICIT NONE
!== Basic parameter setting ===================================================
! 
INTEGER, PARAMETER :: NMAX = 50, NP=800
!
! NMAX : Maximum number of possible contagious particles in DOI
!
!== Dynamic variables =========================================================
!
! 
! NP : Number of particles
REAL*8 :: xm(NP)
REAL*8 :: xx(NP,2), xv(NP,2), xa(NP,2)
!
! xm : mass of particle
! xx, xv, xa : position, velocity, acceleration of particle
! W : kernel function
!
!== Static variables ==========================================================
! 
!INTEGER :: i, j
REAL*8 :: TimeToStop, TimeInter, MatDat(20,1)
!
!== Intrinsic function variables ==============================================
! 
REAL::t0, t1, t2
!
! calculation time count
t0 = 0.0
t1 = secnds(t0)
!== Parsing input data ========================================================
! 
Call Parsor(NP, MatDat, xm, xx, xv, TimeToStop, TimeInter)
!== Set initial variables =====================================================
! 
!========================================================================
!
! Start main loop
Call Solver(NP, MatDat, xm, xx, xv, xa, TimeToStop, TimeInter)
!
! End of program
t2 = secnds(t1)
print *, "Approximate cpu time is", t2, "seconds"
!
STOP
END PROGRAM SPH2D



