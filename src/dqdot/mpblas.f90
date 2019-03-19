!** MPBLAS
SUBROUTINE MPBLAS(I1)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to DQDOTA and DQDOTI
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (MPBLAS-A)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     This subroutine is called to set up Brent's 'mp' package
  !     for use by the extended precision inner products from the BLAS.
  !
  !     In the SLATEC library we require the Extended Precision MP number
  !     to have a mantissa twice as long as Double Precision numbers.
  !     The calculation of MPT (and MPMXR which is the actual array size)
  !     in this routine will give 2x (or slightly more) on the machine
  !     that we are running on.  The INTEGER array size of 30 was chosen
  !     to be slightly longer than the longest INTEGER array needed on
  !     any machine that we are currently aware of.
  !
  !***
  ! **See also:**  DQDOTA, DQDOTI
  !***
  ! **References:**  R. P. Brent, A Fortran multiple-precision arithmetic
  !                 package, ACM Transactions on Mathematical Software 4,
  !                 1 (March 1978), pp. 57-70.
  !               R. P. Brent, MP, a Fortran multiple-precision arithmetic
  !                 package, Algorithm 524, ACM Transactions on Mathema-
  !                 tical Software 4, 1 (March 1978), pp. 71-81.
  !***
  ! **Routines called:**  I1MACH, XERMSG
  !***
  ! COMMON BLOCKS    MPCOM

  !* REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !   930124  Increased Array size in MPCON for SUN -r8, and calculate
  !               size for Quad Precision for 2x DP.  (RWC)
  
  INTEGER I1, I1MACH, MPB, mpbexp, MPLun, MPM, MPMxr, MPR, MPT
  COMMON /MPCOM / MPB, MPT, MPM, MPLun, MPMxr, MPR(30)
  !* FIRST EXECUTABLE STATEMENT  MPBLAS
  I1 = 1
  !
  !     For full extended precision accuracy, MPB should be as large as
  !     possible, subject to the restrictions in Brent's paper.
  !
  !     Statements below are for an integer wordlength of  48, 36, 32,
  !     24, 18, and 16.  Pick one, or generate a new one.
  !       48     MPB = 4194304
  !       36     MPB =   65536
  !       32     MPB =   16384
  !       24     MPB =    1024
  !       18     MPB =     128
  !       16     MPB =      64
  !
  mpbexp = I1MACH(8)/2 - 2
  MPB = 2**mpbexp
  !
  !     Set up remaining parameters
  !                  UNIT FOR ERROR MESSAGES
  MPLun = I1MACH(4)
  !                  NUMBER OF MP DIGITS
  MPT = (2*I1MACH(14)+mpbexp-1)/mpbexp
  !                  DIMENSION OF R
  MPMxr = MPT + 4
  !
  IF ( MPMxr>30 ) THEN
    CALL XERMSG('SLATEC','MPBLAS',&
      'Array space not sufficient for Quad Precision 2x Double Precision, Proceeding.',1,1)
    MPT = 26
    MPMxr = 30
  ENDIF
  !                  EXPONENT RANGE
  MPM = MIN(32767,I1MACH(9)/4-1)
END SUBROUTINE MPBLAS
