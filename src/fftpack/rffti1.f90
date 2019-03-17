!DECK RFFTI1
SUBROUTINE RFFTI1(N,Wa,Ifac)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  RFFTI1
  !***PURPOSE  Initialize a real and an integer work array for RFFTF1 and
  !            RFFTB1.
  !***LIBRARY   SLATEC (FFTPACK)
  !***CATEGORY  J1A1
  !***TYPE      SINGLE PRECISION (RFFTI1-S, CFFTI1-C)
  !***KEYWORDS  FFTPACK, FOURIER TRANSFORM
  !***AUTHOR  Swarztrauber, P. N., (NCAR)
  !***DESCRIPTION
  !
  !   Subroutine RFFTI1 initializes the work arrays WA and IFAC which are
  !   used in both RFFTF1 and RFFTB1.  The prime factorization of N and a
  !   tabulation of the trigonometric functions are computed and stored in
  !   IFAC and WA, respectively.
  !
  !   Input Argument
  !
  !   N       the length of the sequence to be transformed.
  !
  !   Output Arguments
  !
  !   WA      a real work array which must be dimensioned at least N.
  !
  !   IFAC    an integer work array which must be dimensioned at least 15.
  !
  !   The same work arrays can be used for both RFFTF1 and RFFTB1 as long
  !   as N remains unchanged.  Different WA and IFAC arrays are required
  !   for different values of N.  The contents of WA and IFAC must not be
  !   changed between calls of RFFTF1 or RFFTB1.
  !
  !***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
  !                 Computations (G. Rodrigue, ed.), Academic Press,
  !                 1982, pp. 51-83.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   830401  Modified to use SLATEC library source file format.
  !   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
  !           (a) changing dummy array size declarations (1) to (*),
  !           (b) changing references to intrinsic function FLOAT
  !               to REAL, and
  !           (c) changing definition of variable TPI by using
  !               FORTRAN intrinsic functions instead of DATA
  !               statements.
  !   881128  Modified by Dick Valent to meet prologue standards.
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900131  Routine changed from subsidiary to user-callable.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  RFFTI1
  REAL arg, argh, argld, fi, tpi, Wa
  INTEGER i, ib, ido, Ifac, ii, ip, ipm, is, j, k1, l1, l2, ld, &
    N, nf, nfm1, nl, nq, nr, ntry
  INTEGER ntryh
  DIMENSION Wa(*), Ifac(*), ntryh(4)
  SAVE ntryh
  DATA ntryh(1), ntryh(2), ntryh(3), ntryh(4)/4, 2, 3, 5/
  !***FIRST EXECUTABLE STATEMENT  RFFTI1
  nl = N
  nf = 0
  j = 0
  100  j = j + 1
  IF ( j<=4 ) THEN
    ntry = ntryh(j)
  ELSE
    ntry = ntry + 2
  ENDIF
  DO
    nq = nl/ntry
    nr = nl - ntry*nq
    IF ( nr/=0 ) GOTO 100
    nf = nf + 1
    Ifac(nf+2) = ntry
    nl = nq
    IF ( ntry==2 ) THEN
      IF ( nf/=1 ) THEN
        DO i = 2, nf
          ib = nf - i + 2
          Ifac(ib+2) = Ifac(ib+1)
        ENDDO
        Ifac(3) = 2
      ENDIF
    ENDIF
    IF ( nl==1 ) THEN
      Ifac(1) = N
      Ifac(2) = nf
      tpi = 8.*ATAN(1.)
      argh = tpi/N
      is = 0
      nfm1 = nf - 1
      l1 = 1
      IF ( nfm1==0 ) RETURN
      DO k1 = 1, nfm1
        ip = Ifac(k1+2)
        ld = 0
        l2 = l1*ip
        ido = N/l2
        ipm = ip - 1
        DO j = 1, ipm
          ld = ld + l1
          i = is
          argld = ld*argh
          fi = 0.
          DO ii = 3, ido, 2
            i = i + 2
            fi = fi + 1.
            arg = fi*argld
            Wa(i-1) = COS(arg)
            Wa(i) = SIN(arg)
          ENDDO
          is = is + ido
        ENDDO
        l1 = l2
      ENDDO
      EXIT
    ENDIF
  ENDDO
END SUBROUTINE RFFTI1