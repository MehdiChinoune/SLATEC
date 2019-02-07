!*==RFFTF1.f90  processed by SPAG 6.72Dc at 10:55 on  6 Feb 2019
!DECK RFFTF1
SUBROUTINE RFFTF1(N,C,Ch,Wa,Ifac)
  IMPLICIT NONE
  !*--RFFTF15
  !*** Start of declarations inserted by SPAG
  REAL C, Ch, Wa
  INTEGER i, idl1, ido, Ifac, ip, iw, ix2, ix3, ix4, k1, kh, l1, &
    l2, N, na, nf
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  RFFTF1
  !***PURPOSE  Compute the forward transform of a real, periodic sequence.
  !***LIBRARY   SLATEC (FFTPACK)
  !***CATEGORY  J1A1
  !***TYPE      SINGLE PRECISION (RFFTF1-S, CFFTF1-C)
  !***KEYWORDS  FFTPACK, FOURIER TRANSFORM
  !***AUTHOR  Swarztrauber, P. N., (NCAR)
  !***DESCRIPTION
  !
  !   Subroutine RFFTF1 computes the Fourier coefficients of a real
  !   periodic sequence (Fourier analysis).  The transform is defined
  !   below at output parameter C.
  !
  !   The arrays WA and IFAC which are used by subroutine RFFTB1 must be
  !   initialized by calling subroutine RFFTI1.
  !
  !   Input Arguments
  !
  !   N       the length of the array R to be transformed.  The method
  !           is most efficient when N is a product of small primes.
  !           N may change so long as different work arrays are provided.
  !
  !   C       a real array of length N which contains the sequence
  !           to be transformed.
  !
  !   CH      a real work array of length at least N.
  !
  !   WA      a real work array which must be dimensioned at least N.
  !
  !   IFAC    an integer work array which must be dimensioned at least 15.
  !
  !           The WA and IFAC arrays must be initialized by calling
  !           subroutine RFFTI1, and different WA and IFAC arrays must be
  !           used for each different value of N.  This initialization
  !           does not have to be repeated so long as N remains unchanged.
  !           Thus subsequent transforms can be obtained faster than the
  !           first.  The same WA and IFAC arrays can be used by RFFTF1
  !           and RFFTB1.
  !
  !   Output Argument
  !
  !   C       C(1) = the sum from I=1 to I=N of R(I)
  !
  !           If N is even set L = N/2; if N is odd set L = (N+1)/2
  !
  !             then for K = 2,...,L
  !
  !                C(2*K-2) = the sum from I = 1 to I = N of
  !
  !                     C(I)*COS((K-1)*(I-1)*2*PI/N)
  !
  !                C(2*K-1) = the sum from I = 1 to I = N of
  !
  !                    -C(I)*SIN((K-1)*(I-1)*2*PI/N)
  !
  !           If N is even
  !
  !                C(N) = the sum from I = 1 to I = N of
  !
  !                     (-1)**(I-1)*C(I)
  !
  !   Notes:  This transform is unnormalized since a call of RFFTF1
  !           followed by a call of RFFTB1 will multiply the input
  !           sequence by N.
  !
  !           WA and IFAC contain initialization calculations which must
  !           not be destroyed between calls of subroutine RFFTF1 or
  !           RFFTB1.
  !
  !***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
  !                 Computations (G. Rodrigue, ed.), Academic Press,
  !                 1982, pp. 51-83.
  !***ROUTINES CALLED  RADF2, RADF3, RADF4, RADF5, RADFG
  !***REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   830401  Modified to use SLATEC library source file format.
  !   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
  !           changing dummy array size declarations (1) to (*).
  !   881128  Modified by Dick Valent to meet prologue standards.
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900131  Routine changed from subsidiary to user-callable.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  RFFTF1
  DIMENSION Ch(*), C(*), Wa(*), Ifac(*)
  !***FIRST EXECUTABLE STATEMENT  RFFTF1
  nf = Ifac(2)
  na = 1
  l2 = N
  iw = N
  DO k1 = 1, nf
    kh = nf - k1
    ip = Ifac(kh+3)
    l1 = l2/ip
    ido = N/l2
    idl1 = ido*l1
    iw = iw - (ip-1)*ido
    na = 1 - na
    IF ( ip==4 ) THEN
      ix2 = iw + ido
      ix3 = ix2 + ido
      IF ( na/=0 ) THEN
        CALL RADF4(ido,l1,Ch,C,Wa(iw),Wa(ix2),Wa(ix3))
      ELSE
        CALL RADF4(ido,l1,C,Ch,Wa(iw),Wa(ix2),Wa(ix3))
      ENDIF
    ELSEIF ( ip/=2 ) THEN
      IF ( ip==3 ) THEN
        ix2 = iw + ido
        IF ( na/=0 ) THEN
          CALL RADF3(ido,l1,Ch,C,Wa(iw),Wa(ix2))
        ELSE
          CALL RADF3(ido,l1,C,Ch,Wa(iw),Wa(ix2))
        ENDIF
      ELSEIF ( ip/=5 ) THEN
        IF ( ido==1 ) na = 1 - na
        IF ( na/=0 ) THEN
          CALL RADFG(ido,ip,l1,idl1,Ch,Ch,Ch,C,C,Wa(iw))
          na = 0
        ELSE
          CALL RADFG(ido,ip,l1,idl1,C,C,C,Ch,Ch,Wa(iw))
          na = 1
        ENDIF
      ELSE
        ix2 = iw + ido
        ix3 = ix2 + ido
        ix4 = ix3 + ido
        IF ( na/=0 ) THEN
          CALL RADF5(ido,l1,Ch,C,Wa(iw),Wa(ix2),Wa(ix3),Wa(ix4))
        ELSE
          CALL RADF5(ido,l1,C,Ch,Wa(iw),Wa(ix2),Wa(ix3),Wa(ix4))
        ENDIF
      ENDIF
    ELSEIF ( na/=0 ) THEN
      CALL RADF2(ido,l1,Ch,C,Wa(iw))
    ELSE
      CALL RADF2(ido,l1,C,Ch,Wa(iw))
    ENDIF
    l2 = l1
  ENDDO
  IF ( na==1 ) RETURN
  DO i = 1, N
    C(i) = Ch(i)
  ENDDO
END SUBROUTINE RFFTF1
