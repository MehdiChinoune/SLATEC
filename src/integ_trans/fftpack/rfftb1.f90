!** RFFTB1
SUBROUTINE RFFTB1(N,C,Ch,Wa,Ifac)
  !>
  !***
  !  Compute the backward fast Fourier transform of a real
  !            coefficient array.
  !***
  ! **Library:**   SLATEC (FFTPACK)
  !***
  ! **Category:**  J1A1
  !***
  ! **Type:**      SINGLE PRECISION (RFFTB1-S, CFFTB1-C)
  !***
  ! **Keywords:**  FFTPACK, FOURIER TRANSFORM
  !***
  ! **Author:**  Swarztrauber, P. N., (NCAR)
  !***
  ! **Description:**
  !
  !   Subroutine RFFTB1 computes the real periodic sequence from its
  !   Fourier coefficients (Fourier synthesis).  The transform is defined
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
  !   C       For N even and for I = 1,...,N
  !
  !                C(I) = C(1)+(-1)**(I-1)*C(N)
  !
  !                     plus the sum from K=2 to K=N/2 of
  !
  !                      2.*C(2*K-2)*COS((K-1)*(I-1)*2*PI/N)
  !
  !                     -2.*C(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)
  !
  !           For N odd and for I = 1,...,N
  !
  !                C(I) = C(1) plus the sum from K=2 to K=(N+1)/2 of
  !
  !                     2.*C(2*K-2)*COS((K-1)*(I-1)*2*PI/N)
  !
  !                    -2.*C(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)
  !
  !   Notes:  This transform is unnormalized since a call of RFFTF1
  !           followed by a call of RFFTB1 will multiply the input
  !           sequence by N.
  !
  !           WA and IFAC contain initialization calculations which must
  !           not be destroyed between calls of subroutine RFFTF1 or
  !           RFFTB1.
  !
  !***
  ! **References:**  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
  !                 Computations (G. Rodrigue, ed.), Academic Press,
  !                 1982, pp. 51-83.
  !***
  ! **Routines called:**  RADB2, RADB3, RADB4, RADB5, RADBG

  !* REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   830401  Modified to use SLATEC library source file format.
  !   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
  !           changing dummy array size declarations (1) to (*).
  !   881128  Modified by Dick Valent to meet prologue standards.
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900131  Routine changed from subsidiary to user-callable.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  REAL C(*), Ch(*), Wa(*)
  INTEGER i, idl1, ido, Ifac(*), ip, iw, ix2, ix3, ix4, k1, l1, l2, N, na, nf
  !* FIRST EXECUTABLE STATEMENT  RFFTB1
  nf = Ifac(2)
  na = 0
  l1 = 1
  iw = 1
  DO k1 = 1, nf
    ip = Ifac(k1+2)
    l2 = ip*l1
    ido = N/l2
    idl1 = ido*l1
    IF ( ip==4 ) THEN
      ix2 = iw + ido
      ix3 = ix2 + ido
      IF ( na/=0 ) THEN
        CALL RADB4(ido,l1,Ch,C,Wa(iw),Wa(ix2),Wa(ix3))
      ELSE
        CALL RADB4(ido,l1,C,Ch,Wa(iw),Wa(ix2),Wa(ix3))
      END IF
      na = 1 - na
    ELSEIF ( ip==2 ) THEN
      IF ( na/=0 ) THEN
        CALL RADB2(ido,l1,Ch,C,Wa(iw))
      ELSE
        CALL RADB2(ido,l1,C,Ch,Wa(iw))
      END IF
      na = 1 - na
    ELSEIF ( ip==3 ) THEN
      ix2 = iw + ido
      IF ( na/=0 ) THEN
        CALL RADB3(ido,l1,Ch,C,Wa(iw),Wa(ix2))
      ELSE
        CALL RADB3(ido,l1,C,Ch,Wa(iw),Wa(ix2))
      END IF
      na = 1 - na
    ELSEIF ( ip/=5 ) THEN
      IF ( na/=0 ) THEN
        CALL RADBG(ido,ip,l1,idl1,Ch,Ch,Ch,C,C,Wa(iw))
      ELSE
        CALL RADBG(ido,ip,l1,idl1,C,C,C,Ch,Ch,Wa(iw))
      END IF
      IF ( ido==1 ) na = 1 - na
    ELSE
      ix2 = iw + ido
      ix3 = ix2 + ido
      ix4 = ix3 + ido
      IF ( na/=0 ) THEN
        CALL RADB5(ido,l1,Ch,C,Wa(iw),Wa(ix2),Wa(ix3),Wa(ix4))
      ELSE
        CALL RADB5(ido,l1,C,Ch,Wa(iw),Wa(ix2),Wa(ix3),Wa(ix4))
      END IF
      na = 1 - na
    END IF
    l1 = l2
    iw = iw + (ip-1)*ido
  END DO
  IF ( na==0 ) RETURN
  DO i = 1, N
    C(i) = Ch(i)
  END DO
END SUBROUTINE RFFTB1
