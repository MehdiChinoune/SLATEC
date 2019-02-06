!*==CFFTB1.f90  processed by SPAG 6.72Dc at 10:55 on  6 Feb 2019
!DECK CFFTB1
      SUBROUTINE CFFTB1(N,C,Ch,Wa,Ifac)
      IMPLICIT NONE
!*--CFFTB15
!*** Start of declarations inserted by SPAG
      REAL C , Ch , Wa
      INTEGER i , idl1 , ido , idot , Ifac , ip , iw , ix2 , ix3 , ix4 , k1 , 
     &        l1 , l2 , N , n2 , na , nac , nf
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  CFFTB1
!***PURPOSE  Compute the unnormalized inverse of CFFTF1.
!***LIBRARY   SLATEC (FFTPACK)
!***CATEGORY  J1A2
!***TYPE      COMPLEX (RFFTB1-S, CFFTB1-C)
!***KEYWORDS  FFTPACK, FOURIER TRANSFORM
!***AUTHOR  Swarztrauber, P. N., (NCAR)
!***DESCRIPTION
!
!  Subroutine CFFTB1 computes the backward complex discrete Fourier
!  transform (the Fourier synthesis).  Equivalently, CFFTB1 computes
!  a complex periodic sequence from its Fourier coefficients.
!  The transform is defined below at output parameter C.
!
!  A call of CFFTF1 followed by a call of CFFTB1 will multiply the
!  sequence by N.
!
!  The arrays WA and IFAC which are used by subroutine CFFTB1 must be
!  initialized by calling subroutine CFFTI1 (N, WA, IFAC).
!
!  Input Parameters
!
!  N       the length of the complex sequence C.  The method is
!          more efficient when N is the product of small primes.
!
!  C       a complex array of length N which contains the sequence
!
!  CH      a real work array of length at least 2*N
!
!  WA      a real work array which must be dimensioned at least 2*N.
!
!  IFAC    an integer work array which must be dimensioned at least 15.
!
!          The WA and IFAC arrays must be initialized by calling
!          subroutine CFFTI1 (N, WA, IFAC), and different WA and IFAC
!          arrays must be used for each different value of N.  This
!          initialization does not have to be repeated so long as N
!          remains unchanged.  Thus subsequent transforms can be
!          obtained faster than the first.  The same WA and IFAC arrays
!          can be used by CFFTF1 and CFFTB1.
!
!  Output Parameters
!
!  C       For J=1,...,N
!
!              C(J)=the sum from K=1,...,N of
!
!                 C(K)*EXP(I*(J-1)*(K-1)*2*PI/N)
!
!                         where I=SQRT(-1)
!
!  NOTE:   WA and IFAC contain initialization calculations which must
!          not be destroyed between calls of subroutine CFFTF1 or CFFTB1
!
!***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
!                 Computations (G. Rodrigue, ed.), Academic Press,
!                 1982, pp. 51-83.
!***ROUTINES CALLED  PASSB, PASSB2, PASSB3, PASSB4, PASSB5
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   830401  Modified to use SLATEC library source file format.
!   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
!           changing dummy array size declarations (1) to (*).
!   881128  Modified by Dick Valent to meet prologue standards.
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900131  Routine changed from subsidiary to user-callable.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CFFTB1
      DIMENSION Ch(*) , C(*) , Wa(*) , Ifac(*)
!***FIRST EXECUTABLE STATEMENT  CFFTB1
      nf = Ifac(2)
      na = 0
      l1 = 1
      iw = 1
      DO k1 = 1 , nf
        ip = Ifac(k1+2)
        l2 = ip*l1
        ido = N/l2
        idot = ido + ido
        idl1 = idot*l1
        IF ( ip==4 ) THEN
          ix2 = iw + idot
          ix3 = ix2 + idot
          IF ( na/=0 ) THEN
            CALL PASSB4(idot,l1,Ch,C,Wa(iw),Wa(ix2),Wa(ix3))
          ELSE
            CALL PASSB4(idot,l1,C,Ch,Wa(iw),Wa(ix2),Wa(ix3))
          ENDIF
          na = 1 - na
        ELSEIF ( ip==2 ) THEN
          IF ( na/=0 ) THEN
            CALL PASSB2(idot,l1,Ch,C,Wa(iw))
          ELSE
            CALL PASSB2(idot,l1,C,Ch,Wa(iw))
          ENDIF
          na = 1 - na
        ELSEIF ( ip==3 ) THEN
          ix2 = iw + idot
          IF ( na/=0 ) THEN
            CALL PASSB3(idot,l1,Ch,C,Wa(iw),Wa(ix2))
          ELSE
            CALL PASSB3(idot,l1,C,Ch,Wa(iw),Wa(ix2))
          ENDIF
          na = 1 - na
        ELSEIF ( ip/=5 ) THEN
          IF ( na/=0 ) THEN
            CALL PASSB(nac,idot,ip,l1,idl1,Ch,Ch,Ch,C,C,Wa(iw))
          ELSE
            CALL PASSB(nac,idot,ip,l1,idl1,C,C,C,Ch,Ch,Wa(iw))
          ENDIF
          IF ( nac/=0 ) na = 1 - na
        ELSE
          ix2 = iw + idot
          ix3 = ix2 + idot
          ix4 = ix3 + idot
          IF ( na/=0 ) THEN
            CALL PASSB5(idot,l1,Ch,C,Wa(iw),Wa(ix2),Wa(ix3),Wa(ix4))
          ELSE
            CALL PASSB5(idot,l1,C,Ch,Wa(iw),Wa(ix2),Wa(ix3),Wa(ix4))
          ENDIF
          na = 1 - na
        ENDIF
        l1 = l2
        iw = iw + (ip-1)*idot
      ENDDO
      IF ( na==0 ) RETURN
      n2 = N + N
      DO i = 1 , n2
        C(i) = Ch(i)
      ENDDO
      END SUBROUTINE CFFTB1
