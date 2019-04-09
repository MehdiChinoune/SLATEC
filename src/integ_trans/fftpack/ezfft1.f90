!** EZFFT1
SUBROUTINE EZFFT1(N,Wa,Ifac)
  IMPLICIT NONE
  !>
  !***
  !  EZFFTI calls EZFFT1 with appropriate work array
  !            partitioning.
  !***
  ! **Library:**   SLATEC (FFTPACK)
  !***
  ! **Type:**      SINGLE PRECISION (EZFFT1-S)
  !***
  ! **Author:**  Swarztrauber, P. N., (NCAR)
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   830401  Modified to use SLATEC library source file format.
  !   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
  !           (a) changing dummy array size declarations (1) to (*),
  !           (b) changing references to intrinsic function FLOAT
  !               to REAL, and
  !           (c) changing definition of variable TPI by using
  !               FORTRAN intrinsic function ATAN instead of a DATA
  !               statement.
  !   881128  Modified by Dick Valent to meet prologue standards.
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)

  REAL arg1, argh, ch1, ch1h, dch1, dsh1, sh1, tpi, Wa(*)
  INTEGER i, ib, ido, Ifac(*), ii, ip, ipm, is, j, k1, l1, l2, N, &
    nf, nfm1, nl, nq, nr, ntry
  INTEGER, PARAMETER :: ntryh(4) = [ 4, 2, 3, 5 ]
  !* FIRST EXECUTABLE STATEMENT  EZFFT1
  tpi = 8.*ATAN(1.)
  nl = N
  nf = 0
  j = 0
  100  j = j + 1
  IF ( j<=4 ) THEN
    ntry = ntryh(j)
  ELSE
    ntry = ntry + 2
  END IF
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
        END DO
        Ifac(3) = 2
      END IF
    END IF
    IF ( nl==1 ) THEN
      Ifac(1) = N
      Ifac(2) = nf
      argh = tpi/N
      is = 0
      nfm1 = nf - 1
      l1 = 1
      IF ( nfm1==0 ) RETURN
      DO k1 = 1, nfm1
        ip = Ifac(k1+2)
        l2 = l1*ip
        ido = N/l2
        ipm = ip - 1
        arg1 = l1*argh
        ch1 = 1.
        sh1 = 0.
        dch1 = COS(arg1)
        dsh1 = SIN(arg1)
        DO j = 1, ipm
          ch1h = dch1*ch1 - dsh1*sh1
          sh1 = dch1*sh1 + dsh1*ch1
          ch1 = ch1h
          i = is + 2
          Wa(i-1) = ch1
          Wa(i) = sh1
          IF ( ido>=5 ) THEN
            DO ii = 5, ido, 2
              i = i + 2
              Wa(i-1) = ch1*Wa(i-3) - sh1*Wa(i-2)
              Wa(i) = ch1*Wa(i-2) + sh1*Wa(i-3)
            END DO
          END IF
          is = is + ido
        END DO
        l1 = l2
      END DO
      EXIT
    END IF
  END DO
END SUBROUTINE EZFFT1
