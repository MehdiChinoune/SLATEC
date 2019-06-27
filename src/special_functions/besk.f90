!** BESK
PURE SUBROUTINE BESK(X,Fnu,Kode,N,Y,Nz)
  !> Implement forward recursion on the three term recursion relation for
  !  a sequence of non-negative order Bessel functions K_{FNU+I-1}(X),
  !  or scaled Bessel functions
  !  EXP(X)*K_{FNU+I-1}(X), I=1,...,N for REAL(SP), positive X and non-negative orders FNU.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C10B3
  !***
  ! **Type:**      SINGLE PRECISION (BESK-S, DBESK-D)
  !***
  ! **Keywords:**  K BESSEL FUNCTION, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract
  !         BESK implements forward recursion on the three term
  !         recursion relation for a sequence of non-negative order Bessel
  !         functions K_{FNU+I-1}(X), or scaled Bessel functions
  !         EXP(X)*K_{FNU+I-1}(X), I=1,...,N for real X > 0.0E0 and
  !         non-negative orders FNU.  If FNU < NULIM, orders FNU and
  !         FNU+1 are obtained from BESKNU to start the recursion.  If
  !         FNU >= NULIM, the uniform asymptotic expansion is used for
  !         orders FNU and FNU+1 to start the recursion.  NULIM is 35 or
  !         70 depending on whether N=1 or N >= 2.  Under and overflow
  !         tests are made on the leading term of the asymptotic expansion
  !         before any extensive computation is done.
  !
  !     Description of Arguments
  !
  !         Input
  !           X      - X > 0.0E0
  !           FNU    - order of the initial K function, FNU >= 0.0E0
  !           KODE   - a parameter to indicate the scaling option
  !                    KODE=1 returns Y(I)=       K_{FNU+I-1}(X),
  !                                        I=1,...,N
  !                    KODE=2 returns Y(I)=EXP(X)*K_{FNU+I-1}(X),
  !                                        I=1,...,N
  !           N      - number of members in the sequence, N >= 1
  !
  !         Output
  !           y      - a vector whose first n components contain values
  !                    for the sequence
  !                    Y(I)=       K_{FNU+I-1}(X), I=1,...,N  or
  !                    Y(I)=EXP(X)*K_{FNU+I-1}(X), I=1,...,N
  !                    depending on KODE
  !           NZ     - number of components of Y set to zero due to
  !                    underflow with KODE=1,
  !                    NZ=0  , normal return, computation completed
  !                    NZ /= 0, first NZ components of Y set to zero
  !                             due to underflow, Y(I)=0.0E0, I=1,...,NZ
  !
  !     Error Conditions
  !         Improper input arguments - a fatal error
  !         Overflow - a fatal error
  !         Underflow with KODE=1 -  a non-fatal error (NZ /= 0)
  !
  !***
  ! **References:**  F. W. J. Olver, Tables of Bessel Functions of Moderate
  !                 or Large Orders, NPL Mathematical Tables 6, Her
  !                 Majesty's Stationery Office, London, 1962.
  !               N. M. Temme, On the numerical evaluation of the modified
  !                 Bessel function of the third kind, Journal of
  !                 Computational Physics 19, (1975), pp. 324-337.
  !***
  ! **Routines called:**  ASYIK, BESK0, BESK0E, BESK1, BESK1E, BESKNU,
  !                    I1MACH, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   790201  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  USE service, ONLY : R1MACH, I1MACH
  !
  INTEGER, INTENT(IN) :: Kode, N
  INTEGER, INTENT(OUT) :: Nz
  REAL(SP), INTENT(IN) :: Fnu, X
  REAL(SP), INTENT(OUT) :: Y(N)
  INTEGER :: i, j, k, mz, nb, nd, nn, nud
  REAL(SP) :: cn, dnu, elim, etx, flgik, fn, fnn, gln, gnu, rtz, &
    s, s1, s2, t, tm, trx, w(2), xlim, zn
  INTEGER, PARAMETER :: nulim(2) = [ 35, 70 ]
  !* FIRST EXECUTABLE STATEMENT  BESK
  nn = -I1MACH(12)
  elim = 2.303_SP*(nn*R1MACH(5)-3._SP)
  xlim = R1MACH(1)*1.E+3_SP
  IF( Kode<1 .OR. Kode>2 ) THEN
    ERROR STOP 'BESK : SCALING OPTION, KODE, NOT 1 OR 2'
  ELSEIF( Fnu<0._SP ) THEN
    ERROR STOP 'BESK : ORDER, FNU < 0'
  ELSEIF( X<=0._SP ) THEN
    ERROR STOP 'BESK : X <= 0'
  ELSEIF( X<xlim ) THEN
    ERROR STOP 'BESK : OVERFLOW, FNU OR N TOO LARGE OR X TOO SMALL'
  ELSEIF( N<1 ) THEN
    ERROR STOP 'BESK : N < 1'
  ELSE
    etx = Kode - 1
    !
    !     ND IS A DUMMY VARIABLE FOR N
    !     GNU IS A DUMMY VARIABLE FOR FNU
    !     NZ = NUMBER OF UNDERFLOWS ON KODE=1
    !
    nd = N
    Nz = 0
    nud = INT(Fnu)
    dnu = Fnu - nud
    gnu = Fnu
    nn = MIN(2,nd)
    fn = Fnu + N - 1
    fnn = fn
    IF( fn<2._SP ) THEN
      !
      !     UNDERFLOW TEST FOR KODE=1
      IF( Kode==2 ) GOTO 600
      IF( X<=elim ) GOTO 600
      GOTO 700
    ELSE
      !
      !     OVERFLOW TEST  (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION)
      !     FOR THE LAST ORDER, FNU+N-1>=NULIM
      !
      zn = X/fn
      IF( zn==0._SP ) THEN
        ERROR STOP 'BESK : OVERFLOW, FNU OR N TOO LARGE OR X TOO SMALL'
      ELSE
        rtz = SQRT(1._SP+zn*zn)
        gln = LOG((1._SP+rtz)/zn)
        t = rtz*(1._SP-etx) + etx/(zn+rtz)
        cn = -fn*(t-gln)
        IF( cn>elim ) THEN
          ERROR STOP 'BESK : OVERFLOW, FNU OR N TOO LARGE OR X TOO SMALL'
        ELSEIF( nud<nulim(nn) ) THEN
          !
          IF( Kode==2 ) GOTO 300
          !
          !     UNDERFLOW TEST (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION IN X)
          !     FOR ORDER DNU
          !
          IF( X<=elim ) GOTO 300
          GOTO 700
        ELSEIF( nn==1 ) THEN
          GOTO 200
        END IF
      END IF
    END IF
  END IF
  !
  !     UNDERFLOW TEST (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION)
  !     FOR THE FIRST ORDER, FNU>=NULIM
  !
  100  fn = gnu
  zn = X/fn
  rtz = SQRT(1._SP+zn*zn)
  gln = LOG((1._SP+rtz)/zn)
  t = rtz*(1._SP-etx) + etx/(zn+rtz)
  cn = -fn*(t-gln)
  200 CONTINUE
  IF( cn<-elim ) GOTO 700
  !
  !     ASYMPTOTIC EXPANSION FOR ORDERS FNU AND FNU+1>=NULIM
  !
  flgik = -1._SP
  CALL ASYIK(X,gnu,Kode,flgik,rtz,cn,nn,Y)
  IF( nn==1 ) GOTO 800
  trx = 2._SP/X
  tm = (gnu+gnu+2._SP)/X
  GOTO 500
  300 CONTINUE
  IF( dnu/=0._SP ) THEN
    nb = 2
    IF( nud==0 .AND. nd==1 ) nb = 1
    CALL BESKNU(X,dnu,Kode,nb,w,Nz)
    s1 = w(1)
    IF( nb==1 ) GOTO 400
    s2 = w(2)
  ELSE
    IF( Kode==2 ) THEN
      s1 = BESK0E(X)
    ELSE
      s1 = BESK0(X)
    END IF
    IF( nud==0 .AND. nd==1 ) GOTO 400
    IF( Kode==2 ) THEN
      s2 = BESK1E(X)
    ELSE
      s2 = BESK1(X)
    END IF
  END IF
  trx = 2._SP/X
  tm = (dnu+dnu+2._SP)/X
  !     FORWARD RECUR FROM DNU TO FNU+1 TO GET Y(1) AND Y(2)
  IF( nd==1 ) nud = nud - 1
  IF( nud>0 ) THEN
    DO i = 1, nud
      s = s2
      s2 = tm*s2 + s1
      s1 = s
      tm = tm + trx
    END DO
    IF( nd==1 ) s1 = s2
  ELSEIF( nd<=1 ) THEN
    s1 = s2
  END IF
  400  Y(1) = s1
  IF( nd==1 ) GOTO 800
  Y(2) = s2
  500 CONTINUE
  IF( nd/=2 ) THEN
    !     FORWARD RECUR FROM FNU+2 TO FNU+N-1
    DO i = 3, nd
      Y(i) = tm*Y(i-1) + Y(i-2)
      tm = tm + trx
    END DO
  END IF
  GOTO 800
  !     OVERFLOW TEST
  600 CONTINUE
  IF( fn>1._SP ) THEN
    IF( -fn*(LOG(X)-0.693E0_SP)>elim ) THEN
      ERROR STOP 'BESK : OVERFLOW, FNU OR N TOO LARGE OR X TOO SMALL'
    END IF
  END IF
  IF( dnu==0._SP ) THEN
    j = nud
    IF( j/=1 ) THEN
      j = j + 1
      IF( Kode==2 ) THEN
        Y(j) = BESK0E(X)
      ELSE
        Y(j) = BESK0(X)
      END IF
      IF( nd==1 ) GOTO 800
      j = j + 1
    END IF
    IF( Kode==2 ) THEN
      Y(j) = BESK1E(X)
    ELSE
      Y(j) = BESK1(X)
    END IF
  ELSE
    CALL BESKNU(X,Fnu,Kode,nd,Y,mz)
  END IF
  GOTO 800
  700 CONTINUE
  DO
    !
    !     UPDATE PARAMETERS ON UNDERFLOW
    !
    nud = nud + 1
    nd = nd - 1
    IF( nd==0 ) EXIT
    nn = MIN(2,nd)
    gnu = gnu + 1._SP
    IF( fnn>=2._SP ) THEN
      IF( nud>=nulim(nn) ) GOTO 100
    END IF
  END DO
  800  Nz = N - nd
  IF( Nz==0 ) RETURN
  IF( nd/=0 ) THEN
    DO i = 1, nd
      j = N - i + 1
      k = nd - i + 1
      Y(j) = Y(k)
    END DO
  END IF
  DO i = 1, Nz
    Y(i) = 0._SP
  END DO

  RETURN
END SUBROUTINE BESK