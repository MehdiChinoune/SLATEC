!*==CBUNI.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK CBUNI
SUBROUTINE CBUNI(Z,Fnu,Kode,N,Y,Nz,Nui,Nlast,Fnul,Tol,Elim,Alim)
  IMPLICIT NONE
  !*--CBUNI5
  !***BEGIN PROLOGUE  CBUNI
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to CBESI and CBESK
  !***LIBRARY   SLATEC
  !***TYPE      ALL (CBUNI-A, ZBUNI-A)
  !***AUTHOR  Amos, D. E., (SNL)
  !***DESCRIPTION
  !
  !     CBUNI COMPUTES THE I BESSEL FUNCTION FOR LARGE ABS(Z).GT.
  !     FNUL AND FNU+N-1.LT.FNUL. THE ORDER IS INCREASED FROM
  !     FNU+N-1 GREATER THAN FNUL BY ADDING NUI AND COMPUTING
  !     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR I(FNU,Z)
  !     ON IFORM=1 AND THE EXPANSION FOR J(FNU,Z) ON IFORM=2
  !
  !***SEE ALSO  CBESI, CBESK
  !***ROUTINES CALLED  CUNI1, CUNI2, R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  CBUNI
  COMPLEX cscl, cscr, cy, rz, st, s1, s2, Y, Z
  REAL Alim, ax, ay, dfnu, Elim, Fnu, fnui, Fnul, gnu, Tol, xx, &
    yy, ascle, bry, str, sti, stm, R1MACH
  INTEGER i, iflag, iform, k, Kode, N, nl, Nlast, Nui, nw, Nz
  DIMENSION Y(N), cy(2), bry(3)
  !***FIRST EXECUTABLE STATEMENT  CBUNI
  Nz = 0
  xx = REAL(Z)
  yy = AIMAG(Z)
  ax = ABS(xx)*1.7321E0
  ay = ABS(yy)
  iform = 1
  IF ( ay>ax ) iform = 2
  IF ( Nui==0 ) THEN
    IF ( iform==2 ) THEN
      !-----------------------------------------------------------------------
      !     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
      !     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
      !     AND HPI=PI/2
      !-----------------------------------------------------------------------
      CALL CUNI2(Z,Fnu,Kode,N,Y,nw,Nlast,Fnul,Tol,Elim,Alim)
    ELSE
      !-----------------------------------------------------------------------
      !     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
      !     -PI/3.LE.ARG(Z).LE.PI/3
      !-----------------------------------------------------------------------
      CALL CUNI1(Z,Fnu,Kode,N,Y,nw,Nlast,Fnul,Tol,Elim,Alim)
    ENDIF
    IF ( nw>=0 ) GOTO 100
  ELSE
    fnui = Nui
    dfnu = Fnu + (N-1)
    gnu = dfnu + fnui
    IF ( iform==2 ) THEN
      !-----------------------------------------------------------------------
      !     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
      !     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
      !     AND HPI=PI/2
      !-----------------------------------------------------------------------
      CALL CUNI2(Z,gnu,Kode,2,cy,nw,Nlast,Fnul,Tol,Elim,Alim)
    ELSE
      !-----------------------------------------------------------------------
      !     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
      !     -PI/3.LE.ARG(Z).LE.PI/3
      !-----------------------------------------------------------------------
      CALL CUNI1(Z,gnu,Kode,2,cy,nw,Nlast,Fnul,Tol,Elim,Alim)
    ENDIF
    IF ( nw>=0 ) THEN
      IF ( nw/=0 ) THEN
        Nlast = N
        GOTO 99999
      ELSE
        ay = ABS(cy(1))
        !----------------------------------------------------------------------
        !     SCALE BACKWARD RECURRENCE, BRY(3) IS DEFINED BUT NEVER USED
        !----------------------------------------------------------------------
        bry(1) = 1.0E+3*R1MACH(1)/Tol
        bry(2) = 1.0E0/bry(1)
        bry(3) = bry(2)
        iflag = 2
        ascle = bry(2)
        ax = 1.0E0
        cscl = CMPLX(ax,0.0E0)
        IF ( ay<=bry(1) ) THEN
          iflag = 1
          ascle = bry(1)
          ax = 1.0E0/Tol
          cscl = CMPLX(ax,0.0E0)
        ELSEIF ( ay>=bry(2) ) THEN
          iflag = 3
          ascle = bry(3)
          ax = Tol
          cscl = CMPLX(ax,0.0E0)
        ENDIF
        ay = 1.0E0/ax
        cscr = CMPLX(ay,0.0E0)
        s1 = cy(2)*cscl
        s2 = cy(1)*cscl
        rz = CMPLX(2.0E0,0.0E0)/Z
        DO i = 1, Nui
          st = s2
          s2 = CMPLX(dfnu+fnui,0.0E0)*rz*s2 + s1
          s1 = st
          fnui = fnui - 1.0E0
          IF ( iflag<3 ) THEN
            st = s2*cscr
            str = REAL(st)
            sti = AIMAG(st)
            str = ABS(str)
            sti = ABS(sti)
            stm = MAX(str,sti)
            IF ( stm>ascle ) THEN
              iflag = iflag + 1
              ascle = bry(iflag)
              s1 = s1*cscr
              s2 = st
              ax = ax*Tol
              ay = 1.0E0/ax
              cscl = CMPLX(ax,0.0E0)
              cscr = CMPLX(ay,0.0E0)
              s1 = s1*cscl
              s2 = s2*cscl
            ENDIF
          ENDIF
        ENDDO
        Y(N) = s2*cscr
        IF ( N==1 ) RETURN
        nl = N - 1
        fnui = nl
        k = nl
        DO i = 1, nl
          st = s2
          s2 = CMPLX(Fnu+fnui,0.0E0)*rz*s2 + s1
          s1 = st
          st = s2*cscr
          Y(k) = st
          fnui = fnui - 1.0E0
          k = k - 1
          IF ( iflag<3 ) THEN
            str = REAL(st)
            sti = AIMAG(st)
            str = ABS(str)
            sti = ABS(sti)
            stm = MAX(str,sti)
            IF ( stm>ascle ) THEN
              iflag = iflag + 1
              ascle = bry(iflag)
              s1 = s1*cscr
              s2 = st
              ax = ax*Tol
              ay = 1.0E0/ax
              cscl = CMPLX(ax,0.0E0)
              cscr = CMPLX(ay,0.0E0)
              s1 = s1*cscl
              s2 = s2*cscl
            ENDIF
          ENDIF
        ENDDO
        RETURN
      ENDIF
    ENDIF
  ENDIF
  Nz = -1
  IF ( nw==(-2) ) Nz = -2
  RETURN
  100  Nz = nw
  RETURN
  99999 CONTINUE
  END SUBROUTINE CBUNI
