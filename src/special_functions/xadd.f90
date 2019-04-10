!** XADD
SUBROUTINE XADD(X,Ix,Y,Iy,Z,Iz,Ierror)
  USE XBLK ,ONLY: RADixx, RADixl, RAD2l, L
  IMPLICIT NONE
  !>
  !***
  !  To provide single-precision floating-point arithmetic
  !            with an extended exponent range.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  A3D
  !***
  ! **Type:**      SINGLE PRECISION (XADD-S, DXADD-D)
  !***
  ! **Keywords:**  EXTENDED-RANGE SINGLE-PRECISION ARITHMETIC
  !***
  ! **Author:**  Lozier, Daniel W., (National Bureau of Standards)
  !           Smith, John M., (NBS and George Mason University)
  !***
  ! **Description:**
  !     REAL X, Y, Z
  !     INTEGER IX, IY, IZ
  !
  !                  FORMS THE EXTENDED-RANGE SUM  (Z,IZ) =
  !                  (X,IX) + (Y,IY).  (Z,IZ) IS ADJUSTED
  !                  BEFORE RETURNING. THE INPUT OPERANDS
  !                  NEED NOT BE IN ADJUSTED FORM, BUT THEIR
  !                  PRINCIPAL PARTS MUST SATISFY
  !                  RADIX**(-2L).LE.ABS(X).LE.RADIX**(2L),
  !                  RADIX**(-2L).LE.ABS(Y).LE.RADIX**(2L).
  !
  !***
  ! **See also:**  XSET
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  XADJ
  !***
  ! COMMON BLOCKS    XBLK2

  !* REVISION HISTORY  (YYMMDD)
  !   820712  DATE WRITTEN
  !   881020  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
  !   901019  Revisions to prologue.  (DWL and WRB)
  !   901106  Changed all specific intrinsics to generic.  (WRB)
  !           Corrected order of sections in prologue and added TYPE
  !           section.  (WRB)
  !   920127  Revised PURPOSE section of prologue.  (DWL)

  INTEGER i, i1, i2, Ierror, is, j, Ix, Iy, Iz
  REAL s, t
  REAL X, Y, Z
  !
  !
  !   THE CONDITIONS IMPOSED ON L AND KMAX BY THIS SUBROUTINE
  ! ARE
  !     (1) 1 .LT. L .LE. 0.5*LOGR(0.5*DZERO)
  !
  !     (2) NRADPL .LT. L .LE. KMAX/6
  !
  !     (3) KMAX .LE. (2**NBITS - 4*L - 1)/2
  !
  ! THESE CONDITIONS MUST BE MET BY APPROPRIATE CODING
  ! IN SUBROUTINE XSET.
  !
  !* FIRST EXECUTABLE STATEMENT  XADD
  Ierror = 0
  IF ( X==0.0 ) THEN
    Z = Y
    Iz = Iy
    CALL XADJ(Z,Iz,Ierror)
    RETURN
  ELSEIF ( Y/=0.0 ) THEN
    IF ( Ix<0.OR.Iy<0 ) THEN
      IF ( Ix>=0.OR.Iy>=0 ) THEN
        IF ( ABS(Ix)>6*L.OR.ABS(Iy)>6*L ) THEN
          IF ( Ix>=0 ) THEN
            Z = X
            Iz = Ix
          ELSE
            Z = Y
            Iz = Iy
          END IF
          CALL XADJ(Z,Iz,Ierror)
          RETURN
        END IF
      END IF
    END IF
    i = Ix - Iy
    IF ( i<0 ) THEN
      s = Y
      is = Iy
      t = X
    ELSEIF ( i==0 ) THEN
      IF ( ABS(X)>1.0.AND.ABS(Y)>1.0 ) THEN
        s = X/RADixl
        t = Y/RADixl
        Z = s + t
        Iz = Ix + L
      ELSEIF ( ABS(X)<1.0.AND.ABS(Y)<1.0 ) THEN
        s = X*RADixl
        t = Y*RADixl
        Z = s + t
        Iz = Ix - L
      ELSE
        Z = X + Y
        Iz = Ix
      END IF
      CALL XADJ(Z,Iz,Ierror)
      RETURN
    ELSE
      s = X
      is = Ix
      t = Y
    END IF
    !
    !  AT THIS POINT, THE ONE OF (X,IX) OR (Y,IY) THAT HAS THE
    ! LARGER AUXILIARY INDEX IS STORED IN (S,IS). THE PRINCIPAL
    ! PART OF THE OTHER INPUT IS STORED IN T.
    !
    i1 = ABS(i)/L
    i2 = MOD(ABS(i),L)
    IF ( ABS(t)>=RADixl ) THEN
      j = i1 - 2
      IF ( j<0 ) GOTO 200
      t = t*RADixx**(-i2)/RAD2l
      GOTO 300
    ELSE
      IF ( ABS(t)>=1.0 ) GOTO 200
      IF ( RADixl*ABS(t)<1.0 ) THEN
        j = i1 + 1
        t = t*RADixx**(L-i2)
        GOTO 300
      END IF
    END IF
  ELSE
    Z = X
    Iz = Ix
    CALL XADJ(Z,Iz,Ierror)
    RETURN
  END IF
  100  j = i1
  t = t*RADixx**(-i2)
  GOTO 300
  200  j = i1 - 1
  IF ( j<0 ) GOTO 100
  t = t*RADixx**(-i2)/RADixl
  !
  !  AT THIS POINT, SOME OR ALL OF THE DIFFERENCE IN THE
  ! AUXILIARY INDICES HAS BEEN USED TO EFFECT A LEFT SHIFT
  ! OF T.  THE SHIFTED VALUE OF T SATISFIES
  !
  !       RADIX**(-2*L) .LE. ABS(T) .LE. 1.0
  !
  ! AND, IF J=0, NO FURTHER SHIFTING REMAINS TO BE DONE.
  !
  300 CONTINUE
  IF ( j/=0 ) THEN
    IF ( ABS(s)<RADixl.AND.j<=3 ) THEN
      IF ( ABS(s)>=1.0 ) THEN
        SELECT CASE (j)
          CASE (1)
            s = s*RADixl
            GOTO 400
          CASE (2,3)
            GOTO 350
          CASE DEFAULT
        END SELECT
      END IF
      IF ( RADixl*ABS(s)>=1.0 ) THEN
        SELECT CASE (j)
          CASE (1)
            s = s*RADixl
          CASE (2)
            s = s*RADixl
            s = s*RADixl
          CASE (3)
            GOTO 350
          CASE DEFAULT
            GOTO 320
        END SELECT
        GOTO 400
      END IF
      320 CONTINUE
      SELECT CASE (j)
        CASE (1)
          s = s*RADixl
          GOTO 400
        CASE (2)
        CASE (3)
          s = s*RADixl
        CASE DEFAULT
          GOTO 350
      END SELECT
      s = s*RADixl
      s = s*RADixl
      GOTO 400
    END IF
    350  Z = s
    Iz = is
    CALL XADJ(Z,Iz,Ierror)
    RETURN
  END IF
  !
  !   AT THIS POINT, THE REMAINING DIFFERENCE IN THE
  ! AUXILIARY INDICES HAS BEEN USED TO EFFECT A RIGHT SHIFT
  ! OF S.  IF THE SHIFTED VALUE OF S WOULD HAVE EXCEEDED
  ! RADIX**L, THEN (S,IS) IS RETURNED AS THE VALUE OF THE
  ! SUM.
  !
  400 CONTINUE
  IF ( ABS(s)>1.0.AND.ABS(t)>1.0 ) THEN
    s = s/RADixl
    t = t/RADixl
    Z = s + t
    Iz = is - j*L + L
  ELSEIF ( ABS(s)<1.0.AND.ABS(t)<1.0 ) THEN
    s = s*RADixl
    t = t*RADixl
    Z = s + t
    Iz = is - j*L - L
  ELSE
    Z = s + t
    Iz = is - j*L
  END IF
  CALL XADJ(Z,Iz,Ierror)
  RETURN
END SUBROUTINE XADD
