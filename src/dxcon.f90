!DECK DXCON
SUBROUTINE DXCON(X,Ix,Ierror)
  IMPLICIT NONE
  INTEGER i, i1, icase, Ierror, ispace, itemp, j, j1, j2
  !***BEGIN PROLOGUE  DXCON
  !***PURPOSE  To provide double-precision floating-point arithmetic
  !            with an extended exponent range.
  !***LIBRARY   SLATEC
  !***CATEGORY  A3D
  !***TYPE      DOUBLE PRECISION (XCON-S, DXCON-D)
  !***KEYWORDS  EXTENDED-RANGE DOUBLE-PRECISION ARITHMETIC
  !***AUTHOR  Lozier, Daniel W., (National Bureau of Standards)
  !           Smith, John M., (NBS and George Mason University)
  !***DESCRIPTION
  !     DOUBLE PRECISION X
  !     INTEGER IX
  !
  !                  CONVERTS (X,IX) = X*RADIX**IX
  !                  TO DECIMAL FORM IN PREPARATION FOR
  !                  PRINTING, SO THAT (X,IX) = X*10**IX
  !                  WHERE 1/10 .LE. ABS(X) .LT. 1
  !                  IS RETURNED, EXCEPT THAT IF
  !                  (ABS(X),IX) IS BETWEEN RADIX**(-2L)
  !                  AND RADIX**(2L) THEN THE REDUCED
  !                  FORM WITH IX = 0 IS RETURNED.
  !
  !***SEE ALSO  DXSET
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  DXADJ, DXC210, DXRED
  !***COMMON BLOCKS    DXBLK2
  !***REVISION HISTORY  (YYMMDD)
  !   820712  DATE WRITTEN
  !   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
  !   901019  Revisions to prologue.  (DWL and WRB)
  !   901106  Changed all specific intrinsics to generic.  (WRB)
  !           Corrected order of sections in prologue and added TYPE
  !           section.  (WRB)
  !   920127  Revised PURPOSE section of prologue.  (DWL)
  !***END PROLOGUE  DXCON
  REAL(8) :: X
  INTEGER Ix
  !
  !   THE CONDITIONS IMPOSED ON L AND KMAX BY THIS SUBROUTINE
  ! ARE
  !    (1) 4 .LE. L .LE. 2**NBITS - 1 - KMAX
  !
  !    (2) KMAX .LE. ((2**NBITS)-2)/LOG10R - L
  !
  ! THESE CONDITIONS MUST BE MET BY APPROPRIATE CODING
  ! IN SUBROUTINE DXSET.
  !
  REAL(8) :: RADix, RADixl, RAD2l, DLG10r
  INTEGER L, L2, KMAx
  COMMON /DXBLK2/ RADix, RADixl, RAD2l, DLG10r, L, L2, KMAx
  SAVE /DXBLK2/ , ispace
  !
  REAL(8) :: a, b, z
  !
  DATA ispace/1/
  !   THE PARAMETER ISPACE IS THE INCREMENT USED IN FORM-
  ! ING THE AUXILIARY INDEX OF THE DECIMAL EXTENDED-RANGE
  ! FORM. THE RETURNED VALUE OF IX WILL BE AN INTEGER MULT-
  ! IPLE OF ISPACE. ISPACE MUST SATISFY 1 .LE. ISPACE .LE.
  ! L/2. IF A VALUE GREATER THAN 1 IS TAKEN, THE RETURNED
  ! VALUE OF X WILL SATISFY 10**(-ISPACE) .LE. ABS(X) .LE. 1
  ! WHEN (ABS(X),IX) .LT. RADIX**(-2L), AND 1/10 .LE. ABS(X)
  ! .LT. 10**(ISPACE-1) WHEN (ABS(X),IX) .GT. RADIX**(2L).
  !
  !***FIRST EXECUTABLE STATEMENT  DXCON
  Ierror = 0
  CALL DXRED(X,Ix,Ierror)
  IF ( Ierror/=0 ) RETURN
  IF ( Ix/=0 ) THEN
    CALL DXADJ(X,Ix,Ierror)
    IF ( Ierror/=0 ) RETURN
    !
    ! CASE 1 IS WHEN (X,IX) IS LESS THAN RADIX**(-2L) IN MAGNITUDE,
    ! CASE 2 IS WHEN (X,IX) IS GREATER THAN RADIX**(2L) IN MAGNITUDE.
    itemp = 1
    icase = (3+SIGN(itemp,Ix))/2
    IF ( icase==2 ) THEN
      IF ( ABS(X)<1.0D0 ) THEN
        X = X*RADixl
        Ix = Ix - L
      ENDIF
    ELSEIF ( ABS(X)>=1.0D0 ) THEN
      X = X/RADixl
      Ix = Ix + L
    ENDIF
    !
    ! AT THIS POINT, RADIX**(-L) .LE. ABS(X) .LT. 1.0D0     IN CASE 1,
    !                      1.0D0 .LE. ABS(X) .LT. RADIX**L  IN CASE 2.
    i = INT( LOG10(ABS(X))/DLG10r )
    a = RADix**i
    IF ( icase==2 ) THEN
      DO WHILE ( a>ABS(X) )
        i = i - 1
        a = a/RADix
      ENDDO
      DO WHILE ( ABS(X)>=RADix*a )
        i = i + 1
        a = a*RADix
      ENDDO
    ELSE
      DO WHILE ( a>RADix*ABS(X) )
        i = i - 1
        a = a/RADix
      ENDDO
      DO WHILE ( ABS(X)>=a )
        i = i + 1
        a = a*RADix
      ENDDO
    ENDIF
    !
    ! AT THIS POINT I IS SUCH THAT
    ! RADIX**(I-1) .LE. ABS(X) .LT. RADIX**I      IN CASE 1,
    !     RADIX**I .LE. ABS(X) .LT. RADIX**(I+1)  IN CASE 2.
    itemp = INT( ispace/DLG10r )
    a = RADix**itemp
    b = 10.0D0**ispace
    DO WHILE ( a>b )
      itemp = itemp - 1
      a = a/RADix
    ENDDO
    DO WHILE ( b>=a*RADix )
      itemp = itemp + 1
      a = a*RADix
    ENDDO
    !
    ! AT THIS POINT ITEMP IS SUCH THAT
    ! RADIX**ITEMP .LE. 10**ISPACE .LT. RADIX**(ITEMP+1).
    IF ( itemp<=0 ) THEN
      ! ITEMP = 0 IF, AND ONLY IF, ISPACE = 1 AND RADIX = 16.0D0
      X = X*RADix**(-i)
      Ix = Ix + i
      CALL DXC210(Ix,z,j,Ierror)
      IF ( Ierror/=0 ) RETURN
      X = X*z
      Ix = j
      IF ( icase==1 ) GOTO 50
      IF ( icase==2 ) GOTO 100
    ENDIF
    i1 = i/itemp
    X = X*RADix**(-i1*itemp)
    Ix = Ix + i1*itemp
    !
    ! AT THIS POINT,
    ! RADIX**(-ITEMP) .LE. ABS(X) .LT. 1.0D0        IN CASE 1,
    !           1.0D0 .LE. ABS(X) .LT. RADIX**ITEMP IN CASE 2.
    CALL DXC210(Ix,z,j,Ierror)
    IF ( Ierror/=0 ) RETURN
    j1 = j/ispace
    j2 = j - j1*ispace
    X = X*z*10.0D0**j2
    Ix = j1*ispace
    !
    ! AT THIS POINT,
    !  10.0D0**(-2*ISPACE) .LE. ABS(X) .LT. 1.0D0                IN CASE 1,
    !           10.0D0**-1 .LE. ABS(X) .LT. 10.0D0**(2*ISPACE-1) IN CASE 2.
    IF ( icase==2 ) GOTO 100
    50 CONTINUE
    DO WHILE ( b*ABS(X)<1.0D0 )
      X = X*b
      Ix = Ix - ispace
    ENDDO
    RETURN
    100 CONTINUE
    DO WHILE ( 10.0D0*ABS(X)>=b )
      X = X/b
      Ix = Ix + ispace
    ENDDO
  ENDIF
  RETURN
END SUBROUTINE DXCON
