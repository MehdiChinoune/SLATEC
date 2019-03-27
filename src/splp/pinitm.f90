!** PINITM
SUBROUTINE PINITM(M,N,Sx,Ix,Lmx,Ipagef)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to SPLP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (PINITM-S, DPINTM-D)
  !***
  ! **Author:**  Hanson, R. J., (SNLA)
  !           Wisniewski, J. A., (SNLA)
  !***
  ! **Description:**
  !
  !     PINITM LIMITS THE TYPE OF STORAGE TO A SEQUENTIAL SCHEME.
  !     THE MATRIX IS STORED BY COLUMNS.
  !     SPARSE MATRIX INITIALIZATION SUBROUTINE.
  !
  !            M=NUMBER OF ROWS OF THE MATRIX.
  !            N=NUMBER OF COLUMNS OF THE MATRIX.
  !  SX(*),IX(*)=THE WORK ARRAYS WHICH ARE USED TO STORE THE SPARSE
  !              MATRIX.  THESE ARRAYS ARE AUTOMATICALLY MAINTAINED BY
  !              THE PACKAGE FOR THE USER.
  !          LMX=LENGTH OF THE WORK ARRAY SX(*).
  !              LMX MUST BE AT LEAST N+7 WHERE
  !              FOR GREATEST EFFICIENCY LMX SHOULD BE AT LEAST N+NZ+6
  !              WHERE NZ IS THE MAXIMUM NUMBER OF NONZEROES TO BE
  !              STORED IN THE MATRIX.  VALUES OF LMX BETWEEN N+7 AND
  !              N+NZ+6 WILL CAUSE DEMAND PAGING TO OCCUR.
  !              THIS IS IMPLEMENTED BY THE PACKAGE.
  !              IX(*) MUST BE DIMENSIONED AT LEAST LMX
  !      IPAGEF=UNIT NUMBER WHERE DEMAND PAGES WILL BE STORED.
  !
  !     THIS SUBROUTINE IS A MODIFICATION OF THE SUBROUTINE LINITM,
  !     SANDIA LABS. REPT. SAND78-0785.
  !     MODIFICATIONS BY K.L. HIEBERT AND R.J. HANSON
  !     REVISED 811130-1000
  !     REVISED YYMMDD-HHMM
  !
  !***
  ! **See also:**  SPLP
  !***
  ! **Routines called:**  XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900328  Added TYPE section.  (WRB)
  !   910403  Updated AUTHOR and DESCRIPTION sections.  (WRB)
  
  INTEGER i, iopt, Ipagef, Ix(*), Lmx, lp4, M, N, n20008, n20012, nerr
  REAL Sx(Lmx), zero, one
  SAVE zero, one
  DATA zero, one/0.E0, 1.E0/
  !* FIRST EXECUTABLE STATEMENT  PINITM
  iopt = 1
  !
  !     CHECK FOR INPUT ERRORS.
  !
  IF ( M<=0.OR.N<=0 ) THEN
    nerr = 55
    CALL XERMSG('SLATEC','PINITM','MATRIX DIMENSION M OR N .LE. 0.',nerr,iopt)
  ENDIF
  !
  !     VERIFY IF VALUE OF LMX IS LARGE ENOUGH.
  !
  IF ( Lmx<N+7 ) THEN
    nerr = 55
    CALL XERMSG('SLATEC','PINITM','THE VALUE OF LMX IS TOO SMALL.',nerr,iopt)
  ENDIF
  !
  !     INITIALIZE DATA STRUCTURE INDEPENDENT VALUES.
  !
  Sx(1) = zero
  Sx(2) = zero
  Sx(3) = Ipagef
  Ix(1) = Lmx
  Ix(2) = M
  Ix(3) = N
  Ix(4) = 0
  Sx(Lmx-1) = zero
  Sx(Lmx) = -one
  Ix(Lmx-1) = -1
  lp4 = N + 4
  !
  !     INITIALIZE DATA STRUCTURE DEPENDENT VALUES.
  !
  i = 4
  n20008 = lp4
  DO WHILE ( (n20008-i)>=0 )
    Sx(i) = zero
    i = i + 1
  ENDDO
  i = 5
  n20012 = lp4
  DO WHILE ( (n20012-i)>=0 )
    Ix(i) = lp4
    i = i + 1
  ENDDO
  Sx(N+5) = zero
  Ix(N+5) = 0
  Ix(Lmx) = 0
  !
  !     INITIALIZATION COMPLETE.
  !
END SUBROUTINE PINITM
