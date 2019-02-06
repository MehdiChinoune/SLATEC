!*==DEVPCK.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DEVPCK
      SUBROUTINE DEVPCK(Lout,Kprint,X,Y,F,Fx,Fy,Xe,Ye,Fe,De,Fe2,Fail)
      IMPLICIT NONE
!*--DEVPCK5
!***BEGIN PROLOGUE  DEVPCK
!***SUBSIDIARY
!***PURPOSE  Test usage of increment argument in DPCHFD and DPCHFE for
!            DPCHQ1.
!***LIBRARY   SLATEC (PCHIP)
!***TYPE      DOUBLE PRECISION (EVPCCK-S, DEVPCK-D)
!***KEYWORDS  PCHIP EVALUATOR QUICK CHECK
!***AUTHOR  Fritsch, F. N., (LLNL)
!***DESCRIPTION
!
! ---- CODE TO TEST USAGE OF INCREMENT ARGUMENT IN DPCHFD AND DPCHFE ---
!
!     EVALUATES A BICUBIC FUNCTION AND ITS FIRST PARTIAL DERIVATIVES
!     ON A 4X6 MESH CONTAINED IN A 10X10 ARRAY.
!
!     INTERPOLATION OF THESE DATA ALONG MESH LINES IN EITHER DIMENSION
!     SHOULD AGREE WITH CORRECT FUNCTION WITHIN ROUNDOFF ERROR.
!
!     ARRAYS ARE ARGUMENTS ONLY TO ALLOW SHARING STORAGE WITH OTHER
!     TEST ROUTINES.
!
!     NOTE:  RUN WITH KPRINT=4 FOR FULL GORY DETAILS (10 PAGES WORTH).
!
!
!     FORTRAN INTRINSICS USED:  ABS.
!     FORTRAN LIBRARY ROUTINES USED:  (WRITE).
!     SLATEC LIBRARY ROUTINES USED:  DPCHFD, DPCHFE, D1MACH.
!
!***ROUTINES CALLED  D1MACH, DPCHFD, DPCHFE
!***REVISION HISTORY  (YYMMDD)
!   820601  DATE WRITTEN
!   820714  CONVERTED TO QUICK CHECK FOR SLATEC LIBRARY.
!   820715  1. CORRECTED SOME FORMATS.
!           2. ADDED CALL TO D1MACH TO SET MACHEP.
!   890406  1. Modified to make sure final elements of X and XE
!             agree, to avoid possible failure due to roundoff
!             error.
!           2. Added printout of TOL in case of failure.
!           3. Removed unnecessary IMPLICIT declaration.
!           4. Corrected a few S.P. constants to D.P.
!           5. Minor cosmetic changes.
!   890706  Cosmetic changes to prologue.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   891004  Cosmetic changes to prologue.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  Revised prologue and improved some output formats.  (FNF)
!   900316  Additional minor cosmetic changes.  (FNF)
!   900321  Made miscellaneous cosmetic changes.  (FNF)
!   901130  Made many changes to output:  (FNF)
!           1. Reduced amount of output for KPRINT=3.  (Now need to
!              use KPRINT=4 for full output.)
!           2. Added 1P's to formats and revised some to reduce maximum
!              line length.
!   910708  Minor modifications in use of KPRINT.  (WRB)
!   930317  Improved output formats.  (FNF)
!***END PROLOGUE  DEVPCK
!
!  Declare arguments.
!
      INTEGER Lout , Kprint
      LOGICAL Fail
      DOUBLE PRECISION X(10) , Y(10) , F(10,10) , Fx(10,10) , Fy(10,10) , Xe(51)
     &                 , Ye(51) , Fe(51) , De(51) , Fe2(51)
!
!  DECLARATIONS.
!
      INTEGER i , ier2 , ierr , inc , j , k , ne , nerr , nmax , nx , ny
      LOGICAL faild , faile , failoc , skip
      DOUBLE PRECISION dermax , derr , dtrue , dx , fdiff , fdifmx , fermax , 
     &                 ferr , ftrue , machep , tol , pdermx , pdifmx , pfermx , 
     &                 zero
      DOUBLE PRECISION D1MACH
!
!  DEFINE TEST FUNCTION AND DERIVATIVES.
!
      DOUBLE PRECISION ax , ay , FCN , DFDX , DFDY
      FCN(ax,ay) = ax*(ay*ay)*(ax*ax+1.D0)
      DFDX(ax,ay) = (ay*ay)*(3.D0*ax*ax+1.D0)
      DFDY(ax,ay) = 2.D0*ax*ay*(ax*ax+1.D0)
!
      DATA nmax/10/ , nx/4/ , ny/6/
      DATA ne/51/
      DATA zero/0.D0/
!
!  INITIALIZE.
!
!***FIRST EXECUTABLE STATEMENT  DEVPCK
      machep = D1MACH(4)
!       Following tolerance is looser than S.P. version to avoid
!       spurious failures on some systems.
      tol = 25.D0*machep
!
      Fail = .FALSE.
!
!  SET UP 4-BY-6 MESH IN A 10-BY-10 ARRAY:
!     X =  0.25(0.25)1.   ;
!     Y = -0.75(0.5 )1.75 .
!
      DO i = 1 , nx - 1
        X(i) = 0.25D0*i
      ENDDO
      X(nx) = 1.D0
      DO j = 1 , ny
        Y(j) = 0.5D0*j - 1.25D0
        DO i = 1 , nx
          F(i,j) = FCN(X(i),Y(j))
          Fx(i,j) = DFDX(X(i),Y(j))
          Fy(i,j) = DFDY(X(i),Y(j))
        ENDDO
      ENDDO
!
!  SET UP EVALUATION POINTS:
!     XE =  0.(0.02)1. ;
!     YE = -2.(0.08)2. .
!
      dx = 1.D0/(ne-1)
      DO k = 1 , ne - 1
        Xe(k) = dx*(k-1)
        Ye(k) = 4.D0*Xe(k) - 2.D0
      ENDDO
      Xe(ne) = 1.D0
      Ye(ne) = 2.D0
!
      IF ( Kprint>=3 ) WRITE (Lout,99001)
!
!  FORMATS.
!
99001 FORMAT ('1'//10X,'TEST DPCHFE AND DPCHFD')
      IF ( Kprint>=2 ) WRITE (Lout,99002)
99002 FORMAT (//10X,'DEVPCK RESULTS'/10X,'--------------')
!
!  EVALUATE ON HORIZONTAL MESH LINES (Y FIXED, X RUNNING) ..............
!
      nerr = 0
      inc = 1
      skip = .FALSE.
      DO j = 1 , ny
!        --------------------------------------------------------------
        CALL DPCHFD(nx,X,F(1,j),Fx(1,j),inc,skip,ne,Xe,Fe,De,ierr)
!        --------------------------------------------------------------
        IF ( Kprint>=3 ) WRITE (Lout,99003) inc , 'J' , j , 'Y' , Y(j) , ierr
        IF ( ierr<0 ) THEN
!
          failoc = .TRUE.
          IF ( Kprint>=2 ) WRITE (Lout,99011) ierr
        ELSE
          IF ( Kprint>3 ) WRITE (Lout,99004) 'X'
!
!        DPCHFE SHOULD AGREE EXACTLY WITH DPCHFD.
!
!        -----------------------------------------------------------
          CALL DPCHFE(nx,X,F(1,j),Fx(1,j),inc,skip,ne,Xe,Fe2,ier2)
!        -----------------------------------------------------------
!
          DO k = 1 , ne
            ftrue = FCN(Xe(k),Y(j))
            ferr = Fe(k) - ftrue
            dtrue = DFDX(Xe(k),Y(j))
            derr = De(k) - dtrue
            IF ( Kprint>3 ) WRITE (Lout,99005) Xe(k) , ftrue , Fe(k) , ferr , 
     &                             dtrue , De(k) , derr
            IF ( k==1 ) THEN
!              INITIALIZE.
              fermax = ABS(ferr)
              pfermx = Xe(1)
              dermax = ABS(derr)
              pdermx = Xe(1)
              fdifmx = ABS(Fe2(1)-Fe(1))
              pdifmx = Xe(1)
            ELSE
!              SELECT.
              ferr = ABS(ferr)
              IF ( ferr>fermax ) THEN
                fermax = ferr
                pfermx = Xe(k)
              ENDIF
              derr = ABS(derr)
              IF ( derr>dermax ) THEN
                dermax = derr
                pdermx = Xe(k)
              ENDIF
              fdiff = ABS(Fe2(k)-Fe(k))
              IF ( fdiff>fdifmx ) THEN
                fdifmx = fdiff
                pdifmx = Xe(k)
              ENDIF
            ENDIF
          ENDDO
!
          faild = (fermax>tol) .OR. (dermax>tol)
          faile = fdifmx/=zero
          failoc = faild .OR. faile .OR. (ierr/=13) .OR. (ier2/=ierr)
!
          IF ( failoc.AND.(Kprint>=2) ) WRITE (Lout,99006) 'J' , j , 'Y' , Y(j)
!
          IF ( (Kprint>=3).OR.(faild.AND.(Kprint==2)) ) WRITE (Lout,99007)
     &         fermax , pfermx , dermax , pdermx
          IF ( faild.AND.(Kprint>=2) ) WRITE (Lout,99010) tol
!
          IF ( (Kprint>=3).OR.(faile.AND.(Kprint==2)) ) WRITE (Lout,99008)
     &         fdifmx , pdifmx
!
          IF ( (ierr/=13).AND.(Kprint>=2) ) WRITE (Lout,99009) 'D' , ierr , 13
!
          IF ( (ier2/=ierr).AND.(Kprint>=2) ) WRITE (Lout,99009) 'E' , ier2 , 
     &         ierr
        ENDIF
!
        IF ( failoc ) nerr = nerr + 1
        Fail = Fail .OR. failoc
      ENDDO
!
      IF ( Kprint>=2 ) THEN
        IF ( nerr>0 ) THEN
          WRITE (Lout,99012) nerr , 'J'
        ELSE
          WRITE (Lout,99013) 'J'
        ENDIF
      ENDIF
!
!  EVALUATE ON VERTICAL MESH LINES (X FIXED, Y RUNNING) ................
!
      nerr = 0
      inc = nmax
      skip = .FALSE.
      DO i = 1 , nx
!        --------------------------------------------------------------
        CALL DPCHFD(ny,Y,F(i,1),Fy(i,1),inc,skip,ne,Ye,Fe,De,ierr)
!        --------------------------------------------------------------
        IF ( Kprint>=3 ) WRITE (Lout,99003) inc , 'I' , i , 'X' , X(i) , ierr
        IF ( ierr<0 ) THEN
!
          failoc = .TRUE.
          IF ( Kprint>=2 ) WRITE (Lout,99011) ierr
        ELSE
          IF ( Kprint>3 ) WRITE (Lout,99004) 'Y'
!
!        DPCHFE SHOULD AGREE EXACTLY WITH DPCHFD.
!
!        -----------------------------------------------------------
          CALL DPCHFE(ny,Y,F(i,1),Fy(i,1),inc,skip,ne,Ye,Fe2,ier2)
!        -----------------------------------------------------------
!
          DO k = 1 , ne
            ftrue = FCN(X(i),Ye(k))
            ferr = Fe(k) - ftrue
            dtrue = DFDY(X(i),Ye(k))
            derr = De(k) - dtrue
            IF ( Kprint>3 ) WRITE (Lout,99005) Ye(k) , ftrue , Fe(k) , ferr , 
     &                             dtrue , De(k) , derr
            IF ( k==1 ) THEN
!              INITIALIZE.
              fermax = ABS(ferr)
              pfermx = Ye(1)
              dermax = ABS(derr)
              pdermx = Ye(1)
              fdifmx = ABS(Fe2(1)-Fe(1))
              pdifmx = Ye(1)
            ELSE
!              SELECT.
              ferr = ABS(ferr)
              IF ( ferr>fermax ) THEN
                fermax = ferr
                pfermx = Ye(k)
              ENDIF
              derr = ABS(derr)
              IF ( derr>dermax ) THEN
                dermax = derr
                pdermx = Ye(k)
              ENDIF
              fdiff = ABS(Fe2(k)-Fe(k))
              IF ( fdiff>fdifmx ) THEN
                fdifmx = fdiff
                pdifmx = Ye(k)
              ENDIF
            ENDIF
          ENDDO
!
          faild = (fermax>tol) .OR. (dermax>tol)
          faile = fdifmx/=zero
          failoc = faild .OR. faile .OR. (ierr/=20) .OR. (ier2/=ierr)
!
          IF ( failoc.AND.(Kprint>=2) ) WRITE (Lout,99006) 'I' , i , 'X' , X(i)
!
          IF ( (Kprint>=3).OR.(faild.AND.(Kprint==2)) ) WRITE (Lout,99007)
     &         fermax , pfermx , dermax , pdermx
          IF ( faild.AND.(Kprint>=2) ) WRITE (Lout,99010) tol
!
          IF ( (Kprint>=3).OR.(faile.AND.(Kprint==2)) ) WRITE (Lout,99008)
     &         fdifmx , pdifmx
!
          IF ( (ierr/=20).AND.(Kprint>=2) ) WRITE (Lout,99009) 'D' , ierr , 20
!
          IF ( (ier2/=ierr).AND.(Kprint>=2) ) WRITE (Lout,99009) 'E' , ier2 , 
     &         ierr
        ENDIF
!
        IF ( failoc ) nerr = nerr + 1
        Fail = Fail .OR. failoc
      ENDDO
!
      IF ( Kprint>=2 ) THEN
        IF ( nerr>0 ) THEN
          WRITE (Lout,99012) nerr , 'I'
        ELSE
          WRITE (Lout,99013) 'I'
        ENDIF
      ENDIF
!
!  TERMINATE.
!
      RETURN
99003 FORMAT (//20X,'DPCHFD INCREMENT TEST -- INCFD = ',I2/15X,'ON ',A1,
     &        '-LINE ',I2,',  ',A1,' =',F8.4,'  --  IERR =',I3)
99004 FORMAT (/3X,A1,'E',10X,'F',8X,'FE',9X,'DIFF',13X,'D',8X,'DE',9X,'DIFF')
99005 FORMAT (F7.2,2(2X,2F10.5,1P,D15.5,0P))
99006 FORMAT (/' ***** DPCHFD AND/OR DPCHFE FAILED ON ',A1,'-LINE ',I1,',  ',A1,
     &        ' =',F8.4)
99007 FORMAT (/19X,'  MAXIMUM ERROR IN FUNCTION =',1P,1P,D13.5,0P,' (AT',F6.2,
     &        '),'/33X,'IN DERIVATIVE =',1P,D13.5,0P,' (AT',F6.2,').')
99008 FORMAT ('  MAXIMUM DIFFERENCE BETWEEN DPCHFE AND DPCHFD =',1P,D13.5,0P,
     &        ' (AT',F6.2,').')
99009 FORMAT (/'  DPCHF',A1,' RETURNED IERR = ',I2,' INSTEAD OF ',I2)
99010 FORMAT ('  *** BOTH SHOULD BE .LE. TOL =',1P,D12.5,' ***')
99011 FORMAT (//' ***** ERROR ***** DPCHFD RETURNED IERR =',I5//)
99012 FORMAT (//' ***** ERROR ***** DPCHFD AND/OR DPCHFE FAILED ON',I2,1X,A1,
     &        '-LINES.'//)
99013 FORMAT (/' DPCHFD AND DPCHFE OK ON ',A1,'-LINES.')
!------------- LAST LINE OF DEVPCK FOLLOWS -----------------------------
      END SUBROUTINE DEVPCK
