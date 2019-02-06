!*==SPINIT.f90  processed by SPAG 6.72Dc at 11:02 on  6 Feb 2019
!DECK SPINIT
      SUBROUTINE SPINIT(Mrelas,Nvars,Costs,Bl,Bu,Ind,Primal,Info,Amat,Csc,
     &                  Costsc,Colnrm,Xlamda,Anorm,Rhs,Rhsnrm,Ibasis,Ibb,Imat,
     &                  Lopt)
      IMPLICIT NONE
!*--SPINIT7
!*** Start of declarations inserted by SPAG
      INTEGER i , Info , ip , iplace , j , Mrelas , n20007 , n20019 , n20028 , 
     &        n20041 , n20056 , n20066 , n20070 , n20074 , n20078 , Nvars
      REAL SASUM
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  SPINIT
!***SUBSIDIARY
!***PURPOSE  Subsidiary to SPLP
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (SPINIT-S, DPINIT-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO
!     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES.
!
!     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/.
!     /REAL (12 BLANKS)/DOUBLE PRECISION/,/SCOPY/DCOPY/
!     REVISED 810519-0900
!     REVISED YYMMDD-HHMM
!
!     INITIALIZATION SUBROUTINE FOR SPLP(*) PACKAGE.
!
!***SEE ALSO  SPLP
!***ROUTINES CALLED  PNNZRS, SASUM, SCOPY
!***REVISION HISTORY  (YYMMDD)
!   811215  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890605  Removed unreferenced labels.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  SPINIT
      REAL aij , Amat(*) , Anorm , Bl(*) , Bu(*) , cmax , Colnrm(*) , Costs(*) , 
     &     Costsc , Csc(*) , csum , one , Primal(*) , Rhs(*) , Rhsnrm , scalr , 
     &     testsc , Xlamda , zero
      INTEGER Ibasis(*) , Ibb(*) , Imat(*) , Ind(*)
      LOGICAL contin , usrbas , colscp , cstscp , minprb , Lopt(8)
!
!***FIRST EXECUTABLE STATEMENT  SPINIT
      zero = 0.
      one = 1.
      contin = Lopt(1)
      usrbas = Lopt(2)
      colscp = Lopt(5)
      cstscp = Lopt(6)
!
!     SCALE DATA. NORMALIZE BOUNDS. FORM COLUMN CHECK SUMS.
!
!     INITIALIZE ACTIVE BASIS MATRIX.
      minprb = Lopt(7)
!
!     PROCEDURE (SCALE DATA. NORMALIZE BOUNDS. FORM COLUMN CHECK SUMS)
!
!     DO COLUMN SCALING IF NOT PROVIDED BY THE USER.
      IF ( .NOT.colscp ) THEN
        j = 1
        n20007 = Nvars
        DO WHILE ( (n20007-j)>=0 )
          cmax = zero
          i = 0
          DO
            CALL PNNZRS(i,aij,iplace,Amat,Imat,j)
            IF ( i/=0 ) THEN
              cmax = MAX(cmax,ABS(aij))
            ELSE
              IF ( cmax/=zero ) THEN
                Csc(j) = one/cmax
              ELSE
                Csc(j) = one
              ENDIF
              j = j + 1
              EXIT
            ENDIF
          ENDDO
        ENDDO
      ENDIF
!
!     FORM CHECK SUMS OF COLUMNS. COMPUTE MATRIX NORM OF SCALED MATRIX.
      Anorm = zero
      j = 1
      n20019 = Nvars
      DO WHILE ( (n20019-j)>=0 )
        Primal(j) = zero
        csum = zero
        i = 0
        DO
          CALL PNNZRS(i,aij,iplace,Amat,Imat,j)
          IF ( i>0 ) THEN
            Primal(j) = Primal(j) + aij
            csum = csum + ABS(aij)
          ELSE
            IF ( Ind(j)==2 ) Csc(j) = -Csc(j)
            Primal(j) = Primal(j)*Csc(j)
            Colnrm(j) = ABS(Csc(j)*csum)
            Anorm = MAX(Anorm,Colnrm(j))
            j = j + 1
            EXIT
          ENDIF
        ENDDO
      ENDDO
!
!     IF THE USER HAS NOT PROVIDED COST VECTOR SCALING THEN SCALE IT
!     USING THE MAX. NORM OF THE TRANSFORMED COST VECTOR, IF NONZERO.
      testsc = zero
      j = 1
      n20028 = Nvars
      DO WHILE ( (n20028-j)>=0 )
        testsc = MAX(testsc,ABS(Csc(j)*Costs(j)))
        j = j + 1
      ENDDO
      IF ( .NOT.cstscp ) THEN
        IF ( testsc<=zero ) THEN
          Costsc = one
        ELSE
          Costsc = one/testsc
        ENDIF
      ENDIF
      Xlamda = (Costsc+Costsc)*testsc
      IF ( Xlamda==zero ) Xlamda = one
!
!     IF MAXIMIZATION PROBLEM, THEN CHANGE SIGN OF COSTSC AND LAMDA
!     =WEIGHT FOR PENALTY-FEASIBILITY METHOD.
      IF ( .NOT.minprb ) Costsc = -Costsc
!:CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     PROCEDURE (INITIALIZE RHS(*),IBASIS(*), AND IBB(*))
!
!     INITIALLY SET RIGHT-HAND SIDE VECTOR TO ZERO.
      CALL SCOPY(Mrelas,zero,0,Rhs,1)
!
!     TRANSLATE RHS ACCORDING TO CLASSIFICATION OF INDEPENDENT VARIABLES
      j = 1
      n20041 = Nvars
      DO WHILE ( (n20041-j)>=0 )
        IF ( Ind(j)==1 ) THEN
          scalr = -Bl(j)
        ELSEIF ( Ind(j)==2 ) THEN
          scalr = -Bu(j)
        ELSEIF ( Ind(j)==3 ) THEN
          scalr = -Bl(j)
        ELSEIF ( Ind(j)==4 ) THEN
          scalr = zero
        ENDIF
        IF ( scalr==zero ) THEN
          j = j + 1
        ELSE
          i = 0
          DO
            CALL PNNZRS(i,aij,iplace,Amat,Imat,j)
            IF ( i>0 ) THEN
              Rhs(i) = scalr*aij + Rhs(i)
            ELSE
              j = j + 1
              EXIT
            ENDIF
          ENDDO
        ENDIF
      ENDDO
!
!     TRANSLATE RHS ACCORDING TO CLASSIFICATION OF DEPENDENT VARIABLES.
      i = Nvars + 1
      n20056 = Nvars + Mrelas
      DO WHILE ( (n20056-i)>=0 )
        IF ( Ind(i)==1 ) THEN
          scalr = Bl(i)
        ELSEIF ( Ind(i)==2 ) THEN
          scalr = Bu(i)
        ELSEIF ( Ind(i)==3 ) THEN
          scalr = Bl(i)
        ELSEIF ( Ind(i)==4 ) THEN
          scalr = zero
        ENDIF
        Rhs(i-Nvars) = Rhs(i-Nvars) + scalr
        i = i + 1
      ENDDO
      Rhsnrm = SASUM(Mrelas,Rhs,1)
!
!     IF THIS IS NOT A CONTINUATION OR THE USER HAS NOT PROVIDED THE
!     INITIAL BASIS, THEN THE INITIAL BASIS IS COMPRISED OF THE
!     DEPENDENT VARIABLES.
      IF ( .NOT.(contin.OR.usrbas) ) THEN
        j = 1
        n20066 = Mrelas
        DO WHILE ( (n20066-j)>=0 )
          Ibasis(j) = Nvars + j
          j = j + 1
        ENDDO
      ENDIF
!
!     DEFINE THE ARRAY IBB(*)
      j = 1
      n20070 = Nvars + Mrelas
      DO WHILE ( (n20070-j)>=0 )
        Ibb(j) = 1
        j = j + 1
      ENDDO
      j = 1
      n20074 = Mrelas
      DO WHILE ( (n20074-j)>=0 )
        Ibb(Ibasis(j)) = -1
        j = j + 1
      ENDDO
!
!     DEFINE THE REST OF IBASIS(*)
      ip = Mrelas
      j = 1
      n20078 = Nvars + Mrelas
      DO WHILE ( (n20078-j)>=0 )
        IF ( Ibb(j)>0 ) THEN
          ip = ip + 1
          Ibasis(ip) = j
        ENDIF
        j = j + 1
      ENDDO
      RETURN
      END SUBROUTINE SPINIT
