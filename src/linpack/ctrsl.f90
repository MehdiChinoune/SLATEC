!DECK CTRSL
SUBROUTINE CTRSL(T,Ldt,N,B,Job,Info)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CTRSL
  !***PURPOSE  Solve a system of the form  T*X=B or CTRANS(T)*X=B, where
  !            T is a triangular matrix.  Here CTRANS(T) is the conjugate
  !            transpose.
  !***LIBRARY   SLATEC (LINPACK)
  !***CATEGORY  D2C3
  !***TYPE      COMPLEX (STRSL-S, DTRSL-D, CTRSL-C)
  !***KEYWORDS  LINEAR ALGEBRA, LINPACK, TRIANGULAR LINEAR SYSTEM,
  !             TRIANGULAR MATRIX
  !***AUTHOR  Stewart, G. W., (U. of Maryland)
  !***DESCRIPTION
  !
  !     CTRSL solves systems of the form
  !
  !                   T * X = B
  !     or
  !                   CTRANS(T) * X = B
  !
  !     where T is a triangular matrix of order N.  Here CTRANS(T)
  !     denotes the conjugate transpose of the matrix T.
  !
  !     On Entry
  !
  !         T         COMPLEX(LDT,N)
  !                   T contains the matrix of the system.  The zero
  !                   elements of the matrix are not referenced, and
  !                   the corresponding elements of the array can be
  !                   used to store other information.
  !
  !         LDT       INTEGER
  !                   LDT is the leading dimension of the array T.
  !
  !         N         INTEGER
  !                   N is the order of the system.
  !
  !         B         COMPLEX(N).
  !                   B contains the right hand side of the system.
  !
  !         JOB       INTEGER
  !                   JOB specifies what kind of system is to be solved.
  !                   If JOB is
  !
  !                        00   solve T*X = B, T lower triangular,
  !                        01   solve T*X = B, T upper triangular,
  !                        10   solve CTRANS(T)*X = B, T lower triangular,
  !                        11   solve CTRANS(T)*X = B, T upper triangular.
  !
  !     On Return
  !
  !         B         B contains the solution, if INFO .EQ. 0.
  !                   Otherwise B is unaltered.
  !
  !         INFO      INTEGER
  !                   INFO contains zero if the system is nonsingular.
  !                   Otherwise INFO contains the index of
  !                   the first zero diagonal element of T.
  !
  !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***ROUTINES CALLED  CAXPY, CDOTC
  !***REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  CTRSL
  INTEGER Ldt, N, Job, Info
  COMPLEX T(Ldt,*), B(*)
  !
  !
  COMPLEX CDOTC, temp
  INTEGER case, j, jj
  REAL, EXTERNAL :: CABS1
  !***FIRST EXECUTABLE STATEMENT  CTRSL
  !
  !        CHECK FOR ZERO DIAGONAL ELEMENTS.
  !
  DO Info = 1, N
    IF ( CABS1(T(Info,Info))==0.0E0 ) RETURN
  ENDDO
  Info = 0
  !
  !        DETERMINE THE TASK AND GO TO IT.
  !
  case = 1
  IF ( MOD(Job,10)/=0 ) case = 2
  IF ( MOD(Job,100)/10/=0 ) case = case + 2
  SELECT CASE (case)
    CASE (2)
      !
      !        SOLVE T*X=B FOR T UPPER TRIANGULAR.
      !
      B(N) = B(N)/T(N,N)
      IF ( N>=2 ) THEN
        DO jj = 2, N
          j = N - jj + 1
          temp = -B(j+1)
          CALL CAXPY(j,temp,T(1,j+1),1,B(1),1)
          B(j) = B(j)/T(j,j)
        ENDDO
      ENDIF
    CASE (3)
      !
      !        SOLVE CTRANS(T)*X=B FOR T LOWER TRIANGULAR.
      !
      B(N) = B(N)/CONJG(T(N,N))
      IF ( N>=2 ) THEN
        DO jj = 2, N
          j = N - jj + 1
          B(j) = B(j) - CDOTC(jj-1,T(j+1,j),1,B(j+1),1)
          B(j) = B(j)/CONJG(T(j,j))
        ENDDO
      ENDIF
    CASE (4)
      !
      !        SOLVE CTRANS(T)*X=B FOR T UPPER TRIANGULAR.
      !
      B(1) = B(1)/CONJG(T(1,1))
      IF ( N>=2 ) THEN
        DO j = 2, N
          B(j) = B(j) - CDOTC(j-1,T(1,j),1,B(1),1)
          B(j) = B(j)/CONJG(T(j,j))
        ENDDO
      ENDIF
    CASE DEFAULT
      !
      !        SOLVE T*X=B FOR T LOWER TRIANGULAR
      !
      B(1) = B(1)/T(1,1)
      IF ( N>=2 ) THEN
        DO j = 2, N
          temp = -B(j-1)
          CALL CAXPY(N-j+1,temp,T(j,j-1),1,B(j),1)
          B(j) = B(j)/T(j,j)
        ENDDO
      ENDIF
  END SELECT
  RETURN
END SUBROUTINE CTRSL
