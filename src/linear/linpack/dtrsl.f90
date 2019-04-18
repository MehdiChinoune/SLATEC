!** DTRSL
SUBROUTINE DTRSL(T,Ldt,N,B,Job,Info)
  !>
  !  Solve a system of the form  T*X=B or TRANS(T)*X=B, where
  !            T is a triangular matrix.
  !***
  ! **Library:**   SLATEC (LINPACK)
  !***
  ! **Category:**  D2A3
  !***
  ! **Type:**      DOUBLE PRECISION (STRSL-S, DTRSL-D, CTRSL-C)
  !***
  ! **Keywords:**  LINEAR ALGEBRA, LINPACK, TRIANGULAR LINEAR SYSTEM,
  !             TRIANGULAR MATRIX
  !***
  ! **Author:**  Stewart, G. W., (U. of Maryland)
  !***
  ! **Description:**
  !
  !     DTRSL solves systems of the form
  !
  !                   T * X = B
  !     or
  !                   TRANS(T) * X = B
  !
  !     where T is a triangular matrix of order N.  Here TRANS(T)
  !     denotes the transpose of the matrix T.
  !
  !     On Entry
  !
  !         T         DOUBLE PRECISION(LDT,N)
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
  !         B         DOUBLE PRECISION(N).
  !                   B contains the right hand side of the system.
  !
  !         JOB       INTEGER
  !                   JOB specifies what kind of system is to be solved.
  !                   If JOB is
  !
  !                        00   solve T*X=B, T lower triangular,
  !                        01   solve T*X=B, T upper triangular,
  !                        10   solve TRANS(T)*X=B, T lower triangular,
  !                        11   solve TRANS(T)*X=B, T upper triangular.
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
  !***
  ! **References:**  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***
  ! **Routines called:**  DAXPY, DDOT

  !* REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER Ldt, N, Job, Info
  REAL(8) :: T(Ldt,*), B(*)
  !
  !
  REAL(8) :: temp
  INTEGER case, j, jj
  !* FIRST EXECUTABLE STATEMENT  DTRSL
  !
  !        CHECK FOR ZERO DIAGONAL ELEMENTS.
  !
  DO Info = 1, N
    IF ( T(Info,Info)==0.0D0 ) RETURN
  END DO
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
          CALL DAXPY(j,temp,T(1,j+1),1,B(1),1)
          B(j) = B(j)/T(j,j)
        END DO
      END IF
    CASE (3)
      !
      !        SOLVE TRANS(T)*X=B FOR T LOWER TRIANGULAR.
      !
      B(N) = B(N)/T(N,N)
      IF ( N>=2 ) THEN
        DO jj = 2, N
          j = N - jj + 1
          B(j) = B(j) - DDOT(jj-1,T(j+1,j),1,B(j+1),1)
          B(j) = B(j)/T(j,j)
        END DO
      END IF
    CASE (4)
      !
      !        SOLVE TRANS(T)*X=B FOR T UPPER TRIANGULAR.
      !
      B(1) = B(1)/T(1,1)
      IF ( N>=2 ) THEN
        DO j = 2, N
          B(j) = B(j) - DDOT(j-1,T(1,j),1,B(1),1)
          B(j) = B(j)/T(j,j)
        END DO
      END IF
    CASE DEFAULT
      !
      !        SOLVE T*X=B FOR T LOWER TRIANGULAR
      !
      B(1) = B(1)/T(1,1)
      IF ( N>=2 ) THEN
        DO j = 2, N
          temp = -B(j-1)
          CALL DAXPY(N-j+1,temp,T(j,j-1),1,B(j),1)
          B(j) = B(j)/T(j,j)
        END DO
      END IF
  END SELECT
  RETURN
END SUBROUTINE DTRSL
