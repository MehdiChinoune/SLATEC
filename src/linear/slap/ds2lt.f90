!** DS2LT
PURE SUBROUTINE DS2LT(N,Nelt,Ia,Ja,A,Isym,Nel,Iel,Jel,El)
  !> Lower Triangle Preconditioner SLAP Set Up.
  !  Routine to store the lower triangle of a matrix stored in the SLAP Column format.
  !***
  ! **Library:**   SLATEC (SLAP)
  !***
  ! **Category:**  D2E
  !***
  ! **Type:**      DOUBLE PRECISION (SS2LT-S, DS2LT-D)
  !***
  ! **Keywords:**  LINEAR SYSTEM, LOWER TRIANGLE, SLAP SPARSE
  !***
  ! **Author:**  Greenbaum, Anne, (Courant Institute)
  !           Seager, Mark K., (LLNL)
  !             Lawrence Livermore National Laboratory
  !             PO BOX 808, L-60
  !             Livermore, CA 94550 (510) 423-3141
  !             seager@llnl.gov
  !***
  ! **Description:**
  !
  !- Usage:
  !     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
  !     INTEGER NEL, IEL(NEL), JEL(NEL)
  !     DOUBLE PRECISION A(NELT), EL(NEL)
  !
  !     CALL DS2LT( N, NELT, IA, JA, A, ISYM, NEL, IEL, JEL, EL )
  !
  !- Arguments:
  ! N      :IN       Integer
  !         Order of the Matrix.
  ! NELT   :IN       Integer.
  !         Number of non-zeros stored in A.
  ! IA     :IN       Integer IA(NELT).
  ! JA     :IN       Integer JA(NELT).
  ! A      :IN       Double Precision A(NELT).
  !         These arrays should hold the matrix A in the SLAP Column
  !         format.  See "Description", below.
  ! ISYM   :IN       Integer.
  !         Flag to indicate symmetric storage format.
  !         If ISYM=0, all non-zero entries of the matrix are stored.
  !         If ISYM=1, the matrix is symmetric, and only the lower
  !         triangle of the matrix is stored.
  ! NEL    :OUT      Integer.
  !         Number of non-zeros in the lower triangle of A.   Also
  !         corresponds to the length of the IEL, JEL, EL arrays.
  ! IEL    :OUT      Integer IEL(NEL).
  ! JEL    :OUT      Integer JEL(NEL).
  ! EL     :OUT      Double Precision     EL(NEL).
  !         IEL, JEL, EL contain the lower triangle of the A matrix
  !         stored in SLAP Column format.  See "Description", below,
  !         for more details bout the SLAP Column format.
  !
  !- Description
  !       =================== S L A P Column format ==================
  !       This routine  requires that  the matrix A  be stored in  the
  !       SLAP Column format.  In this format the non-zeros are stored
  !       counting down columns (except for  the diagonal entry, which
  !       must appear first in each  "column")  and are stored  in the
  !       double precision array A.   In other words,  for each column
  !       in the matrix put the diagonal entry in  A.  Then put in the
  !       other non-zero  elements going down  the column (except  the
  !       diagonal) in order.   The  IA array holds the  row index for
  !       each non-zero.  The JA array holds the offsets  into the IA,
  !       A arrays  for  the  beginning  of each   column.   That  is,
  !       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the
  !       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1),
  !       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column.
  !       Note that we always have  JA(N+1) = NELT+1,  where N is  the
  !       number of columns in  the matrix and NELT  is the number  of
  !       non-zeros in the matrix.
  !
  !       Here is an example of the  SLAP Column  storage format for a
  !       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a
  !       column):
  !
  !           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
  !                              1  2  3    4  5    6  7    8    9 10 11
  !       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
  !       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
  !       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
  !       | 0  0  0 44  0|
  !       |51  0 53  0 55|
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   890404  DATE WRITTEN
  !   890404  Previous REVISION DATE
  !   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
  !   890922  Numerous changes to prologue to make closer to SLATEC standard.  (FNF)
  !   890929  Numerous changes to reduce SP/DP differences.  (FNF)
  !   910411  Prologue converted to Version 4.0 format.  (BAB)
  !   920511  Added complete declaration section.  (WRB)
  !   930701  Updated CATEGORY section.  (FNF, WRB)

  !     .. Scalar Arguments ..
  INTEGER, INTENT(IN) :: Isym, N, Nelt
  INTEGER, INTENT(INOUT) :: Nel
  !     .. Array Arguments ..
  INTEGER, INTENT(IN) :: Ia(Nelt), Ja(Nelt)
  INTEGER, INTENT(OUT) :: Iel(Nel), Jel(Nel)
  REAL(DP), INTENT(IN) :: A(Nelt)
  REAL(DP), INTENT(OUT) :: El(Nel)
  !     .. Local Scalars ..
  INTEGER :: i, icol, j, jbgn, jend
  !* FIRST EXECUTABLE STATEMENT  DS2LT
  IF( Isym==0 ) THEN
    !
    !         The matrix is stored non-symmetricly.  Pick out the lower
    !         triangle.
    !
    Nel = 0
    DO icol = 1, N
      Jel(icol) = Nel + 1
      jbgn = Ja(icol)
      jend = Ja(icol+1) - 1
      DO j = jbgn, jend
        IF( Ia(j)>=icol ) THEN
          Nel = Nel + 1
          Iel(Nel) = Ia(j)
          El(Nel) = A(j)
        END IF
      END DO
    END DO
    Jel(N+1) = Nel + 1
  ELSE
    !
    !         The matrix is symmetric and only the lower triangle is
    !         stored.  Copy it to IEL, JEL, EL.
    !
    Nel = Nelt
    DO i = 1, Nelt
      Iel(i) = Ia(i)
      El(i) = A(i)
    END DO
    DO i = 1, N + 1
      Jel(i) = Ja(i)
    END DO
  END IF
  !------------- LAST LINE OF DS2LT FOLLOWS ----------------------------
END SUBROUTINE DS2LT