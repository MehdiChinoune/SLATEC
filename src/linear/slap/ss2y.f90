!** SS2Y
SUBROUTINE SS2Y(N,Nelt,Ia,Ja,A,Isym)
  IMPLICIT NONE
  !>
  !***
  !  SLAP Triad to SLAP Column Format Converter.
  !            Routine to convert from the SLAP Triad to SLAP Column
  !            format.
  !***
  ! **Library:**   SLATEC (SLAP)
  !***
  ! **Category:**  D1B9
  !***
  ! **Type:**      SINGLE PRECISION (SS2Y-S, DS2Y-D)
  !***
  ! **Keywords:**  LINEAR SYSTEM, SLAP SPARSE
  !***
  ! **Author:**  Seager, Mark K., (LLNL)
  !             Lawrence Livermore National Laboratory
  !             PO BOX 808, L-60
  !             Livermore, CA 94550 (510) 423-3141
  !             seager@llnl.gov
  !***
  ! **Description:**
  !
  !- Usage:
  !     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
  !     REAL    A(NELT)
  !
  !     CALL SS2Y( N, NELT, IA, JA, A, ISYM )
  !
  !- Arguments:
  ! N      :IN       Integer
  !         Order of the Matrix.
  ! NELT   :IN       Integer.
  !         Number of non-zeros stored in A.
  ! IA     :INOUT    Integer IA(NELT).
  ! JA     :INOUT    Integer JA(NELT).
  ! A      :INOUT    Real A(NELT).
  !         These arrays should hold the matrix A in either the SLAP
  !         Triad format or the SLAP Column format.  See "Description",
  !         below.  If the SLAP Triad format is used, this format is
  !         translated to the SLAP Column format by this routine.
  ! ISYM   :IN       Integer.
  !         Flag to indicate symmetric storage format.
  !         If ISYM=0, all non-zero entries of the matrix are stored.
  !         If ISYM=1, the matrix is symmetric, and only the lower
  !         triangle of the matrix is stored.
  !
  !- Description:
  !       The Sparse Linear Algebra Package (SLAP) utilizes two matrix
  !       data structures: 1) the  SLAP Triad  format or  2)  the SLAP
  !       Column format.  The user can hand this routine either of the
  !       of these data structures.  If the SLAP Triad format is give
  !       as input then this routine transforms it into SLAP Column
  !       format.  The way this routine tells which format is given as
  !       input is to look at JA(N+1).  If JA(N+1) = NELT+1 then we
  !       have the SLAP Column format.  If that equality does not hold
  !       then it is assumed that the IA, JA, A arrays contain the
  !       SLAP Triad format.
  !
  !       =================== S L A P Triad format ===================
  !       This routine requires that the  matrix A be   stored in  the
  !       SLAP  Triad format.  In  this format only the non-zeros  are
  !       stored.  They may appear in  *ANY* order.  The user supplies
  !       three arrays of  length NELT, where  NELT is  the number  of
  !       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For
  !       each non-zero the user puts the row and column index of that
  !       matrix element  in the IA and  JA arrays.  The  value of the
  !       non-zero   matrix  element is  placed  in  the corresponding
  !       location of the A array.   This is  an  extremely  easy data
  !       structure to generate.  On  the  other hand it   is  not too
  !       efficient on vector computers for  the iterative solution of
  !       linear systems.  Hence,   SLAP changes   this  input    data
  !       structure to the SLAP Column format  for  the iteration (but
  !       does not change it back).
  !
  !       Here is an example of the  SLAP Triad   storage format for a
  !       5x5 Matrix.  Recall that the entries may appear in any order.
  !
  !           5x5 Matrix      SLAP Triad format for 5x5 matrix on left.
  !                              1  2  3  4  5  6  7  8  9 10 11
  !       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
  !       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
  !       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
  !       | 0  0  0 44  0|
  !       |51  0 53  0 55|
  !
  !       =================== S L A P Column format ==================
  !
  !       This routine  requires that  the matrix A  be stored in  the
  !       SLAP Column format.  In this format the non-zeros are stored
  !       counting down columns (except for  the diagonal entry, which
  !       must appear first in each  "column")  and are stored  in the
  !       real array A.  In other words, for each column in the matrix
  !       put the diagonal entry in A.  Then put in the other non-zero
  !       elements going down   the  column (except  the diagonal)  in
  !       order.  The IA array holds the row  index for each non-zero.
  !       The JA array holds the offsets into the IA, A arrays for the
  !       beginning of   each    column.    That  is,    IA(JA(ICOL)),
  !       A(JA(ICOL)) points to the beginning of the ICOL-th column in
  !       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the
  !       end  of   the ICOL-th  column.  Note   that  we  always have
  !       JA(N+1) = NELT+1, where  N  is the number of columns in  the
  !       matrix and  NELT   is the number of non-zeros in the matrix.
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
  ! **Routines called:**  QS2I1R

  !* REVISION HISTORY  (YYMMDD)
  !   871119  DATE WRITTEN
  !   881213  Previous REVISION DATE
  !   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
  !   890922  Numerous changes to prologue to make closer to SLATEC
  !           standard.  (FNF)
  !   890929  Numerous changes to reduce SP/DP differences.  (FNF)
  !   910411  Prologue converted to Version 4.0 format.  (BAB)
  !   910502  Corrected C***FIRST EXECUTABLE STATEMENT line.  (FNF)
  !   920511  Added complete declaration section.  (WRB)
  !   930701  Updated CATEGORY section.  (FNF, WRB)

  !     .. Scalar Arguments ..
  INTEGER Isym, N, Nelt
  !     .. Array Arguments ..
  REAL A(Nelt)
  INTEGER Ia(Nelt), Ja(Nelt)
  !     .. Local Scalars ..
  REAL temp
  INTEGER i, ibgn, icol, iend, itemp, j
  !     .. External Subroutines ..
  EXTERNAL :: QS2I1R
  !* FIRST EXECUTABLE STATEMENT  SS2Y
  !
  !         Check to see if the (IA,JA,A) arrays are in SLAP Column
  !         format.  If it's not then transform from SLAP Triad.
  !
  IF ( Ja(N+1)==Nelt+1 ) RETURN
  !
  !         Sort into ascending order by COLUMN (on the ja array).
  !         This will line up the columns.
  !
  CALL QS2I1R(Ja,Ia,A,Nelt,1)
  !
  !         Loop over each column to see where the column indices change
  !         in the column index array ja.  This marks the beginning of the
  !         next column.
  !
  !VD$R NOVECTOR
  Ja(1) = 1
  DO icol = 1, N - 1
    DO j = Ja(icol) + 1, Nelt
      IF ( Ja(j)/=icol ) THEN
        Ja(icol+1) = j
        EXIT
      END IF
    END DO
  END DO
  Ja(N+1) = Nelt + 1
  !
  !         Mark the n+2 element so that future calls to a SLAP routine
  !         utilizing the YSMP-Column storage format will be able to tell.
  !
  Ja(N+2) = 0
  !
  !         Now loop through the IA array making sure that the diagonal
  !         matrix element appears first in the column.  Then sort the
  !         rest of the column in ascending order.
  !
  DO icol = 1, N
    ibgn = Ja(icol)
    iend = Ja(icol+1) - 1
    DO i = ibgn, iend
      IF ( Ia(i)==icol ) THEN
        !
        !              Swap the diagonal element with the first element in the
        !              column.
        !
        itemp = Ia(i)
        Ia(i) = Ia(ibgn)
        Ia(ibgn) = itemp
        temp = A(i)
        A(i) = A(ibgn)
        A(ibgn) = temp
        EXIT
      END IF
    END DO
    ibgn = ibgn + 1
    IF ( ibgn<iend ) THEN
      DO i = ibgn, iend
        DO j = i + 1, iend
          IF ( Ia(i)>Ia(j) ) THEN
            itemp = Ia(i)
            Ia(i) = Ia(j)
            Ia(j) = itemp
            temp = A(i)
            A(i) = A(j)
            A(j) = temp
          END IF
        END DO
      END DO
    END IF
  END DO
  !------------- LAST LINE OF SS2Y FOLLOWS ----------------------------
END SUBROUTINE SS2Y
