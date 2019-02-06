!*==SS2LT.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK SS2LT
      SUBROUTINE SS2LT(N,Nelt,Ia,Ja,A,Isym,Nel,Iel,Jel,El)
      IMPLICIT NONE
!*--SS2LT5
!***BEGIN PROLOGUE  SS2LT
!***PURPOSE  Lower Triangle Preconditioner SLAP Set Up.
!            Routine to store the lower triangle of a matrix stored
!            in the SLAP Column format.
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  D2E
!***TYPE      SINGLE PRECISION (SS2LT-S, DS2LT-D)
!***KEYWORDS  LINEAR SYSTEM, LOWER TRIANGLE, SLAP SPARSE
!***AUTHOR  Greenbaum, Anne, (Courant Institute)
!           Seager, Mark K., (LLNL)
!             Lawrence Livermore National Laboratory
!             PO BOX 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!             seager@llnl.gov
!***DESCRIPTION
!
! *Usage:
!     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
!     INTEGER NEL, IEL(NEL), JEL(NEL)
!     REAL    A(NELT), EL(NEL)
!
!     CALL SS2LT( N, NELT, IA, JA, A, ISYM, NEL, IEL, JEL, EL )
!
! *Arguments:
! N      :IN       Integer
!         Order of the Matrix.
! NELT   :IN       Integer.
!         Number of non-zeros stored in A.
! IA     :IN       Integer IA(NELT).
! JA     :IN       Integer JA(NELT).
! A      :IN       Real A(NELT).
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
! EL     :OUT      Real     EL(NEL).
!         IEL, JEL, EL contain the lower triangle of the A matrix
!         stored in SLAP Column format.  See "Description", below,
!         for more details bout the SLAP Column format.
!
! *Description
!       =================== S L A P Column format ==================
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
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   871119  DATE WRITTEN
!   881213  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   920511  Added complete declaration section.  (WRB)
!   930701  Updated CATEGORY section.  (FNF, WRB)
!***END PROLOGUE  SS2LT
!     .. Scalar Arguments ..
      INTEGER Isym , N , Nel , Nelt
!     .. Array Arguments ..
      REAL A(Nelt) , El(Nelt)
      INTEGER Ia(Nelt) , Iel(Nel) , Ja(Nelt) , Jel(Nel)
!     .. Local Scalars ..
      INTEGER i , icol , j , jbgn , jend
!***FIRST EXECUTABLE STATEMENT  SS2LT
      IF ( Isym==0 ) THEN
!
!         The matrix is stored non-symmetricly.  Pick out the lower
!         triangle.
!
        Nel = 0
        DO icol = 1 , N
          Jel(icol) = Nel + 1
          jbgn = Ja(icol)
          jend = Ja(icol+1) - 1
!VD$ NOVECTOR
          DO j = jbgn , jend
            IF ( Ia(j)>=icol ) THEN
              Nel = Nel + 1
              Iel(Nel) = Ia(j)
              El(Nel) = A(j)
            ENDIF
          ENDDO
        ENDDO
        Jel(N+1) = Nel + 1
      ELSE
!
!         The matrix is symmetric and only the lower triangle is
!         stored.  Copy it to IEL, JEL, EL.
!
        Nel = Nelt
        DO i = 1 , Nelt
          Iel(i) = Ia(i)
          El(i) = A(i)
        ENDDO
        DO i = 1 , N + 1
          Jel(i) = Ja(i)
        ENDDO
      ENDIF
!------------- LAST LINE OF SS2LT FOLLOWS ----------------------------
      END SUBROUTINE SS2LT
