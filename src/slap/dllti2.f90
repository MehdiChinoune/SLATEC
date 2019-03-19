!** DLLTI2
SUBROUTINE DLLTI2(N,B,X,Nel,Iel,Jel,El,Dinv)
  IMPLICIT NONE
  !>
  !***
  !  SLAP Backsolve routine for LDL' Factorization.
  !            Routine to solve a system of the form  L*D*L' X = B,
  !            where L is a unit lower triangular matrix and D is a
  !            diagonal matrix and ' means transpose.
  !***
  ! **Library:**   SLATEC (SLAP)
  !***
  ! **Category:**  D2E
  !***
  ! **Type:**      DOUBLE PRECISION (SLLTI2-S, DLLTI2-D)
  !***
  ! **Keywords:**  INCOMPLETE FACTORIZATION, ITERATIVE PRECONDITION, SLAP,
  !             SPARSE, SYMMETRIC LINEAR SYSTEM SOLVE
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
  !     INTEGER N, NEL, IEL(NEL), JEL(NEL)
  !     DOUBLE PRECISION B(N), X(N), EL(NEL), DINV(N)
  !
  !     CALL DLLTI2( N, B, X, NEL, IEL, JEL, EL, DINV )
  !
  !- Arguments:
  ! N      :IN       Integer
  !         Order of the Matrix.
  ! B      :IN       Double Precision B(N).
  !         Right hand side vector.
  ! X      :OUT      Double Precision X(N).
  !         Solution to L*D*L' x = b.
  ! NEL    :IN       Integer.
  !         Number of non-zeros in the EL array.
  ! IEL    :IN       Integer IEL(NEL).
  ! JEL    :IN       Integer JEL(NEL).
  ! EL     :IN       Double Precision     EL(NEL).
  !         IEL, JEL, EL contain the unit lower triangular factor   of
  !         the incomplete decomposition   of the A  matrix  stored in
  !         SLAP Row format.   The diagonal of ones *IS* stored.  This
  !         structure can be set  up  by  the DS2LT routine.  See  the
  !         "Description", below for more details about the  SLAP  Row
  !         format.
  ! DINV   :IN       Double Precision DINV(N).
  !         Inverse of the diagonal matrix D.
  !
  !- Description:
  !       This routine is supplied with  the SLAP package as a routine
  !       to perform the MSOLVE operation in the SCG iteration routine
  !       for  the driver  routine DSICCG.   It must be called via the
  !       SLAP  MSOLVE calling sequence  convention  interface routine
  !       DSLLI.
  !         **** THIS ROUTINE ITSELF DOES NOT CONFORM TO THE ****
  !               **** SLAP MSOLVE CALLING CONVENTION ****
  !
  !       IEL, JEL, EL should contain the unit lower triangular factor
  !       of  the incomplete decomposition of  the A matrix  stored in
  !       SLAP Row format.   This IC factorization  can be computed by
  !       the  DSICS routine.  The  diagonal  (which is all one's) is
  !       stored.
  !
  !       ==================== S L A P Row format ====================
  !
  !       This routine requires  that the matrix A  be  stored  in the
  !       SLAP  Row format.   In this format  the non-zeros are stored
  !       counting across  rows (except for the diagonal  entry, which
  !       must  appear first  in each  "row")  and  are stored  in the
  !       double precision  array A.  In other words, for each row  in
  !       the matrix  put the diagonal  entry in A.   Then put in  the
  !       other  non-zero elements  going across  the row  (except the
  !       diagonal) in order.  The JA array holds the column index for
  !       each non-zero.  The IA array holds the offsets  into the JA,
  !       A  arrays  for  the   beginning  of  each  row.    That  is,
  !       JA(IA(IROW)),A(IA(IROW)) are the first elements of the IROW-
  !       th row in  JA and A,  and  JA(IA(IROW+1)-1), A(IA(IROW+1)-1)
  !       are  the last elements  of the  IROW-th row.   Note  that we
  !       always have  IA(N+1) = NELT+1, where N is the number of rows
  !       in the matrix  and  NELT is the  number of non-zeros  in the
  !       matrix.
  !
  !       Here is an example of the SLAP Row storage format for a  5x5
  !       Matrix (in the A and JA arrays '|' denotes the end of a row):
  !
  !           5x5 Matrix         SLAP Row format for 5x5 matrix on left.
  !                              1  2  3    4  5    6  7    8    9 10 11
  !       |11 12  0  0 15|   A: 11 12 15 | 22 21 | 33 35 | 44 | 55 51 53
  !       |21 22  0  0  0|  JA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
  !       | 0  0 33  0 35|  IA:  1  4  6    8  9   12
  !       | 0  0  0 44  0|
  !       |51  0 53  0 55|
  !
  !       With  the SLAP  Row format  the "inner loop" of this routine
  !       should vectorize   on machines with   hardware  support  for
  !       vector gather/scatter operations.  Your compiler may require
  !       a  compiler directive  to  convince   it that there  are  no
  !       implicit vector  dependencies.  Compiler directives  for the
  !       Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied
  !       with the standard SLAP distribution.
  !
  !***
  ! **See also:**  DSICCG, DSICS
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   871119  DATE WRITTEN
  !   881213  Previous REVISION DATE
  !   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
  !   890922  Numerous changes to prologue to make closer to SLATEC
  !           standard.  (FNF)
  !   890929  Numerous changes to reduce SP/DP differences.  (FNF)
  !   910411  Prologue converted to Version 4.0 format.  (BAB)
  !   920511  Added complete declaration section.  (WRB)
  !   921113  Corrected C***CATEGORY line.  (FNF)
  !   930701  Updated CATEGORY section.  (FNF, WRB)
  
  !     .. Scalar Arguments ..
  INTEGER N, Nel
  !     .. Array Arguments ..
  REAL(8) :: B(N), Dinv(N), El(Nel), X(N)
  INTEGER Iel(Nel), Jel(Nel)
  !     .. Local Scalars ..
  INTEGER i, ibgn, iend, irow
  !* FIRST EXECUTABLE STATEMENT  DLLTI2
  !
  !         Solve  L*y = b,  storing result in x.
  !
  DO i = 1, N
    X(i) = B(i)
  ENDDO
  DO irow = 1, N
    ibgn = Iel(irow) + 1
    iend = Iel(irow+1) - 1
    IF ( ibgn<=iend ) THEN
      !LLL. OPTION ASSERT (NOHAZARD)
      !DIR$ IVDEP
      !VD$ NOCONCUR
      !VD$ NODEPCHK
      DO i = ibgn, iend
        X(irow) = X(irow) - El(i)*X(Jel(i))
      ENDDO
    ENDIF
  ENDDO
  !
  !         Solve  D*Z = Y,  storing result in X.
  !
  DO i = 1, N
    X(i) = X(i)*Dinv(i)
  ENDDO
  !
  !         Solve  L-trans*X = Z.
  !
  DO irow = N, 2, -1
    ibgn = Iel(irow) + 1
    iend = Iel(irow+1) - 1
    IF ( ibgn<=iend ) THEN
      !LLL. OPTION ASSERT (NOHAZARD)
      !DIR$ IVDEP
      !VD$ NOCONCUR
      !VD$ NODEPCHK
      DO i = ibgn, iend
        X(Jel(i)) = X(Jel(i)) - El(i)*X(irow)
      ENDDO
    ENDIF
  ENDDO
  !
  !------------- LAST LINE OF DLLTI2 FOLLOWS ----------------------------
END SUBROUTINE DLLTI2
