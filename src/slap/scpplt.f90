!** SCPPLT
SUBROUTINE SCPPLT(N,Nelt,Ia,Ja,A,Isym,Iunit)
  IMPLICIT NONE
  !>
  !***
  !  Printer Plot of SLAP Column Format Matrix.
  !            Routine to print out a SLAP Column format matrix in a
  !            "printer plot" graphical representation.
  !***
  ! **Library:**   SLATEC (SLAP)
  !***
  ! **Category:**  N1
  !***
  ! **Type:**      SINGLE PRECISION (SCPPLT-S, DCPPLT-D)
  !***
  ! **Keywords:**  DIAGNOSTICS, LINEAR SYSTEM, SLAP SPARSE
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
  !     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, IUNIT
  !     REAL    A(NELT)
  !
  !     CALL SCPPLT( N, NELT, IA, JA, A, ISYM, IUNIT )
  !
  !- Arguments:
  ! N      :IN       Integer
  !         Order of the Matrix.
  !         If N.gt.MAXORD, only the leading MAXORD x MAXORD
  !         submatrix will be printed.  (Currently MAXORD = 225.)
  ! NELT   :IN       Integer.
  !         Number of non-zeros stored in A.
  ! IA     :IN       Integer IA(NELT).
  ! JA     :IN       Integer JA(NELT).
  ! A      :IN       Real A(NELT).
  !         These arrays should hold the matrix A in the SLAP
  !         Column format.  See "Description", below.
  ! ISYM   :IN       Integer.
  !         Flag to indicate symmetric storage format.
  !         If ISYM=0, all non-zero entries of the matrix are stored.
  !         If ISYM=1, the matrix is symmetric, and only the lower
  !         triangle of the matrix is stored.
  ! IUNIT  :IN       Integer.
  !         Fortran logical I/O device unit number to write the matrix
  !         to.  This unit must be connected in a system dependent fashion
  !         to a file or the console or you will get a nasty message
  !         from the Fortran I/O libraries.
  !
  !- Description:
  !       This routine prints out a SLAP  Column format matrix  to the
  !       Fortran logical I/O unit   number  IUNIT.  The  numbers them
  !       selves  are not printed  out, but   rather  a one  character
  !       representation of the numbers.   Elements of the matrix that
  !       are not represented in the (IA,JA,A)  arrays are  denoted by
  !       ' ' character (a blank).  Elements of A that are *ZERO* (and
  !       hence  should  really not be  stored) are  denoted  by a '0'
  !       character.  Elements of A that are *POSITIVE* are denoted by
  !       'D' if they are Diagonal elements  and '#' if  they are off
  !       Diagonal  elements.  Elements of  A that are *NEGATIVE* are
  !       denoted by 'N'  if they  are Diagonal  elements and  '*' if
  !       they are off Diagonal elements.
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
  !- Cautions:
  !     This routine will attempt to write to the Fortran logical output
  !     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that
  !     this logical unit is attached to a file or terminal before calling
  !     this routine with a non-zero value for IUNIT.  This routine does
  !     not check for the validity of a non-zero IUNIT unit number.
  !
  !- Portability:
  !     This routine, as distributed, can generate lines up to 229
  !     characters long.  Some Fortran systems have more restricted
  !     line lengths.  Change parameter MAXORD and the large number
  !     in FORMAT 1010 to reduce this line length.
  !
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
  !   921007  Replaced hard-wired 225 with parameter MAXORD.  (FNF)
  !   921021  Corrected syntax of CHARACTER declaration.  (FNF)
  !   930701  Updated CATEGORY section.  (FNF, WRB)

  !     .. Scalar Arguments ..
  INTEGER Isym, Iunit, N, Nelt
  !     .. Array Arguments ..
  REAL A(Nelt)
  INTEGER Ia(Nelt), Ja(Nelt)
  !     .. Parameters ..
  INTEGER, PARAMETER :: MAXORD = 225
  !     .. Local Scalars ..
  INTEGER i, icol, irow, j, jbgn, jend, nmax
  !     .. Local Arrays ..
  CHARACTER chmat(MAXORD)*(MAXORD)
  !     .. Intrinsic Functions ..
  INTRINSIC MIN, MOD, REAL
  !* FIRST EXECUTABLE STATEMENT  SCPPLT
  !
  !         Set up the character matrix...
  !
  nmax = MIN(MAXORD,N)
  DO i = 1, nmax
    chmat(i)(1:nmax) = ' '
  ENDDO
  DO icol = 1, nmax
    jbgn = Ja(icol)
    jend = Ja(icol+1) - 1
    DO j = jbgn, jend
      irow = Ia(j)
      IF ( irow<=nmax ) THEN
        IF ( Isym/=0 ) THEN
          !         Put in non-sym part as well...
          IF ( A(j)==0.0E0 ) THEN
            chmat(irow)(icol:icol) = '0'
          ELSEIF ( A(j)>0.0E0 ) THEN
            chmat(irow)(icol:icol) = '#'
          ELSE
            chmat(irow)(icol:icol) = '*'
          ENDIF
        ENDIF
        IF ( irow==icol ) THEN
          !         Diagonal entry.
          IF ( A(j)==0.0E0 ) THEN
            chmat(irow)(icol:icol) = '0'
          ELSEIF ( A(j)>0.0E0 ) THEN
            chmat(irow)(icol:icol) = 'D'
          ELSE
            chmat(irow)(icol:icol) = 'N'
          ENDIF
          !         Off-Diagonal entry
        ELSEIF ( A(j)==0.0E0 ) THEN
          chmat(irow)(icol:icol) = '0'
        ELSEIF ( A(j)>0.0E0 ) THEN
          chmat(irow)(icol:icol) = '#'
        ELSE
          chmat(irow)(icol:icol) = '*'
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  !
  !         Write out the heading.
  WRITE (Iunit,99001) N, Nelt, REAL(Nelt)/(N*N)
  !
  99001 FORMAT (/'**** Picture of Column SLAP matrix follows ****'/&
    ' N, NELT and Density = ',2I10,E16.7)
  WRITE (Iunit,99002) (MOD(i,10),i=1,nmax)
  !      The following assumes MAXORD.le.225.
  99002 FORMAT (4X,225(I1))
  !
  !         Write out the character representations matrix elements.
  DO irow = 1, nmax
    WRITE (Iunit,99003) irow, chmat(irow)(1:nmax)
    99003 FORMAT (1X,I3,A)
  ENDDO
  RETURN
  !------------- LAST LINE OF SCPPLT FOLLOWS ----------------------------
END SUBROUTINE SCPPLT
