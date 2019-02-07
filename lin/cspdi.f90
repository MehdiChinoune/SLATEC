!*==CSPDI.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK CSPDI
SUBROUTINE CSPDI(Ap,N,Kpvt,Det,Work,Job)
  IMPLICIT NONE
  !*--CSPDI5
  !***BEGIN PROLOGUE  CSPDI
  !***PURPOSE  Compute the determinant and inverse of a complex symmetric
  !            matrix stored in packed form using the factors from CSPFA.
  !***LIBRARY   SLATEC (LINPACK)
  !***CATEGORY  D2C1, D3C1
  !***TYPE      COMPLEX (SSPDI-S, DSPDI-D, CHPDI-C, CSPDI-C)
  !***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX,
  !             PACKED, SYMMETRIC
  !***AUTHOR  Bunch, J., (UCSD)
  !***DESCRIPTION
  !
  !     CSPDI computes the determinant and inverse
  !     of a complex symmetric matrix using the factors from CSPFA,
  !     where the matrix is stored in packed form.
  !
  !     On Entry
  !
  !        AP      COMPLEX (N*(N+1)/2)
  !                the output from CSPFA.
  !
  !        N       INTEGER
  !                the order of the matrix A .
  !
  !        KVPT    INTEGER(N)
  !                the pivot vector from CSPFA.
  !
  !        WORK    COMPLEX(N)
  !                work vector.  Contents ignored.
  !
  !        JOB     INTEGER
  !                JOB has the decimal expansion  AB  where
  !                   if  B .NE. 0, the inverse is computed,
  !                   if  A .NE. 0, the determinant is computed.
  !
  !                For example, JOB = 11  gives both.
  !
  !     On Return
  !
  !        Variables not requested by JOB are not used.
  !
  !        AP     contains the upper triangle of the inverse of
  !               the original matrix, stored in packed form.
  !               The columns of the upper triangle are stored
  !               sequentially in a one-dimensional array.
  !
  !        DET    COMPLEX(2)
  !               determinant of original matrix.
  !               Determinant = DET(1) * 10.0**DET(2)
  !               with 1.0 .LE. ABS(DET(1)) .LT. 10.0
  !               or DET(1) = 0.0.
  !
  !     Error Condition
  !
  !        A division by zero will occur if the inverse is requested
  !        and  CSPCO  has set RCOND .EQ. 0.0
  !        or  CSPFA  has set  INFO .NE. 0 .
  !
  !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***ROUTINES CALLED  CAXPY, CCOPY, CDOTU, CSWAP
  !***REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891107  Corrected category and modified routine equivalence
  !           list.  (WRB)
  !   891107  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  CSPDI
  INTEGER N , Job
  COMPLEX Ap(*) , Work(*) , Det(2)
  INTEGER Kpvt(*)
  !
  COMPLEX ak , akkp1 , akp1 , CDOTU , d , t , temp
  REAL ten
  INTEGER ij , ik , ikp1 , iks , j , jb , jk , jkp1
  INTEGER k , kk , kkp1 , km1 , ks , ksj , kskp1 , kstep
  LOGICAL noinv , nodet
  REAL, EXTERNAL :: CABS1
  !
  !***FIRST EXECUTABLE STATEMENT  CSPDI
  noinv = MOD(Job,10)==0
  nodet = MOD(Job,100)/10==0
  !
  IF ( .NOT.(nodet) ) THEN
    Det(1) = (1.0E0,0.0E0)
    Det(2) = (0.0E0,0.0E0)
    ten = 10.0E0
    t = (0.0E0,0.0E0)
    ik = 0
    DO k = 1 , N
      kk = ik + k
      d = Ap(kk)
      !
      !           CHECK IF 1 BY 1
      !
      IF ( Kpvt(k)<=0 ) THEN
        !
        !              2 BY 2 BLOCK
        !              USE DET (D  T)  =  (D/T * C - T) * T
        !                      (T  C)
        !              TO AVOID UNDERFLOW/OVERFLOW TROUBLES.
        !              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG.
        !
        IF ( CABS1(t)/=0.0E0 ) THEN
          d = t
          t = (0.0E0,0.0E0)
        ELSE
          ikp1 = ik + k
          kkp1 = ikp1 + k
          t = Ap(kkp1)
          d = (d/t)*Ap(kkp1+1) - t
        ENDIF
      ENDIF
      !
      IF ( .NOT.(nodet) ) THEN
        Det(1) = d*Det(1)
        IF ( CABS1(Det(1))/=0.0E0 ) THEN
          DO WHILE ( CABS1(Det(1))<1.0E0 )
            Det(1) = CMPLX(ten,0.0E0)*Det(1)
            Det(2) = Det(2) - (1.0E0,0.0E0)
          ENDDO
          DO WHILE ( CABS1(Det(1))>=ten )
            Det(1) = Det(1)/CMPLX(ten,0.0E0)
            Det(2) = Det(2) + (1.0E0,0.0E0)
          ENDDO
        ENDIF
      ENDIF
      ik = ik + k
    ENDDO
  ENDIF
  !
  !     COMPUTE INVERSE(A)
  !
  IF ( .NOT.(noinv) ) THEN
    k = 1
    ik = 0
    DO WHILE ( k<=N )
      km1 = k - 1
      kk = ik + k
      ikp1 = ik + k
      IF ( Kpvt(k)<0 ) THEN
        !
        !              2 BY 2
        !
        kkp1 = ikp1 + k
        t = Ap(kkp1)
        ak = Ap(kk)/t
        akp1 = Ap(kkp1+1)/t
        akkp1 = Ap(kkp1)/t
        d = t*(ak*akp1-(1.0E0,0.0E0))
        Ap(kk) = akp1/d
        Ap(kkp1+1) = ak/d
        Ap(kkp1) = -akkp1/d
        IF ( km1>=1 ) THEN
          CALL CCOPY(km1,Ap(ikp1+1),1,Work,1)
          ij = 0
          DO j = 1 , km1
            jkp1 = ikp1 + j
            Ap(jkp1) = CDOTU(j,Ap(ij+1),1,Work,1)
            CALL CAXPY(j-1,Work(j),Ap(ij+1),1,Ap(ikp1+1),1)
            ij = ij + j
          ENDDO
          Ap(kkp1+1) = Ap(kkp1+1) + CDOTU(km1,Work,1,Ap(ikp1+1),1)
          Ap(kkp1) = Ap(kkp1) + CDOTU(km1,Ap(ik+1),1,Ap(ikp1+1),1)
          CALL CCOPY(km1,Ap(ik+1),1,Work,1)
          ij = 0
          DO j = 1 , km1
            jk = ik + j
            Ap(jk) = CDOTU(j,Ap(ij+1),1,Work,1)
            CALL CAXPY(j-1,Work(j),Ap(ij+1),1,Ap(ik+1),1)
            ij = ij + j
          ENDDO
          Ap(kk) = Ap(kk) + CDOTU(km1,Work,1,Ap(ik+1),1)
        ENDIF
        kstep = 2
      ELSE
        !
        !              1 BY 1
        !
        Ap(kk) = (1.0E0,0.0E0)/Ap(kk)
        IF ( km1>=1 ) THEN
          CALL CCOPY(km1,Ap(ik+1),1,Work,1)
          ij = 0
          DO j = 1 , km1
            jk = ik + j
            Ap(jk) = CDOTU(j,Ap(ij+1),1,Work,1)
            CALL CAXPY(j-1,Work(j),Ap(ij+1),1,Ap(ik+1),1)
            ij = ij + j
          ENDDO
          Ap(kk) = Ap(kk) + CDOTU(km1,Work,1,Ap(ik+1),1)
        ENDIF
        kstep = 1
      ENDIF
      !
      !           SWAP
      !
      ks = ABS(Kpvt(k))
      IF ( ks/=k ) THEN
        iks = (ks*(ks-1))/2
        CALL CSWAP(ks,Ap(iks+1),1,Ap(ik+1),1)
        ksj = ik + ks
        DO jb = ks , k
          j = k + ks - jb
          jk = ik + j
          temp = Ap(jk)
          Ap(jk) = Ap(ksj)
          Ap(ksj) = temp
          ksj = ksj - (j-1)
        ENDDO
        IF ( kstep/=1 ) THEN
          kskp1 = ikp1 + ks
          temp = Ap(kskp1)
          Ap(kskp1) = Ap(kkp1)
          Ap(kkp1) = temp
        ENDIF
      ENDIF
      ik = ik + k
      IF ( kstep==2 ) ik = ik + k + 1
      k = k + kstep
    ENDDO
  ENDIF
END SUBROUTINE CSPDI
