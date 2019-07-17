MODULE lapack
  USE service, ONLY : SP, DP
  IMPLICIT NONE

  INTERFACE
    PURE SUBROUTINE SHSEQR(Job,Compz,N,Ilo,Ihi,H,Ldh,Wr,Wi,Z,Ldz,Work,Lwork,Info)
      IMPORT SP
      INTEGER, INTENT(IN) :: Ihi, Ilo, Ldh, Ldz, Lwork, N
      INTEGER, INTENT(OUT) :: Info
      CHARACTER, INTENT(IN) :: Compz, Job
      REAL(SP), INTENT(INOUT) ::  H(Ldh,N), Z(Ldz,N)
      REAL(SP), INTENT(OUT) :: Wi(N), Wr(N), WORK(Lwork)
    END SUBROUTINE SHSEQR
    PURE SUBROUTINE CHSEQR(Job,Compz,N,Ilo,Ihi,H,Ldh,W,Z,Ldz,Work,Lwork,Info)
      IMPORT SP
      INTEGER, INTENT(IN) :: Ihi, Ilo, Ldh, Ldz, Lwork, N
      INTEGER, INTENT(OUT) :: Info
      CHARACTER, INTENT(IN) :: Compz, Job
      COMPLEX(SP), INTENT(INOUT) ::  H(Ldh,N), Z(Ldz,N)
      COMPLEX(SP), INTENT(OUT) :: W(N), WORK(Lwork)
    END SUBROUTINE CHSEQR
    PURE SUBROUTINE SGTSV(N,Nrhs,Dl,D,Du,B,Ldb,Info )
      IMPORT SP
      INTEGER, INTENT(IN) :: Ldb, N, Nrhs
      INTEGER, INTENT(OUT) :: Info
      REAL(SP), INTENT(INOUT) :: B(Ldb,Nrhs), D(N), Dl(N-1), Du(N-1)
    END SUBROUTINE SGTSV
    PURE SUBROUTINE DGTSV(N,Nrhs,Dl,D,Du,B,Ldb,Info )
      IMPORT DP
      INTEGER, INTENT(IN) :: Ldb, N, Nrhs
      INTEGER, INTENT(OUT) :: Info
      REAL(DP), INTENT(INOUT) :: B(Ldb,Nrhs), D(N), Dl(N-1), Du(N-1)
    END SUBROUTINE DGTSV
    PURE SUBROUTINE SGBTRS(Trans,N,Kl,Ku,Nrhs,Ab,Ldab,Ipiv,B,Ldb,Info)
      IMPORT SP
      CHARACTER, INTENT(IN) :: Trans
      INTEGER, INTENT(IN) :: Kl, Ku, Ldab, Ldb, N, Nrhs
      INTEGER, INTENT(OUT) :: Info
      INTEGER, INTENT(IN) :: Ipiv(N)
      REAL(SP), INTENT(IN) ::  Ab(Ldab,N)
      REAL(SP), INTENT(INOUT) :: B(Ldb,Nrhs)
    END SUBROUTINE SGBTRS
    PURE SUBROUTINE DGBTRS(Trans,N,Kl,Ku,Nrhs,Ab,Ldab,Ipiv,B,Ldb,Info)
      IMPORT DP
      CHARACTER, INTENT(IN) :: Trans
      INTEGER, INTENT(IN) :: Kl, Ku, Ldab, Ldb, N, Nrhs
      INTEGER, INTENT(OUT) :: Info
      INTEGER, INTENT(IN) :: Ipiv(N)
      REAL(DP), INTENT(IN) ::  Ab(Ldab,N)
      REAL(DP), INTENT(INOUT) :: B(Ldb,Nrhs)
    END SUBROUTINE DGBTRS
    PURE SUBROUTINE CGBTRS(Trans,N,Kl,Ku,Nrhs,Ab,Ldab,Ipiv,B,Ldb,Info)
      IMPORT SP
      CHARACTER, INTENT(IN) :: Trans
      INTEGER, INTENT(IN) :: Kl, Ku, Ldab, Ldb, N, Nrhs
      INTEGER, INTENT(OUT) :: Info
      INTEGER, INTENT(IN) :: Ipiv(N)
      COMPLEX(SP), INTENT(IN) ::  Ab(Ldab,N)
      COMPLEX(SP), INTENT(INOUT) :: B(Ldb,Nrhs)
    END SUBROUTINE CGBTRS
    PURE SUBROUTINE SGETRS(Trans,N,Nrhs,A,Lda,Ipiv,B,Ldb,Info)
      IMPORT SP
      CHARACTER, INTENT(IN) :: Trans
      INTEGER, INTENT(IN) :: Lda, Ldb, N, Nrhs
      INTEGER, INTENT(OUT) :: Info
      INTEGER, INTENT(IN) :: Ipiv(N)
      REAL(SP), INTENT(IN) :: A(Lda,N)
      REAL(SP), INTENT(INOUT) :: B(Ldb,Nrhs)
    END SUBROUTINE SGETRS
    PURE SUBROUTINE DGETRS(Trans,N,Nrhs,A,Lda,Ipiv,B,Ldb,Info)
      IMPORT DP
      CHARACTER, INTENT(IN) :: Trans
      INTEGER, INTENT(IN) :: Lda, Ldb, N, Nrhs
      INTEGER, INTENT(OUT) :: Info
      INTEGER, INTENT(IN) :: Ipiv(N)
      REAL(DP), INTENT(IN) :: A(Lda,N)
      REAL(DP), INTENT(INOUT) :: B(Ldb,Nrhs)
    END SUBROUTINE DGETRS
    PURE SUBROUTINE CGETRS(Trans,N,Nrhs,A,Lda,Ipiv,B,Ldb,Info)
      IMPORT SP
      CHARACTER, INTENT(IN) :: Trans
      INTEGER, INTENT(IN) :: Lda, Ldb, N, Nrhs
      INTEGER, INTENT(OUT) :: Info
      INTEGER, INTENT(IN) :: Ipiv(N)
      COMPLEX(SP), INTENT(IN) :: A(Lda,N)
      COMPLEX(SP), INTENT(INOUT) :: B(Ldb,Nrhs)
    END SUBROUTINE CGETRS
  END INTERFACE

END MODULE lapack