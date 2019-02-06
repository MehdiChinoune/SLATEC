!*==LA05BS.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK LA05BS
      SUBROUTINE LA05BS(A,Ind,Ia,N,Ip,Iw,W,G,B,Trans)
      IMPLICIT NONE
!*--LA05BS5
!*** Start of declarations inserted by SPAG
      INTEGER i , Ia , ii , j , k , k2 , kk , kl , kll , kp , kpc , l1 , LCOl , 
     &        LENl , LENu , LP , LROw , N , n1 , NCP
      INTEGER nz
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  LA05BS
!***SUBSIDIARY
!***PURPOSE  Subsidiary to SPLP
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (LA05BS-S, LA05BD-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     THIS SUBPROGRAM IS A SLIGHT MODIFICATION OF A SUBPROGRAM
!     FROM THE C. 1979 AERE HARWELL LIBRARY.  THE NAME OF THE
!     CORRESPONDING HARWELL CODE CAN BE OBTAINED BY DELETING
!     THE FINAL LETTER =S= IN THE NAMES USED HERE.
!     REVISED SEP. 13, 1979.
!
!     ROYALTIES HAVE BEEN PAID TO AERE-UK FOR USE OF THEIR CODES
!     IN THE PACKAGE GIVEN HERE.  ANY PRIMARY USAGE OF THE HARWELL
!     SUBROUTINES REQUIRES A ROYALTY AGREEMENT AND PAYMENT BETWEEN
!     THE USER AND AERE-UK.  ANY USAGE OF THE SANDIA WRITTEN CODES
!     SPLP( ) (WHICH USES THE HARWELL SUBROUTINES) IS PERMITTED.
!
! IP(I,1),IP(I,2) POINT TO START OF ROW/COLUMN I OF U.
! IW(I,1),IW(I,2) ARE LENGTHS OF ROW/COL I OF U.
! IW(.,3),IW(.,4) HOLD ROW/COL NUMBERS IN PIVOTAL ORDER.
!
!***SEE ALSO  SPLP
!***ROUTINES CALLED  XERMSG, XSETUN
!***COMMON BLOCKS    LA05DS
!***REVISION HISTORY  (YYMMDD)
!   811215  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900402  Added TYPE section.  (WRB)
!   920410  Corrected second dimension on IW declaration.  (WRB)
!***END PROLOGUE  LA05BS
      REAL A(Ia) , B(*) , am , W(*) , G , SMAll
      LOGICAL Trans
      INTEGER Ind(Ia,2) , Iw(N,8)
      INTEGER Ip(N,2)
      COMMON /LA05DS/ SMAll , LP , LENl , LENu , NCP , LROw , LCOl
!***FIRST EXECUTABLE STATEMENT  LA05BS
      IF ( G<0. ) THEN
!
        CALL XSETUN(LP)
        IF ( LP>0 ) CALL XERMSG('SLATEC','LA05BS',
     &                          'EARLIER ENTRY GAVE ERROR RETURN.',-8,2)
      ELSE
        kll = Ia - LENl + 1
        IF ( Trans ) THEN
!
!     MULTIPLY VECTOR BY INVERSE OF TRANSPOSE OF U
          DO i = 1 , N
            W(i) = B(i)
            B(i) = 0.
          ENDDO
          DO ii = 1 , N
            i = Iw(ii,4)
            am = W(i)
            IF ( am/=0. ) THEN
              j = Iw(ii,3)
              kp = Ip(j,1)
              am = am/A(kp)
              B(j) = am
              kl = Iw(j,1) + kp - 1
              IF ( kp/=kl ) THEN
                k2 = kp + 1
                DO k = k2 , kl
                  i = Ind(k,2)
                  W(i) = W(i) - am*A(k)
                ENDDO
              ENDIF
            ENDIF
          ENDDO
!
!     MULTIPLY VECTOR BY INVERSE OF TRANSPOSE OF L
          IF ( kll>Ia ) RETURN
          DO k = kll , Ia
            j = Ind(k,2)
            IF ( B(j)/=0. ) THEN
              i = Ind(k,1)
              B(i) = B(i) + A(k)*B(j)
            ENDIF
          ENDDO
        ELSE
!
!     MULTIPLY VECTOR BY INVERSE OF L
          IF ( LENl>0 ) THEN
            l1 = Ia + 1
            DO kk = 1 , LENl
              k = l1 - kk
              i = Ind(k,1)
              IF ( B(i)/=0. ) THEN
                j = Ind(k,2)
                B(j) = B(j) + A(k)*B(i)
              ENDIF
            ENDDO
          ENDIF
          DO i = 1 , N
            W(i) = B(i)
            B(i) = 0.
          ENDDO
!
!     MULTIPLY VECTOR BY INVERSE OF U
          n1 = N + 1
          DO ii = 1 , N
            i = n1 - ii
            i = Iw(i,3)
            am = W(i)
            kp = Ip(i,1)
            IF ( kp<=0 ) THEN
              kp = -kp
              Ip(i,1) = kp
              nz = Iw(i,1)
              kl = kp - 1 + nz
              k2 = kp + 1
              DO k = k2 , kl
                j = Ind(k,2)
                am = am - A(k)*B(j)
              ENDDO
            ENDIF
            IF ( am/=0. ) THEN
              j = Ind(kp,2)
              B(j) = am/A(kp)
              kpc = Ip(j,2)
              kl = Iw(j,2) + kpc - 1
              IF ( kl/=kpc ) THEN
                k2 = kpc + 1
                DO k = k2 , kl
                  i = Ind(k,1)
                  Ip(i,1) = -ABS(Ip(i,1))
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDIF
      END SUBROUTINE LA05BS
