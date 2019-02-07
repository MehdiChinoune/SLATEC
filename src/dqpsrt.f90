!*==DQPSRT.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DQPSRT
SUBROUTINE DQPSRT(Limit,Last,Maxerr,Ermax,Elist,Iord,Nrmax)
  IMPLICIT NONE
  !*--DQPSRT5
  !***BEGIN PROLOGUE  DQPSRT
  !***SUBSIDIARY
  !***PURPOSE  This routine maintains the descending ordering in the
  !            list of the local error estimated resulting from the
  !            interval subdivision process. At each call two error
  !            estimates are inserted using the sequential search
  !            method, top-down for the largest error estimate and
  !            bottom-up for the smallest error estimate.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (QPSRT-S, DQPSRT-D)
  !***KEYWORDS  SEQUENTIAL SORTING
  !***AUTHOR  Piessens, Robert
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !           de Doncker, Elise
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !***DESCRIPTION
  !
  !           Ordering routine
  !           Standard fortran subroutine
  !           Double precision version
  !
  !           PARAMETERS (MEANING AT OUTPUT)
  !              LIMIT  - Integer
  !                       Maximum number of error estimates the list
  !                       can contain
  !
  !              LAST   - Integer
  !                       Number of error estimates currently in the list
  !
  !              MAXERR - Integer
  !                       MAXERR points to the NRMAX-th largest error
  !                       estimate currently in the list
  !
  !              ERMAX  - Double precision
  !                       NRMAX-th largest error estimate
  !                       ERMAX = ELIST(MAXERR)
  !
  !              ELIST  - Double precision
  !                       Vector of dimension LAST containing
  !                       the error estimates
  !
  !              IORD   - Integer
  !                       Vector of dimension LAST, the first K elements
  !                       of which contain pointers to the error
  !                       estimates, such that
  !                       ELIST(IORD(1)),...,  ELIST(IORD(K))
  !                       form a decreasing sequence, with
  !                       K = LAST if LAST.LE.(LIMIT/2+2), and
  !                       K = LIMIT+1-LAST otherwise
  !
  !              NRMAX  - Integer
  !                       MAXERR = IORD(NRMAX)
  !
  !***SEE ALSO  DQAGE, DQAGIE, DQAGPE, DQAWSE
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   800101  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !***END PROLOGUE  DQPSRT
  !
  REAL(8) :: Elist , Ermax , errmax , errmin
  INTEGER i , ibeg , ido , Iord , isucc , j , jbnd , jupbn , k , Last , &
    Limit , Maxerr , Nrmax
  DIMENSION Elist(*) , Iord(*)
  !
  !           CHECK WHETHER THE LIST CONTAINS MORE THAN
  !           TWO ERROR ESTIMATES.
  !
  !***FIRST EXECUTABLE STATEMENT  DQPSRT
  IF ( Last>2 ) THEN
    !
    !           THIS PART OF THE ROUTINE IS ONLY EXECUTED IF, DUE TO A
    !           DIFFICULT INTEGRAND, SUBDIVISION INCREASED THE ERROR
    !           ESTIMATE. IN THE NORMAL CASE THE INSERT PROCEDURE SHOULD
    !           START AFTER THE NRMAX-TH LARGEST ERROR ESTIMATE.
    !
    errmax = Elist(Maxerr)
    IF ( Nrmax/=1 ) THEN
      ido = Nrmax - 1
      DO i = 1 , ido
        isucc = Iord(Nrmax-1)
        ! ***JUMP OUT OF DO-LOOP
        IF ( errmax<=Elist(isucc) ) EXIT
        Iord(Nrmax) = isucc
        Nrmax = Nrmax - 1
      ENDDO
    ENDIF
    !
    !           COMPUTE THE NUMBER OF ELEMENTS IN THE LIST TO BE MAINTAINED
    !           IN DESCENDING ORDER. THIS NUMBER DEPENDS ON THE NUMBER OF
    !           SUBDIVISIONS STILL ALLOWED.
    !
    jupbn = Last
    IF ( Last>(Limit/2+2) ) jupbn = Limit + 3 - Last
    errmin = Elist(Last)
    !
    !           INSERT ERRMAX BY TRAVERSING THE LIST TOP-DOWN,
    !           STARTING COMPARISON FROM THE ELEMENT ELIST(IORD(NRMAX+1)).
    !
    jbnd = jupbn - 1
    ibeg = Nrmax + 1
    IF ( ibeg<=jbnd ) THEN
      DO i = ibeg , jbnd
        isucc = Iord(i)
        ! ***JUMP OUT OF DO-LOOP
        IF ( errmax>=Elist(isucc) ) GOTO 100
        Iord(i-1) = isucc
      ENDDO
    ENDIF
    Iord(jbnd) = Maxerr
    Iord(jupbn) = Last
  ELSE
    Iord(1) = 1
    Iord(2) = 2
  ENDIF
  GOTO 300
  !
  !           INSERT ERRMIN BY TRAVERSING THE LIST BOTTOM-UP.
  !
  100  Iord(i-1) = Maxerr
  k = jbnd
  DO j = i , jbnd
    isucc = Iord(k)
    ! ***JUMP OUT OF DO-LOOP
    IF ( errmin<Elist(isucc) ) GOTO 200
    Iord(k+1) = isucc
    k = k - 1
  ENDDO
  Iord(i) = Last
  GOTO 300
  200  Iord(k+1) = Last
  !
  !           SET MAXERR AND ERMAX.
  !
  300  Maxerr = Iord(Nrmax)
  Ermax = Elist(Maxerr)
END SUBROUTINE DQPSRT
