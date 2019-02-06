*DECK OUTERR
      SUBROUTINE OUTERR (METHOD, IERR, IOUT, NFAIL, ISTDO, ITER, ERR)
C***BEGIN PROLOGUE  OUTERR
C***SUBSIDIARY
C***PURPOSE  Output error messages for the SLAP Quick Check.
C***LIBRARY   SLATEC (SLAP)
C***TYPE      SINGLE PRECISION (OUTERR-S, DUTERR-D)
C***AUTHOR  Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (510) 423-3141
C             seager@llnl.gov
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   881010  DATE WRITTEN
C   881213  Previous REVISION DATE
C   890920  Converted prologue to SLATEC 4.0 format.  (FNF)
C   920511  Added complete declaration section.  (WRB)
C   921021  Added 1P's to output formats.  (FNF)
C***END PROLOGUE  OUTERR
C     .. Scalar Arguments ..
      REAL ERR
      INTEGER IERR, IOUT, ISTDO, ITER, NFAIL
      CHARACTER METHOD*6
C***FIRST EXECUTABLE STATEMENT  OUTERR
      IF( IERR.NE.0 ) NFAIL = NFAIL + 1
      IF( IOUT.EQ.1 .AND. IERR.NE.0 ) THEN
         WRITE(ISTDO,1000) METHOD
      ENDIF
      IF( IOUT.EQ.2 ) THEN
         IF( IERR.EQ.0 ) THEN
            WRITE(ISTDO,1010) METHOD
         ELSE
            WRITE(ISTDO,1020) METHOD,IERR,ITER,ERR
         ENDIF
      ENDIF
      IF( IOUT.GE.3 ) THEN
         IF( IERR.EQ.0 ) THEN
            WRITE(ISTDO,1030) METHOD,IERR,ITER,ERR
         ELSE
            WRITE(ISTDO,1020) METHOD,IERR,ITER,ERR
         ENDIF
      ENDIF
      RETURN
 1000 FORMAT( 1X,A6,' : **** FAILURE ****')
 1010 FORMAT( 1X,A6,' : **** PASSED  ****')
 1020 FORMAT(' **************** WARNING ***********************'/
     $       ' **** ',A6,' Quick Test FAILED: IERR = ',I5,' ****'/
     $       ' **************** WARNING ***********************'/
     $       ' Iteration Count = ',I3,' Stop Test = ',1P,E12.6)
 1030 FORMAT(' ***************** PASSED ***********************'/
     $       ' **** ',A6,' Quick Test PASSED: IERR = ',I5,' ****'/
     $       ' ***************** PASSED ***********************'/
     $       ' Iteration Count = ',I3,' Stop Test = ',1P,E12.6)
C------------- LAST LINE OF OUTERR FOLLOWS ----------------------------
      END
