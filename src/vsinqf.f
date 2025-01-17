      SUBROUTINE VSINQF(M,N,X,XT,MDIMX,WSAVE)
C***BEGIN PROLOGUE  VSINQF
C***DATE WRITTEN   860701   (YYMMDD)
C***REVISION DATE  900509   (YYMMDD)
C***CATEGORY NO.  J1A3
C***KEYWORDS  FAST FOURIER TRANSFORM, SINE TRANSFORM, ODD WAVE
C             NUMBERS, MULTIPLE SEQUENCES
C***AUTHOR  BOISVERT, R. F. (NIST)
C***PURPOSE  Forward sine transform, odd wave numbers, M sequences.
C***DESCRIPTION
C
C  Subroutine VSINQF computes the forward fast Fourier sine transform
C  of M quarter wave sequences.  That is, sine series representations
C  with only odd wave numbers.  The transform is defined below at output
C  parameter X.
C
C  The array WSAVE which is used by subroutine VSINQF must be
C  initialized by calling subroutine VSINQI(N,WSAVE).
C
C
C  Input Parameters
C
C  M       the number of sequences to be transformed.
C
C  N       the length of the sequences to be transformed.  The method
C          is most efficient when N is a product of small primes.
C
C  X       an array of size at least X(MDIMX,N) which contains the
C          the sequences to be transformed.  The sequences are stored
C          in the ROWS of X.  Thus, the Jth sequence is stored in
C          X(J,I), I=1,..,N.
C
C  XT      a work array of size at least XT(MDIMX,N).
C
C  MDIMX   the first dimension of the array X exactly as it appears in
C          the calling program.
C
C  WSAVE   a work array which must be dimensioned at least 2*N+15
C          in the program that calls VSINQF.  The WSAVE array must be
C          initialized by calling subroutine VSINQI(N,WSAVE), and a
C          different WSAVE array must be used for each different
C          value of N.  This initialization does not have to be
C          repeated so long as N remains unchanged.
C
C  Output Parameters
C
C  X       For I=1,...,N and J=1,...,M
C
C               X(I) = ( (-1)**(I-1)*X(N)
C
C                        + the sum from K=1 to K=N-1 of
C
C                        2*X(K)*SIN((2*I-1)*K*PI/(2*N)) )/SQRT(4*N)
C
C  WSAVE   contains initialization calculations which must not
C          be destroyed between calls of VSINQF or VSINQB.
C
C  -----------------------------------------------------------------
C
C  NOTE  -  A call of VSINQF followed immediately by a call of
C           of VSINQB will return the original sequences X.  Thus,
C           VSINQB is the correctly normalized inverse VSINQF.
C
C  -----------------------------------------------------------------
C
C  VSINQF is a straightforward extension of the subprogram SINQF to
C  handle M simultaneous sequences.  SINQF was originally developed
C  by P. N. Swarztrauber of NCAR.
C
 
C***ROUTINES CALLED  VCOSQF
C***END PROLOGUE  VSINQF
      DIMENSION       X(MDIMX,*), XT(MDIMX,*), WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  VSINQF
      IF (M .LE. 0)  GO TO 900
      IF (N .LE. 1)  GO TO 900
C
C     ... PREPROCESSING
C
      NS2 = N/2
      DO 110 K=1,NS2
         KC = N-K
         DO 110 J=1,M
            XHOLD = X(J,K)
            X(J,K) = X(J,KC+1)
            X(J,KC+1) = XHOLD
  110 CONTINUE
C
C     ... COSINE QUARTER WAVE TRANSFORM
C
      CALL VCOSQF (M,N,X,XT,MDIMX,WSAVE)
C
C     ... POSTPROCESSING
C
      DO 120 K=2,N,2
         DO 120 J=1,M
            X(J,K) = -X(J,K)
  120 CONTINUE
C
C  EXIT
C
  900 CONTINUE
      RETURN
      END
