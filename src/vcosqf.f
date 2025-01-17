      SUBROUTINE VCOSQF(M,N,X,XT,MDIMX,WSAVE)
C***BEGIN PROLOGUE  VCOSQF
C***DATE WRITTEN   860701   (YYMMDD)
C***REVISION DATE  900509   (YYMMDD)
C***CATEGORY NO.  J1A3
C***KEYWORDS  FAST FOURIER TRANSFORM, COSINE TRANSFORM, ODD WAVE
C             NUMBERS, MULTIPLE SEQUENCES
C***AUTHOR  BOISVERT, R. F. (NIST)
C***PURPOSE  Forward cosine transform, odd wave numbers, M sequences.
C***DESCRIPTION
C
C  Subroutine VCOSQF computes the forward fast Fourier cosine transform
C  of M quarter wave sequences.  That is, cosine series representations
C  with only odd wave numbers.  The transform is defined below at output
C  parameter X.
C
C  The array WSAVE which is used by subroutine VCOSQF must be
C  initialized by calling subroutine VCOSQI(N,WSAVE).
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
C          in the program that calls VCOSQF.  The WSAVE array must be
C          initialized by calling subroutine VCOSQI(N,WSAVE), and a
C          different WSAVE array must be used for each different
C          value of N.  This initialization does not have to be
C          repeated so long as N remains unchanged.
C
C  Output Parameters
C
C  X       For I=1,...,N and J=1,...,M
C
C               X(I) = ( X(1) + the sum from K=2 to K=N of
C
C                  2*X(K)*COS((2*I-1)*(K-1)*PI/(2*N)) )/SQRT(4*N)
C
C  WSAVE   contains initialization calculations which must not
C          be destroyed between calls of VCOSQF or VCOSQB.
C
C  -----------------------------------------------------------------
C
C  NOTE  -  A call of VCOSQF followed immediately by a call of
C           of VCOSQB will return the original sequences X.  Thus,
C           VCOSQB is the correctly normalized inverse VCOSQF.
C
C  -----------------------------------------------------------------
C
C  VCOSQF is a straightforward extension of the subprogram COSQF to
C  handle M simultaneous sequences.  COSQF was originally developed
C  by P. N. Swarztrauber of NCAR.
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
C               pp. 51-83.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  VCOSQF
      DIMENSION       X(MDIMX,*), XT(MDIMX,*), WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  VCOSQF
      IF (M .LE. 0)  GO TO 900
      IF (N .GT. 2)  GO TO 300
      IF (N .LT. 2)  GO TO 900
C
C  CASE  N = 2
C
      SQRT2 = SQRT(2.0)
      SCALE = 0.50E0/SQRT2
      DO 210 J=1,M
         TSQX = SQRT2*X(J,2)
         X(J,2) = SCALE*(X(J,1)-TSQX)
         X(J,1) = SCALE*(X(J,1)+TSQX)
  210 CONTINUE
      GO TO 900
C
C  CASE N .GT. 2
C
  300 CONTINUE
C
C     ... PREPROCESSING
C
      NS2 = (N+1)/2
      NP2 = N+2
      DO 310 K=2,NS2
         KC = NP2-K
         DO 310 J=1,M
            XT(J,K) = X(J,K)+X(J,KC)
            XT(J,KC) = X(J,K)-X(J,KC)
  310 CONTINUE
      MODN = MOD(N,2)
      IF (MODN .EQ. 0) THEN
         DO 320 J=1,M
            XT(J,NS2+1) = X(J,NS2+1)+X(J,NS2+1)
  320    CONTINUE
      ENDIF
      DO 330 K=2,NS2
         KC = NP2-K
         DO 330 J=1,M
            X(J,K) = WSAVE(K-1)*XT(J,KC)+WSAVE(KC-1)*XT(J,K)
            X(J,KC) = WSAVE(K-1)*XT(J,K)-WSAVE(KC-1)*XT(J,KC)
  330 CONTINUE
      IF (MODN .EQ. 0) THEN
         DO 340 J=1,M
            X(J,NS2+1) = WSAVE(NS2)*XT(J,NS2+1)
  340    CONTINUE
      ENDIF
C
C     ... REAL, PERIODIC TRANSFORM
C
      CALL VRFFTF (M,N,X,XT,MDIMX,WSAVE(N+1))
C
C     ... POSTPROCESSING
C
      DO 350 I=3,N,2
         DO 350 J=1,M
            XIM1 = X(J,I-1)-X(J,I)
            X(J,I) = X(J,I-1)+X(J,I)
            X(J,I-1) = XIM1
  350 CONTINUE
C
C     ... NORMALIZATION
C
      SCALE = 0.5
      DO 360 I=1,N
         DO 360 J=1,M
            X(J,I) = SCALE*X(J,I)
  360 CONTINUE
C
C  EXIT
C
  900 CONTINUE
      RETURN
      END
