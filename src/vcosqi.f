      SUBROUTINE VCOSQI(N,WSAVE)
C***BEGIN PROLOGUE  VCOSQI
C***DATE WRITTEN   860701   (YYMMDD)
C***REVISION DATE  900509   (YYMMDD)
C***CATEGORY NO.  J1A3
C***KEYWORDS  FAST FOURIER TRANSFORM, COSINE TRANSFORM, ODD WAVE
C             NUMBERS, MULTIPLE SEQUENCES
C***AUTHOR  BOISVERT, R. F. (NIST)
C***PURPOSE  Initialize for VCOSQF and VCOSQB.
C***DESCRIPTION
C
C  Subroutine VCOSQI initializes the array WSAVE which is used in
C  both VCOSQF and VCOSQB.  The prime factorization of N together with
C  a tabulation of the trigonometric functions are computed and
C  stored in WSAVE.
C
C  Input Parameter
C
C  N       the length of the array to be transformed.  The method
C          is most efficient when N is a product of small primes.
C
C  Output Parameter
C
C  WSAVE   a work array which must be dimensioned at least 2*N+15.
C          The same work array can be used for both VCOSQF and VCOSQB
C          as long as N remains unchanged.  Different WSAVE arrays
C          are required for different values of N.  The contents of
C          WSAVE must not be changed between calls of VCOSQF or VCOSQB.
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
C               pp. 51-83.
C***ROUTINES CALLED  VRFFTI
C***END PROLOGUE  VCOSQI
      DIMENSION       WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  VCOSQI
      PIH = 0.5*PIMACH(1.0)
      DT = PIH/REAL(N)
      FK = 0.
      DO 101 K=1,N
         FK = FK+1.
         WSAVE(K) = COS(FK*DT)
  101 CONTINUE
      CALL VRFFTI (N,WSAVE(N+1))
      RETURN
      END
