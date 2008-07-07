      Module CSTEED
C***  This module is for information & storage only.
C     (it is not essential to working of the code)

      REAL(kind=selected_real_kind(12,100)) RERR
      INTEGER NFP,N11,NPQ(2),N20,KAS(2)
      End Module CSTEED
      
      
      SUBROUTINE COULCC(XX,ETA1,ZLMIN,NL, FC,GC,FCP,GCP, SIG,
     X                  MODE1,KFN,IFAIL)
C  COULCC36 for f90
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  COMPLEX COULOMB WAVEFUNCTION PROGRAM USING STEED'S METHOD           C
C                                                                      C
C  I.J. Thompson and A.R. Barnett, Sept. 1983                          C
C                                                                      C
C  original program  RCWFN       in    CPC  8 (1974) 377-395           C
C                 +  RCWFF       in    CPC 11 (1976) 141-142           C
C                 +  COULFG      in    CPC 27 (1982) 147-166           C
C  description of real algorithm in    CPC 21 (1981) 297-314           C
C  description of complex algorithm    JCP 64 (1986) 490-509           C
C  this version written up       in    CPC 36 (1985) 363-372 ACDP      C
C          revised for f90       in    CPC XX (YYYY) ZZZ-ZZZ           C
C                                                                      C
C  COULCC returns F,G,G',G',SIG for complex XX, ETA1, and ZLMIN,       C
C   for NL integer-spaced lambda values ZLMIN to ZLMIN+NL-1 inclusive, C
C   thus giving  complex-energy solutions to the Coulomb Schrodinger   C
C   equation,to the Klein-Gordon equation and to suitable forms of     C
C   the Dirac equation ,also spherical & cylindrical Bessel equations  C
C                                                                      C
C  if |MODE1|= 1  get F,G,F',G'   for integer-spaced lambda values     C
C            = 2      F,G      unused arrays must be dimensioned in    C
C            = 3      F,  F'          call to at least length (1)      C
C            = 4      F                                                C
C            = 11 get F,H+,F',H+' ) if KFN=0, H+ = G + i.F        )    C
C            = 12     F,H+        )       >0, H+ = J + i.Y = H(1) ) in C
C            = 21 get F,H-,F',H-' ) if KFN=0, H- = G - i.F        ) GC C
C            = 22     F,H-        )       >0, H- = J - i.Y = H(2) )    C
C                                                                      C
C     if MODE1<0 then the values returned are scaled by an exponential C
C                factor (dependent only on XX) to bring nearer unity   C
C                the functions for large |XX|, small ETA & |ZL| < |XX| C
C        Define SCALE = (  0        if MODE1 > 0                       C
C                       (  AIMAG(XX) if MODE1 < 0  &  KFN < 3          C
C                       (  REAL(XX) if MODE1 < 0  &  KFN = 3           C
C        then FC = EXP(-ABS(SCALE)) * ( F, j, J, or I)                 C
C         and GC = EXP(-ABS(SCALE)) * ( G, y, or Y )                   C
C               or EXP(SCALE)       * ( H+, H(1), or K)                C
C               or EXP(-SCALE)      * ( H- or H(2) )                   C
C                                                                      C
C  if  KFN  =  0,-1  complex Coulomb functions are returned   F & G    C
C           =  1   spherical Bessel      "      "     "       j & y    C
C           =  2 cylindrical Bessel      "      "     "       J & Y    C
C           =  3 modified cyl. Bessel    "      "     "       I & K    C
C                                                                      C
C          and where Coulomb phase shifts put in SIG if KFN=0 (not -1) C
C                                                                      C
C  The use of MODE and KFN is independent                              C
C    (except that for KFN=3,  H(1) & H(2) are not given)               C
C                                                                      C
C  With negative orders lambda, COULCC can still be used but with      C
C    reduced accuracy as CF1 becomes unstable. The user is thus        C
C    strongly advised to use reflection formulae based on              C
C    H+-(ZL,,) = H+-(-ZL-1,,) * exp +-i(sig(ZL)-sig(-ZL-1)-(ZL+1/2)pi) C
C                                                                      C
C  Precision:  results to within 2-3 decimals of 'machine accuracy',   C
C              except in the following cases:                          C
C              (1) if CF1A fails because X too small or ETA too large  C
C               the F solution  is less accurate if it decreases with  C
C               decreasing lambda (e.g. for lambda.LE.-1 & ETA.NE.0)   C
C              (2) if ETA is large (e.g. >> 50) and X inside the       C
C                turning point, then progressively less accuracy       C
C              (3) errors grow near (but not at) the logarithmic       C
C               singularity in the irregular solution for abs(X) << 1. C
C               If D = distance of ZLMIN from the logarithmic soln.,   C
C                then relative errors in GC will be up to min(D,1/D).  C
C                                                                      C
C              RERR in module CSTEED traces the main roundoff errors.  C
C                                                                      C
C   COULCC is coded for:                                               C
C          double precision at least 12-digits and D+-100 exponents    C
C          extended precision at least 24-digits and D+-100 exponents  C
C          Change these values to your requirements.                   C
C                                                                      C
C   IFAIL  on input   = 0 : no printing of error messages              C
C                    ne 0 : print error messages on file 6             C
C   IFAIL  in output = -2 : argument out of range                      C
C                    = -1 : one of the continued fractions failed,     C
C                           or arithmetic check before final recursion C
C                    =  0 : All Calculations satisfactory              C
C                    ge 0 : results available for orders up to & at    C
C                             position NL-IFAIL in the output arrays.  C
C                    = -3 : values at ZLMIN not found as over/underflowC
C                    = -4 : roundoff errors make results meaningless   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     Machine dependent constants :                                    C
C                                                                      C
C     ACCUR    target bound on relative error (except near 0 crossings)C
C               (ACCUR should be at least 100 * ACC8)                  C
C     ACC8     smallest number with 1+ACC8 .ne.1 in REAL(kind=kreal)   C
C     ACC16    smallest number with 1+ACC16.ne.1 in REAL*16 arithmetic C
C                                                                      C
C     ROUTINES CALLED :       LOGAM/CLOGAM/CDIGAM,                     C
C                             F20, CF1A, RCF, CF1C, CF2, F11, CF1R     C
C     Intrinsic functions :   MIN, MAX, SQRT, REAL, AIMAG, ABS, LOG, EXP,
C      (Generic names)        NINT, MOD, ATAN, ATAN2, COS, SIN, CMPLX,
C                             SIGN, CONJG, INT, TANH                   C
C     Note: Statement fntn.   NINTC = integer nearest to a complex no. C
C                                                                      C
C     Parameters determining region of calculations :                  C
C                                                                      C
C        R20      estimate of (2F0 iterations)/(CF2 iterations)        C
C        ASYM     minimum X/(ETA**2+L) for CF1A to converge easily     C
C        XNEAR    minimum ABS(X) for CF2 to converge accurately        C
C        LIMIT    maximum no. iterations for CF1, CF2, and 1F1 series  C
C        JMAX     size of work arrays for Pade accelerations           C
C        NDROP    number of successive decrements to define instabilityC
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      Use CSTEED
      parameter(kreal=selected_real_kind(12,100))
      
      IMPLICIT COMPLEX(kind=kreal) (A-H,O-Z)
      PARAMETER(JMAX=50)
      DIMENSION FC(NL),GC(NL),FCP(NL),GCP(NL),SIG(NL),XRCF(JMAX,4)
      LOGICAL PR,ETANE0,IFCP,RLEL,DONEM,UNSTAB,ZLNEG,AXIAL,NOCF2,NPINT
      REAL(kind=kreal) ERR,ABSC,ACCUR,ACCT,ACC8,ACCH,ACC16,ACCB, 
     X       XNEAR,CF1R,DRL,RX,RETA,XL,
     X       ZERO,ONE,TWO,HALF,HPI,TLOG,FPMAX,FPMIN,FPLMIN,FPLMAX,
     X       PACCQ,EPS,OFF,SCALE,SF,SFSH,TA,RK,OMEGA,R20,ASYM,ABSX
C
      intent(in) :: XX,ETA1,ZLMIN,NL,MODE1,KFN
      intent(out):: FC,GC,FCP,GCP,SIG
      intent(inout):: IFAIL
C
      DATA ZERO,ONE,TWO,LIMIT /0.0D+0, 1.0D+0, 2.0D+0, 20000 /,
     X     HALF, CI / 0.5D+0, (0D+0, 1D+0) /,
     X     R20,ASYM,XNEAR,NDROP / 3., 3., .5, 5 /, ACCUR / 1d-14 /

      DRL(X) = REAL(X,kreal)
!      DRL(X) = REAL(X)
      NINTC(W) = NINT(REAL(W))
      ABSC(W) = ABS(REAL(W)) + ABS(AIMAG(W))
      NPINT(W,ACCB) = ABSC(NINTC(W)-W).LT.ACCB .AND. REAL(W).LT.HALF
C
      MODE = MOD(ABS(MODE1),10)
      IFCP = MOD(MODE,2).EQ.1
      PR = IFAIL.NE.0
      IFAIL = -2
      N11   = 0
      NFP   = 0
      KAS(1)   = 0
      KAS(2)   = 0
      NPQ(1)   = 0
      NPQ(2)   = 0
      N20 = 0
      HPI = TWO*ATAN(ONE)
      TLOG = LOG(TWO)
!				Machine-dependent constants
      ACC8 = epsilon(ONE)
      FPMIN = TINY(ONE) / ACC8
      FPMAX = HUGE(ONE) * ACC8
      FPLMAX= log(FPMAX)
      FPLMIN = log(FPMIN)
      FPHMIN = sqrt(FPMIN)
      ACCUR = MAX(ACCUR, 50*ACC8)
      ACCT = ACCUR * .5
C                       initialise the log-gamma function :
      CALL GAMMA0(ACC8,FPLMIN) 
      ACCH  = SQRT(ACCUR)
      ACCB  = SQRT(ACCH)
      RERR = ACCT
C
      CIK = ONE
         IF(KFN.GE.3) CIK = CI * SIGN(ONE,FPMIN-AIMAG(XX))
      X     = XX * CIK
      ETA   = ETA1
      IF(KFN .GT. 0) ETA = ZERO
         ETANE0  = ABSC(ETA).GT.ACC8
         ETAI = ETA*CI
      DELL  = ZERO
      IF(KFN .GE. 2)  DELL = HALF
      ZM1   = ZLMIN - DELL
      SCALE = ZERO
      IF(MODE1.LT.0) SCALE = AIMAG(X)
C
      M1 = 1
      L1  = M1 + NL - 1
      RLEL = ABS(AIMAG(ETA)) + ABS(AIMAG(ZM1)) .LT. ACC8
      ABSX = ABS(X)
      AXIAL = RLEL .AND. ABS(AIMAG(X)) .LT. ACC8 * ABSX
      IF(MODE.LE.2 .AND. ABSX.LT.FPMIN) GO TO 310
      XI  = ONE/X
      XLOG = LOG(X)
C            log with cut along the negative real axis! see also OMEGA
      ID = 1
      DONEM = .FALSE.
         UNSTAB = .FALSE.
      LF = M1
      IFAIL = -1
   10    ZLM = ZM1 + LF - M1
         ZLL = ZM1 + L1 - M1
C
C ***       ZLL  is final lambda value, or 0.5 smaller for J,Y Bessels
C
              Z11 = ZLL
              IF(ID.LT.0) Z11 = ZLM
              P11 = CI*SIGN(ONE,ACC8-AIMAG(ETA))
      LAST = L1
C
C ***       Find phase shifts and Gamow factor at lambda = ZLL
C
      PK = ZLL + ONE
      AA = PK - ETAI
      AB = PK + ETAI
      BB = TWO*PK
         ZLNEG = NPINT(BB,ACCB)
                     CLGAA = CLOGAM(AA)
                     CLGAB = CLGAA
         IF(ETANE0.AND..NOT.RLEL)  CLGAB = CLOGAM(AB)
         IF(ETANE0.AND.     RLEL)  CLGAB = CONJG(CLGAA)
         SIGMA = (CLGAA - CLGAB) * CI*HALF
         IF(KFN.EQ.0) SIG(L1) = SIGMA
         IF(.NOT.ZLNEG) CLL = ZLL*TLOG- HPI*ETA - CLOGAM(BB)
     X                                          + (CLGAA+CLGAB)*HALF
              THETA  = X - ETA*(XLOG+TLOG) - ZLL*HPI + SIGMA
C
      TA = (AIMAG(AA)**2+AIMAG(AB)**2+ABS(REAL(AA))+ABS(REAL(AB)))*HALF
      IF(ID.GT.0 .AND. ABSX .LT. TA*ASYM .AND. .NOT.ZLNEG) GO TO 20
C
C ***         use CF1 instead of CF1A, if predicted to converge faster,
C                 (otherwise using CF1A as it treats negative lambda &
C                  recurrence-unstable cases properly)
C
           RK = SIGN(ONE, REAL(X) + ACC8)
           P =  THETA
           IF(RK.LT.0) P = -X + ETA*(LOG(-X)+TLOG)-ZLL*HPI-SIGMA
      F = RK * CF1A(X*RK,ETA*RK,ZLL,P,ACCT,JMAX/2,NFP,FEST,ERR,
     X                              FPMAX,XRCF,XRCF(1,3), XRCF(1,4))
      FESL = LOG(FEST) + ABS(AIMAG(X))
         NFP = - NFP
      IF(NFP.LT.0   .OR.(UNSTAB.AND.ERR.LT.ACCB)) GO TO 40
      IF(.NOT.ZLNEG .OR. UNSTAB.AND.ERR.GT.ACCB)  GO TO 20
         IF(PR) WRITE(6,1060) '-L',ERR
         IF(ERR.GT.ACCB) GO TO 280
         GO TO 40
C
C ***    evaluate CF1  =  f   =  F'(ZLL,ETA,X)/F(ZLL,ETA,X)
C
   20 IF(AXIAL) THEN
C                                                        REAL VERSION
      RX = X; RETA = ETA; XL = ZLL
      F = CF1R(RX,RETA,XL,ACC8,SF ,RK,  ETANE0,LIMIT,ERR,NFP,
     X         ACCH,FPMIN,FPMAX,PR,'COULCC')
          FCL = SF
          TPK1= RK
         ELSE
C                                                        COMPLEX VERSION
      F = CF1C(X,ETA,ZLL,ACC8,FCL,TPK1,ETANE0,LIMIT,ERR,NFP,
     X         ACCH,FPMIN,FPMAX,PR,'COULCC')
         ENDIF
      IF(ERR.GT.ONE) GO TO 390
C
C ***  Make a simple check for CF1 being badly unstable:
C
      IF(ID.LT.0) GO TO 30
      UNSTAB = REAL((ONE-ETA*XI)*CI*AIMAG(THETA)/F).GT.ZERO
     X .AND..NOT.AXIAL .AND. ABS(AIMAG(THETA)).GT.-LOG(ACC8)*.5
     X .AND. ABSC(ETA)+ABSC(ZLL).LT.ABSC(X)
      IF(UNSTAB) GO TO 60
C
C *** compare accumulated phase FCL with asymptotic phase for G(k+1) :
C     to determine estimate of F(ZLL) (with correct sign) to start recur
C
   30 W   =  X*X  *(HALF/TPK1 + ONE/TPK1**2) + ETA*(ETA-TWO*X)/TPK1
      FESL   = (ZLL+ONE) * XLOG + CLL - W - LOG(FCL)
   40 FESL = FESL - ABS(SCALE)
          RK   =        MAX(REAL(FESL), FPLMIN*HALF)
          FESL = CMPLX(MIN(RK,   FPLMAX*HALF ) , AIMAG(FESL),kreal)
      FEST= EXP(FESL)
C
           RERR = MAX(RERR, ERR, ACC8 * ABS(REAL(THETA)) )
C
      FCL = FEST
      FPL = FCL*F
      IF(IFCP) FCP(L1) = FPL
               FC (L1) = FCL
C
C *** downward recurrence to lambda = ZLM. array GC,if present,stores RL
C
      I  = MAX(-ID, 0)
      ZL  = ZLL + I
         MONO = 0
        OFF = ABS(FCL)
         TA = ABSC(SIGMA)
      DO 70  L  = L1-ID,LF,-ID
         IF(ETANE0) THEN
               IF(RLEL) THEN
                    DSIG = ATAN2(DRL(ETA),DRL(ZL))
                    RL = SQRT(DRL(ZL)**2 + DRL(ETA)**2)
                  ELSE
                    AA = ZL - ETAI
                    BB = ZL + ETAI
                    IF(ABSC(AA).LT.ACCH.OR.ABSC(BB).LT.ACCH) GOTO 50
                    DSIG = (LOG(AA) - LOG(BB)) * CI*HALF
                    RL = AA * EXP(CI*DSIG)
                 ENDIF
             IF(ABSC(SIGMA).LT.TA*HALF) THEN
C               re-calculate SIGMA because of accumulating roundoffs:
                SL =(CLOGAM(ZL+I-ETAI)-CLOGAM(ZL+I+ETAI))*CI*HALF
                RL = (ZL - ETAI) * EXP(CI*ID*(SIGMA - SL))
                SIGMA = SL
                TA = ZERO
              ELSE
                SIGMA = SIGMA - DSIG * ID
              ENDIF
                TA = MAX(TA, ABSC(SIGMA))
             SL    =  ETA  + ZL*ZL*XI
                PL = ZERO
                IF(ABSC(ZL).GT.ACCH) PL = (SL*SL - RL*RL)/ZL
             FCL1  = (FCL *SL + ID*ZL*FPL)/RL
              SF = ABS(FCL1)
                       IF(SF.GT.FPMAX) GO TO 350
             FPL   = (FPL *SL + ID*PL*FCL)/RL
             IF(MODE .LE. 1) GCP(L+ID)= PL * ID
        ELSE
C                               ETA = 0, including Bessels.  NB RL==SL
           RL = ZL* XI
           FCL1 = FCL * RL + FPL*ID
              SF = ABS(FCL1)
                      IF(SF.GT.FPMAX) GO TO 350
           FPL  =(FCL1* RL - FCL) * ID
        ENDIF
C             IF(ABSC(FCL1).LT.ABSC(FCL)) THEN
              IF(SF.LT.OFF) THEN
                 MONO = MONO + 1
                ELSE
                 MONO = 0
                ENDIF
         FCL   =  FCL1
           OFF = SF
         FC(L) =  FCL
         IF(IFCP) FCP(L)  = FPL
           IF(KFN.EQ.0) SIG(L) = SIGMA
           IF(MODE .LE. 2) GC(L+ID) = RL
      ZL = ZL - ID
      IF(MONO.LT.NDROP) GO TO 70
      IF(AXIAL .OR. REAL(ZLM)*ID.GT.-NDROP.AND..NOT.ETANE0) GO TO 70
         UNSTAB = .TRUE.
C
C ***    take action if cannot or should not recur below this ZL:
   50    ZLM = ZL
         LF = L
            IF(ID.LT.0) GO TO 380
         IF(.NOT.UNSTAB) LF = L + 1
         IF(L+MONO.LT.L1-2 .OR. ID.LT.0 .OR. .NOT.UNSTAB) GO TO 80
C             otherwise, all L values (for stability) should be done
C                        in the reverse direction:
         GO TO 60
   70 CONTINUE
      GO TO 80
   60       ID = -1
            LF = L1
            L1 = M1
            RERR = ACCT
            GO TO 10
   80 IF(FCL .EQ. ZERO) FCL = + ACC8
      F  = FPL/FCL
C
C *** Check, if second time around, that the 'f' values agree!
C
      IF(ID.GT.0) FIRST = F
      IF(DONEM) RERR = MAX(RERR, ABSC(F-FIRST)/ABSC(F))
      IF(DONEM) GO TO 90
C
       NOCF2 = .FALSE.
      THETAM  = X - ETA*(XLOG+TLOG) - ZLM*HPI + SIGMA
C
C *** on left x-plane, determine OMEGA by requiring cut on -x axis
C     on right x-plane, choose OMEGA (using estimate based on THETAM)
C       so H(omega) is smaller and recurs upwards accurately.
C     (x-plane boundary is shifted to give CF2(LH) a chance to converge)
C
                           OMEGA = SIGN(ONE,AIMAG(X)+ACC8)
      IF(REAL(X).GE.XNEAR) OMEGA = SIGN(ONE,AIMAG(THETAM)+ACC8)
      IF(AXIAL)             OMEGA = ONE

C
         SFSH = EXP(OMEGA*SCALE - ABS(SCALE))
         OFF=EXP(MIN(TWO * MAX(ABS(AIMAG(X)),ABS(AIMAG(THETAM)),
     X                         ABS(AIMAG(ZLM))*3 ) , FPLMAX) )
          EPS = MAX(ACC8 , ACCT * HALF / OFF)
C
C ***    Try first estimated omega, then its opposite,
C        to find the H(omega) linearly independent of F
C        i.e. maximise  CF1-CF2 = 1/(F H(omega)) , to minimise H(omega)
C
   90 DO 100 L=1,2
         LH = 1
         IF(OMEGA.LT.ZERO) LH = 2
      PM = CI*OMEGA
      ETAP = ETA * PM
         IF(DONEM) GO TO 130
         PQ1 = ZERO
         PACCQ = ONE
         KASE = 0
C
C ***            Check for small X, i.e. whether to avoid CF2 :
C
      IF(MODE.GE.3 .AND. ABSX.LT.ONE ) GO TO 190
      IF(MODE.LT.3 .AND. (NOCF2 .OR. ABSX.LT.XNEAR .AND.
     X   ABSC(ETA)*ABSX .LT. 5 .AND. ABSC(ZLM).LT.4)) THEN
        KASE = 5
        GO TO 120
        ENDIF
C
C ***  Evaluate   CF2 : PQ1 = p + i.omega.q  at lambda = ZLM
C
         PQ1 = CF2(X,ETA,ZLM,PM,EPS,LIMIT,ERR,NPQ(LH),ACC8,ACCH,
     X             PR,ACCUR,DELL,'COULCC')
C
       ERR = ERR * MAX(ONE,ABSC(PQ1)/MAX(ABSC(F-PQ1),ACC8) )
       IF(ERR.LT.ACCH)       GO TO 110
C
C *** check if impossible to get F-PQ accurately because of cancellation
             NOCF2 = REAL(X).LT.XNEAR .AND. ABS(AIMAG(X)).LT.-LOG(ACC8)
C                original guess for OMEGA (based on THETAM) was wrong
C                Use KASE 5 or 6 if necessary if Re(X) < XNEAR
  100            OMEGA = - OMEGA
                IF(UNSTAB) GO TO 360
                IF(REAL(X).LT.-XNEAR .AND. PR) WRITE(6,1060) '-X',ERR
  110     RERR = MAX(RERR,ERR)
C
C ***  establish case of calculation required for irregular solution
C
  120 IF(KASE.GE.5) GO TO 130
      IF(REAL(X) .GT. XNEAR) THEN
C          estimate errors if KASE 2 or 3 were to be used:
         PACCQ = EPS * OFF * ABSC(PQ1) / MAX(ABS(AIMAG(PQ1)),ACC8)
        ENDIF
      IF(PACCQ .LT. ACCUR) THEN
          KASE = 2
          IF(AXIAL) KASE = 3
      ELSE
          KASE = 1
          IF(NPQ(1) * R20 .LT. JMAX)     KASE = 4
C             i.e. change to kase=4 if the 2F0 predicted to converge
      ENDIF
  130 GO TO (190,140,150,170,190,190),  ABS(KASE)
  140    IF(.NOT.DONEM)
C
C ***  Evaluate   CF2 : PQ2 = p - i.omega.q  at lambda = ZLM   (Kase 2)
C
     X  PQ2 = CF2(X,ETA,ZLM,-PM,EPS,LIMIT,ERR,NPQ(3-LH),ACC8,ACCH,
     X             PR,ACCUR,DELL,'COULCC')
C
        P     = (PQ2 + PQ1) * HALF
        Q     = (PQ2 - PQ1) * HALF*PM
      GO TO 160
  150   P     = REAL(PQ1,kreal)
        Q     = AIMAG(PQ1)
C
C ***   With Kase = 3 on the real axes, P and Q are real & PQ2 = PQ1*
C
        PQ2 = CONJG(PQ1)
C
C *** solve for FCM = F at lambda = ZLM,then find norm factor W=FCM/FCL
C
  160 W   = (PQ1 - F) * (PQ2 - F)
         SF = EXP(-ABS(SCALE))
      FCM = SQRT(Q / W) * SF
C                  any SQRT given here is corrected by
C                  using sign for FCM nearest to phase of FCL
      IF(REAL(FCM/FCL).LT.ZERO) FCM  = - FCM
      GAM = (F - P)/Q
         TA = ABSC(GAM + PM)
         PACCQ= EPS * MAX(TA,ONE/TA)
      HCL = FCM * (GAM + PM) * (SFSH/(SF*SF))
C
      IF(PACCQ.GT.ACCUR .AND. KASE.GT.0) THEN
C                                    Consider a KASE = 1 Calculation
          F11V= F11(X,ETA,Z11,P11,ACCT,LIMIT,0,ERR,N11,FPMAX,ACC8,ACC16)
          IF(ERR.LT.PACCQ) GO TO 200
          ENDIF
      RERR=MAX(RERR,PACCQ)
      GO TO 230
C
C *** Arrive here if KASE = 4
C     to evaluate the exponentially decreasing H(LH) directly.
C
  170  IF(DONEM) GO TO 180
      AA = ETAP - ZLM
      BB = ETAP + ZLM + ONE
      F20V = F20(AA,BB,-HALF*PM*XI, ACCT,JMAX,ERR,FPMAX,N20,XRCF)
        IF(N20.LE.0) GO TO 190
        RERR = MAX(RERR,ERR)
         HCL = FPMIN
         IF(ABS(REAL(PM*THETAM)+OMEGA*SCALE).GT.FPLMAX) GO TO 330
  180 HCL = F20V * EXP(PM * THETAM + OMEGA*SCALE)
      FCM = SFSH / ((F - PQ1) * HCL )
      GO TO 230
C
C *** Arrive here if KASE=1   (or if 2F0 tried mistakenly & failed)
C
C           for small values of X, calculate F(X,SL) directly from 1F1
C               using REAL*16 arithmetic if possible.
C           where Z11 = ZLL if ID>0, or = ZLM if ID<0
C
  190 F11V = F11(X,ETA,Z11,P11,ACCT,LIMIT,0,ERR,N11,FPMAX,ACC8,ACC16)
C
  200       IF(N11.LT.0) THEN
C                               F11 failed from BB = negative integer
               WRITE(6,1060) '-L',ONE
               GO TO 390
               ENDIF
            IF(ERR.GT.PACCQ .AND. PACCQ.LT.ACCB) THEN
C                               Consider a KASE 2 or 3 calculation :
                KASE = -2
                IF(AXIAL) KASE = -3
                GO TO 130
                ENDIF
         RERR = MAX(RERR, ERR)
         IF(ERR.GT.FPMAX) GO TO 370
         IF(ID.LT.0) CLL = Z11*TLOG- HPI*ETA - CLOGAM(BB)
     X                       + CLOGAM(Z11 + ONE + P11*ETA) - P11*SIGMA
      EK   = (Z11+ONE)*XLOG - P11*X + CLL  - ABS(SCALE)
      IF(ID.GT.0) EK = EK - FESL + LOG(FCL)
         IF(REAL(EK).GT.FPLMAX) GO TO 350
         IF(REAL(EK).LT.FPLMIN) GO TO 340
      FCM = F11V * EXP(EK)
C
      IF(KASE.GE.5) THEN
        IF(ABSC(ZLM+ZLM-NINTC(ZLM+ZLM)).LT.ACCH) KASE = 6
C
C ***  For abs(X) < XNEAR, then CF2 may not converge accurately, so
C ***      use an expansion for irregular soln from origin :
C
         SL = ZLM
            ZLNEG = REAL(ZLM) .LT. -ONE + ACCB
         IF(KASE.EQ.5 .OR. ZLNEG) SL = - ZLM - ONE
         PK = SL + ONE
            AA = PK - ETAP
            AB = PK + ETAP
            BB = TWO*PK
                     CLGAA = CLOGAM(AA)
                     CLGAB = CLGAA
         IF(ETANE0)  CLGAB = CLOGAM(AB)
                     CLGBB = CLOGAM(BB)
           IF(KASE.EQ.6 .AND. .NOT.ZLNEG) THEN
              IF(NPINT(AA,ACCUR)) CLGAA = CLGAB - TWO*PM*SIGMA
              IF(NPINT(AB,ACCUR)) CLGAB = CLGAA + TWO*PM*SIGMA
             ENDIF
          CLL = SL*TLOG- HPI*ETA - CLGBB + (CLGAA + CLGAB) * HALF
          DSIG = (CLGAA - CLGAB) * PM*HALF
             IF(KASE.EQ.6) P11 = - PM
          EK  = PK * XLOG - P11*X + CLL  - ABS(SCALE)
                     SF = EXP(-ABS(SCALE))
                     CHI = ZERO
       IF(.NOT.( KASE.EQ.5 .OR. ZLNEG ) ) GO TO 210
C
C *** Use  G(l)  =  (cos(CHI) * F(l) - F(-l-1)) /  sin(CHI)
C
C      where CHI = sig(l) - sig(-l-1) - (2l+1)*pi/2
C
         CHI = SIGMA - DSIG - (ZLM-SL) * HPI
         F11V=F11(X,ETA,SL,P11,ACCT,LIMIT,0,ERR,NPQ(1),FPMAX,ACC8,ACC16)
                    RERR = MAX(RERR,ERR)
            IF(KASE.EQ.6) GO TO 210
         FESL = F11V * EXP( EK )
         FCL1 = EXP(PM*CHI) * FCM
         HCL = FCL1 - FESL
               RERR=MAX(RERR,ACCT*MAX(ABSC(FCL1),ABSC(FESL))/ABSC(HCL))
         HCL = HCL / SIN(CHI) * (SFSH/(SF*SF))
       GO TO 220
C
C *** Use the logarithmic expansion for the irregular solution (KASE 6)
C        for the case that BB is integral so sin(CHI) would be zero.
C
  210    RL = BB - ONE
         N  = NINTC(RL)
         ZLOG = XLOG + TLOG - PM*HPI
         CHI = CHI + PM * THETAM + OMEGA * SCALE + AB * ZLOG
            AA  = ONE - AA
         IF(NPINT(AA,ACCUR)) THEN
            HCL = ZERO
         ELSE
               IF(ID.GT.0 .AND. .NOT.ZLNEG) F11V = FCM * EXP(-EK)
            HCL = EXP(CHI - CLGBB - CLOGAM(AA)) * (-1)**(N+1)
     X              * ( F11V * ZLOG +
     X      F11(X,ETA,SL,-PM,ACCT,LIMIT,2,ERR,NPQ(2),FPMAX,ACC8,ACC16))
                RERR = MAX(RERR,ERR)
            ENDIF
         IF(N.GT.0) THEN
             EK = CHI + CLOGAM(RL) - CLGAB - RL*ZLOG
             DF =F11(X,ETA,-SL-ONE,-PM,ZERO,N,0,ERR,L,FPMAX,ACC8,ACC16)
             HCL = HCL + EXP(EK) * DF
            ENDIF
         RERR = MAX(RERR,TWO*ABS(BB-NINTC(BB)))
C
  220    PQ1 = F - SFSH/(FCM * HCL)
      ELSE
           IF(MODE.LE.2) HCL = SFSH/((F - PQ1) * FCM)
           KASE = 1
      ENDIF
C
C ***  Now have absolute normalisations for Coulomb Functions
C          FCM & HCL  at lambda = ZLM
C      so determine linear transformations for Functions required :
C
  230 IH = ABS(MODE1) / 10
        IF(KFN.EQ.3) IH = (3-AIMAG(CIK))/2  + HALF
      P11 = ONE
      IF(IH.EQ.1) P11 = CI
      IF(IH.EQ.2) P11 = -CI
                  DF = - PM
      IF(IH.GE.1) DF = - PM + P11
          IF(ABSC(DF).LT.ACCH) DF = ZERO
C
C *** Normalisations for spherical or cylindrical Bessel functions
C
                          ALPHA = ZERO
          IF(KFN  .EQ. 1) ALPHA = XI
          IF(KFN  .GE. 2) ALPHA = XI*HALF
                          BETA  = ONE
          IF(KFN  .EQ. 1) BETA  = XI
          IF(KFN  .GE. 2) BETA  = SQRT(XI/HPI)
          IF(KFN  .GE. 2 .AND. REAL(BETA).LT.ZERO) BETA  = - BETA
C
      AA = ONE
      IF(KFN.GT.0) AA = -P11 * BETA
      IF(KFN.GE.3) THEN
C                        Calculate rescaling factors for I & K output
         P = EXP((ZLM+DELL) * HPI * CIK)
         AA= BETA * HPI * P
         BETA = BETA / P
         Q = CIK * ID
        ENDIF
C                        Calculate rescaling factors for GC output
      IF(IH.EQ.0) THEN
         TA = ABS(SCALE) + AIMAG(PM)*SCALE
         RK = ZERO
         IF(TA.LT.FPLMAX) RK = EXP(-TA)
       ELSE
         TA = ABS(SCALE) + AIMAG(P11)*SCALE
C
         IF(ABSC(DF).GT.ACCH .AND. TA.GT.FPLMAX) GO TO 320
         IF(ABSC(DF).GT.ACCH) DF = DF * EXP(TA)
         SF = TWO * (LH-IH) * SCALE
         RK = ZERO
         IF(SF.GT.FPLMAX) GO TO 320
         IF(SF.GT.FPLMIN) RK = EXP(SF)
      ENDIF
C
         KAS((3-ID)/2) = KASE
      W = FCM / FCL
         IF(LOG(ABSC(W))+LOG(ABSC(FC(LF))) .LT. FPLMIN) GO TO 340
         IF(MODE.GE.3) GO TO 240
            IF(ABSC(F-PQ1) .LT. ACCH*ABSC(F) .AND. PR)
     X                             WRITE(6,1020) LH,ZLM+DELL
      HPL = HCL * PQ1
         IF(ABSC(HPL).LT.FPMIN.OR.ABSC(HCL).LT.FPMIN) GO TO 330
C
C *** IDward recurrence from HCL,HPL(LF) (stored GC(L) is RL if reqd)
C *** renormalise FC,FCP at each lambda
C ***    ZL   = ZLM - MIN(ID,0) here
C
  240 DO 270 L = LF,L1,ID
                     FCL = W* FC(L)
                      IF(ABSC(FCL).LT.FPMIN) GO TO 340
            IF(IFCP) FPL = W*FCP(L)
                     FC(L)  = BETA * FCL
            IF(IFCP) FCP(L) = BETA * (FPL - ALPHA * FCL) * CIK
                     FC(L)  = TIDY(FC(L),ACCUR)
            IF(IFCP) FCP(L) = TIDY(FCP(L),ACCUR)
       IF(MODE .GE. 3) GO TO 260
       IF(L.EQ.LF)  GO TO 250
                      ZL = ZL + ID
                      ZID= ZL * ID
                      RL = GC(L)
         IF(ETANE0)   THEN
                      SL = ETA + ZL*ZL*XI
            IF(MODE.EQ.1) THEN
              PL = GCP(L)
            ELSE
              PL = ZERO
              IF(ABSC(ZL).GT.ACCH) PL = (SL*SL - RL*RL)/ZID
            ENDIF
           HCL1     = (SL*HCL - ZID*HPL) / RL
           HPL      = (SL*HPL - PL *HCL) / RL
         ELSE
           HCL1 = RL * HCL - HPL * ID
           HPL  = (HCL - RL * HCL1) * ID
         ENDIF
         HCL      = HCL1
         IF(ABSC(HCL).GT.FPMAX) GO TO 320
  250    GC(L) = AA * (RK * HCL + DF * FCL)
      IF(MODE.EQ.1) GCP(L) = (AA *(RK*HPL +DF*FPL) - ALPHA * GC(L)) *CIK
         GC(L) = TIDY(GC(L),ACCUR)
      IF(MODE.EQ.1) GCP(L) = TIDY(GCP(L),ACCUR)
         IF(KFN.GE.3) AA = AA * Q
  260    IF(KFN.GE.3) BETA = - BETA * Q
  270  LAST = MIN(LAST,(L1 - L)*ID)
C
C *** Come here after all soft errors to determine how many L values ok
C
  280  IF(ID.GT.0 .OR.  LAST.EQ.0) IFAIL = LAST
       IF(ID.LT.0 .AND. LAST.NE.0) IFAIL = -3
C
C *** Come here after ALL errors for this L range (ZLM,ZLL)
C
  290 IF(ID.GT.0 .AND. LF.NE.M1) GO TO 300
         IF(IFAIL.LT.0) RETURN
         IF(RERR.GT.ACCB) WRITE(6,1070) RERR
         IF(RERR.GT.0.1) IFAIL = -4
         RETURN
C
C *** so on first block, 'F' started decreasing monotonically,
C                        or hit bound states for low ZL.
C     thus redo M1 to LF-1 in reverse direction
C      i.e. do CF1A at ZLMIN & CF2 at ZLM (midway between ZLMIN & ZLMAX)
C
  300 ID = -1
      IF(.NOT.UNSTAB) LF = LF - 1
      DONEM = UNSTAB
      LF = MIN(LF,L1)
      L1 = M1
      GO TO 10
C
C ***    error messages
C
  310 IF(PR) WRITE (6,1000) XX
 1000 FORMAT(/' COULCC: CANNOT CALCULATE IRREGULAR SOLUTIONS FOR X =',
     X 1P,2D10.2,', AS ABS(X) IS TOO SMALL'/)
      RETURN
  320 IF(PR) WRITE(6,1010) ZL+DELL,'IR',HCL,'MORE',FPMAX
 1010 FORMAT(' COULCC: AT ZL =',2F8.3,' ',A2,'REGULAR SOLUTION (',1P,
     X 2E10.1,') WILL BE ',A4,' THAN',E10.1)
      GO TO 280
  330 IF(PR) WRITE(6,1010) ZL+DELL,'IR',HCL,'LESS',FPMIN
      GO TO 280
  340 IF(PR) WRITE(6,1010) ZL+DELL,'  ',FCL,'LESS',FPMIN
      GO TO 280
  350 IF(PR) WRITE(6,1010) ZL+DELL,'  ',FCL,'MORE',FPMAX
      GO TO 280
 1020 FORMAT('0COULCC WARNING: LINEAR INDEPENDENCE BETWEEN ''F'' AND ''H
     X(',I1,')'' IS LOST AT ZL =',2F7.2,' (EG. COULOMB EIGENSTATE, OR CF
     X1 UNSTABLE)'/)
  360 IF(PR) WRITE(6,1030) ZLL+DELL
 1030 FORMAT(' COULCC: (ETA&L)/X TOO LARGE FOR CF1A, AND CF1 UNSTABLE AT
     X L =',2F8.2)
      GO TO 280
  370 IF(PR) WRITE(6,1040) Z11,I
 1040 FORMAT(' COULCC: OVERFLOW IN 1F1 SERIES AT ZL =',2F8.3,' AT TERM',
     X I5)
      GO TO 390
  380 IF(PR) WRITE(6,1050) ZLMIN,ZLM,ZLM+ONE,ZLMIN+NL-ONE
 1050 FORMAT(' COULCC: BOTH BOUND-STATE POLES AND F-INSTABILITIES OCCUR'
     X ,', OR MULTIPLE INSTABILITIES PRESENT.'
     X,/,' TRY CALLING TWICE,  FIRST FOR ZL FROM',2F8.3,' TO',2F8.3,
     X ' (INCL.)',/,20X,     'SECOND FOR ZL FROM',2F8.3,' TO',2F8.3)
C     GO TO 390
  390 IFAIL = -1
      GO TO 290
 1060 FORMAT('0COULCC WARNING: AS ''',A2,''' REFLECTION RULES NOT USED,
     XERRORS CAN BE UP TO',1P,D12.2/)
 1070 FORMAT('0COULCC WARNING: OVERALL ROUNDOFF ERROR APPROX.',1P,E11.1)
      END
      FUNCTION CF2(X,ETA,ZL,PM,EPS,LIMIT,ERR,NPQ,ACC8,ACCH,
     X             PR,ACCUR,DELL,CALLER)
      parameter(kreal=selected_real_kind(12,100))
      IMPLICIT COMPLEX(kind=kreal)(A-H,O-Z)
      LOGICAL PR
      REAL(kind=kreal) EPS,ERR,ACC8,ACCH,ACCUR,TA,RK,
     X       ABSC,ZERO,HALF,ONE,TWO
      CHARACTER*6 CALLER
      DATA ZERO,HALF,ONE,TWO / 0D+0, .5D+0, 1D+0, 2D+0 /
      ABSC(W) = ABS(REAL(W)) + ABS(AIMAG(W))
C
C                                    (omega)        (omega)
C *** Evaluate  CF2  = p + PM.q  =  H   (ETA,X)' / H   (ETA,X)
C                                    ZL             ZL
C     where PM = omega.i
C
      TA = TWO*LIMIT
      E2MM1 = ETA*ETA + ZL*ZL + ZL
      ETAP = ETA * PM
      XI = ONE/X
      WI = TWO*ETAP
      RK = ZERO
      PQ = (ONE - ETA*XI) * PM
      AA = -E2MM1 + ETAP
      BB = TWO*(X - ETA + PM)
         RL = XI * PM
      IF(ABSC(BB).LT.ACCH) THEN
         RL = RL * AA / (AA + RK + WI)
         PQ = PQ + RL * (BB + TWO*PM)
            AA = AA + TWO*(RK+ONE+WI)
            BB = BB + (TWO+TWO)*PM
            RK = RK + (TWO+TWO)
         ENDIF
      DD = ONE/BB
      DL = AA*DD* RL
   10 PQ    = PQ + DL
         RK = RK + TWO
         AA = AA + RK + WI
         BB = BB + TWO*PM
         DD = ONE/(AA*DD + BB)
         DL = DL*(BB*DD - ONE)
            ERR = ABSC(DL)/ABSC(PQ)
         IF(ERR.GE.MAX(EPS,ACC8*RK*HALF) .AND. RK.LE.TA) GO TO 10
C
         NPQ   = RK/TWO
         PQ    = PQ + DL
           IF(PR.AND.NPQ.GE.LIMIT-1 .AND. ERR.GT.ACCUR)
     X             WRITE(6,1000) CALLER,INT(AIMAG(PM)),NPQ,ERR,ZL+DELL
 1000 FORMAT(' ',A6,': CF2(',I2,') NOT CONVERGED FULLY IN ',I7,
     X' ITERATIONS, SO ERROR IN IRREGULAR SOLUTION =',1P,D11.2,' AT ZL
     X=', 0P,2F8.3)
      CF2 = PQ
      RETURN
      END
      FUNCTION F11(X,ETA,ZL,P,EPS,LIMIT,KIN,ERR,NITS,FPMAX,ACC8,ACC16)
      parameter(kreal=selected_real_kind(12,100))
!      parameter(kqreal=selected_real_kind(24,100))
      parameter(kqreal=selected_real_kind(12,100))
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      COMPLEX(kind=kreal) X,ETA,ZL,P,AA,BB,Z,F11,CDIGAM,CI
      COMPLEX(kind=kreal) DD,G,F,AI,BI,T
      LOGICAL ZLLIN
      REAL(kind=kqreal) AR,BR,GR,GI,DR,DI,TR,TI,UR,UI,FI,FI1,DEN,RZ
      DATA ZERO,ONE,TWO / 0D+0, 1D+0, 2D+0 /, CI / (0D+0, 1D+0) /
      ABSC(AA) = ABS(REAL(AA)) + ABS(AIMAG(AA))
      NINTC(AA) = NINT(REAL(AA))
   5  FORMAT(A, I12)
   7  FORMAT(A, 4G12.2)
C      PRINT 7, "X=     ", X, ETA
C      PRINT 7, "ETA=   ", ETA
C      PRINT 7, "ZL=    ", ZL
C      PRINT 7, "P=     ", P
C      PRINT 7, "EPS=   ", EPS
C      PRINT 5, "LIMIT= ", LIMIT
C      PRINT 5, "KIN=   ", KIN
C      PRINT 7, "ERR=   ", ERR
C      PRINT 5, "NITS=  ", NITS
C      PRINT 7, "FPMAX= ", FPMAX
C      PRINT 7, "ACC8=  ", ACC8
C      PRINT 7, "ACC16= ", ACC16
C      PRINT 7, "ONE=   ", ONE
C      PRINT 7, "TWO=   ", TWO
C      PRINT 7, "************"

C
C *** evaluate the HYPERGEOMETRIC FUNCTION 1F1
C                                        i
C            F (AA;BB; Z) = SUM  (AA)   Z / ( (BB)  i! )
C           1 1              i       i            i
C
C     to accuracy EPS with at most LIMIT terms.
C  If KIND = 0 : using extended precision but real arithmetic only,
C            1 : using normal precision in complex arithmetic,
C   or       2 : using normal complex arithmetic, but with CDIGAM factor
C
C  where
         AA = ZL+ONE - ETA*P
         BB = TWO*(ZL+ONE)
C  and
         Z  = TWO*P*X

         ACC16 = epsilon(AR)
	 KIND = KIN
	 if(KIND==0.and.(kreal==kqreal.or.ACC16>ACC8*.1)) KIND=1
C
         ZLLIN = REAL(BB).LE.ZERO .AND. ABS(BB-NINTC(BB)).LT.ACC8**0.75
             IF(.NOT.ZLLIN.OR.REAL(BB)+LIMIT.LT.1.5) GO TO 10
                NITS = -1
                RETURN
   10 IF(LIMIT.LE.0) THEN
         F11 = ZERO
         ERR = ZERO
         NITS= 1
         RETURN
         ENDIF
      TA = ONE
      RK = ONE
      IF(KIND.LE.0.AND.ABSC(Z)*ABSC(AA).GT.ABSC(BB) * 1.0) THEN
         DR = ONE
         DI = ZERO
         GR = ONE
         GI = ZERO
         AR = AA   ! real part
         BR = BB   ! real part
         RZ = Z    ! real part
         FI = ZERO
      DO 20 I=2,LIMIT
         FI1 = FI + ONE
         TR = BR * FI1
         TI = AIMAG(BB) * FI1
         DEN= ONE / (TR*TR + TI*TI)
         UR = (AR*TR + AIMAG(AA)*TI) * DEN
         UI = (AIMAG(AA)*TR - AR*TI) * DEN
         TR = UR*GR - UI*GI
         TI = UR*GI + UI*GR
         GR = RZ * TR - AIMAG(Z)*TI
         GI = RZ * TI + AIMAG(Z)*TR
         DR = DR + GR
         DI = DI + GI
            ERR = ABS(GR) + ABS(GI)
               IF(ERR.GT.FPMAX) GO TO 60
            RK  = ABS(DR) + ABS(DI)
            TA = MAX(TA,RK)
         IF(ERR.LT.RK*EPS .OR. I.GE.4.AND.ERR.LT.ACC16) GO TO 30
         FI = FI1
         AR = AR + ONE
   20    BR = BR + ONE
C
   30    F11 = DR + CI * DI
         ERR = ACC16 * TA / RK
C
      ELSE
C* ---------------------------------- alternative code
C*    If REAL*16 arithmetic is not available, (or already using it!),
C*    then use KIND > 0
         G = ONE
          F = ONE
          IF(KIND.GE.2) F = CDIGAM(AA) - CDIGAM(BB) - CDIGAM(G)
         DD = F
         DO 40 I=2,LIMIT
            AI = AA + (I-2)
            BI = BB + (I-2)
            R  = I-ONE
         G = G * Z * AI / (BI * R)
         IF(KIND.GE.2)
C                              multiply by (psi(a+r)-psi(b+r)-psi(1+r))
     X        F = F + ONE/AI - ONE/BI - ONE/R
         T  = G * F
         DD = DD + T
            ERR = ABSC(T)
               IF(ERR.GT.FPMAX) GO TO 60
            RK = ABSC(DD)
         TA = MAX(TA,RK)
         IF(ERR.LT.RK*EPS.OR.ERR.LT.ACC8.AND.I.GE.4) GO TO 50
   40    CONTINUE
 
   50    ERR = ACC8 * TA / RK
         F11 = DD
C* ------------------------------------------- end of alternative code
      ENDIF
   60    NITS = I
      RETURN
      END
      FUNCTION CF1R(X,ETA,ZL,EPS,FCL,TPK1,ETANE0,LIMIT,ERR,NFP,
     X              ACCH,FPMIN,FPMAX,PR,CALLER)
      parameter(kreal=selected_real_kind(12,100))
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      PARAMETER(ZERO=0D0, ONE=1D0, TWO=2D0)
      LOGICAL PR,ETANE0
      CHARACTER*6 CALLER
C
C
C ***    Evaluate CF1  =  F   =  F'(ZL,ETA,X)/F(ZL,ETA,X)
C
C        using real arithmetic
C
      FCL = ONE
      XI = ONE/X
      PK  = ZL + ONE
      PX  = PK  + LIMIT
      EK  = ETA/PK
      F   =  EK + PK*XI
         IF(ABS(F).LT.FPMIN) F = FPMIN
         D = ZERO
         C = F
           SMALL    = SQRT(FPMIN)
           RK2 = ONE + EK*EK
C
C ***   begin CF1 loop on PK = k = lambda + 1
C
   10 PK1 = PK + ONE
         TPK1 = PK + PK1
         IF(ETANE0) THEN
           EK  = ETA / PK
           RK2 = ONE + EK*EK
           TK  = TPK1*(XI + EK/PK1)
         ELSE
           TK  = TPK1*XI
         ENDIF
        C  =  TK - RK2 / C
        D  =  TK - RK2 * D
         IF(ABS(C).LT.FPMIN) C = FPMIN
         IF(ABS(D).LT.FPMIN) D = FPMIN
         D = ONE/D
         DF = D * C
         F  = F * DF
            FCL = FCL * D * TPK1*XI
            IF(ABS(FCL).LT.SMALL) FCL = FCL / SMALL
            IF(ABS(FCL).GT.FPMAX) FCL = FCL * FPMIN
         PK = PK1
                          IF( PK .GT. PX ) GO TO 50
      IF(ABS(DF-ONE) .GE. EPS)             GO TO 10
                NFP = PK - ZL - 1
                ERR = EPS * SQRT(REAL(NFP))
      CF1R = F
      RETURN
   50 IF(PR) WRITE (6,1010) CALLER,LIMIT,ABS(X)
 1010 FORMAT(' ',A6,': CF1 HAS FAILED TO CONVERGE AFTER ',I10  ,' ITERAT
     XIONS AS ABS(X) =',F15.0)
      ERR = TWO
      RETURN
      END
      FUNCTION CF1C(X,ETA,ZL,EPS,FCL,TPK1,ETANE0,LIMIT,ERR,NFP,
     X              ACCH,FPMIN,FPMAX,PR,CALLER)
      parameter(kreal=selected_real_kind(12,100))
      IMPLICIT COMPLEX(kind=kreal)(A-H,O-Z)
      LOGICAL PR,ETANE0
      REAL(kind=kreal) EPS,ERR,ACCH,FPMIN,FPMAX,ABSC,SMALL,PX, 
     X                 ZERO,ONE,TWO
      PARAMETER(ZERO=0D0, ONE=1D0, TWO=2D0)
      CHARACTER*6 CALLER
      ABSC(W) = ABS(REAL(W)) + ABS(AIMAG(W))
C
C
C ***    Evaluate CF1  =  F   =  F'(ZL,ETA,X)/F(ZL,ETA,X)
C
C        using complex arithmetic
C
      FCL = ONE
      XI = ONE/X
      PK  = ZL + ONE
      PX  = PK  + LIMIT
      EK  = ETA/PK
      F   =  EK + PK*XI
         IF(ABSC(F).LT.FPMIN) F = FPMIN
         D = ZERO
         C = F
           SMALL    = SQRT(FPMIN)
           RK2 = ONE + EK*EK
C
C ***   begin CF1 loop on PK = k = lambda + 1
C
   10 PK1 = PK + ONE
         TPK1 = PK + PK1
         IF(ETANE0) THEN
           EK  = ETA / PK
           RK2 = ONE + EK*EK
           TK  = TPK1*(XI + EK/PK1)
         ELSE
           TK  = TPK1*XI
         ENDIF
        C  =  TK - RK2 / C
        D  =  TK - RK2 * D
         IF(ABSC(C).LT.FPMIN) C = FPMIN
         IF(ABSC(D).LT.FPMIN) D = FPMIN
         D = ONE/D
         DF = D * C
         F  = F * DF
            FCL = FCL * D * TPK1*XI
            IF(ABSC(FCL).LT.SMALL) FCL = FCL / SMALL
            IF(ABSC(FCL).GT.FPMAX) FCL = FCL * FPMIN
         PK = PK1
                       IF(REAL(PK).GT. PX ) GO TO 50
      IF(ABSC(DF-ONE) .GE. EPS)             GO TO 10
                NFP = PK - ZL - 1
                ERR = EPS * SQRT(REAL(NFP))
      CF1C = F
      RETURN
   50 IF(PR) WRITE (6,1010) CALLER,LIMIT,ABS(X)
 1010 FORMAT(' ',A6,': CF1 HAS FAILED TO CONVERGE AFTER ',I10  ,' ITERAT
     XIONS AS ABS(X) =',F15.0)
      ERR = TWO
      RETURN
      END
      FUNCTION F20(AA,BB,Z,EPS,JMAX,RE,FPMAX,N,X)
C
C     evaluate the HYPERGEOMETRIC FUNCTION 2F0
C                                             i
C            F (AA,BB;;Z) = SUM  (AA)  (BB)  Z / i!
C           2 0              i       i     i
C
C     to accuracy EPS with at most JMAX terms.
C
C     if the terms start diverging,
C     the corresponding continued fraction is found by RCF
C     & evaluated progressively by Steed's method to obtain convergence.
C
C      useful number also input:  FPMAX = near-largest f.p. number
C
      parameter(kreal=selected_real_kind(12,100))
      IMPLICIT COMPLEX(kind=kreal)(A-H,O-Z)
      DIMENSION X(JMAX,4)
      LOGICAL FINITE
      REAL(kind=kreal) EP,EPS,AT,ATL,ABSC,RE,FPMAX
      DATA ONE,ZERO / (1D+0,0D+0), (0D+0,0D+0) /
      ABSC(W) = ABS(REAL(W)) + ABS(AIMAG(W))
      NINTC(W) = NINT(REAL(REAL(W)))
C
      RE = 0.0
      X(1,1) = ONE
      SUM = X(1,1)
      ATL = ABSC(X(1,1))
         F    = SUM
         D = ONE
         DF   = SUM
      J = 0
      EP = EPS * JMAX *10.
      MA = - NINTC(AA)
      MB = - NINTC(BB)
      FINITE = ABS(ABS(REAL(AA))-MA).LT.EP .AND. ABS(AIMAG(AA)).LT.EP
     X    .OR. ABS(ABS(REAL(BB))-MB).LT.EP .AND. ABS(AIMAG(BB)).LT.EP
      IMAX = JMAX
      IF(FINITE.AND.MA.GE.0) IMAX = MIN(MA+1,IMAX)
      IF(FINITE.AND.MB.GE.0) IMAX = MIN(MB+1,IMAX)
      DO 10 I=2,IMAX
      X(I,1) = X(I-1,1) * Z * (AA+I-2) * (BB+I-2) / (I-1)
         IF(ABSC(X(I,1)).GT.FPMAX) GO TO 40
      AT = ABSC(X(I,1))
         IF(J.EQ.0) THEN
                 SUM = SUM + X(I,1)
                 IF(AT .LT. ABSC(SUM)*EPS) GO TO 20
               ENDIF
      IF(FINITE) GO TO 10
      IF(J.GT.0 .OR. AT.GT.ATL .OR. I.GE.JMAX-2) J = J + 1
         IF(J.EQ.0) GO TO 10
	 II = I
         CALL RCF(X(1,1),X(1,2),J,II,X(1,3),EPS)
              IF(II.LT.0) GO TO 40
            DO 50 K=MAX(J,2),I
            D = ONE/(D*X(K,2) + ONE)
            DF = DF*(D - ONE)
            F = F + DF
            IF(ABSC(DF) .LT. ABSC(F)*EPS) GO TO 30
            IF(DF.EQ.ZERO.AND.F.EQ.ZERO.AND.I.GE.4) GO TO 30
   50       CONTINUE
         J = I
   10 ATL = AT
      IF(.NOT.FINITE) I = -JMAX
   20 N = I
       F20 = SUM
       IF(.NOT.FINITE) RE  = AT / ABSC(SUM)
       RETURN
   30 F20 = F
      RE = ABSC(DF) / ABSC(F)
      N = K
      RETURN
   40 I = 0
      GO TO 20
      END
      FUNCTION CF1A(RHO,ETA,XL,PSI,EPS,NMAX,NUSED,FCL,RE,FPMAX,XX,G,C)
C
C     evaluate the ASYMPTOTIC EXPANSION for the
C            LOGARITHMIC DERIVATIVE OF THE REGULAR SOLUTION
C
C ***        CF1A  =  f   =  F'(XL,ETA,RHO)/F(XL,ETA,RHO)
C
C      that is valid for REAL(RHO)>0, and best for RHO >> ETA**2, XL,
C      and is derived from the 2F0 expansions for H+ and H-
C      e.g. by Froeberg (Rev. Mod. Physics Vol 27, p399 , 1955)
C      Some lines of this subprogram are for convenience copied from
C           Takemasa, Tamura & Wolter CPC 17 (1979) 351.
C
C     Evaluate to accuracy EPS with at most NMAX terms.
C
C     If the terms start diverging,
C     the corresponding continued fraction is found by RCF
C     & evaluated progressively by Steed's method to obtain convergence.
C
C      useful number also input:  FPMAX = near-largest f.p. number
C
      parameter(kreal=selected_real_kind(12,100))
      IMPLICIT COMPLEX(kind=kreal)(A-H,O-Z)
      DIMENSION XX(2,NMAX),G(NMAX),C(NMAX)
      REAL(kind=kreal) RE,EPS,T1,T2,T3,ZERO,ONE,TWO,AT,ATL,ABSC,FPMAX
      DATA ZERO,ONE,TWO,CI / 0D+0, 1D+0, 2D+0, (0D+0,1D+0) /
      ABSC(W) = ABS(REAL(W)) + ABS(AIMAG(W))
C
      HPI = TWO*ATAN(ONE)
      T1 = SIN(REAL(PSI,kreal))
      T2 = COS(REAL(PSI,kreal))
      ATL= TANH(AIMAG(PSI))
C             GIVE COS(PSI)/COSH(IM(PSI)), WHICH ALWAYS HAS CORRECT SIGN
          COSL = CMPLX( T2 , -T1 * ATL ,kreal)
      TANL = CMPLX(T1,T2*ATL,kreal) / COSL
      RE = ZERO
      XLL1= XL*(XL+ONE)
      ETASQ = ETA*ETA
      SL1=ONE
      SL=SL1
      SC1=ZERO
      SC=SC1
      TL1=SC
      TL=TL1
      TC1=ONE-ETA/RHO
      TC=TC1
      FCL  = TL + SL*TANL
      G(1) = (TC + SC*TANL) / FCL
      GLAST = G(1)
      ATL = ABSC(GLAST)
         F    = GLAST
         D = ONE
         DF   = GLAST
      J = 0
      DO 10 N=2,NMAX
      T1=N-1
      T2=TWO*T1-ONE
      T3=T1*(T1-ONE)
      DENOM=TWO*RHO*T1
      C1=(ETA*T2)/DENOM
      C2=(ETASQ+XLL1-T3)/DENOM
      SL2=C1*SL1-C2*TL1
      TL2=C1*TL1+C2*SL1
      SC2=C1*SC1-C2*TC1-SL2/RHO
      TC2=C1*TC1+C2*SC1-TL2/RHO
      SL=SL+SL2
      TL=TL+TL2
      SC=SC+SC2
      TC=TC+TC2
      SL1=SL2
      TL1=TL2
      SC1=SC2
      TC1=TC2
      FCL  =  TL + SL*TANL
         IF(ABSC(FCL).GT.FPMAX .OR. ABSC(FCL).LT.1./FPMAX) GO TO 40
      GSUM = (TC + SC*TANL) / FCL
      G(N) = GSUM - GLAST
      GLAST = GSUM
         AT = ABSC(G(N))
         IF(AT.LT.ABSC(GSUM)*EPS) GO TO 20
      IF(J.GT.0 .OR. AT.GT.ATL .OR. N.GE.NMAX-2) J = J + 1
         IF(J.EQ.0) GO TO 10
              NN = N
            CALL RCF(G,C,J,NN,XX,EPS)
              IF(NN.LT.0) GO TO 40
            DO 60 K=MAX(J,2),N
               D = ONE/(D*C(K) + ONE)
               DF = DF*(D - ONE)
               F = F + DF
         IF(ABSC(DF) .LT. ABSC(F)*EPS) GO TO 30
         IF(DF.EQ.ZERO.AND.F.EQ.ZERO.AND.N.GE.4) GO TO 30
   60         CONTINUE
         J = N
   10    ATL = AT
      K = -NMAX
      GO TO 30
   20 FCL = FCL * COSL
         CF1A = GSUM
         RE = AT / ABSC(GSUM)
         NUSED = N
         RETURN
   30 CF1A = F
      FCL = FCL * COSL
         RE = ABSC(DF) / ABSC(F)
         NUSED = K
      RETURN
   40 CF1A = G(1)
      FCL = 1.0
      RE = 1.0
      NUSED = 0
      RETURN
      END

      MODULE RCFCM
      COMPLEX(kind=selected_real_kind(12,100)) X0,X1
      LOGICAL EVEN
      INTEGER M2M1,MP12,M
      END MODULE RCFCM
      
      SUBROUTINE RCF(A,B,IBEG,INUM,XX,EPS)
C
C*******************************************************************
C
C  RCF converts polynomial A to the corresponding continued
C         fraction, in 'normal'  form with coefficients B
C         by the 'P algorithmn' of Patry & Gupta
C
C   A(z) = A1/z + A2/z**3 + A3/z**5 + ... + An/z**(2n-1)
C
C   B(z) = B1/z+ B2/z+ B3/z+ .../(z+ Bn/z)
C
C  data:
C   A     vector A(k), k=1,INUM         input
C   B     vector B(k), k=IBEG,INUM      output
C   IBEG  order of first coef. calc.    input
C   INUM  order of A, even or odd       input
C   XX    auxiliary vector of length .ge. length of vector B
C         caller provides space for A,B,XX
C     Note that neither of the first two terms A(1) A(2) should be zero
C             & the user can start the calculation with any value of
C                IBEG provided the c.f. coefs have been already
C                calculated up to INUM = IBEG-1
C             & the method breaks down as soon as the absolute value
C                of a c.f. coef. is less than EPS.    At the time of the
C                break up XX(1) has been replaced by 1E-50, and INUM has
C                been replaced by minus times the number of this coef.
C   algorithm: J.Patry & S.Gupta,
C              EIR-bericht nr. 247,
C              Eidg. Institut fur Reaktorforschung Wuerenlingen
C              Wueringlingen, Schweiz.
C              November 1973
C   see also:  Haenggi,Roesel & Trautmann,
C              Jnl. Computational Physics, vol 137, pp242-258 (1980)
C   note:      restart procedure modified by I.J.Thompson
C
C*******************************************************************
C
      Use RCFCM
      parameter(kreal=selected_real_kind(12,100))
      IMPLICIT COMPLEX(kind=kreal)(A-H,O-Z)
      DIMENSION A(100),B(100),XX(2,100)
      REAL(kind=kreal) EPS

C     ibn = ibeg + inum - 1
      IBN = INUM
C                             B(IBN) is last value set on this call
      IF(IBEG.GT.4 .AND. M .NE. IBEG-1) GO TO 90
C                             B(M) is last value set in previous call
      IF(IBEG.GT.4) GO TO 50
      IF(IBEG.EQ.4) GO TO 20
      B(1) = A(1)
      IF(IBN.GE.2) B(2) = - A(2)/A(1)
      IF(IBN.LT.3) GO TO 10
      X0 = A(3) / A(2)
      XX(2,1) = B(2)
      XX(1,1) = - X0
      XX(1,2) = 0.
      B(3) = -X0 - B(2)
      X0 = -B(3) * A(2)
      M = 3
      MP12 = 2
      EVEN = .TRUE.
      IF(IBN.GT.3) GO TO 20
   10 RETURN
   20 IF(ABS(B(3)) .LT. EPS*ABS(X0)) GOTO 80
      M = 4
   30 X1 = A(M)
      M2M1 = MP12
      MP12 = M2M1 + 1
      IF(EVEN) MP12 = M2M1
      DO 40 K=2,MP12
   40 X1 = X1 + A(M-K+1) * XX(1,K-1)
      B(M) = - X1/X0
      IF(M.GE.IBN) RETURN
   50 IF(ABS(B(M)).LT.EPS*ABS(X0)) GO TO 80
      K = M2M1
   60 XX(2,K) = XX(1,K) + B(M) * XX(2,K-1)
      K = K-1
      IF(K.GT.1) GO TO 60
      XX(2,1) = XX(1,1) + B(M)
      DO 70 K=1,M2M1
      X0 = XX(2,K)
      XX(2,K) = XX(1,K)
   70 XX(1,K) = X0
      X0 = X1
      XX(1,M2M1+1) = 0.
      M = M+1
      EVEN = .NOT.EVEN
      GO TO 30
   80 INUM = -M
C     XX(1,1) = 1.E-50
C     PRINT 1000,M
C1000 FORMAT('0RCF: ZERO CF COEFFICIENT AT POSITION ',I4/)
      RETURN
   90 PRINT 1000,M,IBEG-1
 1000 FORMAT('0RCF: LAST CALL SET M =',I4,', BUT RESTART REQUIRES',I4)
      STOP
      END
      
      MODULE LOGAMMA
      parameter(kreal=selected_real_kind(12,100))
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      REAL(kind=kreal) X0,PI,ALPI,HL2P,ACCUR,B(15)
      REAL(kind=kreal) FPLMIN
      INTEGER NX0,NT
      DATA LERR /6/, NX0 /6/, NB /15/,
     X  ZERO,ONE,TWO,FOUR,HALF,QUART /0D0,1D0,2D0,4D0,.5D0,.25D0/
      END MODULE LOGAMMA
      
      
      
      FUNCTION CLOGAM(Z)
C
C     this routine computes the logarithm of the gamma function gamma(z)
C     for any complex argument 'Z' 
C     to any accuracy preset by a GAMMA0 function call.
C
      use LOGAMMA
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      COMPLEX(kind=kreal) Z,U,V,H,R,CLOGAM,SER
C
      SAVE
      X=REAL(Z,kreal)
      T=AIMAG(Z)
      MX = INT(REAL(ACCUR*100 - X))
      IF(ABS(ABS(X)-MX) + ABS(T).LT.ACCUR*50) GO TO 60
      F=ABS(T)
      V=CMPLX(X,F,kreal)
      IF(X .LT. ZERO) V=ONE-V
      H=ZERO
      C=REAL(V,kreal)
      N=NX0-INT(C)
      IF(N .LT. 0) GO TO 30
      H=V
      D=AIMAG(V)
      A=ATAN2(D,C)
      IF(N .EQ. 0) GO TO 20
      DO 10 I = 1,N
      C=C+ONE
      V=CMPLX(C,D,kreal)
      H=H*V
   10 A=A+ATAN2(D,C)
   20 H=CMPLX(HALF*LOG(REAL(H,kreal)**2+AIMAG(H)**2),A,kreal)
      V=V+ONE
   30 R=ONE/V**2
      SER = B(NT)
      DO 40 J=2,NT
        K = NT+1 - J
   40 SER = B(K) + R*SER
      CLOGAM = HL2P+(V-HALF)*LOG(V)-V + SER/V - H
      IF(X .GE. ZERO) GO TO 50
C
      A= INT(X)-ONE
      C=PI*(X-A)
      D=PI*F
C     E=EXP(-TWO*D)
        E = ZERO
        F = -TWO*D
        IF(F.GT.FPLMIN) E = EXP(F)
      F=SIN(C)
      E= D + HALF*LOG(E*F**2+QUART*(ONE-E)**2)
      F=ATAN2(COS(C)*TANH(D),F)-A*PI
      CLOGAM=ALPI-CMPLX(E,F,kreal)-CLOGAM
C
   50 IF(SIGN(ONE,T) .LT. -HALF) CLOGAM=CONJG(CLOGAM)
      RETURN
C
   60 WRITE(LERR,1000) X
 1000 FORMAT(1X,'CLOGAM: ARGUMENT IS NON POSITIVE INTEGER = ',F20.2)
      CLOGAM = ZERO
      RETURN
      END
C
      FUNCTION CDIGAM(Z)
C
C     this routine computes the logarithmic derivative of the gamma
C     function  psi(Z) = digamma(Z) = d (ln gamma(Z))/dZ  for any
C     complex argument Z, to any accuracy preset by CALL LOGAM(ACC)
C
      use LOGAMMA
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      COMPLEX(kind=kreal) Z,U,V,H,R,CDIGAM,SER
      U=Z
      X=REAL(U,kreal)
      A=ABS(X)
      IF(ABS(AIMAG(U)) + ABS(A + INT(X)) .LT. ACCUR) GO TO 110
      IF(X .LT. ZERO) U=-U
      V=U
      H=ZERO
      N=NX0-INT(A)
      IF(N .LT. 0) GO TO 90
      H=ONE/V
      IF(N .EQ. 0) GO TO 80
      DO 70 I = 1,N
      V=V+ONE
   70 H=H+ONE/V
   80 V=V+ONE
   90 R=ONE/V**2
      SER = B(NT) * (2*NT-1)
      DO 100 J=2,NT
        K = NT+1 - J
  100 SER = B(K)*(2*K-1) + R*SER
      CDIGAM = LOG(V) - HALF/V - R*SER - H
      IF(X .GE. ZERO) RETURN
      H=PI*U
      CDIGAM = CDIGAM + ONE/U + PI*COS(H)/SIN(H)
      RETURN
C
  110 WRITE(LERR,1000) X
 1000 FORMAT(1X,'CDIGAM: ARGUMENT IS NON POSITIVE INTEGER = ',F20.2)
      CDIGAM=ZERO
      RETURN
      END
C
      SUBROUTINE GAMMA0(ACC,FPLMINI)
C
C      initialisation call for calculations to accuracy 'ACC'
C
      use LOGAMMA
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      REAL(kind=kreal) BN(15),BD(15)
      DATA BN(1),BD(1)    / +1D+0,   6D+0 /,
     X     BN(2),BD(2)    / -1D+0,  30D+0 /,
     X     BN(3),BD(3)    / +1D+0,  42D+0 /,
     X     BN(4),BD(4)    / -1D+0,  30D+0 /,
     X     BN(5),BD(5)    / +5D+0,  66D+0 /,
     X     BN(6),BD(6)    /          -691D+0,  2730D+0/,
     X     BN(7),BD(7)    /          +  7D+0,     6D+0/,
     X     BN(8),BD(8)    /         -3617D+0,   510D+0/,
     X     BN(9),BD(9)    /         43867D+0,   798D+0/,
     X     BN(10),BD(10)  /       -174611D+0,   330D+0/,
     X     BN(11),BD(11)  /        854513D+0,   138D+0/,
     X     BN(12),BD(12)  /    -236364091D+0,  2730D+0/,
     X     BN(13),BD(13)  /     + 8553103D+0,     6D+0/,
     X     BN(14),BD(14)  /  -23749461029D+0,   870D+0/,
     X     BN(15),BD(15)  / 8615841276005D+0, 14322D+0/

      FPLMIN = FPLMINI
      NX0 = 6
      X0  = NX0 + ONE
      PI = FOUR*ATAN(ONE)
      ALPI = LOG(PI)
      HL2P = LOG(TWO*PI) * HALF
      ACCUR = ACC
      DO 120 K=1,NB
       F21 = K*2 - ONE
       B(K) = BN(K) / (BD(K) * K*TWO * F21)
       ERR = ABS(B(K)) * K*TWO / X0**F21
  120 IF(ERR.LT.ACC) GO TO 130
       NX0 = INT((ERR/ACC)**(ONE/F21) * X0)
       K = NB
  130 NT = K
C     print *,' logam requires k = ',k ,' with cutoff at x =',nx0+1
      RETURN
      END
      FUNCTION TIDY(Z,ACC)
C                     TIDY A COMPLEX NUMBER
      parameter(kreal=selected_real_kind(12,100))
      REAL(kind=kreal) X,Y,ACC,AZ
      COMPLEX(kind=kreal) Z,TIDY
C
      X = REAL(Z,kreal)
      Y = AIMAG(Z)
      AZ= (ABS(X) + ABS(Y)) * ACC * 5
      IF(ABS(X) .LT. AZ) X = 0D+0
      IF(ABS(Y) .LT. AZ) Y = 0D+0
      TIDY = CMPLX(X,Y,kreal)
      RETURN
      END

C     The unadulterated confluent hypergeometric function
C     
      FUNCTION HYPERGF(AA,BB,X,EPS,LIMIT,KIN,ERR,NITS,FPMAX,ACC8,ACC16)
      parameter(kreal=selected_real_kind(12,100))
!      parameter(kqreal=selected_real_kind(24,100))
      parameter(kqreal=selected_real_kind(12,100))
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      COMPLEX(kind=kreal) X,AA,BB,Z,HYPERGF,CDIGAM,CI
      COMPLEX(kind=kreal) DD,G,F,AI,BI,T
      LOGICAL ZLLIN
      REAL(kind=kqreal) AR,BR,GR,GI,DR,DI,TR,TI,UR,UI,FI,FI1,DEN,RZ
      DATA ZERO,ONE,TWO / 0D+0, 1D+0, 2D+0 /, CI / (0D+0, 1D+0) /
      ABSC(AA) = ABS(REAL(AA)) + ABS(AIMAG(AA))
      NINTC(AA) = NINT(REAL(AA))
   5  FORMAT(A, I12)
   7  FORMAT(A, 4ES12.2)
C      PRINT 7, "X=     ", X
C      PRINT 7, "EPS=   ", EPS
C      PRINT 5, "LIMIT= ", LIMIT
C      PRINT 5, "KIN=   ", KIN
C      PRINT 7, "ERR=   ", ERR
C      PRINT 5, "NITS=  ", NITS
C      PRINT 7, "FPMAX= ", FPMAX
C      PRINT 7, "ACC8=  ", ACC8
C      PRINT 7, "ACC16= ", ACC16
C      PRINT 7, "ONE=   ", ONE
C      PRINT 7, "TWO=   ", TWO
C      PRINT 7, "************"

C
C *** evaluate the HYPERGEOMETRIC FUNCTION 1F1
C                                        i
C            F (AA;BB; Z) = SUM  (AA)   Z / ( (BB)  i! )
C           1 1              i       i            i
C
C     to accuracy EPS with at most LIMIT terms.
C  If KIND = 0 : using extended precision but real arithmetic only,
C            1 : using normal precision in complex arithmetic,
C   or       2 : using normal complex arithmetic, but with CDIGAM factor
C
C  where
C  and
         Z  = X

         ACC16 = epsilon(AR)
	 KIND = KIN
	 if(KIND==0.and.(kreal==kqreal.or.ACC16>ACC8*.1)) KIND=1
C
         ZLLIN = REAL(BB).LE.ZERO .AND. ABS(BB-NINTC(BB)).LT.ACC8**0.75
             IF(.NOT.ZLLIN.OR.REAL(BB)+LIMIT.LT.1.5) GO TO 10
                NITS = -1
                RETURN
   10 IF(LIMIT.LE.0) THEN
         HYPERGF = ZERO
         ERR = ZERO
         NITS= 1
         RETURN
         ENDIF
      TA = ONE
      RK = ONE
      IF(KIND.LE.0.AND.ABSC(Z)*ABSC(AA).GT.ABSC(BB) * 1.0) THEN
         DR = ONE
         DI = ZERO
         GR = ONE
         GI = ZERO
         AR = AA   ! real part
         BR = BB   ! real part
         RZ = Z    ! real part
         FI = ZERO
      DO 20 I=2,LIMIT
         FI1 = FI + ONE
         TR = BR * FI1
         TI = AIMAG(BB) * FI1
         DEN= ONE / (TR*TR + TI*TI)
         UR = (AR*TR + AIMAG(AA)*TI) * DEN
         UI = (AIMAG(AA)*TR - AR*TI) * DEN
         TR = UR*GR - UI*GI
         TI = UR*GI + UI*GR
         GR = RZ * TR - AIMAG(Z)*TI
         GI = RZ * TI + AIMAG(Z)*TR
         DR = DR + GR
         DI = DI + GI
            ERR = ABS(GR) + ABS(GI)
               IF(ERR.GT.FPMAX) GO TO 60
            RK  = ABS(DR) + ABS(DI)
            TA = MAX(TA,RK)
         IF(ERR.LT.RK*EPS .OR. I.GE.4.AND.ERR.LT.ACC16) GO TO 30
         FI = FI1
         AR = AR + ONE
   20    BR = BR + ONE
C
   30    HYPERGF = DR + CI * DI
         ERR = ACC16 * TA / RK
C
      ELSE
C* ---------------------------------- alternative code
C*    If REAL*16 arithmetic is not available, (or already using it!),
C*    then use KIND > 0
         G = ONE
          F = ONE
          IF(KIND.GE.2) F = CDIGAM(AA) - CDIGAM(BB) - CDIGAM(G)
         DD = F
         DO 40 I=2,LIMIT
            AI = AA + (I-2)
            BI = BB + (I-2)
            R  = I-ONE
         G = G * Z * AI / (BI * R)
         IF(KIND.GE.2)
C                              multiply by (psi(a+r)-psi(b+r)-psi(1+r))
     X        F = F + ONE/AI - ONE/BI - ONE/R
         T  = G * F
         DD = DD + T
            ERR = ABSC(T)
               IF(ERR.GT.FPMAX) GO TO 60
            RK = ABSC(DD)
         TA = MAX(TA,RK)
         IF(ERR.LT.RK*EPS.OR.ERR.LT.ACC8.AND.I.GE.4) GO TO 50
   40    CONTINUE
 
   50    ERR = ACC8 * TA / RK
         HYPERGF = DD
C* ------------------------------------------- end of alternative code
      ENDIF
   60    NITS = I
      RETURN
      END
