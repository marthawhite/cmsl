* SUBROUTINE PLQDF1             ALL SYSTEMS                   91/12/01
* PORTABILITY : ALL SYSTEMS
* 91/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* DUAL RANGE SPACE QUADRATIC PROGRAMMING METHOD FOR MINIMAX
* APPROXIMATION WITH LINEAR CONSTRAINTS.
*
* PARAMETERS :
*  II  NF  NUMBER OF VARIABLES.
*  II  NA  NUMBER OF LINEAR APPROXIMATED FUNCTIONS.
*  II  NC  NUMBER OF LINEAR CONSTRAINTS.
*  RI  X(NF)  VECTOR OF VARIABLES.
*  II  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS.
*  RI  XL(NF)  VECTOR CONTAINING LOWER BOUNDS FOR VARIABLES.
*  RI  XU(NF)  VECTOR CONTAINING UPPER BOUNDS FOR VARIABLES.
*  RI  AF(NA)  VECTOR CONTAINING VALUES OF THE APPROXIMATED
*         FUNCTIONS.
*  RO  AFD(NA)  VECTOR CONTAINING INCREMENTS OF THE APPROXIMATED
*         FUNCTIONS.
*  II  IA(NA)  VECTOR CONTAINING TYPES OF DEVIATIONS.
*  IO  IAA(NF+1)  VECTOR CONTAINING INDICES OF ACTIVE FUNCTIONS.
*  RI  AG(NF*NA)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
*         APPROXIMATED FUNCTIONS.
*  RO  AR((NF+1)*(NF+2)/2)  TRIANGULAR DECOMPOSITION OF KERNEL OF THE
*         ORTHOGONAL PROJECTION.
*  RO  AZ(NF+1)  VECTOR OF LAGRANGE MULTIPLIERS.
*  RI  CF(NF)  VECTOR CONTAINING VALUES OF THE CONSTRAINT FUNCTIONS.
*  II  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
*  RI  CL(NC)  VECTOR CONTAINING LOWER BOUNDS FOR CONSTRAINT FUNCTIONS.
*  RI  CU(NC)  VECTOR CONTAINING UPPER BOUNDS FOR CONSTRAINT FUNCTIONS.
*  RI  CG(NF*NC)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
*         CONSTRAINTS.
*  RO  G(NF+1)  GRADIENT OF THE LAGRANGIAN FUNCTION.
*  RU  H(NF*(NF+1)/2)  TRIANGULAR DECOMPOSITION OR INVERSION OF THE
*         HESSIAN MATRIX APPROXIMATION.
*  RO  S(NF+1)  DIRECTION VECTOR.
*  II  MFP  TYPE OF FEASIBLE POINT. MFP=1-ARBITRARY FEASIBLE POINT.
*         MFP=2-OPTIMUM FEASIBLE POINT. MFP=3-REPEATED SOLUTION.
*  II  KBF  SPECIFICATION OF SIMPLE BOUNDS. KBF=0-NO SIMPLE BOUNDS.
*         KBF=1-ONE SIDED SIMPLE BOUNDS. KBF=2=TWO SIDED SIMPLE BOUNDS.
*  II  KBC  SPECIFICATION OF LINEAR CONSTRAINTS. KBC=0-NO LINEAR
*         CONSTRAINTS. KBC=1-ONE SIDED LINEAR CONSTRAINTS. KBC=2=TWO
*         SIDED LINEAR CONSTRAINTS.
*  IU  IDECF  DECOMPOSITION INDICATOR. IDECF=0-NO DECOMPOSITION.
*         IDECF=1-GILL-MURRAY DECOMPOSITION. IDECF=9-INVERSION.
*         IDECF=10-DIAGONAL MATRIX.
*  RI  ETA0  MACHINE PRECISION.
*  RI  ETA2  TOLERANCE FOR POSITIVE DEFINITENESS OF THE HESSIAN MATRIX.
*  RI  ETA9  MAXIMUM FOR REAL NUMBERS.
*  RI  EPS7  TOLERANCE FOR LINEAR INDEPENDENCE OF CONSTRAINTS.
*  RI  EPS9  TOLERANCE FOR ACTIVITY OF CONSTRAINTS.
*  RO  XNORM  VALUE OF LINEARIZED MINIMAX FUNCTION.
*  RO  UMAX  MAXIMUM ABSOLUTE VALUE OF A NEGATIVE LAGRANGE MULTIPLIER.
*  RO  GMAX  MAXIMUM ABSOLUTE VALUE OF A PARTIAL DERIVATIVE.
*  IO  N  DIMENSION OF MANIFOLD DEFINED BY ACTIVE CONSTRAINTS.
*  IO  ITERQ  TYPE OF FEASIBLE POINT. ITERQ=1-ARBITRARY FEASIBLE POINT.
*         ITERQ=2-OPTIMUM FEASIBLE POINT. ITERQ=-1 FEASIBLE POINT DOES
*         NOT EXISTS. ITERQ=-2 OPTIMUM FEASIBLE POINT DOES NOT EXISTS.
*
* SUBPROGRAMS USED :
*  S   PLMINA  DETERMINATION OF THE NEW ACTIVE FUNCTION.
*  S   PLMINL  DETERMINATION OF THE NEW ACTIVE LINEAR CONSTRAINT.
*  S   PLMINS  DETERMINATION OF THE NEW ACTIVE SIMPLE BOUND.
*  S   PLMINT  DETERMINATION OF THE NEW ACTIVE TRUST REGION BOUND.
*  S   PLADF1  CONSTRAINT ADDITION.
*  S   PLRMF0  CONSTRAINT DELETION.
*  S   MXPDGF  GILL-MURRAY DECOMPOSITION OF A DENSE SYMMETRIC MATRIX.
*  S   MXPDGB  BACK SUBSTITUTION AFTER GILL-MURRAY DECOMPOSITION.
*  S   MXDPRB  BACK SUBSTITUTION AFTER CHOLESKI DECOMPOSITION.
*  S   MXDSMM  MATRIX VECTOR PRODUCT.
*  S   MXVDIR  VECTOR AUGMENTED BY THE SCALED VECTOR.
*  RF  MXVDOT  DOT PRODUCT OF TWO VECTORS.
*  S   MXVCOP  COPYING OF A VECTOR.
*  S   MXVINA  ABSOLUTE VALUES OF ELEMENTS OF INTEGER VECTOR.
*  S   MXVINV  CHANGE OF INTEGER VECTOR AFTER CONSTRAINT ADDITION.
*  S   MXVNEG  COPYING OF A VECTOR WITH CHANGE OF THE SIGN.
*  S   MXVSET  INITIATION OF A VECTOR.
*
* L.LUKSAN: DUAL METHOD FOR SOLVING A SPECIAL PROBLEM OF QUADRATIC
* PROGRAMMING AS A SUBPROBLEM AT LINEARLY CONSTRAINED NONLINEAR MINIMAX
* APPROXIMATION. KYBERNETIKA 20 (1984) 445-457.
*
      SUBROUTINE PLQDF1(NF,NA,NC,X,IX,XL,XU,AF,AFD,IA,IAA,AG,AR,AZ,CF,
     +                  IC,CL,CU,CG,G,H,S,MFP,KBF,KBC,IDECF,ETA0,ETA2,
     +                  ETA9,EPS7,EPS9,XNORM,UMAX,GMAX,N,ITERQ)
C     .. Scalar Arguments ..
      DOUBLE PRECISION EPS7,EPS9,ETA0,ETA2,ETA9,GMAX,UMAX,XNORM
      INTEGER IDECF,ITERQ,KBC,KBF,MFP,N,NA,NC,NF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION AF(*),AFD(*),AG(*),AR(*),AZ(*),CF(*),CG(*),CL(*),
     +                 CU(*),G(*),H(*),S(*),X(*),XL(*),XU(*)
      INTEGER IA(*),IAA(*),IC(*),IX(*)
C     ..
C     .. Scalars in Common ..
      INTEGER NADD,NDECF,NFG,NFH,NFV,NIT,NRED,NREM,NRES
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION BET,CON,E,GAM,PAR,Q,QO,SNORM,STEP,STEP1,STEP2,T,
     +                 TEMP
      INTEGER I,IER,INEW,INF,IOLD,J,K,KA,KC,KNEW,KOLD,KREM,L,NAA,NAR
      INTEGER NCYC,MCYC
C     ..
C     .. External Functions ..
      DOUBLE PRECISION MXVDOT,MXVMAX
      EXTERNAL MXVDOT,MXVMAX
C     ..
C     .. External Subroutines ..
      EXTERNAL MXDPGB,MXDPGF,MXDPRB,MXDSMM,MXVCOP,MXVDIR,MXVINA,MXVINV,
     +         MXVMUL,MXVNEG,MXVSET,PLADF1,PLDLAG,PLMINA,PLMINL,PLMINS,
     +         PLRMF0
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MIN,SIGN
C     ..
C     .. Common blocks ..
      COMMON /STAT/NDECF,NRES,NRED,NREM,NADD,NIT,NFV,NFG,NFH
C     ..
      NADD = 0
      NREM = 0
      T = 1.0D0
      CON = ETA9
      IF (IDECF.LT.0) IDECF = 1
      IF (IDECF.EQ.0) THEN
*
*     GILL-MURRAY DECOMPOSITION
*
          TEMP = ETA2
          CALL MXDPGF(NF,H,INF,TEMP,STEP)
          NDECF = NDECF + 1
          IDECF = 1
      END IF

      IF (IDECF.GE.2 .AND. IDECF.LE.8) THEN
          ITERQ = -10
          RETURN

      END IF
*
*     INITIATION
*
      NRED = 0
      NCYC = 0
      MCYC = MAX(1000,100*(NA+NC))
      IF (MFP.EQ.3) GO TO 10
      N = NF
      NAA = 0
      NAR = 0
C      XNORM = -ETA9
      XNORM = -1.0D20
      Q = 2.0D0 * XNORM
      CALL MXVINA(NA,IA)
      IF (KBF.GT.0) CALL MXVINA(NF,IX)
      IF (KBC.GT.0) CALL MXVINA(NC,IC)
*
*     DIRECTION DETERMINATION
*
   10 CALL MXVSET(NF,0.0D0,S)
      NCYC = NCYC + 1
      DO 20 J = 1,NAA
          L = IAA(J)
          IF (L.GT.NC) THEN
              L = L - NC
              CALL MXVDIR(NF,AZ(J),AG((L-1)*NF+1),S,S)

          ELSE IF (L.GT.0) THEN
              CALL MXVDIR(NF,AZ(J),CG((L-1)*NF+1),S,S)

          ELSE
              L = -L
              S(L) = S(L) + AZ(J)
          END IF

   20 CONTINUE
      CALL MXVCOP(NF,S,G)
      IF (NAA.GT.0) THEN
          IF (IDECF.EQ.1) THEN
              CALL MXDPGB(NF,H,S,0)

          ELSE IF (IDECF.EQ.9) THEN
              CALL MXDSMM(NF,H,G,S)

          ELSE
              CALL MXVMUL(NF,H,S,S,-1)
          END IF

      END IF
*
*     INITIAL MINIMAX VARIABLE
*
      IF (NAA.EQ.1) THEN
          IF (INEW.LE.NC) THEN
          ELSE
              TEMP = AF(INEW-NC) + MXVDOT(NF,AG((INEW-NC-1)*NF+1),S)
              XNORM = -SIGN(1,KNEW)*TEMP
          ENDIF
      END IF
*
*     CHECK OF FEASIBILITY
*
      INEW = 0
      PAR = 0.0D0
      CALL PLMINA(NF,NA,NC,AF,AFD,IA,AG,S,INEW,KNEW,EPS9,XNORM,PAR)
      IF (NAA.GT.0) THEN
          CALL PLMINL(NF,NC,CF,IC,CL,CU,CG,S,KBC,INEW,KNEW,EPS9,PAR)
          CALL PLMINS(NF,IX,X,XL,XU,S,KBF,INEW,KNEW,EPS9,PAR)
      END IF

      IF (INEW.EQ.0) THEN
*
*     SOLUTION ACHIEVED
*
          CALL MXVNEG(NF,G,G)
          ITERQ = 2
          RETURN

      ELSE IF (NCYC.GE.MCYC) THEN
*
*     CYCLE OCCURRED
*
          CALL MXVNEG(NF,G,G)
          ITERQ = 4
          RETURN

      ELSE
*
*     IRREGULAR CASE
*
          QO = Q
          Q=0.5D0 * MXVDOT(NF,G,S) + XNORM
          IF (Q.LE.QO) THEN
              CALL MXVNEG(NF,G,G)
              ITERQ = 3
              RETURN
          END IF

          SNORM = 0.0D0
      END IF

   30 IER = 0
*
*     STEPSIZE DETERMINATION
*
      CALL PLADF1(NF,NC,IA,IAA,AG,AR,CG,H,S,G,IDECF,N,INEW,KNEW,IER,
     +            EPS7,GMAX,UMAX,E,T)
      CALL PLDLAG(NF,NC,IA,IAA,S,N,KOLD)
      IF (KOLD.EQ.0) THEN
*
*     ZERO STEPSIZE
*
          STEP1 = 0.0D0
          STEP = STEP1
          SNORM = SIGN(1,KNEW)
          XNORM = XNORM - PAR

      ELSE
*
*     PRIMAL STEPSIZE
*
          CALL MXDPRB(NAA,AR,S,1)
          BET = E - MXVDOT(NAA,S,G)
          GAM = BET/MXVDOT(NAA,S,S)
          UMAX = BET*GAM + UMAX
          IF (UMAX.LE.EPS7*GMAX) THEN
              STEP1 = CON

          ELSE
              STEP1 = -PAR/UMAX
          END IF
*
*     DUAL STEPSIZE
*
          CALL MXDPRB(NAA,AR,S,-1)
          CALL MXDPRB(NAA,AR,G,-1)
          CALL MXVDIR(NAA,GAM,S,G,G)
          IF (KNEW.LT.0) CALL MXVNEG(NAA,G,G)
          STEP = MXVMAX(NAA,G)
          IOLD = 0
          STEP2 = CON
          DO 40 J = 1,NAA
              L = IAA(J)
              IF (L.GT.NC) THEN
                  L = L - NC
                  K = IA(L)

              ELSE IF (L.GT.0) THEN
                  K = IC(L)

              ELSE
                  L = -L
                  K = IX(L)
              END IF

              IF (K.LE.-5) THEN

              ELSE IF ((K.EQ.-1.OR.K.EQ.-3.) .AND. G(J).LE.0.0D0) THEN

              ELSE IF ((K.EQ.-2.OR.K.EQ.-4.) .AND. G(J).GE.0.0D0) THEN

              ELSE IF (ABS(G(J)).LE.ETA0*STEP) THEN

              ELSE
                  TEMP = AZ(J)/G(J)
                  IF (STEP2.GT.TEMP) THEN
                      IOLD = J
                      STEP2 = TEMP
                  END IF

              END IF

   40     CONTINUE
*
*     FINAL STEPSIZE
*
          STEP = MIN(STEP1,STEP2)
          IF (STEP.GE.CON) THEN
*
*     FEASIBLE SOLUTION DOES NOT EXIST
*
              ITERQ = -1
              RETURN

          END IF
*
*     NEW LAGRANGE MULTIPLIERS
*
          CALL MXVDIR(NAA,-STEP,G,AZ,AZ)
          SNORM = SNORM + SIGN(1,KNEW)*STEP
          XNORM = XNORM + SIGN(1,KNEW)*STEP*GAM
          PAR = PAR - (STEP/STEP1)*PAR
      END IF

      IF (STEP.EQ.STEP1) THEN
          IF (N.LT.0) THEN
*
*     IMPOSSIBLE SITUATION
*
              ITERQ = -5
              RETURN

          END IF
*
*     CONSTRAINT ADDITION
*
          IF (IER.EQ.0) THEN
              N = N - 1
              NAA = NAA + 1
              NAR = NAR + NAA
              AZ(NAA) = SNORM
          END IF

          IF (INEW.GT.NC) THEN
              KA = INEW - NC
              CALL MXVINV(IA,KA,KNEW)

          ELSE IF (INEW.GT.0) THEN
              KC = INEW
              CALL MXVINV(IC,KC,KNEW)

          ELSE IF (ABS(KNEW).EQ.1) THEN
              I = -INEW
              CALL MXVINV(IX,I,KNEW)

          ELSE
              I = -INEW
              IF (KNEW.GT.0) IX(I) = -3
              IF (KNEW.LT.0) IX(I) = -4
          END IF

          NRED = NRED + 1
          NADD = NADD + 1
          GO TO 10

      ELSE
*
*     CONSTRAINT DELETION
*
          DO 50 J = IOLD,NAA - 1
              AZ(J) = AZ(J+1)
   50     CONTINUE
          CALL PLRMF0(NF,NC,IX,IA,IAA,AR,IC,G,N,IOLD,KREM,IER)
          NAR = NAR - NAA
          NAA = NAA - 1
          CALL MXVINA(NA,IA)
          IF (KBC.GT.0) CALL MXVINA(NC,IC)
          IF (KBF.GT.0) CALL MXVINA(NF,IX)
          DO 60 J = 1,NAA
              L = IAA(J)
              IF (L.GT.NC) THEN
                  L = L - NC
                  IA(L) = -IA(L)

              ELSE IF (L.GT.0) THEN
                  IC(L) = -IC(L)

              ELSE
                  L = -L
                  IX(L) = -IX(L)
              END IF

   60     CONTINUE
          GO TO 30

      END IF

      END
* SUBROUTINE PLMINA             ALL SYSTEMS                   90/12/01
* PORTABILITY : ALL SYSTEMS
* 90/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* DETERMINATION OF THE NEW ACTIVE FUNCTION.
*
* PARAMETERS :
*  II  NF  DECLARED NUMBER OF VARIABLES.
*  II  NA  NUMBER OF CURRENT LINEAR APPROXIMATED FUNCTIONS.
*  II  NC  NUMBER OF CURRENT LINEAR CONSTRAINTS.
*  RI  AF(NA)  VECTOR CONTAINING VALUES OF THE APPROXIMATED
*         FUNCTIONS.
*  RO  AFD(NA)  VECTOR CONTAINING INCREMENTS OF THE APPROXIMATED
*         FUNCTIONS.
*  II  IA(NA)  VECTOR CONTAINING TYPES OF DEVIATIONS.
*  RI  AG(NF*NA)  VECTOR CONTAINING SCALING PARAMETERS.
*  RI  S(NF)  DIRECTION VECTOR.
*  IO  INEW  INDEX OF THE NEW ACTIVE FUNCTION.
*  IO  KNEW  SIGNUM OF THE NEW ACTIVE GRADIENT.
*  RI  EPS9  TOLERANCE FOR ACTIVE FUNCTIONS.
*  RO  XNORM  VALUE OF LINEARIZED MINIMAX FUNCTION.
*  RA  PAR  AUXILIARY VARIABLE.
*
* SUBPROGRAMS USED :
*  RF  MXVDOT  DOT PRODUCT OF TWO VECTORS.
*
      SUBROUTINE PLMINA(NF,NA,NC,AF,AFD,IA,AG,S,INEW,KNEW,EPS9,XNORM,
     +                  PAR)
C     .. Scalar Arguments ..
      DOUBLE PRECISION EPS9,PAR,XNORM
      INTEGER INEW,KNEW,NA,NC,NF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION AF(*),AFD(*),AG(*),S(*)
      INTEGER IA(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION POM,TEMP
      INTEGER JCG,KA
C     ..
C     .. External Functions ..
      DOUBLE PRECISION MXVDOT
      EXTERNAL MXVDOT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN
C     ..
      JCG = 1
      DO 10 KA = 1,NA
          IF (IA(KA).GT.0) THEN
              TEMP = MXVDOT(NF,AG(JCG),S)
              AFD(KA) = TEMP
              TEMP = AF(KA) + TEMP
              IF (IA(KA).EQ.1 .OR. IA(KA).GE.3) THEN
                  POM = XNORM + TEMP
                  IF (POM.LT.MIN(PAR,-EPS9*MAX(ABS(XNORM),1.0D0))) THEN
                      INEW = KA + NC
                      KNEW = 1
                      PAR = POM
                  END IF

              END IF

              IF (IA(KA).EQ.2 .OR. IA(KA).GE.3) THEN
                  POM = XNORM - TEMP
                  IF (POM.LT.MIN(PAR,-EPS9*MAX(ABS(XNORM),1.0D0))) THEN
                      INEW = KA + NC
                      KNEW = -1
                      PAR = POM
                  END IF

              END IF

          END IF

          JCG = JCG + NF
   10 CONTINUE
      RETURN

      END
* SUBROUTINE PLMINL             ALL SYSTEMS                   90/12/01
* PORTABILITY : ALL SYSTEMS
* 90/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* DETERMINATION OF THE NEW ACTIVE LINEAR CONSTRAINT.
*
* PARAMETERS :
*  II  NF  DECLARED NUMBER OF VARIABLES.
*  II  NC  NUMBER OF CONSTRAINTS.
*  RI  CF(NC)  VECTOR CONTAINING VALUES OF THE CONSTRAINT FUNCTIONS.
*  II  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
*  RI  CL(NC)  VECTOR CONTAINING LOWER BOUNDS FOR CONSTRAINT FUNCTIONS.
*  RI  CU(NC)  VECTOR CONTAINING UPPER BOUNDS FOR CONSTRAINT FUNCTIONS.
*  RI  CG(NF*NC)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
*         CONSTRAINTS.
*  RI  S(NF)  DIRECTION VECTOR.
*  II  KBC  SPECIFICATION OF LINEAR CONSTRAINTS. KBC=0-NO LINEAR
*         CONSTRAINTS. KBC=1-ONE SIDED LINEAR CONSTRAINTS. KBC=2=TWO
*         SIDED LINEAR CONSTRAINTS.
*  IO  INEW  INDEX OF THE NEW ACTIVE CONSTRAINT.
*  IO  KNEW  SIGNUM OF THE NEW ACTIVE NORMAL.
*  RI  EPS9  TOLERANCE FOR ACTIVE CONSTRAINTS.
*  RA  PAR  AUXILIARY VARIABLE.
*
* SUBPROGRAMS USED :
*  RF  MXVDOT  DOT PRODUCT OF TWO VECTORS.
*
      SUBROUTINE PLMINL(NF,NC,CF,IC,CL,CU,CG,S,KBC,INEW,KNEW,EPS9,PAR)
C     .. Scalar Arguments ..
      DOUBLE PRECISION EPS9,PAR
      INTEGER INEW,KBC,KNEW,NC,NF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CF(*),CG(*),CL(*),CU(*),S(*)
      INTEGER IC(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION POM,TEMP
      INTEGER JCG,KC
C     ..
C     .. External Functions ..
      DOUBLE PRECISION MXVDOT
      EXTERNAL MXVDOT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN
C     ..
      IF (KBC.GT.0) THEN
          JCG = 1
          DO 10 KC = 1,NC
              IF (IC(KC).GT.0) THEN
                  TEMP = CF(KC) + MXVDOT(NF,CG(JCG),S)
                  IF (IC(KC).EQ.1 .OR. IC(KC).GE.3) THEN
                      POM = TEMP - CL(KC)
                      IF (POM.LT.MIN(PAR,-EPS9*MAX(ABS(CL(KC)),
     +                    1.0D0))) THEN
                          INEW = KC
                          KNEW = 1
                          PAR = POM
                      END IF

                  END IF

                  IF (IC(KC).EQ.2 .OR. IC(KC).GE.3) THEN
                      POM = CU(KC) - TEMP
                      IF (POM.LT.MIN(PAR,-EPS9*MAX(ABS(CU(KC)),
     +                    1.0D0))) THEN
                          INEW = KC
                          KNEW = -1
                          PAR = POM
                      END IF

                  END IF

              END IF

              JCG = JCG + NF
   10     CONTINUE
      END IF

      RETURN

      END
* SUBROUTINE PLMINS             ALL SYSTEMS                   91/12/01
* PORTABILITY : ALL SYSTEMS
* 91/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* DETERMINATION OF THE NEW ACTIVE SIMPLE BOUND.
*
* PARAMETERS :
*  II  NF DECLARED NUMBER OF VARIABLES.
*  II  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS.
*  RI  XO(NF)  SAVED VECTOR OF VARIABLES.
*  RI  XL(NF)  VECTOR CONTAINING LOWER BOUNDS FOR VARIABLES.
*  RI  XU(NF)  VECTOR CONTAINING UPPER BOUNDS FOR VARIABLES.
*  RI  S(NF)  DIRECTION VECTOR.
*  II  KBF  SPECIFICATION OF SIMPLE BOUNDS. KBF=0-NO SIMPLE BOUNDS.
*         KBF=1-ONE SIDED SIMPLE BOUNDS. KBF=2=TWO SIDED SIMPLE BOUNDS.
*  IO  INEW  INDEX OF THE NEW ACTIVE CONSTRAINT.
*  IO  KNEW  SIGNUM OF THE NEW NORMAL.
*  RI  EPS9  TOLERANCE FOR ACTIVE CONSTRAINTS.
*  RA  PAR  AUXILIARY VARIABLE.
*
      SUBROUTINE PLMINS(NF,IX,XO,XL,XU,S,KBF,INEW,KNEW,EPS9,PAR)
C     .. Scalar Arguments ..
      DOUBLE PRECISION EPS9,PAR
      INTEGER INEW,KBF,KNEW,NF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION S(*),XL(*),XO(*),XU(*)
      INTEGER IX(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION POM,TEMP
      INTEGER I
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN
C     ..
      IF (KBF.GT.0) THEN
          DO 10 I = 1,NF
              IF (IX(I).GT.0) THEN
                  TEMP = 1.0D0
                  IF (IX(I).EQ.1 .OR. IX(I).GE.3) THEN
                      POM = XO(I) + S(I)*TEMP - XL(I)
                      IF (POM.LT.MIN(PAR,-EPS9*MAX(ABS(XL(I)),
     +                    TEMP))) THEN
                          INEW = -I
                          KNEW = 1
                          PAR = POM
                      END IF

                  END IF

                  IF (IX(I).EQ.2 .OR. IX(I).GE.3) THEN
                      POM = XU(I) - S(I)*TEMP - XO(I)
                      IF (POM.LT.MIN(PAR,-EPS9*MAX(ABS(XU(I)),
     +                    TEMP))) THEN
                          INEW = -I
                          KNEW = -1
                          PAR = POM
                      END IF

                  END IF

              END IF

   10     CONTINUE
      END IF

      RETURN

      END
* SUBROUTINE PLDLAG               ALL SYSTEMS                91/12/01
* PORTABILITY : ALL SYSTEMS
* 91/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* VECTOR OF LAGRANGE MULTIPLIERS FOR DUAL QP METHOD IS DETERMINED.
*
* PARAMETERS :
*  II  NF  DECLARED NUMBER OF VARIABLES.
*  II  NC  NUMBER OF LINEARIZED CONSTRAINTS.
*  II  IA(NA)  VECTOR CONTAINING TYPES OF DEVIATIONS.
*  II  IAA(NF+1)  VECTOR CONTAINING INDICES OF ACTIVE FUNCTIONS.
*  RO  AZ(NF+1)  OUTPUT VECTOR.
*  II  N  ACTUAL NUMBER OF VARIABLES.
*  IA  KOLD  AUXILIARY VARIABLE.
*
      SUBROUTINE PLDLAG(NF,NC,IA,IAA,AZ,N,KOLD)
C     .. Scalar Arguments ..
      INTEGER KOLD,N,NC,NF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION AZ(*)
      INTEGER IA(*),IAA(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER J,L,NAA
C     ..
      NAA = NF - N
      KOLD = 0
      DO 10 J = 1,NAA
          L = IAA(J)
          IF (L.GT.NC) THEN
              L = L - NC
              TEMP = 1.0D0
              IF (IA(L).EQ.-2 .OR. IA(L).EQ.-4) TEMP = -TEMP
              AZ(J) = TEMP
              KOLD = 1

          ELSE
              AZ(J) = 0.0D0
          END IF

   10 CONTINUE
      RETURN

      END
* SUBROUTINE PLADF1               ALL SYSTEMS                91/12/01
* PORTABILITY : ALL SYSTEMS
* 91/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* TRIANGULAR DECOMPOSITION OF KERNEL OF THE GENERAL PROJECTION
* IS UPDATED AFTER FUNCTION OR CONSTRAINT ADDITION.
*
* PARAMETERS :
*  II  NF  DECLARED NUMBER OF VARIABLES.
*  II  NC  NUMBER OF LINEARIZED CONSTRAINTS.
*  II  IA(NA)  VECTOR CONTAINING TYPES OF DEVIATIONS.
*  IU  IAA(NF+1)  VECTOR CONTAINING INDICES OF ACTIVE FUNCTIONS.
*  RI  AG(NF*NA)  MATRIX WHOSE COLUMNS ARE GRADIENTS OF THE LINEAR
*          APPROXIMATED FUNCTIONS.
*  RU  AR((NF+1)*(NF+2)/2)  TRIANGULAR DECOMPOSITION OF KERNEL OF
*        THE ORTHOGONAL PROJECTION.
*  RI  CG(NF*NC)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
*         CONSTRAINTS.
*  RI  H(NF*(NF+1)/2)  TRIANGULAR DECOMPOSITION OR INVERSION OF THE
*         HESSIAN MATRIX APPROXIMATION.
*  RA  S(NF+1)  AUXILIARY VECTOR.
*  RO  G(NF+1)  VECTOR USED IN THE DUAL RANGE SPACE QUADRATIC
*        PROGRAMMING METHOD.
*  IU  IDECF  DECOMPOSITION INDICATOR. IDECF=0-NO DECOMPOSITION.
*         IDECF=1-GILL-MURRAY DECOMPOSITION. IDECF=9-INVERSION.
*         IDECF=10-DIAGONAL MATRIX.
*  IU  N  ACTUAL NUMBER OF VARIABLES.
*  II  INEW  INDEX OF THE NEW ACTIVE CONSTRAINT.
*  II  KNEW  SIGNUM OF THE NEW ACTIVE GRADIENT.
*  IO  IER  ERROR INDICATOR.
*  RI  EPS7  TOLERANCE FOR LINEAR INDEPENDENCE OF CONSTRAINTS.
*  RO  GMAX  MAXIMUM ABSOLUTE VALUE OF A PARTIAL DERIVATIVE.
*  RO  UMAX  MAXIMUM ABSOLUTE VALUE OF A NEGATIVE LAGRANGE MULTIPLIER.
*  RO  E  AUXILIARY VARIABLE.
*  RI  T  AUXILIARY VARIABLE.
*
* SUBPROGRAMS USED :
*  S   MXPDGB  BACK SUBSTITUTION AFTER GILL-MURRAY DECOMPOSITION.
*  S   MXDPRB  BACK SUBSTITUTION AFTER CHOLESKI DECOMPOSITION.
*  S   MXDSMM  MATRIX-VECTOR PRODUCT.
*  S   MXDSMV  COPYING OF A ROW OF DENSE SYMMETRIC MATRIX.
*  S   MXVCOP  COPYING OF A VECTOR.
*  RF  MXVDOT  DOT PRODUCT OF TWO VECTORS.
*
      SUBROUTINE PLADF1(NF,NC,IA,IAA,AG,AR,CG,H,S,G,IDECF,N,INEW,KNEW,
     +                  IER,EPS7,GMAX,UMAX,E,T)
C     .. Scalar Arguments ..
      DOUBLE PRECISION E,EPS7,GMAX,T,UMAX
      INTEGER IDECF,IER,INEW,KNEW,N,NC,NF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION AG(*),AR(*),CG(*),G(*),H(*),S(*)
      INTEGER IA(*),IAA(*)
C     ..
C     .. Scalars in Common ..
      INTEGER NADD,NDECF,NFG,NFH,NFV,NIT,NRED,NREM,NRES
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION POM,TEMP
      INTEGER J,JAG,JOB,K,L,NAA,NAR
C     ..
C     .. External Functions ..
      DOUBLE PRECISION MXVDOT
      EXTERNAL MXVDOT
C     ..
C     .. External Subroutines ..
      EXTERNAL MXDPGB,MXDPRB,MXDSMM,MXDSMV,MXVCOP,MXVMUL,MXVSET
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SIGN,SQRT
C     ..
C     .. Common blocks ..
      COMMON /STAT/NDECF,NRES,NRED,NREM,NADD,NIT,NFV,NFG,NFH
C     ..
      JOB = 1
      E = 0.0D0
      IF (INEW.GT.NC) E = SIGN(1,KNEW)
      IER = 0
      IF (JOB.EQ.0 .AND. N.LT.0) IER = 2
      IF (INEW.EQ.0) IER = 3
      IF (IDECF.GE.2 .AND. IDECF.LE.8) IER = -2
      IF (IER.NE.0) RETURN
      NAA = NF - N
      NAR = NAA* (NAA+1)/2
      IF (INEW.GT.NC) THEN
          JAG = (INEW-NC-1)*NF + 1
          IF (IDECF.EQ.1) THEN
              CALL MXVCOP(NF,AG(JAG),S)
              CALL MXDPGB(NF,H,S,0)

          ELSE IF (IDECF.EQ.9) THEN
              CALL MXDSMM(NF,H,AG(JAG),S)

          ELSE
              CALL MXVCOP(NF,AG(JAG),S)
              CALL MXVMUL(NF,H,S,S,-1)
          END IF

          GMAX = MXVDOT(NF,AG(JAG),S) + T

      ELSE IF (INEW.GT.0) THEN
          JAG = (INEW-1)*NF + 1
          IF (IDECF.EQ.1) THEN
              CALL MXVCOP(NF,CG(JAG),S)
              CALL MXDPGB(NF,H,S,0)

          ELSE IF (IDECF.EQ.9) THEN
              CALL MXDSMM(NF,H,CG(JAG),S)

          ELSE
              CALL MXVCOP(NF,CG(JAG),S)
              CALL MXVMUL(NF,H,S,S,-1)
          END IF

          GMAX = MXVDOT(NF,CG(JAG),S)

      ELSE
          K = -INEW
          IF (IDECF.EQ.1) THEN
              CALL MXVSET(NF,0.0D0,S)
              S(K) = 1.0D0
              CALL MXDPGB(NF,H,S,0)

          ELSE IF (IDECF.EQ.9) THEN
              CALL MXDSMV(NF,H,S,K)

          ELSE
              CALL MXVSET(NF,0.0D0,S)
              S(K) = 1.0D0/H(K)
          END IF

          GMAX = S(K)
      END IF

      IF (NAA.GT.0) THEN
          POM = T*E
          DO 10 J = 1,NAA
              L = IAA(J)
              IF (L.GT.NC) THEN
                  L = L - NC
                  G(J) = MXVDOT(NF,AG((L-1)*NF+1),S)
                  IF (INEW.GT.NC) THEN
                      TEMP = POM
                      IF (IA(L).EQ.-2 .OR. IA(L).EQ.-4) TEMP = -TEMP
                      G(J) = G(J) + TEMP
                  END IF

              ELSE IF (L.GT.0) THEN
                  G(J) = MXVDOT(NF,CG((L-1)*NF+1),S)

              ELSE
                  L = -L
                  G(J) = S(L)
              END IF

   10     CONTINUE
      END IF

      IF (N.LT.0) THEN
          CALL MXDPRB(NAA,AR,G,1)
          UMAX = 0.0D0
          IER = 2
          RETURN

      ELSE IF (NAA.EQ.0) THEN
          UMAX = GMAX

      ELSE
          CALL MXDPRB(NAA,AR,G,1)
          UMAX = GMAX - MXVDOT(NAA,G,G)
          CALL MXVCOP(NAA,G,AR(NAR+1))
      END IF

      IF (UMAX.LE.EPS7*GMAX) THEN
          IER = 1
          RETURN

      ELSE
          NAA = NAA + 1
          NAR = NAR + NAA
          IAA(NAA) = INEW
          AR(NAR) = SQRT(UMAX)
          IF (JOB.EQ.0) THEN
              N = N - 1
              NADD = NADD + 1
          END IF

      END IF

      RETURN

      END
* SUBROUTINE PLRMF0             ALL SYSTEMS                   91/12/01
* PORTABILITY : ALL SYSTEMS
* 91/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* OPERATIONS AFTER CONSTRAINT DELETION.
*
* PARAMETERS :
*  II  NF  DECLARED NUMBER OF VARIABLES.
*  II  NC  NUMBER OF CONSTRAINTS.
*  II  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS.
*  II  IA(NA)  VECTOR CONTAINING TYPES OF DEVIATIONS.
*  IU  IAA(NF+1)  VECTOR CONTAINING INDICES OF ACTIVE FUNCTIONS.
*  RU  AR((NF+1)*(NF+2)/2)  TRIANGULAR DECOMPOSITION OF KERNEL OF THE
*         ORTHOGONAL PROJECTION.
*  II  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
*  RA  S(NF+1)  AUXILIARY VECTOR.
*  II  N  ACTUAL NUMBER OF VARIABLES.
*  II  IOLD  INDEX OF THE OLD ACTIVE CONSTRAINT.
*  IO  KREM  AUXILIARY VARIABLE.
*  IO  IER  ERROR INDICATOR.
*
* SUBPROGRAMS USED :
*  S   PLRMR0  CORRECTION OF KERNEL OF THE ORTHOGONAL PROJECTION
*         AFTER CONSTRAINT DELETION.
*
      SUBROUTINE PLRMF0(NF,NC,IX,IA,IAA,AR,IC,S,N,IOLD,KREM,IER)
C     .. Scalar Arguments ..
      INTEGER IER,IOLD,KREM,N,NC,NF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION AR(*),S(*)
      INTEGER IA(*),IAA(*),IC(*),IX(*)
C     ..
C     .. Scalars in Common ..
      INTEGER NADD,NDECF,NFG,NFH,NFV,NIT,NRED,NREM,NRES
C     ..
C     .. Local Scalars ..
      INTEGER L
C     ..
C     .. External Subroutines ..
      EXTERNAL PLRMR0
C     ..
C     .. Common blocks ..
      COMMON /STAT/NDECF,NRES,NRED,NREM,NADD,NIT,NFV,NFG,NFH
C     ..
      CALL PLRMR0(NF,IAA,AR,S,N,IOLD,KREM,IER)
      N = N + 1
      NREM = NREM + 1
      L = IAA(NF-N+1)
      IF (L.GT.NC) THEN
          L = L - NC
          IA(L) = -IA(L)

      ELSE IF (L.GT.0) THEN
          IC(L) = -IC(L)

      ELSE
          L = -L
          IX(L) = -IX(L)
      END IF

      RETURN

      END
* SUBROUTINE PLRMR0               ALL SYSTEMS                91/12/01
* PORTABILITY : ALL SYSTEMS
* 91/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* TRIANGULAR DECOMPOSITION OF KERNEL OF THE ORTHOGONAL PROJECTION IS
* UPDATED AFTER CONSTRAINT DELETION.
*
* PARAMETERS :
*  II  NF  DECLARED NUMBER OF VARIABLES.
*  IU  ICA(NF)  VECTOR CONTAINING INDICES OF ACTIVE CONSTRAINTS.
*  RU  CR(NF*(NF+1)/2)  TRIANGULAR DECOMPOSITION OF KERNEL OF THE
*         ORTHOGONAL PROJECTION.
*  RA  G(NF)  AUXILIARY VECTOR.
*  II  N  ACTUAL NUMBER OF VARIABLES.
*  II  IOLD  INDEX OF THE OLD ACTIVE CONSTRAINT.
*  IO  KREM  AUXILIARY VARIABLE.
*  IO  IER  ERROR INDICATOR.
*
* SUBPROGRAMS USED :
*  S   MXVCOP  COPYING OF A VECTOR.
*  S   MXVORT  DETERMINATION OF AN ELEMENTARY ORTHOGONAL MATRIX FOR
*         PLANE ROTATION.
*  S   MXVROT  PLANE ROTATION OF A VECTOR.
*  S   MXVSET  INITIATION OF A VECTOR.
*
      SUBROUTINE PLRMR0(NF,ICA,CR,G,N,IOLD,KREM,IER)
C     .. Scalar Arguments ..
      INTEGER IER,IOLD,KREM,N,NF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CR(*),G(*)
      INTEGER ICA(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CK,CL
      INTEGER I,J,K,KC,L,NCA
C     ..
C     .. External Subroutines ..
      EXTERNAL MXVCOP,MXVORT,MXVROT,MXVSET
C     ..
      NCA = NF - N
      IF (IOLD.LT.NCA) THEN
          K = IOLD* (IOLD-1)/2
          KC = ICA(IOLD)
          CALL MXVCOP(IOLD,CR(K+1),G)
          CALL MXVSET(NCA-IOLD,0.0D0,G(IOLD+1))
          K = K + IOLD
          DO 20 I = IOLD + 1,NCA
              K = K + I
              CALL MXVORT(CR(K-1),CR(K),CK,CL,IER)
              CALL MXVROT(G(I-1),G(I),CK,CL,IER)
              L = K
              DO 10 J = I,NCA - 1
                  L = L + J
                  CALL MXVROT(CR(L-1),CR(L),CK,CL,IER)
   10         CONTINUE
   20     CONTINUE
          K = IOLD* (IOLD-1)/2
          DO 30 I = IOLD,NCA - 1
              L = K + I
              ICA(I) = ICA(I+1)
              CALL MXVCOP(I,CR(L+1),CR(K+1))
              K = L
   30     CONTINUE
          ICA(NCA) = KC
          CALL MXVCOP(NCA,G,CR(K+1))
          KREM = 1
      END IF

      RETURN

      END
* SUBROUTINE MXDPGB                ALL SYSTEMS                91/12/01
* PORTABILITY : ALL SYSTEMS
* 91/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* SOLUTION OF A SYSTEM OF LINEAR EQUATIONS WITH A DENSE SYMMETRIC
* POSITIVE DEFINITE MATRIX A+E USING THE FACTORIZATION A+E=L*D*TRANS(L)
* OBTAINED BY THE SUBROUTINE MXDPGF.
*
* PARAMETERS :
*  II  N ORDER OF THE MATRIX A.
*  RI  A(N*(N+1)/2) FACTORIZATION A+E=L*D*TRANS(L) OBTAINED BY THE
*         SUBROUTINE MXDPGF.
*  RU  X(N)  ON INPUT THE RIGHT HAND SIDE OF A SYSTEM OF LINEAR
*         EQUATIONS. ON OUTPUT THE SOLUTION OF A SYSTEM OF LINEAR
*         EQUATIONS.
*  II  JOB  OPTION. IF JOB=0 THEN X:=(A+E)**(-1)*X. IF JOB>0 THEN
*         X:=L**(-1)*X. IF JOB<0 THEN X:=TRANS(L)**(-1)*X.
*
* METHOD :
* BACK SUBSTITUTION
*
      SUBROUTINE MXDPGB(N,A,X,JOB)
*     .. Scalar Arguments ..
      INTEGER          JOB,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(*),X(*)
*     ..
*     .. Local Scalars ..
      INTEGER          I,II,IJ,J
*     ..
      IF (JOB.GE.0) THEN
*
*     PHASE 1 : X:=L**(-1)*X
*
          IJ = 0
          DO 20 I = 1,N
              DO 10 J = 1,I - 1
                  IJ = IJ + 1
                  X(I) = X(I) - A(IJ)*X(J)
   10         CONTINUE
              IJ = IJ + 1
   20     CONTINUE
      END IF

      IF (JOB.EQ.0) THEN
*
*     PHASE 2 : X:=D**(-1)*X
*
          II = 0
          DO 30 I = 1,N
              II = II + I
              X(I) = X(I)/A(II)
   30     CONTINUE
      END IF

      IF (JOB.LE.0) THEN
*
*     PHASE 3 : X:=TRANS(L)**(-1)*X
*
          II = N* (N-1)/2
          DO 50 I = N - 1,1,-1
              IJ = II
              DO 40 J = I + 1,N
                  IJ = IJ + J - 1
                  X(I) = X(I) - A(IJ)*X(J)
   40         CONTINUE
              II = II - I
   50     CONTINUE
      END IF

      RETURN

      END
* SUBROUTINE MXDPGF                ALL SYSTEMS                89/12/01
* PORTABILITY : ALL SYSTEMS
* 89/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* FACTORIZATION A+E=L*D*TRANS(L) OF A DENSE SYMMETRIC POSITIVE DEFINITE
* MATRIX A+E WHERE D AND E ARE DIAGONAL POSITIVE DEFINITE MATRICES AND
* L IS A LOWER TRIANGULAR MATRIX. IF A IS SUFFICIENTLY POSITIVE
* DEFINITE THEN E=0.
*
* PARAMETERS :
*  II  N ORDER OF THE MATRIX A.
*  RU  A(N*(N+1)/2)  ON INPUT A GIVEN DENSE SYMMETRIC (USUALLY POSITIVE
*         DEFINITE) MATRIX A STORED IN THE PACKED FORM. ON OUTPUT THE
*         COMPUTED FACTORIZATION A+E=L*D*TRANS(L).
*  IO  INF  AN INFORMATION OBTAINED IN THE FACTORIZATION PROCESS. IF
*         INF=0 THEN A IS SUFFICIENTLY POSITIVE DEFINITE AND E=0. IF
*         INF<0 THEN A IS NOT SUFFICIENTLY POSITIVE DEFINITE AND E>0. IF
*         INF>0 THEN A IS INDEFINITE AND INF IS AN INDEX OF THE
*         MOST NEGATIVE DIAGONAL ELEMENT USED IN THE FACTORIZATION
*         PROCESS.
*  RU  ALF  ON INPUT A DESIRED TOLERANCE FOR POSITIVE DEFINITENESS. ON
*         OUTPUT THE MOST NEGATIVE DIAGONAL ELEMENT USED IN THE
*         FACTORIZATION PROCESS (IF INF>0).
*  RO  TAU  MAXIMUM DIAGONAL ELEMENT OF THE MATRIX E.
*
* METHOD :
* P.E.GILL, W.MURRAY : NEWTON TYPE METHODS FOR UNCONSTRAINED AND
* LINEARLY CONSTRAINED OPTIMIZATION, MATH. PROGRAMMING 28 (1974)
* PP. 311-350.
*
      SUBROUTINE MXDPGF(N,A,INF,ALF,TAU)
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALF,TAU
      INTEGER          INF,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(*)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION BET,DEL,GAM,RHO,SIG,TOL
      INTEGER          I,IJ,IK,J,K,KJ,KK,L
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC        ABS,MAX
*     ..
      L = 0
      INF = 0
      TOL = ALF
*
*     ESTIMATION OF THE MATRIX NORM
*
      ALF = 0.0D0
      BET = 0.0D0
      GAM = 0.0D0
      TAU = 0.0D0
      KK = 0
      DO 20 K = 1,N
          KK = KK + K
          BET = MAX(BET,ABS(A(KK)))
          KJ = KK
          DO 10 J = K + 1,N
              KJ = KJ + J - 1
              GAM = MAX(GAM,ABS(A(KJ)))
   10     CONTINUE
   20 CONTINUE
      BET = MAX(TOL,BET,GAM/N)
*      DEL = TOL*BET
      DEL = TOL*MAX(BET,1.0D0)
      KK = 0
      DO 60 K = 1,N
          KK = KK + K
*
*     DETERMINATION OF A DIAGONAL CORRECTION
*
          SIG = A(KK)
          IF (ALF.GT.SIG) THEN
              ALF = SIG
              L = K
          END IF

          GAM = 0.0D0
          KJ = KK
          DO 30 J = K + 1,N
              KJ = KJ + J - 1
              GAM = MAX(GAM,ABS(A(KJ)))
   30     CONTINUE
          GAM = GAM*GAM
          RHO = MAX(ABS(SIG),GAM/BET,DEL)
          IF (TAU.LT.RHO-SIG) THEN
              TAU = RHO - SIG
              INF = -1
          END IF
*
*     GAUSSIAN ELIMINATION
*
          A(KK) = RHO
          KJ = KK
          DO 50 J = K + 1,N
              KJ = KJ + J - 1
              GAM = A(KJ)
              A(KJ) = GAM/RHO
              IK = KK
              IJ = KJ
              DO 40 I = K + 1,J
                  IK = IK + I - 1
                  IJ = IJ + 1
                  A(IJ) = A(IJ) - A(IK)*GAM
   40         CONTINUE
   50     CONTINUE
   60 CONTINUE
      IF (L.GT.0 .AND. ABS(ALF).GT.DEL) INF = L
      RETURN

      END
* SUBROUTINE MXDPRB                ALL SYSTEMS                89/12/01
* PORTABILITY : ALL SYSTEMS
* 89/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* SOLUTION OF A SYSTEM OF LINEAR EQUATIONS WITH A DENSE SYMMETRIC
* POSITIVE DEFINITE MATRIX A USING THE FACTORIZATION A=TRANS(R)*R.
*
* PARAMETERS :
*  II  N ORDER OF THE MATRIX A.
*  RI  A(N*(N+1)/2) FACTORIZATION A=TRANS(R)*R.
*  RU  X(N)  ON INPUT THE RIGHT HAND SIDE OF A SYSTEM OF LINEAR
*         EQUATIONS. ON OUTPUT THE SOLUTION OF A SYSTEM OF LINEAR
*         EQUATIONS.
*  II  JOB  OPTION. IF JOB=0 THEN X:=A**(-1)*X. IF JOB>0 THEN
*         X:=TRANS(R)**(-1)*X. IF JOB<0 THEN X:=R**(-1)*X.
*
* METHOD :
* BACK SUBSTITUTION
*
      SUBROUTINE MXDPRB(N,A,X,JOB)
*     .. Scalar Arguments ..
      INTEGER          JOB,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(*),X(*)
*     ..
*     .. Local Scalars ..
      INTEGER          I,II,IJ,J
*     ..
      IF (JOB.GE.0) THEN
*
*     PHASE 1 : X:=TRANS(R)**(-1)*X
*
          IJ = 0
          DO 20 I = 1,N
              DO 10 J = 1,I - 1
                  IJ = IJ + 1
                  X(I) = X(I) - A(IJ)*X(J)
   10         CONTINUE
              IJ = IJ + 1
              X(I) = X(I)/A(IJ)
   20     CONTINUE
      END IF

      IF (JOB.LE.0) THEN
*
*     PHASE 2 : X:=R**(-1)*X
*
          II = N* (N+1)/2
          DO 40 I = N,1,-1
              IJ = II
              DO 30 J = I + 1,N
                  IJ = IJ + J - 1
                  X(I) = X(I) - A(IJ)*X(J)
   30         CONTINUE
              X(I) = X(I)/A(II)
              II = II - I
   40     CONTINUE
      END IF

      RETURN

      END
* SUBROUTINE MXDSMM                ALL SYSTEMS                89/12/01
* PORTABILITY : ALL SYSTEMS
* 89/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* MULTIPLICATION OF A DENSE SYMMETRIC MATRIX A BY A VECTOR X.
*
* PARAMETERS :
*  II  N  ORDER OF THE MATRIX A.
*  RI  A(N*(N+1)/2)  DENSE SYMMETRIC MATRIX STORED IN THE PACKED FORM.
*  RI  X(N)  INPUT VECTOR.
*  RO  Y(N)  OUTPUT VECTOR EQUAL TO  A*X.
*
      SUBROUTINE MXDSMM(N,A,X,Y)
*     .. Scalar Arguments ..
      INTEGER          N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(*),X(*),Y(*)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER          I,J,K,L
*     ..
      K = 0
      DO 30 I = 1,N
          TEMP = 0.0D0
          L = K
          DO 10 J = 1,I
              L = L + 1
              TEMP = TEMP + A(L)*X(J)
   10     CONTINUE
          DO 20 J = I + 1,N
              L = L + J - 1
              TEMP = TEMP + A(L)*X(J)
   20     CONTINUE
          Y(I) = TEMP
          K = K + I
   30 CONTINUE
      RETURN

      END
* SUBROUTINE MXDSMV                ALL SYSTEMS                91/12/01
* PORTABILITY : ALL SYSTEMS
* 91/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* K-TH ROW OF A DENSE SYMMETRIC MATRIX A IS COPIED TO THE VECTOR X.
*
* PARAMETERS :
*  II  N  ORDER OF THE MATRIX A.
*  RI  A(N*(N+1)/2)  DENSE SYMMETRIC MATRIX STORED IN THE PACKED FORM.
*  RO  X(N)  OUTPUT VECTOR.
*  II  K  INDEX OF COPIED ROW.
*
      SUBROUTINE MXDSMV(N,A,X,K)
*     .. Scalar Arguments ..
      INTEGER          K,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(*),X(*)
*     ..
*     .. Local Scalars ..
      INTEGER          I,L
*     ..
      L = K* (K-1)/2
      DO 10 I = 1,N
          IF (I.LE.K) THEN
              L = L + 1

          ELSE
              L = L + I - 1
          END IF

          X(I) = A(L)
   10 CONTINUE
      RETURN

      END
* SUBROUTINE MXVCOP                ALL SYSTEMS                88/12/01
* PORTABILITY : ALL SYSTEMS
* 88/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* COPYING OF A VECTOR.
*
* PARAMETERS :
*  II  N  VECTOR DIMENSION.
*  RI  X(N)  INPUT VECTOR.
*  RO  Y(N)  OUTPUT VECTOR WHERE Y:= X.
*
      SUBROUTINE MXVCOP(N,X,Y)
*     .. Scalar Arguments ..
      INTEGER          N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION X(*),Y(*)
*     ..
*     .. Local Scalars ..
      INTEGER          I
*     ..
      DO 10 I = 1,N
          Y(I) = X(I)
   10 CONTINUE
      RETURN

      END
* SUBROUTINE MXVDIR                ALL SYSTEMS                91/12/01
* PORTABILITY : ALL SYSTEMS
* 91/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* VECTOR AUGMENTED BY THE SCALED VECTOR.
*
* PARAMETERS :
*  II  N  VECTOR DIMENSION.
*  RI  A  SCALING FACTOR.
*  RI  X(N)  INPUT VECTOR.
*  RI  Y(N)  INPUT VECTOR.
*  RO  Z(N)  OUTPUT VECTOR WHERE Z:= Y + A*X.
*
      SUBROUTINE MXVDIR(N,A,X,Y,Z)
*     .. Scalar Arguments ..
      DOUBLE PRECISION A
      INTEGER          N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION X(*),Y(*),Z(*)
*     ..
*     .. Local Scalars ..
      INTEGER          I
*     ..
      DO 10 I = 1,N
          Z(I) = Y(I) + A*X(I)
   10 CONTINUE
      RETURN

      END
* FUNCTION MXVDOT                  ALL SYSTEMS                91/12/01
* PORTABILITY : ALL SYSTEMS
* 91/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* DOT PRODUCT OF TWO VECTORS.
*
* PARAMETERS :
*  II  N  VECTOR DIMENSION.
*  RI  X(N)  INPUT VECTOR.
*  RI  Y(N)  INPUT VECTOR.
*  RR  MXVDOT  VALUE OF DOT PRODUCT MXVDOT=TRANS(X)*Y.
*
      DOUBLE PRECISION
     +  FUNCTION MXVDOT(N,X,Y)
*     .. Scalar Arguments ..
      INTEGER          N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION X(*),Y(*)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER          I
*     ..
      TEMP = 0.0D0
      DO 10 I = 1,N
          TEMP = TEMP + X(I)*Y(I)
   10 CONTINUE
      MXVDOT = TEMP
      RETURN

      END
* SUBROUTINE MXVINA             ALL SYSTEMS                   90/12/01
* PORTABILITY : ALL SYSTEMS
* 90/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* ELEMENTS OF THE INTEGER VECTOR ARE REPLACED BY THEIR ABSOLUTE VALUES.
*
* PARAMETERS :
*  II  N DIMENSION OF THE INTEGER VECTOR.
*  IU  IX(N)  INTEGER VECTOR WHICH IS UPDATED SO THAT IX(I):=ABS(IX(I))
*         FOR ALL I.
*
      SUBROUTINE MXVINA(N,IX)
*     .. Scalar Arguments ..
      INTEGER          N
*     ..
*     .. Array Arguments ..
      INTEGER          IX(*)
*     ..
*     .. Local Scalars ..
      INTEGER          I
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC        ABS
*     ..
      DO 10 I = 1,N
          IX(I) = ABS(IX(I))
          IF (IX(I).GT.10) IX(I) = IX(I) - 10
   10 CONTINUE
      RETURN

      END
* SUBROUTINE MXVINV               ALL SYSTEMS                91/12/01
* PORTABILITY : ALL SYSTEMS
* 91/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* CHANGE OF THE INTEGER VECTOR ELEMENT FOR THE CONSTRAINT ADDITION.
*
* PARAMETERS :
*  II  N  VECTOR DIMENSION.
*  IU  IX(N)  INTEGER VECTOR.
*  II  I  INDEX OF THE CHANGED ELEMENT.
*  II  JOB  CHANGE SPECIFICATION
*
      SUBROUTINE MXVINV(IX,I,JOB)
*     .. Scalar Arguments ..
      INTEGER          I,JOB
*     ..
*     .. Array Arguments ..
      INTEGER          IX(*)
*     ..
      IF ((IX(I).EQ.3.OR.IX(I).EQ.5) .AND. JOB.LT.0) IX(I) = IX(I) + 1
      IF ((IX(I).EQ.4.OR.IX(I).EQ.6) .AND. JOB.GT.0) IX(I) = IX(I) - 1
      IX(I) = -IX(I)
      RETURN

      END
* FUNCTION MXVMAX               ALL SYSTEMS                   91/12/01
* PORTABILITY : ALL SYSTEMS
* 91/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* L-INFINITY NORM OF A VECTOR.
*
* PARAMETERS :
*  II  N  VECTOR DIMENSION.
*  RI  X(N)  INPUT VECTOR.
*  RR  MXVMAX  L-INFINITY NORM OF THE VECTOR X.
*
      DOUBLE PRECISION
     +  FUNCTION MXVMAX(N,X)
*     .. Scalar Arguments ..
      INTEGER          N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION X(*)
*     ..
*     .. Local Scalars ..
      INTEGER          I
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC        ABS,MAX
*     ..
      MXVMAX = 0.0D0
      DO 10 I = 1,N
          MXVMAX = MAX(MXVMAX,ABS(X(I)))
   10 CONTINUE
      RETURN

      END
* SUBROUTINE MXVMUL             ALL SYSTEMS                   89/12/01
* PORTABILITY : ALL SYSTEMS
* 89/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* VECTOR IS PREMULTIPLIED BY THE POWER OF A DIAGONAL MATRIX.
*
* PARAMETERS :
*  II  N  VECTOR DIMENSION.
*  RI  D(N)  DIAGONAL MATRIX STORED AS A VECTOR WITH N ELEMENTS.
*  RI  X(N)  INPUT VECTOR.
*  RO  Y(N)  OUTPUT VECTOR WHERE Y:=(D**K)*X.
*  II  K  INTEGER EXPONENT.
*
      SUBROUTINE MXVMUL(N,D,X,Y,K)
*     .. Scalar Arguments ..
      INTEGER          K,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION D(*),X(*),Y(*)
*     ..
*     .. Local Scalars ..
      INTEGER          I
*     ..
*     .. External Subroutines ..
      EXTERNAL         MXVCOP
*     ..
      IF (K.EQ.0) THEN
          CALL MXVCOP(N,X,Y)

      ELSE IF (K.EQ.1) THEN
          DO 10 I = 1,N
              Y(I) = X(I)*D(I)
   10     CONTINUE

      ELSE IF (K.EQ.-1) THEN
          DO 20 I = 1,N
              Y(I) = X(I)/D(I)
   20     CONTINUE

      ELSE
          DO 30 I = 1,N
              Y(I) = X(I)*D(I)**K
   30     CONTINUE
      END IF

      RETURN

      END
* SUBROUTINE MXVNEG                ALL SYSTEMS                88/12/01
* PORTABILITY : ALL SYSTEMS
* 88/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* CHANGE THE SIGNS OF VECTOR ELEMENTS.
*
* PARAMETERS :
*  II  N  VECTOR DIMENSION.
*  RI  X(N)  INPUT VECTOR.
*  RO  Y(N)  OUTPUT VECTOR WHERE Y:= - X.
*
      SUBROUTINE MXVNEG(N,X,Y)
*     .. Scalar Arguments ..
      INTEGER          N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION X(*),Y(*)
*     ..
*     .. Local Scalars ..
      INTEGER          I
*     ..
      DO 10 I = 1,N
          Y(I) = -X(I)
   10 CONTINUE
      RETURN

      END
* SUBROUTINE MXVORT               ALL SYSTEMS                91/12/01
* PORTABILITY : ALL SYSTEMS
* 91/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* DETERMINATION OF AN ELEMENTARY ORTHOGONAL MATRIX FOR PLANE ROTATION.
*
* PARAMETERS :
*  RU  XK  FIRST VALUE FOR PLANE ROTATION (XK IS TRANSFORMED TO
*         SQRT(XK**2+XL**2))
*  RU  XL  SECOND VALUE FOR PLANE ROTATION (XL IS TRANSFORMED TO
*         ZERO)
*  RO  CK  DIAGONAL ELEMENT OF THE ELEMENTARY ORTHOGONAL MATRIX.
*  RO  CL  OFF-DIAGONAL ELEMENT OF THE ELEMENTARY ORTHOGONAL MATRIX.
*  IO  IER  INFORMATION ON THE TRANSFORMATION. IER=0-GENERAL PLANE
*         ROTATION. IER=1-PERMUTATION. IER=2-TRANSFORMATION SUPPRESSED.
*
      SUBROUTINE MXVORT(XK,XL,CK,CL,IER)
*     .. Scalar Arguments ..
      DOUBLE PRECISION CK,CL,XK,XL
      INTEGER          IER
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION DEN,POM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC        ABS,SQRT
*     ..
      IF (XL.EQ.0.0D0) THEN
          IER = 2

      ELSE IF (XK.EQ.0.0D0) THEN
          XK = XL
          XL = 0.0D0
          IER = 1

      ELSE
          IF (ABS(XK).GE.ABS(XL)) THEN
              POM = XL/XK
              DEN = SQRT(1.0D0+POM*POM)
              CK = 1.0D0/DEN
              CL = POM/DEN
              XK = XK*DEN

          ELSE
              POM = XK/XL
              DEN = SQRT(1.0D0+POM*POM)
              CL = 1.0D0/DEN
              CK = POM/DEN
              XK = XL*DEN
          END IF

          XL = 0.0D0
          IER = 0
      END IF

      RETURN

      END
* SUBROUTINE MXVROT               ALL SYSTEMS                91/12/01
* PORTABILITY : ALL SYSTEMS
* 91/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* PLANE ROTATION IS APPLIED TO TWO VALUES.
*
* PARAMETERS :
*  RU  XK  FIRST VALUE FOR PLANE ROTATION.
*  RU  XL  SECOND VALUE FOR PLANE ROTATION.
*  RI  CK  DIAGONAL ELEMENT OF THE ELEMENTARY ORTHOGONAL MATRIX.
*  RI  CL  OFF-DIAGONAL ELEMENT OF THE ELEMENTARY ORTHOGONAL MATRIX.
*  II  IER  INFORMATION ON THE TRANSFORMATION. IER=0-GENERAL PLANE
*         ROTATION. IER=1-PERMUTATION. IER=2-TRANSFORMATION SUPPRESSED.
*
      SUBROUTINE MXVROT(XK,XL,CK,CL,IER)
*     .. Scalar Arguments ..
      DOUBLE PRECISION CK,CL,XK,XL
      INTEGER          IER
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION YK,YL
*     ..
      IF (IER.EQ.0) THEN
          YK = XK
          YL = XL
          XK = CK*YK + CL*YL
          XL = CL*YK - CK*YL

      ELSE IF (IER.EQ.1) THEN
          YK = XK
          XK = XL
          XL = YK
      END IF

      RETURN

      END
* SUBROUTINE MXVSET                ALL SYSTEMS                88/12/01
* PORTABILITY : ALL SYSTEMS
* 88/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* A SCALAR IS SET TO ALL THE ELEMENTS OF A VECTOR.
*
* PARAMETERS :
*  II  N  VECTOR DIMENSION.
*  RI  A  INITIAL VALUE.
*  RO  X(N)  OUTPUT VECTOR SUCH THAT X(I)=A FOR ALL I.
*
      SUBROUTINE MXVSET(N,A,X)
*     .. Scalar Arguments ..
      DOUBLE PRECISION A
      INTEGER          N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION X(*)
*     ..
*     .. Local Scalars ..
      INTEGER          I
*     ..
      DO 10 I = 1,N
          X(I) = A
   10 CONTINUE
      RETURN

      END








