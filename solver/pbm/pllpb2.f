* SUBROUTINE PLLPB2             ALL SYSTEMS                   91/12/01
* PORTABILITY : ALL SYSTEMS
* 91/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* DETERMINATION OF THE INITIAL FEASIBLE POINT AND THE LINEAR PROGRAMMING
* SUBROUTINE.
*
* PARAMETERS :
*  II  NF  NUMBER OF VARIABLES.
*  II  NC  NUMBER OF LINEAR CONSTRAINTS.
*  RI  X(NF)  VECTOR OF VARIABLES.
*  II  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS.
*  RO  XO(NF)  SAVED VECTOR OF VARIABLES.
*  RI  XL(NF)  VECTOR CONTAINING LOWER BOUNDS FOR VARIABLES.
*  RI  XU(NF)  VECTOR CONTAINING UPPER BOUNDS FOR VARIABLES.
*  RI  CF(NF)  VECTOR CONTAINING VALUES OF THE CONSTRAINT FUNCYIONS.
*  RO  CFD(NF)  VECTOR CONTAINING INCREMENTS OF THE CONSTRAINT
*         FUNCTIONS.
*  II  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
*  II  ICA(NF)  VECTOR CONTAINING INDICES OF ACTIVE CONSTRAINTS.
*  RI  CL(NC)  VECTOR CONTAINING LOWER BOUNDS FOR CONSTRAINT FUNCTIONS.
*  RI  CU(NC)  VECTOR CONTAINING UPPER BOUNDS FOR CONSTRAINT FUNCTIONS.
*  RI  CG(NF*NC)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
*         CONSTRAINTS.
*  RI  CR(NF*(NF+1)/2)  TRIANGULAR DECOMPOSITION OF KERNEL OF THE
*         ORTHOGONAL PROJECTION.
*  RO  CZ(NF)  VECTOR OF LAGRANGE MULTIPLIERS.
*  RI  G(NF)  GRADIENT OF THE OBJECTIVE FUNCTION.
*  RO  GO(NF)  SAVED GRADIENT OF THE OBJECTIVE FUNCTION.
*  RI  S(NF)  DIRECTION VECTOR.
*  II  MFP  TYPE OF FEASIBLE POINT. MFP=1-ARBITRARY FEASIBLE POINT.
*         MFP=2-OPTIMUM FEASIBLE POINT. MFP=3-REPEATED SOLUTION.
*  II  KBF  SPECIFICATION OF SIMPLE BOUNDS. KBF=0-NO SIMPLE BOUNDS.
*         KBF=1-ONE SIDED SIMPLE BOUNDS. KBF=2=TWO SIDED SIMPLE BOUNDS.
*  II  KBC  SPECIFICATION OF LINEAR CONSTRAINTS. KBC=0-NO LINEAR
*         CONSTRAINTS. KBC=1-ONE SIDED LINEAR CONSTRAINTS. KBC=2=TWO
*         SIDED LINEAR CONSTRAINTS.
*  RI  ETA9  MAXIMUM FOR REAL NUMBERS.
*  RI  EPS7  TOLERANCE FOR LINEAR INDEPENDENCE OF CONSTRAINTS.
*  RI  EPS9  TOLERANCE FOR ACTIVITY OF CONSTRAINTS.
*  RO  UMAX  MAXIMUM ABSOLUTE VALUE OF A NEGATIVE LAGRANGE MULTIPLIER.
*  RO  GMAX  MAXIMUM ABSOLUTE VALUE OF A PARTIAL DERIVATIVE.
*  IO  N  DIMENSION OF THE MANIFOLD DEFINED BY ACTIVE CONSTRAINTS.
*  IO  ITERL  TYPE OF FEASIBLE POINT. ITERL=1-ARBITRARY FEASIBLE POINT.
*         ITERL=2-OPTIMUM FEASIBLE POINT. ITERL=-1 FEASIBLE POINT DOES
*         NOT EXISTS. ITERL=-2 OPTIMUM FEASIBLE POINT DOES NOT EXISTS.
*
* SUBPROGRAMS USED :
*  S   PLINIT  DETERMINATION OF INITIAL POINT SATISFYING SIMPLE BOUNDS.
*  S   PLMAXL  MAXIMUM STEPSIZE USING LINEAR CONSTRAINTS.
*  S   PLMAXS  MAXIMUM STEPSIZE USING SIMPLE BOUNDS.
*  S   PLMAXT  MAXIMUM STEPSIZE USING TRUST REGION BOUNDS.
*  S   PLNEWL  IDENTIFICATION OF ACTIVE LINEAR CONSTRAINTS.
*  S   PLNEWS  IDENTIFICATION OF ACTIVE SIMPLE BOUNDS.
*  S   PLNEWT  IDENTIFICATION OF ACTIVE TRUST REGION BOUNDS.
*  S   PLDIRL  NEW VALUES OF CONSTRAINT FUNCTIONS.
*  S   PLDIRS  NEW VALUES OF VARIABLES.
*  S   PLSETC  INITIAL VALUES OF CONSTRAINT FUNCTIONS.
*  S   PLSETG  DETERMINATION OF THE FIRST PHASE GRADIENT VECTOR.
*  S   PLGLAG  GRADIENT OF THE LAGRANGIAN FUNCTION IS DETERMINED.
*  S   PLSLAG  NEGATIVE PROJECTED GRADIENT IS DETERMINED.
*  S   PLTLAG  THE OPTIMUM LAGRANGE MULTIPLIER IS DETERMINED.
*  S   PLVLAG  AN AUXILIARY VECTOR IS DETERMINED.
*  S   PLADR0  CONSTRAINT ADDITION.
*  S   PLRMF0  CONSTRAINT DELETION.
*  S   MXDPRB  BACK SUBSTITUTION AFTER CHOLESKI DECOMPOSITION.
*  S   MXDSMI  DETERMINATION OF THE INITIAL UNIT DENSE SYMMETRIC
*         MATRIX.
*  S   MXVCOP  COPYING OF A VECTOR.
*  S   MXVDIF  DIFFERENCE OF TWO VECTORS.
*  S   MXVINA  ABSOLUTE VALUES OF ELEMENTS OF AN INTEGER VECTOR.
*  S   MXVINC  UPDATE OF AN INTEGER VECTOR.
*  S   MXVIND  CHANGE OF THE INTEGER VECTOR FOR CONSTRAINT ADDITION.
*  S   MXVINT  CHANGE OF THE INTEGER VECTOR FOR TRUST REGION BOUND
*         ADDITION.
*  S   MXVMUL  DIAGONAL PREMULTIPLICATION OF A VECTOR.
*  S   MXVNEG  COPYING OF A VECTOR WITH CHANGE OF THE SIGN.
*  S   MXVSET  INITIATION OF A VECTOR.
*
      SUBROUTINE PLLPB2(NF,NC,X,IX,XO,XL,XU,CF,CFD,IC,ICA,CL,CU,CG,CR,
     +                  CZ,G,GO,S,MFP,KBF,KBC,ETA9,EPS7,EPS9,UMAX,GMAX,
     +                  N,ITERL)
C     .. Scalar Arguments ..
      DOUBLE PRECISION EPS7,EPS9,ETA9,GMAX,UMAX
      INTEGER ITERL,KBC,KBF,MFP,N,NC,NF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CF(*),CFD(*),CG(*),CL(*),CR(*),CU(*),CZ(*),G(*),
     +                 GO(*),S(*),X(*),XL(*),XO(*),XU(*)
      INTEGER IC(*),ICA(*),IX(*)
C     ..
C     .. Scalars in Common ..
      INTEGER NADD,NDECF,NFG,NFH,NFV,NIT,NRED,NREM,NRES
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CON,DMAX,POM
      INTEGER I,IER,INEW,IOLD,IPOM,KC,KREM,MODE
C     ..
C     .. External Subroutines ..
      EXTERNAL MXDPRB,MXDSMI,MXVCOP,MXVINA,MXVIND,MXVNEG,PLADR0,PLDIRL,
     +         PLDIRS,PLINIT,PLMAXL,PLMAXS,PLNEWL,PLNEWS,PLRMF0,PLSETC,
     +         PLSETG,PLSLAG,PLTLAG,PLVLAG
C     ..
C     .. Common blocks ..
      COMMON /STAT/NDECF,NRES,NRED,NREM,NADD,NIT,NFV,NFG,NFH
C     ..
      CON = ETA9
*
*     INITIATION
*
      CALL MXVCOP(NF,X,XO)
      CALL MXVCOP(NF,G,GO)
      IPOM = 0
      NRED = 0
      KREM = 0
      ITERL = 1
      DMAX = 0.0D0
      IF (MFP.EQ.3) GO TO 40
      IF (KBF.GT.0) CALL MXVINA(NF,IX)
*
*     SHIFT OF VARIABLES FOR SATISFYING SIMPLE BOUNDS
*
      CALL PLINIT(NF,X,IX,XL,XU,EPS9,KBF,INEW,ITERL)
      IF (ITERL.LT.0) THEN
          GO TO 60

      END IF

      N = NF
      DO 10 I = 1,NF
          IF (KBF.GT.0 .AND. IX(I).LT.0) THEN
              N = N - 1
              ICA(NF-N) = -I
          END IF

   10 CONTINUE
      CALL MXDSMI(NF-N,CR)
      IF (NC.GT.0) THEN
*
*     ADDITION OF ACTIVE CONSTRAINTS AND INITIAL CHECK OF FEASIBILITY
*
          CALL MXVINA(NC,IC)
          IF (NF.GT.N) CALL PLSETC(NF,NC,X,XO,CF,IC,CG,S)
          DO 20 KC = 1,NC
              IF (IC(KC).NE.0) THEN
                  INEW = 0
                  CALL PLNEWL(KC,CF,IC,CL,CU,EPS9,INEW)
                  CALL PLADR0(NF,N,ICA,CG,CR,S,EPS7,GMAX,UMAX,INEW,NADD,
     +                        IER)
                  CALL MXVIND(IC,KC,IER)
                  IF (IC(KC).LT.-10) IPOM = 1
              END IF

   20     CONTINUE
      END IF

   30 IF (IPOM.EQ.1) THEN
*
*     CHECK OF FEASIBILITY AND UPDATE OF THE FIRST PHASE OBJECTIVE
*     FUNCTION
*
          CALL PLSETG(NF,NC,IC,CG,G,INEW)
          IF (INEW.EQ.0) IPOM = 0
      END IF

      IF (IPOM.EQ.0 .AND. ITERL.EQ.0) THEN
*
*     FEASIBILITY ACHIEVED
*
          ITERL = 1
          CALL MXVCOP(NF,GO,G)
          IF (MFP.EQ.1) GO TO 60
      END IF
*
*     LAGRANGE MULTIPLIERS DETERMINATION
*
   40 IF (NF.GT.N) THEN
          CALL PLVLAG(NF,N,NC,ICA,CG,CG,G,CZ)
          CALL MXDPRB(NF-N,CR,CZ,0)
          CALL PLTLAG(NF,N,NC,IX,IC,ICA,CZ,IC,EPS7,UMAX,IOLD)

      ELSE
          IOLD = 0
          UMAX = 0.0D0
      END IF
*
*     PROJECTED GRADIENT DETERMINATION
*
      IF (N.GT.0) THEN
          CALL MXVNEG(NF,G,S)
          CALL PLSLAG(NF,N,NC,ICA,CG,CZ,CG,S,EPS7,GMAX)

      ELSE
          GMAX = 0.0D0
      END IF

      MODE = 1 - IPOM
      INEW = 0
      IF (GMAX.EQ.0.0D0) THEN
*
*     OPTIMUM ON A LINEAR MANIFOLD OBTAINED
*
          IF (IOLD.EQ.0) THEN
              IF (IPOM.EQ.0) THEN
*
*     OPTIMAL SOLUTION ACHIEVED
*
                  ITERL = 2
                  GO TO 60

              ELSE
                  IPOM = 0
                  DO 50 KC = 1,NC
                      IF (IC(KC).LT.-10) THEN
                          INEW = 0
                          CALL PLNEWL(KC,CF,IC,CL,CU,EPS9,INEW)
                          IF (IC(KC).LT.-10) IPOM = 1
                      END IF

   50             CONTINUE
                  IF (IPOM.EQ.0) THEN
*
*     OPTIMAL SOLUTION ACHIEVED
*
                      CALL MXVCOP(NF,GO,G)
                      ITERL = 2
                      GO TO 60

                  ELSE
*
*     FEASIBLE SOLUTION DOES NOT EXIST
*
                      CALL MXVCOP(NF,GO,G)
                      ITERL = -1
                      GO TO 60

                  END IF

              END IF

          ELSE
*
*     CONSTRAINT DELETION
*
              CALL PLRMF0(NF,NC,IX,IC,ICA,CR,IC,S,N,IOLD,KREM,IER)
              DMAX = 0.0D0
              GO TO 40

          END IF

      ELSE
*
*     STEPSIZE SELECTION
*
          POM = CON
          CALL PLMAXL(NF,NC,CF,CFD,IC,CL,CU,CG,S,POM,KBC,KREM,INEW)
          CALL PLMAXS(NF,X,IX,XL,XU,S,POM,KBF,KREM,INEW)
          IF (INEW.EQ.0) THEN
              IF (IPOM.EQ.0) THEN
*
*     BOUNDED SOLUTION DOES NOT EXIST
*
                  ITERL = -2

              ELSE
*
*     FEASIBLE SOLUTION DOES NOT EXIST
*
                  ITERL = -3
              END IF

              GO TO 60

          ELSE
*
*     STEP REALIZATION
*
              CALL PLDIRS(NF,X,IX,S,POM,KBF)
              CALL PLDIRL(NC,CF,CFD,IC,POM,KBC)
*
*     CONSTRAINT ADDITION
*
              IF (INEW.GT.0) THEN
                  KC = INEW
                  INEW = 0
                  CALL PLNEWL(KC,CF,IC,CL,CU,EPS9,INEW)
                  CALL PLADR0(NF,N,ICA,CG,CR,S,EPS7,GMAX,UMAX,INEW,NADD,
     +                        IER)
                  CALL MXVIND(IC,KC,IER)

              ELSE IF (INEW+NF.GE.0) THEN
                  I = -INEW
                  INEW = 0
                  CALL PLNEWS(X,IX,XL,XU,EPS9,I,INEW)
                  CALL PLADR0(NF,N,ICA,CG,CR,S,EPS7,GMAX,UMAX,INEW,NADD,
     +                        IER)
                  CALL MXVIND(IX,I,IER)
              END IF

              DMAX = POM
              NRED = NRED + 1
              GO TO 30

          END IF

      END IF

   60 CONTINUE
      RETURN

      END
* SUBROUTINE PLDIRL               ALL SYSTEMS                   97/12/01
* PORTABILITY : ALL SYSTEMS
* 97/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* DETERMINATION OF THE NEW VALUES OF THE CONSTRAINT FUNCTIONS.
*
* PARAMETERS :
*  II  NC  NUMBER OF CONSTRAINTS.
*  RU  CF(NF)  VECTOR CONTAINING VALUES OF THE CONSTRAINT FUNCTIONS.
*  RI  CFD(NF)  VECTOR CONTAINING INCREMENTS OF THE CONSTRAINT
*         FUNCTIONS.
*  II  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
*  RI  STEP  CURRENT STEPSIZE.
*  II  KBC  SPECIFICATION OF LINEAR CONSTRAINTS. KBC=0-NO LINEAR
*         CONSTRAINTS. KBC=1-ONE SIDED LINEAR CONSTRAINTS. KBC=2=TWO
*         SIDED LINEAR CONSTRAINTS.
*
      SUBROUTINE PLDIRL(NC,CF,CFD,IC,STEP,KBC)
C     .. Scalar Arguments ..
      DOUBLE PRECISION STEP
      INTEGER KBC,NC
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CF(*),CFD(*)
      INTEGER IC(*)
C     ..
C     .. Local Scalars ..
      INTEGER KC
C     ..
      IF (KBC.GT.0) THEN
          DO 10 KC = 1,NC
              IF (IC(KC).GE.0 .AND. IC(KC).LE.10) THEN
                  CF(KC) = CF(KC) + STEP*CFD(KC)

              ELSE IF (IC(KC).LT.-10) THEN
                  CF(KC) = CF(KC) + STEP*CFD(KC)
              END IF

   10     CONTINUE
      END IF

      RETURN

      END
* SUBROUTINE PLDIRS               ALL SYSTEMS                   97/12/01
* PORTABILITY : ALL SYSTEMS
* 97/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* DETERMINATION OF THE NEW VECTOR OF VARIABLES.
*
* PARAMETERS :
*  II  NF  NUMBER OF VARIABLES.
*  RU  X(NF)  VECTOR OF VARIABLES.
*  II  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS.
*  RI  S(NF)  DIRECTION VECTOR.
*  RI  STEP  CURRENT STEPSIZE.
*  II  KBF  SPECIFICATION OF SIMPLE BOUNDS. KBF=0-NO SIMPLE BOUNDS.
*         KBF=1-ONE SIDED SIMPLE BOUNDS. KBF=2=TWO SIDED SIMPLE BOUNDS.
*
      SUBROUTINE PLDIRS(NF,X,IX,S,STEP,KBF)
C     .. Scalar Arguments ..
      DOUBLE PRECISION STEP
      INTEGER KBF,NF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION S(*),X(*)
      INTEGER IX(*)
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
      DO 10 I = 1,NF
          IF (KBF.LE.0) THEN
              X(I) = X(I) + STEP*S(I)

          ELSE IF (IX(I).GE.0 .AND. IX(I).LE.10) THEN
              X(I) = X(I) + STEP*S(I)

          ELSE IF (IX(I).LT.-10) THEN
              X(I) = X(I) + STEP*S(I)
          END IF

   10 CONTINUE
      RETURN

      END
* SUBROUTINE PLINIT             ALL SYSTEMS                   97/12/01
* PORTABILITY : ALL SYSTEMS
* 97/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* DETERMINATION OF THE INITIAL POINT WHICH SATISFIES SIMPLE BOUNDS.
*
* PARAMETERS :
*  II  NF  NUMBER OF VARIABLES.
*  RU  X(NF)  VECTOR OF VARIABLES.
*  II  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS.
*  RI  XL(NF)  VECTOR CONTAINING LOWER BOUNDS FOR VARIABLES.
*  RI  XU(NF)  VECTOR CONTAINING UPPER BOUNDS FOR VARIABLES.
*  RI  EPS9  TOLERANCE FOR ACTIVE CONSTRAINTS.
*  IO  INEW  INDEX OF THE NEW ACTIVE CONSTRAINT.
*  IO  IND  INDICATOR. IF IND.NE.0 THEN TRUST REGION BOUNDS CANNOT
*         BE SATISFIED.
*
* SUBPROGRAMS USED :
*  S   PLNEWS  TEST ON ACTIVITY OF A GIVEN SIMPLE BOUND.
*
      SUBROUTINE PLINIT(NF,X,IX,XL,XU,EPS9,KBF,INEW,IND)
C     .. Scalar Arguments ..
      DOUBLE PRECISION EPS9
      INTEGER IND,INEW,KBF,NF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION X(*),XL(*),XU(*)
      INTEGER IX(*)
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
C     .. External Subroutines ..
      EXTERNAL PLNEWS
C     ..
      IND = 0
      IF (KBF.GT.0) THEN
          DO 10 I = 1,NF
              CALL PLNEWS(X,IX,XL,XU,EPS9,I,INEW)
              IF (IX(I).LT.5) THEN

              ELSE IF (IX(I).EQ.5) THEN
                  IX(I) = -5

              ELSE IF (IX(I).EQ.11 .OR. IX(I).EQ.13) THEN
                  X(I) = XL(I)
                  IX(I) = 10 - IX(I)

              ELSE IF (IX(I).EQ.12 .OR. IX(I).EQ.14) THEN
                  X(I) = XU(I)
                  IX(I) = 10 - IX(I)
              END IF

   10     CONTINUE
      END IF

      RETURN

      END
* SUBROUTINE PLMAXL               ALL SYSTEMS                   97/12/01
* PORTABILITY : ALL SYSTEMS
* 97/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* DETERMINATION OF THE MAXIMUM STEPSIZE USING LINEAR CONSTRAINTS.
*
* PARAMETERS :
*  II  NF  DECLARED NUMBER OF VARIABLES.
*  II  NC  NUMBER OF CURRENT LINEAR CONSTRAINTS.
*  RI  CF(NF)  VECTOR CONTAINING VALUES OF THE CONSTRAINT FUNCYIONS.
*  RO  CFD(NF)  VECTOR CONTAINING INCREMENTS OF THE CONSTRAINT
*         FUNCTIONS.
*  II  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
*  RI  CL(NC)  VECTOR CONTAINING LOWER BOUNDS FOR CONSTRAINT FUNCTIONS.
*  RI  CU(NC)  VECTOR CONTAINING UPPER BOUNDS FOR CONSTRAINT FUNCTIONS.
*  RI  CG(NF*NC)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
*         CONSTRAINTS.
*  RI  S(NF)  DIRECTION VECTOR.
*  RO  STEP  MAXIMUM STEPSIZE.
*  II  KBC  SPECIFICATION OF LINEAR CONSTRAINTS. KBC=0-NO LINEAR
*         CONSTRAINTS. KBC=1-ONE SIDED LINEAR CONSTRAINTS. KBC=2=TWO
*         SIDED LINEAR CONSTRAINTS.
*  II  KREM  INDICATION OF LINEARLY DEPENDENT GRADIENTS.
*  IO  INEW  INDEX OF THE NEW ACTIVE FUNCTION.
*
* SUBPROGRAMS USED :
*  RF  MXVDOT  DOT PRODUCT OF TWO VECTORS.
*
      SUBROUTINE PLMAXL(NF,NC,CF,CFD,IC,CL,CU,CG,S,STEP,KBC,KREM,INEW)
C     .. Scalar Arguments ..
      DOUBLE PRECISION STEP
      INTEGER INEW,KBC,KREM,NC,NF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CF(*),CFD(*),CG(*),CL(*),CU(*),S(*)
      INTEGER IC(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER JCG,KC
C     ..
C     .. External Functions ..
      DOUBLE PRECISION MXVDOT
      EXTERNAL MXVDOT
C     ..
      IF (KBC.GT.0) THEN
          JCG = 1
          DO 10 KC = 1,NC
              IF (KREM.GT.0 .AND. IC(KC).GT.10) IC(KC) = IC(KC) - 10
              IF (IC(KC).GT.0 .AND. IC(KC).LE.10) THEN
                  TEMP = MXVDOT(NF,CG(JCG),S)
                  CFD(KC) = TEMP
                  IF (TEMP.LT.0.0D0) THEN
                      IF (IC(KC).EQ.1 .OR. IC(KC).GE.3) THEN
                          TEMP = (CL(KC)-CF(KC))/TEMP
                          IF (TEMP.LE.STEP) THEN
                              INEW = KC
                              STEP = TEMP
                          END IF

                      END IF

                  ELSE IF (TEMP.GT.0.0D0) THEN
                      IF (IC(KC).EQ.2 .OR. IC(KC).GE.3) THEN
                          TEMP = (CU(KC)-CF(KC))/TEMP
                          IF (TEMP.LE.STEP) THEN
                              INEW = KC
                              STEP = TEMP
                          END IF

                      END IF

                  END IF

              ELSE IF (IC(KC).LT.-10) THEN
                  TEMP = MXVDOT(NF,CG(JCG),S)
                  CFD(KC) = TEMP
                  IF (TEMP.GT.0.0D0) THEN
                      IF (IC(KC).EQ.-11 .OR. IC(KC).EQ.-13 .OR.
     +                    IC(KC).EQ.-15) THEN
                          TEMP = (CL(KC)-CF(KC))/TEMP
                          IF (TEMP.LE.STEP) THEN
                              INEW = KC
                              STEP = TEMP
                          END IF

                      END IF

                  ELSE IF (TEMP.LT.0.0D0) THEN
                      IF (IC(KC).EQ.-12 .OR. IC(KC).EQ.-14 .OR.
     +                    IC(KC).EQ.-16) THEN
                          TEMP = (CU(KC)-CF(KC))/TEMP
                          IF (TEMP.LE.STEP) THEN
                              INEW = KC
                              STEP = TEMP
                          END IF

                      END IF

                  END IF

              END IF

              JCG = JCG + NF
   10     CONTINUE
      END IF

      RETURN

      END
* SUBROUTINE PLMAXS               ALL SYSTEMS                   97/12/01
* PORTABILITY : ALL SYSTEMS
* 97/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* DETERMINATION OF THE MAXIMUM STEPSIZE USING THE SIMPLE BOUNDS
* FOR VARIABLES.
*
* PARAMETERS :
*  II  NF  NUMBER OF VARIABLES.
*  RI  X(NF)  VECTOR OF VARIABLES.
*  II  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS.
*  RI  XL(NF)  VECTOR CONTAINING LOWER BOUNDS FOR VARIABLES.
*  RI  XU(NF)  VECTOR CONTAINING UPPER BOUNDS FOR VARIABLES.
*  RI  S(NF)  DIRECTION VECTOR.
*  RO  STEP  MAXIMUM STEPSIZE.
*  II  KBF  SPECIFICATION OF SIMPLE BOUNDS. KBF=0-NO SIMPLE BOUNDS.
*         KBF=1-ONE SIDED SIMPLE BOUNDS. KBF=2=TWO SIDED SIMPLE BOUNDS.
*  IO  KREM  INDICATION OF LINEARLY DEPENDENT GRADIENTS.
*  IO  INEW  INDEX OF THE NEW ACTIVE CONSTRAINT.
*
      SUBROUTINE PLMAXS(NF,X,IX,XL,XU,S,STEP,KBF,KREM,INEW)
C     .. Scalar Arguments ..
      DOUBLE PRECISION STEP
      INTEGER INEW,KBF,KREM,NF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION S(*),X(*),XL(*),XU(*)
      INTEGER IX(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I
C     ..
      IF (KBF.GT.0) THEN
          DO 10 I = 1,NF
              IF (KREM.GT.0 .AND. IX(I).GT.10) IX(I) = IX(I) - 10
              IF (IX(I).GT.0 .AND. IX(I).LE.10) THEN
                  IF (S(I).LT.0.0D0) THEN
                      IF (IX(I).EQ.1 .OR. IX(I).GE.3) THEN
                          TEMP = (XL(I)-X(I))/S(I)
                          IF (TEMP.LE.STEP) THEN
                              INEW = -I
                              STEP = TEMP
                          END IF

                      END IF

                  ELSE IF (S(I).GT.0.0D0) THEN
                      IF (IX(I).EQ.2 .OR. IX(I).GE.3) THEN
                          TEMP = (XU(I)-X(I))/S(I)
                          IF (TEMP.LE.STEP) THEN
                              INEW = -I
                              STEP = TEMP
                          END IF

                      END IF

                  END IF

              END IF

   10     CONTINUE
      END IF

      KREM = 0
      RETURN

      END
* SUBROUTINE PLNEWL             ALL SYSTEMS                   97/12/01
* PORTABILITY : ALL SYSTEMS
* 97/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* TEST ON ACTIVITY OF A GIVEN LINEAR CONSTRAINT.
*
* PARAMETERS :
*  II  KC  INDEX OF A GIVEN CONSTRAINT.
*  RI  CF(NC)  VECTOR CONTAINING VALUES OF THE CONSTRAINT FUNCTIONS.
*  IU  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
*  RI  CL(NC)  VECTOR CONTAINING LOWER BOUNDS FOR CONSTRAINT FUNCTIONS.
*  RI  CU(NC)  VECTOR CONTAINING UPPER BOUNDS FOR CONSTRAINT FUNCTIONS.
*  RI  EPS9  TOLERANCE FOR ACTIVE CONSTRAINTS.
*  IO  INEW  INDEX OF THE NEW ACTIVE CONSTRAINT.
*
      SUBROUTINE PLNEWL(KC,CF,IC,CL,CU,EPS9,INEW)
C     .. Scalar Arguments ..
      DOUBLE PRECISION EPS9
      INTEGER INEW,KC
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CF(*),CL(*),CU(*)
      INTEGER IC(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION TEMP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX
C     ..
      IF (IC(KC).LT.-10) IC(KC) = -IC(KC) - 10
      IF (IC(KC).LE.0) THEN

      ELSE IF (IC(KC).EQ.1) THEN
          TEMP = EPS9*MAX(ABS(CL(KC)),1.0D0)
          IF (CF(KC).GT.CL(KC)+TEMP) THEN

          ELSE IF (CF(KC).GE.CL(KC)-TEMP) THEN
              IC(KC) = 11
              INEW = KC

          ELSE
              IC(KC) = -11
          END IF

      ELSE IF (IC(KC).EQ.2) THEN
          TEMP = EPS9*MAX(ABS(CU(KC)),1.0D0)
          IF (CF(KC).LT.CU(KC)-TEMP) THEN

          ELSE IF (CF(KC).LE.CU(KC)+TEMP) THEN
              IC(KC) = 12
              INEW = KC

          ELSE
              IC(KC) = -12
          END IF

      ELSE IF (IC(KC).EQ.3 .OR. IC(KC).EQ.4) THEN
          TEMP = EPS9*MAX(ABS(CL(KC)),1.0D0)
          IF (CF(KC).GT.CL(KC)+TEMP) THEN
              TEMP = EPS9*MAX(ABS(CU(KC)),1.0D0)
              IF (CF(KC).LT.CU(KC)-TEMP) THEN

              ELSE IF (CF(KC).LE.CU(KC)+TEMP) THEN
                  IC(KC) = 14
                  INEW = KC

              ELSE
                  IC(KC) = -14
              END IF

          ELSE IF (CF(KC).GE.CL(KC)-TEMP) THEN
              IC(KC) = 13
              INEW = KC

          ELSE
              IC(KC) = -13
          END IF

      ELSE IF (IC(KC).EQ.5 .OR. IC(KC).EQ.6) THEN
          TEMP = EPS9*MAX(ABS(CL(KC)),1.0D0)
          IF (CF(KC).GT.CL(KC)+TEMP) THEN
              TEMP = EPS9*MAX(ABS(CU(KC)),1.0D0)
              IF (CF(KC).LT.CU(KC)-TEMP) THEN

              ELSE IF (CF(KC).LE.CU(KC)+TEMP) THEN
                  IC(KC) = 16
                  INEW = KC

              ELSE
                  IC(KC) = -16
              END IF

          ELSE IF (CF(KC).GE.CL(KC)-TEMP) THEN
              IC(KC) = 15
              INEW = KC

          ELSE
              IC(KC) = -15
          END IF

      END IF

      RETURN

      END
* SUBROUTINE PLNEWS             ALL SYSTEMS                   97/12/01
* PORTABILITY : ALL SYSTEMS
* 97/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* TEST ON ACTIVITY OF A GIVEN SIMPLE BOUND.
*
* PARAMETERS :
*  RI  X(NF)  VECTOR OF VARIABLES.
*  IU  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS.
*  RI  XL(NF)  VECTOR CONTAINING LOWER BOUNDS FOR VARIABLES.
*  RI  XU(NF)  VECTOR CONTAINING UPPER BOUNDS FOR VARIABLES.
*  RI  EPS9  TOLERANCE FOR ACTIVE CONSTRAINTS.
*  II  I  INDEX OF TESTED SIMPLE BOUND.
*  IO  INEW  INDEX OF THE NEW ACTIVE CONSTRAINT.
*
      SUBROUTINE PLNEWS(X,IX,XL,XU,EPS9,I,INEW)
C     .. Scalar Arguments ..
      DOUBLE PRECISION EPS9
      INTEGER I,INEW
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION X(*),XL(*),XU(*)
      INTEGER IX(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION TEMP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX
C     ..
      TEMP = 1.0D0
      IF (IX(I).LE.0) THEN

      ELSE IF (IX(I).EQ.1) THEN
          IF (X(I).LE.XL(I)+EPS9*MAX(ABS(XL(I)),TEMP)) THEN
              IX(I) = 11
              INEW = -I
          END IF

      ELSE IF (IX(I).EQ.2) THEN
          IF (X(I).GE.XU(I)-EPS9*MAX(ABS(XU(I)),TEMP)) THEN
              IX(I) = 12
              INEW = -I
          END IF

      ELSE IF (IX(I).EQ.3 .OR. IX(I).EQ.4) THEN
          IF (X(I).LE.XL(I)+EPS9*MAX(ABS(XL(I)),TEMP)) THEN
              IX(I) = 13
              INEW = -I
          END IF

          IF (X(I).GE.XU(I)-EPS9*MAX(ABS(XU(I)),TEMP)) THEN
              IX(I) = 14
              INEW = -I
          END IF

      END IF

      RETURN

      END
* SUBROUTINE PLSETC             ALL SYSTEMS                   97/12/01
* PORTABILITY : ALL SYSTEMS
* 97/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* DETERMINATION OF INITIAL VALUES OF THE CONSTRAINT FUNCTIONS.
*
* PARAMETERS :
*  II  NF  NUMBER OF VARIABLES.
*  II  NC  NUMBER OF CURRENT LINEAR CONSTRAINTS.
*  RI  X(NF)  VECTOR OF VARIABLES.
*  RI  XO(NF)  SAVED VECTOR OF VARIABLES.
*  RU  CF(NF)  VECTOR CONTAINING VALUES OF THE CONSTRAINT FUNCYIONS.
*  II  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
*  RI  CG(NF*MCL)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
*         CONSTRAINTS.
*  RA  S(NF)  AUXILIARY VECTOR.
*
* SUBPROGRAMS USED :
*  S   MXVDIF  DIFFERENCE OF TWO VECTORS.
*  RF  MXVDOT  DOT PRODUCT OF TWO VECTORS.
*  S   MXVMUL  DIAGONAL PREMULTIPLICATION OF A VECTOR.
*
      SUBROUTINE PLSETC(NF,NC,X,XO,CF,IC,CG,S)
C     .. Scalar Arguments ..
      INTEGER NC,NF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CF(*),CG(*),S(*),X(*),XO(*)
      INTEGER IC(*)
C     ..
C     .. Local Scalars ..
      INTEGER JCG,KC
C     ..
C     .. External Functions ..
      DOUBLE PRECISION MXVDOT
      EXTERNAL MXVDOT
C     ..
C     .. External Subroutines ..
      EXTERNAL MXVDIF
C     ..
      CALL MXVDIF(NF,X,XO,S)
      JCG = 0
      DO 10 KC = 1,NC
          IF (IC(KC).NE.0) CF(KC) = CF(KC) + MXVDOT(NF,S,CG(JCG+1))
          JCG = JCG + NF
   10 CONTINUE
      RETURN

      END
* SUBROUTINE PLSETG             ALL SYSTEMS                   97/12/01
* PORTABILITY : ALL SYSTEMS
* 97/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* GRADIENT DETERMINATION IN THE FIRST PHASE OF LP SUBROUTINE.
*
* PARAMETERS :
*  II  NF  DECLARED NUMBER OF VARIABLES.
*  II  NC  NUMBER OF CONSTRAINTS.
*  II  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
*  RI  CG(NF*NC)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
*         CONSTRAINTS.
*  RO  G(NF)  GRADIENT OF THE OBJECTIVE FUNCTION.
*  IO  INEW  INDEX OF THE NEW ACTIVE CONSTRAINT.
*
* SUBPROGRAMS USED :
*  S   MXVDIR  VECTOR AUGMENTED BY THE SCALED VECTOR.
*  S   MXVSET  INITIATION OF A VECTOR.
*
      SUBROUTINE PLSETG(NF,NC,IC,CG,G,INEW)
C     .. Scalar Arguments ..
      INTEGER INEW,NC,NF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CG(*),G(*)
      INTEGER IC(*)
C     ..
C     .. Local Scalars ..
      INTEGER KC
C     ..
C     .. External Subroutines ..
      EXTERNAL MXVDIR,MXVSET
C     ..
      CALL MXVSET(NF,0.0D0,G)
      INEW = 0
      DO 10 KC = 1,NC
          IF (IC(KC).GE.-10) THEN

          ELSE IF (IC(KC).EQ.-11 .OR. IC(KC).EQ.-13 .OR.
     +             IC(KC).EQ.-15) THEN
              CALL MXVDIR(NF,-1.0D0,CG((KC-1)*NF+1),G,G)
              INEW = 1

          ELSE IF (IC(KC).EQ.-12 .OR. IC(KC).EQ.-14 .OR.
     +             IC(KC).EQ.-16) THEN
              CALL MXVDIR(NF,1.0D0,CG((KC-1)*NF+1),G,G)
              INEW = 1
          END IF

   10 CONTINUE
      RETURN

      END
* SUBROUTINE PLSLAG               ALL SYSTEMS                97/12/01
* PORTABILITY : ALL SYSTEMS
* 97/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* NEGATIVE PROJECTED GRADIENT IS DETERMINED USING LAGRANGE MULTIPLIERS.
*
* PARAMETERS :
*  II  NF  DECLARED NUMBER OF VARIABLES.
*  II  N  ACTUAL NUMBER OF VARIABLES.
*  II  NC  NUMBER OF LINEARIZED CONSTRAINTS.
*  II  IAA(NF+1)  VECTOR CONTAINING INDICES OF ACTIVE FUNCTIONS.
*  RI  AG(NF*NC)  MATRIX WHOSE COLUMNS ARE GRADIENTS OF THE LINEAR
*         APPROXIMATED FUNCTIONS.
*  RO  AZ(NF+1)  VECTOR OF LAGRANGE MULTIPLIERS.
*  RI  CG(NF*NC)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
*         CONSTRAINTS.
*  RO  S(NF)  NEGATIVE PROJECTED GRADIENT OF THE QUADRATIC FUNCTION.
*  RI  EPS7  TOLERANCE FOR LINEAR AND QUADRATIC PROGRAMMING.
*  RO  GMAX  NORM OF THE TRANSFORMED GRADIENT.
*
* SUBPROGRAMS USED :
*  S   UXVDIR  VECTOR AUGMENTED BY THE SCALED VECTOR.
*  RF  UXVMAX  L-INFINITY NORM OF A VECTOR.
*
      SUBROUTINE PLSLAG(NF,N,NC,IAA,AG,AZ,CG,S,EPS7,GMAX)
C     .. Scalar Arguments ..
      DOUBLE PRECISION EPS7,GMAX
      INTEGER N,NC,NF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION AG(*),AZ(*),CG(*),S(*)
      INTEGER IAA(*)
C     ..
C     .. Local Scalars ..
      INTEGER J,L,NAA
C     ..
C     .. External Functions ..
      DOUBLE PRECISION MXVMAX
      EXTERNAL MXVMAX
C     ..
C     .. External Subroutines ..
      EXTERNAL MXVDIR
C     ..
      NAA = NF - N
      DO 10 J = 1,NAA
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

   10 CONTINUE
      GMAX = MXVMAX(NF,S)
      IF (GMAX.LE.EPS7) GMAX = 0.0D0
      RETURN

      END
* SUBROUTINE PLTLAG               ALL SYSTEMS                97/12/01
* PORTABILITY : ALL SYSTEMS
* 97/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* MAXIMUM ABSOLUTE VALUE OF THE NEGATIVE LAGRANGE MULTIPLIER IS
* COMPUTED.
*
* PARAMETERS :
*  II  NF  DECLARED NUMBER OF VARIABLES.
*  II  N  ACTUAL NUMBER OF VARIABLES.
*  II  NC  NUMBER OF LINEARIZED CONSTRAINTS.
*  II  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS.
*  II  IA(NA)  VECTOR CONTAINING TYPES OF DEVIATIONS.
*  II  IAA(NF+1)  VECTOR CONTAINING INDICES OF ACTIVE FUNCTIONS.
*  RI  AZ(NF+1)  VECTOR OF LAGRANGE MULTIPLIERS.
*  II  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
*  RI  EPS7  TOLERANCE FOR LINEAR AND QUADRATIC PROGRAMMING.
*  RO  UMAX  MAXIMUM ABSOLUTE VALUE OF THE NEGATIVE LAGRANGE MULTIPLIER.
*  IO  IOLD  INDEX OF THE REMOVED CONSTRAINT.
*
      SUBROUTINE PLTLAG(NF,N,NC,IX,IA,IAA,AZ,IC,EPS7,UMAX,IOLD)
C     .. Scalar Arguments ..
      DOUBLE PRECISION EPS7,UMAX
      INTEGER IOLD,N,NC,NF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION AZ(*)
      INTEGER IA(*),IAA(*),IC(*),IX(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER J,K,L,NAA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      IOLD = 0
      UMAX = 0.0D0
      NAA = NF - N
      DO 10 J = 1,NAA
          TEMP = AZ(J)
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

          ELSE IF ((K.EQ.-1.OR.K.EQ.-3) .AND. UMAX+TEMP.GE.0.0D0) THEN

          ELSE IF ((K.EQ.-2.OR.K.EQ.-4) .AND. UMAX-TEMP.GE.0.0D0) THEN

          ELSE
              IOLD = J
              UMAX = ABS(TEMP)
          END IF

   10 CONTINUE
      IF (UMAX.LE.EPS7) IOLD = 0
      RETURN

      END
* SUBROUTINE PLVLAG               ALL SYSTEMS                97/12/01
* PORTABILITY : ALL SYSTEMS
* 97/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* GRADIENT OF THE OBJECTIVE FUNCTION IS PREMULTIPLIED BY TRANSPOSE
* OF THE MATRIX WHOSE COLUMNS ARE NORMALS OF CURRENT ACTIVE CONSTRAINTS
* AND GRADIENTS OF CURRENT ACTIVE FUNCTIONS.
*
* PARAMETERS :
*  II  NF  DECLARED NUMBER OF VARIABLES.
*  II  N  ACTUAL NUMBER OF VARIABLES.
*  II  NC  NUMBER OF LINEARIZED CONSTRAINTS.
*  II  IAA(NF+1)  VECTOR CONTAINING INDICES OF ACTIVE FUNCTIONS.
*  RI  AG(NF*NA)  VECTOR CONTAINING SCALING PARAMETERS.
*  RI  CG(NF*NC)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
*         CONSTRAINTS.
*  RI  G(NF)  GRADIENT OF THE OBJECTIVE FUNCTION.
*  RO  GN(NF+1)  OUTPUT VECTOR.
*
* SUBPROGRAMS USED :
*  RF  MXVDOT  DOT PRODUCT OF TWO VECTORS.
*
      SUBROUTINE PLVLAG(NF,N,NC,IAA,AG,CG,G,GN)
C     .. Scalar Arguments ..
      INTEGER N,NC,NF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION AG(*),CG(*),G(*),GN(*)
      INTEGER IAA(*)
C     ..
C     .. Local Scalars ..
      INTEGER J,L,NAA
C     ..
C     .. External Functions ..
      DOUBLE PRECISION MXVDOT
      EXTERNAL MXVDOT
C     ..
      NAA = NF - N
      DO 10 J = 1,NAA
          L = IAA(J)
          IF (L.GT.NC) THEN
              L = L - NC
              GN(J) = MXVDOT(NF,AG((L-1)*NF+1),G)

          ELSE IF (L.GT.0) THEN
              GN(J) = MXVDOT(NF,CG((L-1)*NF+1),G)

          ELSE
              L = -L
              GN(J) = G(L)
          END IF

   10 CONTINUE
      RETURN

      END
* SUBROUTINE PLADR0               ALL SYSTEMS                97/12/01
* PORTABILITY : ALL SYSTEMS
* 97/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* TRIANGULAR DECOMPOSITION OF KERNEL OF THE ORTHOGONAL PROJECTION
* IS UPDATED AFTER CONSTRAINT ADDITION.
*
* PARAMETERS :
*  II  NF  DECLARED NUMBER OF VARIABLES.
*  IU  N  ACTUAL NUMBER OF VARIABLES.
*  IU  ICA(NF)  VECTOR CONTAINING INDICES OF ACTIVE CONSTRAINTS.
*  RI  CG(NF*NC)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
*         CONSTRAINTS.
*  RU  CR(NF*(NF+1)/2)  TRIANGULAR DECOMPOSITION OF KERNEL OF THE
*         ORTHOGONAL PROJECTION.
*  RA  S(NF)  AUXILIARY VECTOR.
*  RI  EPS7  TOLERANCE FOR LINEAR INDEPENDENCE OF CONSTRAINTS.
*  RO  GMAX  MAXIMUM ABSOLUTE VALUE OF A PARTIAL DERIVATIVE.
*  RO  UMAX  MAXIMUM ABSOLUTE VALUE OF A NEGATIVE LAGRANGE MULTIPLIER.
*  II  INEW  INDEX OF THE NEW ACTIVE CONSTRAINT.
*  IU  NADD  NUMBER OF CONSTRAINT ADDITIONS.
*  IO  IER  ERROR INDICATOR.
*
* SUBPROGRAMS USED :
*  S   MXSPRB  SPARSE BACK SUBSTITUTION.
*  S   MXVCOP  COPYING OF A VECTOR.
*  RF  MXVDOT  DOT PRODUCT OF TWO VECTORS.
*
      SUBROUTINE PLADR0(NF,N,ICA,CG,CR,S,EPS7,GMAX,UMAX,INEW,NADD,IER)
C     .. Scalar Arguments ..
      DOUBLE PRECISION EPS7,GMAX,UMAX
      INTEGER IER,INEW,N,NADD,NF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CG(*),CR(*),S(*)
      INTEGER ICA(*)
C     ..
C     .. Local Scalars ..
      INTEGER I,J,K,L,NCA,NCR
C     ..
C     .. External Functions ..
      DOUBLE PRECISION MXVDOT
      EXTERNAL MXVDOT
C     ..
C     .. External Subroutines ..
      EXTERNAL MXDPRB,MXVCOP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C     ..
      IER = 0
      IF (N.LE.0) IER = 2
      IF (INEW.EQ.0) IER = 3
      IF (IER.NE.0) RETURN
      NCA = NF - N
      NCR = NCA* (NCA+1)/2
      IF (INEW.GT.0) THEN
          CALL MXVCOP(NF,CG((INEW-1)*NF+1),S)
          GMAX = MXVDOT(NF,CG((INEW-1)*NF+1),S)
          DO 10 J = 1,NCA
              L = ICA(J)
              IF (L.GT.0) THEN
                  CR(NCR+J) = MXVDOT(NF,CG((L-1)*NF+1),S)

              ELSE
                  I = -L
                  CR(NCR+J) = S(I)
              END IF

   10     CONTINUE

      ELSE
          K = -INEW
          GMAX = 1.0D0
          DO 20 J = 1,NCA
              L = ICA(J)
              IF (L.GT.0) THEN
                  CR(NCR+J) = CG((L-1)*NF+K)*GMAX

              ELSE
                  CR(NCR+J) = 0.0D0
              END IF

   20     CONTINUE
      END IF

      IF (NCA.EQ.0) THEN
          UMAX = GMAX

      ELSE
          CALL MXDPRB(NCA,CR,CR(NCR+1),1)
          UMAX = GMAX - MXVDOT(NCA,CR(NCR+1),CR(NCR+1))
      END IF

      IF (UMAX.LE.EPS7*GMAX) THEN
          IER = 1
          RETURN

      ELSE
          N = N - 1
          NCA = NCA + 1
          NCR = NCR + NCA
          ICA(NCA) = INEW
          CR(NCR) = SQRT(UMAX)
          NADD = NADD + 1
      END IF

      RETURN

      END
* SUBROUTINE MXDSMI                ALL SYSTEMS                88/12/01
* PORTABILITY : ALL SYSTEMS
* 88/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* DENSE SYMMETRIC MATRIX A IS SET TO THE UNIT MATRIX WITH THE SAME
* ORDER.
*
* PARAMETERS :
*  II  N  ORDER OF THE MATRIX A.
*  RO  A(N*(N+1)/2)  DENSE SYMMETRIC MATRIX STORED IN THE PACKED FORM
*         WHICH IS SET TO THE UNIT MATRIX (I.E. A:=I).
*
      SUBROUTINE MXDSMI(N,A)
*     .. Scalar Arguments ..
      INTEGER          N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(*)
*     ..
*     .. Local Scalars ..
      INTEGER          I,M
*     ..
      M = N* (N+1)/2
      DO 10 I = 1,M
          A(I) = 0.0D0
   10 CONTINUE
      M = 0
      DO 20 I = 1,N
          M = M + I
          A(M) = 1.0D0
   20 CONTINUE
      RETURN

      END
* SUBROUTINE MXVDIF                ALL SYSTEMS                88/12/01
* PORTABILITY : ALL SYSTEMS
* 88/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* VECTOR DIFFERENCE.
*
* PARAMETERS :
*  RI  X(N)  INPUT VECTOR.
*  RI  Y(N)  INPUT VECTOR.
*  RO  Z(N)  OUTPUT VECTOR WHERE Z:= X - Y.
*
      SUBROUTINE MXVDIF(N,X,Y,Z)
*     .. Scalar Arguments ..
      INTEGER          N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION X(*),Y(*),Z(*)
*     ..
*     .. Local Scalars ..
      INTEGER          I
*     ..
      DO 10 I = 1,N
          Z(I) = X(I) - Y(I)
   10 CONTINUE
      RETURN

      END
* SUBROUTINE MXVIND               ALL SYSTEMS                91/12/01
* PORTABILITY : ALL SYSTEMS
* 91/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* CHANGE OF THE INTEGER VECTOR ELEMENT FOR THE CONSTRAINT ADDITION.
*
* PARAMETERS :
*  IU  IX(N)  INTEGER VECTOR.
*  II  I  INDEX OF THE CHANGED ELEMENT.
*  II JOB  CHANGE SPECIFICATION. IS JOB.EQ.0 THEN IX(I)=10-IX(I).
*
      SUBROUTINE MXVIND(IX,I,JOB)
      INTEGER IX(*),I,JOB
      IF (JOB.EQ.0) IX(I)=10-IX(I)
      RETURN
      END

