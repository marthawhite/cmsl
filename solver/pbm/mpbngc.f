
C******************************************************************************
      SUBROUTINE MPBNGC(N,X,IX,BL,BU,M,MG,MC,IC,CL,CU,CG,F,RL,LMAX,
     &                  GAM,EPS,FEAS,JMAX,NITER,NFASG,NOUT,IPRINT,IERR,
     &                  IWORK,LIWORK,WORK,LWORK,IUSER,USER)
C******************************************************************************
C                                                                             *
C  The Multiobjective Proximity control Bundle method for nonsmooth,          *
C  Nonconvex and Generally Constrained optimization:                          *
C                                                                             *
C             Minimize       f (x),...,f (x),                                 *
C                    n        1         m                                     *
C              x in R                                                         *
C                                                                             *
C             subject to                                                      *
C                                   T                                         *
C                            c  <= c x <= c  ,                                *
C                             l     g      u                                  *
C                                                                             *
C                            b <= x <= b  and                                 *
C                             l         u                                     *
C                                                                             *
C                            f (x) <= 0,   j = m + 1,...,m + mg.              *
C                             j                                               *
C                                                                             *
C  utilizing new version of PLQDF1 by Ladislav Luksan as a quadratic solver.  *
C                                                                             *
C-----------------------------------------------------------------------------*
C  Version  3.1  of  MPBNGC  dated  September 9, 2004.                        *
C  Last modified by Ladislav Luksan 2009.                                     *
C-----------------------------------------------------------------------------*
C  Copyright by:          ********************************                    *
C---------------          *       Marko M. Mäkelä        *                    *
C                         ********************************                    *
C                         *     University of Turku      *                    *
C                         *  Department of Mathematics   *                    *
C                         *        FI-20014 Turku        *                    *
C                         *           Finland            *                    *
C                         ********************************                    *
C-----------------------------------------------------------------------------*
C  Parameters:                                                                *
C-------------                                                                *
C     N      - The number of variables                         (1 <= N).      *
C     X      - The vector of dimension N.                                     *
C                Input:  The feasible initial approximation to the solution.  *
C                Output: The best approximation to the solution.              *
C     IX     - Types if box constraints for individual variables.             *
C                0 : Free variable.                                           *
C                1 : Lower bound    - BL(I) <= X(I).                          *
C                2 : Upper bound    - X(I)  <= BU(I).                         *
C                3 : Two-side bound - BL(I) <= X(I) <= BU(I).                 *
C                5 : Fixed variable.                                          *
C     BL     - Lower bounds of X.                                             *
C     BU     - Upper bounds of X.                                             *
C     M      - The number of objective functions.              (1 <= M).      *
C     MG     - The number of general constraint functions.     (0 <= MG).     *
C     MC     - The number of linear constraints.               (0 <= MC).     *
C     IC     - Types of individual linear constraints.                        *
C                0 : Constraint is not used.                                  *
C                1 : Lower bound    - CL(I) <= CG(.,I)*X.                     *
C                2 : Upper bound    - CG(.,I)*X <= CU(I).                     *
C                3 : Two-side bound - CL(I) <= CG(.,I)*X <= CU(I).            *
C                5 : Equality       - CG(.,I)*X = CL(I).                      *
C     CL     - Lower bounds of CG(.,I)*X.                                     *
C     CU     - Upper bounds of CG(.,I)*X.                                     *
C     CG     - Coefficients of linear constraints: i-th column of CG          *
C              corresponds to i-th linear constraint.                         *
C     F      - The MM (= M + MG) vector of function values at the solution.   *
C     FASG   - SUBROUTINE FASG(N,X,MM,F,G,IERR,IUSER,USER),                   *
C              supplied by the user. Calculates the function values F(X) and  *
C              the transposed generalized Jacobian (N times MM) matrix G(X).  *
C              If there exists problems in function or subgradient            *
C              calculations, set IERR = 8.                                    *
C     RL     - Line search parameter                           (0 < RL < 0.5).*
C     LMAX   - The maximum number of FASG calls in line search (0 < LMAX).    *
C     GAM    - The MP1 vector of distance measure parameters   (0 <=GAM(I)).  *
C              (GAM(I) = 0  if f is convex).                                  *
C              Here                                                           *
C                       MP1 =  M + 1 if MG > 0,                               *
C                       MP1 =  M     if MG = 0.                               *
C     EPS    - The final objective function accuracy parameter (0 < EPS).     *
C     FEAS   - The tolerance for constraint feasibility        (0 < FEAS).    *
C     JMAX   - The maximum number of stored subgradients       (2 <=JMAX).    *
C     NITER  - Input : The maximum number of iterations        (0 < NITER).   *
C              Output: Number of used iterations.                             *
C     NFASG  - Input : The maximum number of FASG calls        (1 < NFASG).   *
C     NOUT   - The output file number.                                        *
C              Output: Number of used function and subgradient calls.         *
C     IPRINT - Printout control parameter.                                    *
C               -1 : No printout.                                             *
C                0 : Only the error messages.                                 *
C                1 : The final values of the objective functions.             *
C                2 : The whole final solution.                                *
C                3 : At each iteration values of the objective functions.     *
C                4 : At each iteration the whole solution.                    *
C     IERR   - Failure parameter.                                             *
C                0 : Everything is OK.                                        *
C                1 : Number of calls of FASG = NFASG.                         *
C                2 : Number of iterations = NITER.                            *
C                3 : Invalid input parameters.                                *
C                4 : Not enough working space.                                *
C                5 : Failure in quadratic program.                            *
C                6 : The starting point is not feasible.                      *
C                7 : Failure in attaining the demanded accuracy.              *
C                8 : Failure in function or subgradient calculations          *
C                    (assigned by the user).                                  *
C     IWORK  - Working array.                                                 *
C     LIWORK - Dimension of IWORK                                             *
C           LIWORK >= 2*(MP1*(JMAX+1)+MC+N).                                  *
C     WORK   - Working array of dimension LWORK.                              *
C     LWORK  - Dimension of WORK                                              *
C           LWORK >= N*(N+2*MP1*JMAX+2*MP1+2*MC+2*MG+2*M+31)/2                *
C                    +MP1*(6*JMAX+10)+JMAX+5*MC+MG+M+18.                      *
C     IUSER  - Working array of the user.                                     *
C     USER   - Working array of the user.                                     *
C                                                                             *
C-----------------------------------------------------------------------------*
C  Machine-Depended Constants:                                                *
C-----------------------------                                                *
C     DRELPR - Relative precision.                                            *
C              The smallest positive number such that  1.0 + DRELPR > 1.0     *
C                                                                             *
      DOUBLE PRECISION  DRELPR
      PARAMETER        (DRELPR = 2.775576D-17)
C                                                                             *
C******************************************************************************
C
      INTEGER           N,M,MG,MC,IX(N),IC(MC),LMAX,JMAX,NITER,NFASG,
     &                  NOUT,IPRINT,IERR,LIWORK,LWORK,IWORK(LIWORK),
     &                  IUSER(*)
      DOUBLE PRECISION
     &                  X(N),BL(N),BU(N),CL(MC),CU(MC),CG(N,MC),
     &                  F(M+MG),RL,GAM(M+1),EPS,FEAS,WORK(LWORK),
     &                  USER(*)
      EXTERNAL          FASG
C
C********************* Local Variables ********************
C
      INTEGER           JKP,MM,MP1,MJKP,MTOT,LSG,LFK,LFK1,LSK,
     &                  LSKJKP,LALFA,LAL1,LXUD,LFXUD,LAMDA,LSUML,
     &                  LH,LCF,LCFD,LAF,LAFD,LS,LG,LAZ,LAR,
     &                  LIA,LIAA
      DOUBLE PRECISION  ZERO,ZER0
      PARAMETER        (ZERO = 1.0D+05*DRELPR,
     &                  ZER0 = 1.0D+02*DRELPR)
           
C
C************************** Start *************************
C
            
      JKP    = JMAX+1
      MM     = M+MG
      MP1    = M
      IF (MG.GT.0) MP1=MP1+1
      MJKP   = MP1*JKP
      MTOT   = MJKP
      LSG    = 1
      LFK    = LSG+N*(MTOT+MM)
      LFK1   = LFK+MJKP
      LSK    = LFK1+MP1
      LSKJKP = LSK+JMAX
      LALFA  = LSKJKP+MP1
      LAL1   = LALFA+MTOT
      LXUD   = LAL1+MP1
      LFXUD  = LXUD+N
      LAMDA  = LFXUD+MM
      LSUML  = LAMDA+MTOT
      LH     = LSUML+MP1
      LCF    = LH+N
      LCFD   = LCF+MC
      LAF    = LCFD+MC
      LAFD   = LAF+MTOT
      LS     = LAFD+MTOT
      LG     = LS+N+1
      LAZ    = LG+N+1
      LAR    = LAZ+N+1
      LIA    = 1
      LIAA   = LIA+MTOT
      IERR   = -1
C
C******************* Parameter Checking *******************
C

            
      CALL CHECK(N,M,MP1,MG,MC,RL,LMAX,GAM,EPS,FEAS,JMAX,NITER,
     &           NFASG,IPRINT,LIWORK,LWORK,IERR)
C
C************************ Iteration ***********************
C
      IF (IERR.LT.0) THEN
         CALL ITERA(N,X,IX,BL,BU,M,MG,F,FASG,MC,IC,CL,CU,CG,RL,LMAX,
     &              GAM,EPS,FEAS,JMAX,JKP,MJKP,MM,MP1,NITER,NFASG,NOUT,
     &              IPRINT,IERR,ZERO,ZER0,WORK(LSG),WORK(LFK),
     &              WORK(LFK1),WORK(LSK),WORK(LSKJKP),WORK(LALFA),
     &              WORK(LAL1),WORK(LXUD),WORK(LFXUD),WORK(LAMDA),
     &              WORK(LSUML),WORK(LH),WORK(LCF),WORK(LCFD),
     &              WORK(LAF),WORK(LAFD),WORK(LS),WORK(LG),WORK(LAZ),
     &              WORK(LAR),IWORK(LIA),IWORK(LIAA),IUSER,USER)
      END IF
      
      RETURN
      END
C
C******************************************************************************

C******************************************************************************
      SUBROUTINE ITERA(N,XK,IX,XL,XU,M,MG,FXK,FASG,NC,IC,CL,CU,CG,RL,
     &                 LMAX,GAM,EPS,FEAS,JMAX,JKP,MJKP,MM,MP1,NITER,
     &                 NFASG,NOUT,IPRINT,IERR,ZERO,ZER0,SG,FK,FK1,SK,
     &                 SKJKP,ALFA,AL1,XUD,FXUD,DLAMDA,SUML,H,CF,CFD,AF,
     &                 AFD,S,G,AZ,AR,IA,IAA,IUSER,USER)
C
C******************************************************************************
C
      INTEGER           N,M,MG,NC,IX(N),IC(NC),LMAX,JMAX,JKP,MJKP,MM,
     &                  MP1,NITER,NFASG,NOUT,IPRINT,IERR,IA(MJKP),
     &                  IAA(N+1),IUSER(*)
      DOUBLE PRECISION  XK(N),XL(N),XU(N),FXK(MM),CL(NC),CU(NC),
     &                  CG(N,NC),RL,GAM(MP1),EPS,FEAS,ZERO,ZER0,
     &                  SG(N,MJKP+MM),FK(MP1,JKP),FK1(MP1),
     &                  SK(JMAX),SKJKP(MP1),ALFA(MJKP),AL1(MP1),
     &                  XUD(N),FXUD(MM),DLAMDA(MJKP),SUML(MP1),
     &                  H(N),CF(NC),CFD(NC),AF(MJKP),AFD(MJKP),S(N+1),
     &                  G(N+1),AZ(N+1),AR((N+1)*(N+2)),USER(*)
      EXTERNAL          FASG
C
C********************* Local Variables ********************
C
      INTEGER           I,J,K,L,JK,JK1,JK2,NFUN,NM,NZ,IUK,IND1,MSAMEE,
     &                  MAXSME,KBF,KBC,ITERL,KC
      DOUBLE PRECISION  U,UK,PKPK,DKDK,DKNORM,TL,VK,WK,EXWK,SK1,UMIN,
     &                  EPSV,ALPKT,ETA9,EPS7,EPS9,GMAX,UMAX,PRO,MXVDOT
     &
      INTRINSIC         DSQRT,DMIN1,DMAX1,DABS
      PARAMETER        (MAXSME = 10,
     &                  ETA9   = 1.0D+60,
     &                  EPS7   = 1.0D-12,
     &                  EPS9   = 1.0D-10)
*
*     INITIAL OPERATIONS WITH SIMPLE BOUNDS
*
      KBF=0
      DO 10 I = 1,N
          IF (IX(I).GT.0) THEN
              KBF=2
              IF ((IX(I).EQ.3.OR.IX(I).EQ.4).AND.XU(I).LE.XL(I)) THEN
                  XU(I) = XL(I)
                  XK(I) = XL(I)
                  IX(I) = 5
              ELSE IF (IX(I).EQ.5 .OR. IX(I).EQ.6) THEN
                  XL(I) = XK(I)
                  XU(I) = XK(I)
                  IX(I) = 5
              END IF
              IF (IX(I).EQ.1.OR.IX(I).EQ.3) XK(I) = DMAX1(XK(I),XL(I))
              IF (IX(I).EQ.2.OR.IX(I).EQ.3) XK(I) = DMIN1(XK(I),XU(I))
          ENDIF
   10 CONTINUE
*
*     INITIAL OPERATIONS WITH GENERAL LINEAR CONSTRAINTS
*
      KBC=0
      IF (NC.GT.0) THEN
          KBC=2
          DO 20 KC = 1,NC
              IF ((IC(KC).EQ.3.OR.IC(KC).EQ.4) .AND.
     +            CU(KC).LE.CL(KC)) THEN
                  CU(KC) = CL(KC)
                  IC(KC) = 5
              ELSE IF (IC(KC).EQ.5 .OR. IC(KC).EQ.6) THEN
                  CU(KC) = CL(KC)
                  IC(KC) = 5
              END IF
              CF(KC) = MXVDOT(N,XK,CG(1,KC))
   20     CONTINUE
*
*     DETERMINATION OF AN INITIAL FEASIBLE POINT
*
          CALL MXVSET(N,0.0D0,G)
          CALL PLLPB2(N,NC,XK,IX,H,XL,XU,CF,CFD,IC,IAA,CL,CU,CG,AR,AZ,
     +                G,G,S,1,KBF,KBC,ETA9,EPS7,EPS9,UMAX,GMAX,NZ,
     +                ITERL)
      ELSE IF (KBF.GT.0) THEN
          DO 30 I = 1,N
              IF (IX(I).GE.5) IX(I) = -IX(I)
              IF (IX(I).LE.0) THEN
              ELSE IF ((IX(I).EQ.1.OR.IX(I).EQ.3) .AND.
     +                 XK(I).LE.XL(I)) THEN
                  XK(I) = XL(I)
              ELSE IF ((IX(I).EQ.2.OR.IX(I).EQ.3) .AND.
     +                 XK(I).GE.XU(I)) THEN
                  XK(I) = XU(I)
              END IF
              CALL PLNEWS(XK,IX,XL,XU,EPS9,I,ITERL)
              IF (IX(I).GT.10) IX(I) = 10 - IX(I)
   30     CONTINUE
      END IF
C
C********************* Initialization *********************
C
      K=1
      JK=1
      JK1=2
      JK2=3
      IND1=1
      MSAMEE=0
      NFUN=1
      CALL FASG(N,XK,MM,FXK,SG(1,IND1),IERR,IUSER,USER)
      IF (IERR.EQ.8) GOTO 300
      CALL CMAX(N,M,MG,FXK,SG(1,IND1))
      UK=0.0D+00
      DO 40 L=1,MP1
         PKPK=PRO(N,SG(1,L),SG(1,L))
         UK=UK+DSQRT(PKPK)
         FK(L,1)=FXK(L)
         FK(L,JKP)=FXK(L)
         SKJKP(L)=0.0D+00
         ALFA(L)=0.0D+00
         ALFA(MP1+L)=0.0D+00
         DO 35 I=1,N
            SG(I,MP1+L)=SG(I,L)
 35      CONTINUE
 40   CONTINUE
      SK(1)=0.0D+00
      IF (MG.GT.0) THEN
         ALFA(MP1)=-FXK(MP1)
         ALFA(2*MP1)=-FXK(MP1)
         IF (FXK(MP1).GT.FEAS) IERR=6
         FXK(MP1)=0.0D+00
         UK=UK-DSQRT(PKPK)
      END IF
      UK=DMAX1(UK/M,ZERO)
      U=UK
      NM=MP1
      EPSV=1.0D+15
      UMIN=1.0D-10*UK
      IUK=0
C
C************************ Iteration ***********************
C
60    CONTINUE
C
C***************** Quadratic Program Call *****************
C
      CALL QUADNW(N,NM,NC,XK,IX,XL,XU,CF,IC,CL,CU,CG,AF,AFD,SG,AZ,AR,
     &            IA,IAA,S,G,H,UK,ALFA,DLAMDA,VK,KBF,KBC,IERR)
C
C**********************************************************
C        Varmistetaan, ettei menna ulos laatikkorajoitteista.
C        3.7.1996
C        Merkkivirhe jälkimmäisestä korjattu.
C        27.11.2002
C**********************************************************
C
      WK=0.0D+00
      IF (IERR.LT.0) THEN
         ALPKT=0.0D+00
         DO 140 L=1,MP1
            IF (K.EQ.1) THEN
               DLAMDA(L+MP1)=1.0D+00
               DLAMDA(L)=0.0D+00
               SUML(L)=1.0D+00
            ELSE
               SUML(L)=0.0D+00
               DO 90 J=1,JK1
                  SUML(L)=SUML(L)+DLAMDA(L+(J-1)*MP1)
90             CONTINUE
               DO 100 J=1,JK1
                  IF (SUML(L).LE.ZERO) THEN
                     DLAMDA(L+(J-1)*MP1)=1.0D+00/JK1
                     SUML(L)=0.0D+00
                  ELSE
                     DLAMDA(L+(J-1)*MP1)=DLAMDA(L+(J-1)*MP1)/SUML(L)
                  END IF
100            CONTINUE
            END IF
            FK(L,JKP)=DLAMDA(L)*FK(L,JKP)
            SKJKP(L)=DLAMDA(L)*SKJKP(L)
            DO 110 I=1,N
               SG(I,L)=DLAMDA(L)*SG(I,L)
110         CONTINUE
            DO 130 J=1,JK
               FK(L,JKP)=FK(L,JKP)+DLAMDA(L+J*MP1)*FK(L,J)
               SKJKP(L)=SKJKP(L)+DLAMDA(L+J*MP1)*SK(J)
               DO 120 I=1,N
                  SG(I,L)=SG(I,L)+DLAMDA(L+J*MP1)*SG(I,L+J*MP1)
120            CONTINUE
130         CONTINUE
            ALFA(L)=DMAX1(DABS(FXK(L)-FK(L,JKP)),
     &                    GAM(L)*SKJKP(L)*SKJKP(L))
            ALPKT=ALPKT+SUML(L)*ALFA(L)
140      CONTINUE
         PKPK=PRO(N,S,S)
         DKDK=PKPK
         DKNORM=DSQRT(DKDK)
C
C******************* Stopping Criterion *******************
C
         EXWK=WK
         WK=0.5D+00*UK*DKDK+ALPKT
         IF (DABS(WK-EXWK).LT.ZERO) THEN
            MSAMEE=MSAMEE+1
         ELSE
            MSAMEE=0
         END IF
         IF ((WK.LE.EPS).OR.(MSAMEE.GE.MAXSME)) THEN
            IF (WK.LE.EPS) THEN
               IERR=0
            ELSE
               IERR=7
            END IF
            K=K-1
            IF (IPRINT.GT.0)
     &         CALL PRINT(N,XK,M,FXK,NFUN,WK,K,IPRINT,NOUT)
         ELSE
            IF (IPRINT.GT.2)
     &         CALL PRINT(N,XK,M,FXK,NFUN,WK,K-1,IPRINT,NOUT)
C
C************* Line Search and Weight Updating ************
C
         CALL LSAWUD(N,M,U,IUK,UMIN,EPSV,XK,FXK,FASG,RL,LMAX,NC,MG,
     &               GAM,NFASG,NFUN,IERR,SG,S,JK1,MJKP,MM,MP1,FK1,
     &               SK1,AL1,ALPKT,DKNORM,TL,VK,XUD,FXUD,ZER0,FEAS,
     &               IUSER,USER)
C
C**********************************************************
C
            IF (K.EQ.NITER) THEN
               IERR=2
               IF ((IPRINT.GT.0).AND.(IERR.LT.8))
     &            CALL PRINT(N,XK,M,FXK,NFUN,WK,K,IPRINT,NOUT)
            END IF
         END IF
         IF (IERR.LT.0) THEN
C
C****************** Subgradient Deletion ******************
C
            IF (JK.EQ.JMAX) THEN
               CALL DEL(N,M,MP1,MG,NC,JK,JKP,FK,SK,SG,ALFA,AFD)
               JK1=JK+1
               JK2=JK+2
            END IF
C
C***************** Linearization Updating *****************
C
            IF (TL.GT.ZER0) THEN
               DO 160 L=1,MP1
                  FK(L,JKP)=FK(L,JKP)+TL*PRO(N,SG(1,L),S)
                  SKJKP(L)=SKJKP(L)+TL*DKNORM
                  ALFA(L)=DMAX1(DABS(FXK(L)-FK(L,JKP)),
     &                          GAM(L)*SKJKP(L)*SKJKP(L))
160            CONTINUE
            END IF
            NM=MP1*JK2
            DO 170 L=1,MP1
               ALFA(L+MP1*JK1)=AL1(L)
               FK(L,JK1)=FK1(L)
170         CONTINUE
            SK(JK1)=SK1
C
C********************** Serious Step **********************
C
            IF (TL.GT.ZER0) THEN
               DO 190 J=1,JK
                  SK(J)=SK(J)+TL*DKNORM
                  DO 180 L=1,MP1
                     IF (K.EQ.1) THEN
                        FK(L,J)=FK(L,J)+TL*AFD(L+(J-1)*MP1)
                     ELSE
                        FK(L,J)=FK(L,J)+TL*AFD(L+J*MP1)
                     END IF
                     ALFA(L+J*MP1)=DMAX1(DABS(FXK(L)-FK(L,J)),
     &                                   GAM(L)*SK(J)*SK(J))
180               CONTINUE
190            CONTINUE
               DO 200 J=1,NC
                  CF(J)=PRO(N,CG(1,J),XK)
200            CONTINUE
            END IF
C
C************************ Updating ************************
C
            K=K+1
            JK=JK+1
            JK1=JK+1
            JK2=JK+2
            UK=U
            GOTO 60
         END IF
      END IF
300   CONTINUE
      NITER=K
      NFASG=NFUN
      RETURN
      END
C
C******************************************************************************

C******************************************************************************
      SUBROUTINE CHECK(N,M,MP1,MG,MC,RL,LMAX,GAM,EPS,FEAS,JMAX,NITER,
     &                 NFASG,IPRINT,LIWORK,LWORK,IERR)
C******************************************************************************
C                                                                             *
C  Procedure, which checks the input parameters.                              *
C                                                                             *
C******************************************************************************
C
      INTEGER           N,M,MP1,MG,MC,LMAX,JMAX,NITER,NFASG,IPRINT,
     &                  LIWORK,LWORK,IERR
      DOUBLE PRECISION  RL,EPS,FEAS,
     &                  GAM(MP1)
C
C********************* Local Variables ********************
C
      INTEGER           LIW,LW,L
      DOUBLE PRECISION  ZERO,HALF
      PARAMETER (ZERO = 0.0D+00,
     &           HALF = 0.5D+00)
C
C************************** Start *************************
C
      LIW = 2*(MP1*(JMAX+1)+N)
      LW  = N*(N+2*MP1*JMAX+2*MP1+2*MG+2*M+31)/2
     &      +MP1*(6*JMAX+10)+JMAX+2*MC+MG+M+18

      IF ((N.LT.1).OR.(M.LT.0).OR.(MG.LT.0).OR.(MC.LT.0).OR.
     &   (RL.LE.ZERO).OR.(RL.GE.HALF).OR.(LMAX.LE.0).OR.
     &   (EPS.LE.ZERO).OR.(FEAS.LE.ZERO).OR.(JMAX.LT.2).OR.
     &   (NITER.LE.0).OR.(NFASG.LE.1).OR.(IPRINT.LT.-1).OR.
     &   (IPRINT.GT.4)) THEN
         IERR=3
      ELSE
         IF ((LIWORK.LT.LIW).OR.(LWORK.LT.LW)) IERR=4
      END IF
      DO 10 L=1,MP1
         IF (GAM(L).LT.ZERO) IERR=3
10    CONTINUE
      RETURN
      END
C
C******************************************************************************

C******************************************************************************
      SUBROUTINE DEL(N,M,MP1,MG,MC,JK,JKP,FK,SK,SG,ALFA,AFD)
C******************************************************************************
C                                                                             *
C  Procedure, which deletes the smallest indexis from the vectors FK, SK,     *
C  SG and ALFA.                                                               *
C                                                                             *
C******************************************************************************
C
      INTEGER           N,M,MP1,MG,MC,JK,JKP
      DOUBLE PRECISION  FK(MP1,JKP),SK(JK),SG(N,MP1*JKP+M+MG),
     &                  ALFA(MP1*JKP),AFD(MP1*JKP)
C
C********************* Local Variables ********************
C
      INTEGER           I,J,L
C
C*********************** Iteration ************************
C
      DO 30 J=1,JK-1
         SK(J)=SK(J+1)
         DO 20 L=1,MP1
            FK(L,J)=FK(L,J+1)
            ALFA(L+J*MP1)=ALFA(L+(J+1)*MP1)
            AFD(L+J*MP1)=AFD(L+(J+1)*MP1)
            DO 10 I=1,N
               SG(I,L+J*MP1)=SG(I,L+(J+1)*MP1)
10          CONTINUE
20       CONTINUE
30    CONTINUE
      DO 50 L=1,MP1
         DO 40 I=1,N
            SG(I,L+MP1*JK)=SG(I,L+MP1*(JK+1))
40       CONTINUE
50    CONTINUE
      JK=JK-1
      RETURN
      END
C
C******************************************************************************

C******************************************************************************
      SUBROUTINE H(M,MM,HXK,LM,FXK,FXUD,NFEA,FEAS)
C******************************************************************************
C                                                                             *
C  Procedure, which calculates the current value of H(y;x).                   *
C                                                                             *
C******************************************************************************
C
      INTEGER           M,MM,LM,NFEA
      DOUBLE PRECISION  HXK,FXK(MM),FXUD(MM),FEAS
C
C********************* Local Variables ********************
C
      INTEGER           L
C
C*********************** Iteration ************************
C
      NFEA=0
      HXK=FXUD(1)-FXK(1)
      LM=1
      DO 10 L=1,M
         IF (HXK.LT.(FXUD(L)-FXK(L))) THEN
            HXK=FXUD(L)-FXK(L)
            LM=L
         END IF
10    CONTINUE
      IF ((MM-M).GT.0) THEN
         IF (FXUD(M+1).GT.FEAS) NFEA=1
      END IF
      RETURN
      END
C
C*****************************************************************************

C******************************************************************************
      SUBROUTINE CMAX(N,M,MG,FXK,SG)
C******************************************************************************
C                                                                             *
C  Procedure, which calculates the maximum of the constraint functions.       *
C                                                                             *
C******************************************************************************
C
      INTEGER           N,M,MG
      DOUBLE PRECISION  FXK(M+MG),SG(N,M+MG)
C
C********************* Local Variables ********************
C
      INTEGER           I,L,IMAX
      DOUBLE PRECISION  FMAX
C
C*********************** Iteration ************************
C
      IF (MG.GT.1) THEN
         FMAX=FXK(M+1)
         IMAX=1
         DO 10 L=2,MG
            IF (FMAX.LT.FXK(M+L)) THEN
               FMAX=FXK(M+L)
               IMAX=L
            END IF
10       CONTINUE
         IF (IMAX.GT.1) THEN
            FXK(M+1)=FMAX
            DO 20 I=1,N
               SG(I,M+1)=SG(I,M+IMAX)
20          CONTINUE
         END IF
      END IF
      RETURN
      END
C
C******************************************************************************

C******************************************************************************
      SUBROUTINE LSAWUD(N,M,UK,IUK,UMIN,EPSV,XK,FXK,FASG,RL,LMAX,MC,MG,
     &                  GAM,NFASG,NFUN,IERR,SG,S,JK1,MJKP,MM,MP1,FK1,
     &                  SK1,AL1,ALPKT,DKNORM,TL,VK,XUD,FXUD,ZER0,FEAS,
     &                  IUSER,USER)
C******************************************************************************
C                                                                             *
C  Line Search Algorithm and Weight Updating Algorithm, which minimizes       *
C  the function                                                               *
C                t  ----> f ( x + t d )                                       *
C                              k     k                                        *
C                                                                             *
C  by using the quadratic interpolation  g : R --> R  such that               *
C                                                                             *
C      g ( 0 ) = f ( x )  ,  g ( 1 ) = f ( x + d )  and  g' ( 0 ) = v         *
C                     k                     k   k                    k        *
C                                                                             *
C  and chooses the weight UK by using the safequarded quadratic interpolation.*
C                                                                             *
C******************************************************************************
C
      INTEGER           N,M,IUK,LMAX,MC,MG,NFASG,NFUN,IERR,JK1,
     &                  MJKP,MM,MP1,IUSER(*)
      DOUBLE PRECISION  UK,UMIN,EPSV,RL,SK1,ALPKT,DKNORM,TL,VK,
     &                  ZER0,FEAS,
     &                  XK(N),FXK(MM),GAM(MP1),SG(N,MJKP+MM),S(N),
     &                  FK1(MP1),AL1(MP1),XUD(N),FXUD(MM),
     &                  USER(*)
C
C********************* Local Variables ********************
C
      INTEGER           I,L,LM,LFUN,NFEA
      DOUBLE PRECISION  U,UINT,PRO,T,TU,HXK,HXTU,RR,TT,ZET
      INTRINSIC         DMAX1,DMIN1,DABS,MAX0,MIN0
      PARAMETER  (RR  = 0.5D+00,
     &            TT  = 0.1D+00)
      EXTERNAL FASG
C
C********************* Initialization *********************
C
      ZET=1.0D+00-1.0D+00/(2.0D+00*(1.0D+00-RL))
      U=UK
      TL=0.0D+00
      TU=1.0D+00
      T=1.0D+00
      DO 10 I=1,N
         XUD(I)=XK(I)+S(I)
10    CONTINUE
      NFUN=NFUN+1
      CALL FASG(N,XUD,MM,FXUD,SG(1,MP1*JK1+1),IERR,IUSER,USER)
      IF (IERR.EQ.8) GOTO 80
      CALL CMAX(N,M,MG,FXUD,SG(1,MP1*JK1+1))
      CALL H(M,MM,HXK,LM,FXK,FXUD,NFEA,FEAS)
      LFUN=1
      IF (NFUN.EQ.NFASG) IERR=1
      HXTU=HXK
      IF (DABS(VK).GE.ZER0) THEN
         UINT=2.0D+00*UK*(1.0D+00-HXK/VK)
      ELSE
         UINT=UK
      END IF
C
C*********************** Iteration ************************
C
20    CONTINUE
      IF ((HXK.LE.(RL*T*VK)).AND.(NFEA.EQ.0)) THEN
         TL=T
         DO 30 I=1,N
            XK(I)=XUD(I)
30       CONTINUE
         DO 40 L=1,M
            FXK(L)=FXUD(L)
40       CONTINUE
      ELSE
         TU=T
         HXTU=HXK
      END IF
      IF (IERR.LT.0) THEN
C
C******************** Long Serious Step *******************
C
         IF (TL.GT.TT) THEN
            DO 50 L=1,MP1
               FK1(L)=FXUD(L)
               AL1(L)=FXK(L)-FXUD(L)
50          CONTINUE
            SK1=0.0D+00
            IF (TL.GT.(1.0D+00-ZER0)) THEN
               IF ((HXK.LE.(RR*VK)).AND.(IUK.GT.0)) THEN
                  U=UINT
               ELSE
                  IF (IUK.GT.3) U=0.5D+00*UK
               END IF
               U=DMAX1(U,DMAX1(0.1D+00*UK,UMIN))
               EPSV=DMAX1(EPSV,-2.0D+00*VK)
               IUK=MAX0(IUK+1,1)
               IF (UK.GT.U) THEN
                  UK=U
                  IUK=1
               END IF
            END IF
         ELSE
C
C*************** Short Serious or Null Step ***************
C
            SK1=(T-TL)*DKNORM
            DO 60 L=1,MP1
               FK1(L)=FXUD(L)-(T-TL)*PRO(N,SG(1,MP1*JK1+L),S)
               AL1(L)=DMAX1(DABS(FXK(L)-FK1(L)),GAM(L)*SK1*SK1)
60          CONTINUE
            IF ((-AL1(LM)+PRO(N,SG(1,MP1*JK1+LM),S)*UK).LT.(RR*VK)
     &         .AND.(LFUN.LT.LMAX)) THEN
C
C************************ Updating ************************
C
               IF (TL.LT.ZER0) THEN
                  T=ZET*TU
                  IF ((TU*VK-HXTU).GT.ZER0)
     &               T=DMAX1(T,(0.5D+00*VK*TU*TU)/(TU*VK-HXTU))
                  T=DMIN1(T,(1.0D+00-ZET)*TU)
               ELSE
                  T=0.5D+00*(TL+TU)
               END IF
               DO 70 I=1,N
                  XUD(I)=XK(I)+T*S(I)
70             CONTINUE
               NFUN=NFUN+1
               LFUN=LFUN+1
               CALL FASG(N,XUD,MM,FXUD,SG(1,MP1*JK1+1),IERR,
     &                   IUSER,USER)
               IF (IERR.EQ.8) GOTO 80
               CALL CMAX(N,M,MG,FXUD,SG(1,MP1*JK1+1))
               CALL H(M,MM,HXK,LM,FXK,FXUD,NFEA,FEAS)
               IF (NFUN.EQ.NFASG) IERR=1
               GOTO 20
            END IF
            IF (TL.LT.ZER0) THEN
               EPSV=DMIN1(EPSV,UK*DKNORM+ALPKT)
               IF ((AL1(LM).GT.DMAX1(EPSV,-10.0D+00*VK)).AND.
     &            (IUK.LT.-3)) U=UINT
               U=DMIN1(U,10.0D+00*UK)
               IUK=MIN0(IUK-1,-1)
               IF (UK.LT.U) THEN
                  UK=U
                  IUK=-1
               END IF
            END IF
         END IF
      END IF
80    CONTINUE
      RETURN
      END
C
C******************************************************************************

C******************************************************************************
      SUBROUTINE PRINT(N,X,M,FX,NFUN,WK,IT,IPRINT,NOUT)
C******************************************************************************
C
C******************************************************************************
C
      INTEGER           N,M,NFUN,IT,IPRINT,NOUT
      DOUBLE PRECISION  WK,
     &                  X(N),FX(M)
C
C********************* Local Variables ********************
C
      INTEGER           I,ITEMP,L
      CHARACTER         FORM1*70 ,FORM2*33, FORM3*33
      INTRINSIC         MOD
C
C**********************************************************
C
      DATA FORM1 /'(1X,''Iter:'',I4,''  Nfun:'',I4,''  f'',IZ,''(x) ='',
     &G14.7,''     Eps ='',G14.7)'/
      DATA FORM2 /'(16X,''       f'',IZ,''(x) ='',G14.7)'/
      DATA FORM3 /'(16X,''        x('',IZ,'') ='',G14.7)'/
C
      IF (IPRINT.GT.0) THEN
         IF (IT.GE.10000) FORM1(14:14)='5'
         ITEMP=INT(LOG10(FLOAT(M)))+1
         WRITE (FORM1(36:36),'(I1)') ITEMP
         WRITE (FORM2(18:18),'(I1)') ITEMP
         WRITE (NOUT,FORM1) IT,NFUN,1,FX(1),WK
         DO 10 L=2,M
            WRITE (NOUT,FORM2) L,FX(L)
10       CONTINUE
         IF (MOD(IPRINT,2).EQ.0) THEN
            WRITE (NOUT,999)
            ITEMP=INT(LOG10(FLOAT(N)))+1
            WRITE (FORM3(20:20),'(I1)') ITEMP
            DO 20 I=1,N
               WRITE (NOUT,FORM3) I,X(I)
20          CONTINUE
            WRITE (NOUT,999)
         END IF
      END IF
      RETURN
999   FORMAT(' ')
      END
C
C******************************************************************************

C******************************************************************************
      DOUBLE PRECISION FUNCTION PRO(N,X,Y)
C******************************************************************************
C                                                                             *
C  Algorithm, which calculates the inner product of two N-dimensional         *
C  vectors X and Y, e. g.                                                     *
C                                                                             *
C                        N                                                    *
C                 PRO = SUM X(I)*Y(I)                                         *
C                       I=1                                                   *
C                                                                             *
C******************************************************************************
C
      INTEGER           N
      DOUBLE PRECISION  X(N),Y(N)
C
C********************* Local Variables ********************
C
      INTEGER           I
      DOUBLE PRECISION  SUM
C
C*********************** Iteration ************************
C
      SUM=0.0D+00
      DO 10 I=1,N
         SUM=SUM+X(I)*Y(I)
10    CONTINUE
      PRO=SUM
      RETURN
      END
C
C******************************************************************************

C******************************************************************************
      SUBROUTINE QUADNW(N,NA,NC,X,IX,XL,XU,CF,IC,CL,CU,CG,AF,AFD,AG,AZ,
     &                  AR,IA,IAA,S,G,H,U,ALFA,DLAMDA,V,KBF,KBC,IERR)
C******************************************************************************
C                                                                             *
C  Procedure, which calls PLQDF1 by Ladislav Luksan.                          *
C                                                                             *
C******************************************************************************
C
      INTEGER           N,NA,NC,IX(N),IC(NC),IA(NA),IAA(N+1),KBF,KBC,
     &                  IERR
      DOUBLE PRECISION  X(N),XL(N),XU(N),CF(NC),CL(NC),CU(NC),CG(NC*N),
     &                  AF(NA),AFD(NA),AG(NA*N),AZ(N+1),
     &                  AR((N+1)*(N+2)/2),S(N+1),G(N+1),H(N),U,ALFA(NA),
     &                  DLAMDA(NA),V
C
C********************* Local Variables ********************
C
      INTEGER           NZ,ITERQ,I,K
      DOUBLE PRECISION  ETA0,ETA2,ETA9,EPS7,EPS9,UMAX,GMAX,MXVDOT
      PARAMETER (ETA0  = 1.0D-15,
     &           ETA2  = 1.0D-12,
     &           ETA9  = 1.0D+60,
     &           EPS7  = 1.0D-14,
     &           EPS9  = 1.0D-12)
      DO 2 I=1,NA
      IA(I)=2
      AF(I)=-ALFA(I)
    2 CONTINUE
      DO 3 I=1,N
      H(I)=U
    3 CONTINUE
C
C***************** Quadratic Program Call *****************
C
      CALL PLQDF1(N,NA,NC,X,IX,XL,XU,AF,AFD,IA,IAA,AG,AR,AZ,CF,IC,CL,
     &           CU,CG,G,H,S,2,KBF,KBC,10,ETA0,ETA2,ETA9,EPS7,
     &           EPS9,V,UMAX,GMAX,NZ,ITERQ)
C
C**********************************************************
C
      IF (ITERQ.LE.0) THEN
          IERR=5
          RETURN
      ENDIF
      K=1
      DO 4 I=1,NA
      DLAMDA(I)=0.0D+00
      IF (IA(I).LT.0) AFD(I)=MXVDOT(N,AG(K),S)
      K=K+N
    4 CONTINUE
      DO 5 I=1,N-NZ
      K=IAA(I)
      IF(K.GT.NC) DLAMDA(K-NC)=-AZ(I)
    5 CONTINUE
      RETURN
      END
