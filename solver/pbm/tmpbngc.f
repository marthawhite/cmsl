C  Test main program for proximal bundle method containing the                *
C  Rosenbrock-function                                                        *
C                                                                             *
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
C*******************************************************************************
C
      PROGRAM TMPBNGC
C
C*******************************************************************************
C
      INTEGER           N,JMAX,M,MC,MG,MP1,NOUT,IPRINT,LIWORK,LWORK
      PARAMETER        (N      = 2,
     &                  JMAX   = 10,
     &                  M      = 1,
     &                  MC     = 0,
     &                  MG     = 0,
c     &                  MG     = 1,
     &                  MP1    = M,
c     &                  MP1    = M+1,
     &                  NOUT   = 6,
     &                  IPRINT = 3,
     &                  LIWORK = 2*(MP1*(JMAX+1)+MC+N),
     &                  LWORK  = N*(N+2*MP1*JMAX+2*MP1+2*MC+2*MG+2*M+31)/2
     &                           +MP1*(6*JMAX+10)+JMAX+5*MC+MG+M+18)
      INTEGER           NITER,NFASG,IERR,I,LMAX,
     &                  IX(N),IC(MC+1),IWORK(LIWORK),IUSER(1)
      DOUBLE PRECISION  RL,EPS,FEAS,
     &                  X(N),F(M+MG),GAM(MP1),BL(N),BU(N),
     &                  CG(N,MC+1),CL(MC+1),CU(MC+1),WORK(LWORK),USER(1)
      EXTERNAL          FASG
C
      NITER =  1000
      NFASG =  1000

      IX(1) =  0
      IX(2) =  0

      IC(1) =  0

      X(1)  = -1.2D+00
      X(2)  =  1.0D+00

      RL=0.01D+00
      LMAX=100
      GAM(1)=0.5D+00

      EPS=1.0D-05
      FEAS=1.0D-09
C
      CALL MPBNGC(N,X,IX,BL,BU,M,MG,MC,IC,CL,CU,CG,F,FASG,RL,LMAX,
     &            GAM,EPS,FEAS,JMAX,NITER,NFASG,NOUT,IPRINT,IERR,
     &            IWORK,LIWORK,WORK,LWORK,IUSER,USER)
C
      PRINT*
      PRINT*,'IERR = ',IERR
      PRINT*
      DO 70 I=1,N
          PRINT*,'X(',I,')=',X(I)
70    CONTINUE
      PRINT*
      PRINT*,'F(X)  = ',F
      PRINT*,'NITER = ',NITER
      PRINT*,'NFASG = ',NFASG
      PRINT*

      STOP
      END

C
C******************************************************************************
C
      SUBROUTINE FASG(N,X,M,F,G,IERR,IUSER,USER)
C     FASG   - SUBROUTINE FASG(N,X,MM,F,G,IERR,IUSER,USER),                   *
C              supplied by the user. Calculates the function values F(X) and  *
C              the transposed generalized Jacobian (N times MM) matrix G(X).  *
C              If there exists problems in function or subgradient            *
C              calculations, set IERR = 8.                                    *
C                                                                             * 
C******************************************************************************
C
      INTEGER           N,M,IUSER(1),IERR
      DOUBLE PRECISION  F(M),X(N),G(N,M),USER(1)

C     Add your functions here.

      F(1)=100*(X(2)-X(1)**2)**2+(1-X(1))**2
      G(1,1)=-400*X(1)*(X(2)-X(1)**2)-2*(1-X(1))
      G(2,1)=200*(X(2)-X(1)**2)

      RETURN
      END

