MODULE solver_eos
!
!
USE def_constants, ONLY: pr
USE properties
!USE peng_robinson
!USE interp_table
!USE var_const, ONLY: guessP, guessp_BC, i_bc
!
IMPLICIT NONE
PRIVATE F_Brent, BrentRoots2

CONTAINS


!===============================================================================
          SUBROUTINE eos_1d(MODE, out_1, out_2, resnorm, Niter,&
     &                             exitflag, GGG, X0, in_1, out3)
!===============================================================================
!         1D        inverse the internal energy to find the temperature
!------------------------------------------------------------------------------
!
!  Input :
!  -------
!     MODE        =  2 Peng-robinson, 4 span-wagner, 5 (P,T)--> v and e
!     GGG         =  e
!     in_1        =  v
!     X0          =  initial guess of T
!
!
!  Output :
!  -------
!
!     out_1       =  final T
!     out_2       =  idem
!     out3       =  idem
!
!-------------------------------------------------------------------
!                                         M. De Lorenzo      03/2016
!-------------------------------------------------------------------
!
!
      IMPLICIT NONE
!
      INTEGER             ::  MAX_ITER
      INTEGER, INTENT(IN) ::  MODE
      INTEGER, INTENT(OUT) :: exitflag, Niter
      REAL(8), INTENT(IN)   :: GGG, X0, in_1
      REAL(8), INTENT(OUT)  :: out_1, out_2, resnorm, out3
      REAL(8)  ::  TOL_X, TOL_FUN, ALPHA, MIN_LAMBDA, MAX_LAMBDA, Jstar,&
     &  lambda, slope, fffold, lambda_min, fff, fff2, A_l, aa, bb, Dx_J,&
     &  discriminant, lambda2, lambda_OLD, XXX, F, delta, XXX_forw, J0,&
     &  XXX_back, F_forw, F_back, dF, TYPX, dx_star, dx, XXX_old, g, J
      REAL(8), DIMENSION (2,2)  ::  B_l
      REAL(8), DIMENSION (2,1)  ::  C_l, aabb
!
!
      XXX = X0                        ! Initial value
!
!   Options
!
      TOL_X      = 1d-20            ! relative max step tolerance
      TOL_FUN    = 1d-12            ! function tolerance
      MAX_ITER   = 500              ! max number of iterations
      ALPHA      = 1d-4             ! criteria for decrease
      MIN_LAMBDA = 1d-1             ! min lambda
      MAX_LAMBDA = 5d-1
      TYPX       = abs(XXX)           ! x scaling value, remove zeros
!
! ---- Function F definition
!      F       = GGG_input - GGG
!
      call F_New_Rap1D(MODE, F, out_2, XXX, GGG, in_1, out3)
!
!
! --------- Jacobian estimation-----------------
!
      Dx_J       = 1d-6         ! finite difference delta
      delta      = Dx_J
      XXX_forw   = XXX + delta
      XXX_back   = XXX - delta
      call F_New_Rap1D(MODE, F_forw, out_2, XXX_forw, GGG, in_1, out3)
      call F_New_Rap1D(MODE, F_back, out_2, XXX_back, GGG, in_1, out3)
      dF       = F_forw - F_back    ! delta F
      J        = dF/Dx_J/2d0        ! derivatives dF/d_n
      J0       = 1d0 / TYPX         ! Jacobian scaling factor
      Jstar    = J/J0;              ! scale Jacobian
      resnorm  = abs(F)             ! Norm of residues
      exitflag = 1
!
!
!--- SOLVER -------------------
!
      Niter  = 0
      lambda = 1d0                ! backtracking
      lambda_OLD = lambda
!
!
!
1       DO WHILE ( ((resnorm>TOL_FUN) .OR. (lambda<1d0)) .AND. &
     &             (exitflag>=0) .AND. (Niter<=MAX_ITER))
             Niter = Niter + 1
!
!--------- Newton-Raphson solver
!
             IF (lambda .eq. 1.0d0) THEN
               dx_star = -F/Jstar        ! calculate Newton step
               dx      = dx_star*TYPX
               g       = F*Jstar
               slope   = g*dx_star
               fffold  = F*F
               XXX_old = XXX
               lambda_min = TOL_X/(ABS(dx)/ABS(XXX_old))
            ENDIF
!
!--------- Check about proximity of XXX and XXX_OLD
!
             IF (lambda < lambda_min) THEN
                exitflag = 2;
                !print*,'XXX is too close to XXX_OLD'
                !EXIT ! OUT NewRap
             ENDIF
!
!--------- Eventually, backtracking of New-Rap step
!
             XXX = XXX_old + dx*lambda
             call F_New_Rap1D(MODE, F, out_2, XXX, GGG, in_1, out3)
!
!--------- Jacobian estimation
!
             XXX_forw   = XXX + delta
             XXX_back   = XXX - delta
             call F_New_Rap1D(MODE, F_forw, out_2, XXX_forw, GGG, in_1,&
     & out3)
             call F_New_Rap1D(MODE, F_back, out_2, XXX_back, GGG, in_1,&
     & out3)
             dF     = F_forw - F_back    ! delta F
             J      = dF/Dx_J/2d0        ! derivatives dF/d_n
             Jstar  = J/J0               ! scale Jacobian
             fff    = F*F                ! f to be minimized
!
!--------- Optimization technique (minimizing fff = SUM(F*F) )
!
             IF (fff > fffold + ALPHA*lambda*slope) THEN
                 IF (lambda .eq. 1d0) THEN
                     lambda = -slope/2d0/(fff-fffold-slope) ! calculate lambda
                 ELSE
                     A_l      =  1d0/(lambda_OLD - lambda2)
                     B_l(1,1) =  1d0/lambda_OLD/lambda_OLD
                     B_l(1,2) = -1d0/lambda2/lambda2
                     B_l(2,1) = -lambda2/lambda_OLD/lambda_OLD
                     B_l(2,2) =  lambda_OLD/lambda2/lambda2
                     C_l(1,1) =  fff-fffold-lambda_OLD*slope
                     C_l(2,1) =  fff2-fffold-lambda2*slope
                     !
                     aabb = A_l* MATMUL(B_l,C_l)
                     aa   = aabb(1,1)
                     bb   = aabb(2,1)
                     IF (aa .EQ. 0.0d0) THEN
                         lambda = -slope/2d0/bb;
                     ELSE
                         discriminant = bb**2d0 - 3d0*aa*slope;
                         IF (discriminant < 0d0) THEN
                             lambda = MAX_LAMBDA*lambda_OLD;
                         ELSEIF (bb <= 0d0) THEN
                             lambda = (-bb + sqrt(discriminant))/3d0/aa
                         ELSE
                             lambda = -slope/(bb + sqrt(discriminant))
                         ENDIF
                     ENDIF
                 lambda = min(lambda,MAX_LAMBDA*lambda_OLD)     ! minimum step length
                 ENDIF
             ELSE
                 lambda = 1d0; ! fraction of Newton step
             ENDIF
             IF (lambda < 1d0) THEN
                 lambda2 = lambda_OLD;
                 fff2    = fff;       ! save 2nd most previous value
                 lambda  = max(lambda,MIN_LAMBDA*lambda_OLD) ! minimum step length
                 GO TO 1
             ENDIF
!
             resnorm = abs(F)
!
      END DO
!
      out_1 = XXX
!

      END SUBROUTINE eos_1d
!
!
!===============================================================================
      SUBROUTINE F_New_Rap1D(MODE, F, out_2, XXX, GGG, in_1, out_3)
!===============================================================================
          IMPLICIT NONE
          INTEGER, INTENT(IN)  :: MODE
          REAL(8), INTENT(IN)  :: in_1, XXX, GGG
          REAL(8), INTENT(OUT) :: F, out_2,out_3
!
!LOCAL
          INTEGER :: i,j
          REAL(8) :: temp, e, v,sound,a_out,res,p_guess
          REAL(8) :: delta_p,press,e_l,e_v,V_v,V_l,qual
          
         

! (v,e) --> T: v=in_1, e=GGG, T=XXX
!          v = in_1
!          temp = XXX
         
          IF (MODE .EQ. 4) THEN
                v = in_1
                temp = XXX
!                       
                CALL inter_energy(temp, v, e)
                F = (e - GGG) * 1d-3
                out_2 = 0.0d0
                out_3 = 0.0d0
!               
          ELSEIF (MODE .EQ. 2) THEN
!               
                v = in_1
                temp = XXX              
 
!                CALL interenergy_pr(temp, v, e)
                F = (e - GGG) * 1d-3
                out_2 = 0.0d0
                out_3 = 0.0d0
!
          ELSEIF (MODE .EQ. 5) THEN
! (T,p) --> rho,e, temp=in_1,
                temp = in_1     
                v = XXX
                CALL pressure(temp,v,press)
                F = (press - GGG) * 10d-6 
               
                CALL inter_energy(temp,v,e)  
                out_2 = e
          ELSEIF (MODE .EQ. 3) THEN
! (v,p) --> e  for outletBC using LOOKUP-table
                v = in_1
                e = XXX
!                       
!                CALL CO2BLLT_EQUI(press,temp,sound,&
!      &              qual, a_out, res, e,v,guessp_BC(i_bc))
                F = (press - GGG) * 1d-6
!                guessp_BC(i_bc) = press
!
                out_2 = 0.0d0
                out_3 = 0.0d0                
!
          ELSEIF (MODE .EQ. 6) THEN
!(e,p)--> v  for outlet BC
                e = in_1
                v = XXX
!                       
!                guessp_BC(i_bc) = (guessp_BC(i_bc)*2.0+guessP(100,i_bc))/3.0
!                print*,'guessp_BC', guessp_BC(i_bc)
!                CALL CO2BLLT_EQUI(press,temp,sound,&
!      &              qual, a_out, res, e,v,guessp_BC(i_bc))
                F = (press - GGG) * 1d-6
!                print*, 'press_outlet',press
!                guessp_BC(i_bc) = press
!
                out_2 = 0.0d0
                out_3 = 0.0d0
          ELSE

                print*,'MODE of EOS unknown'
                STOP
          ENDIF          
!
      END SUBROUTINE F_New_Rap1D
!!
!!
!=======================================================================================
!
SUBROUTINE BrentRoots2(MODE, out_1, out_2, residue, Niter, GGG, lower, upper, in_1, in_2)
!
!=======================================================================================
!
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: MODE
      INTEGER, INTENT(OUT) :: Niter
      REAL(8), INTENT(IN)  :: lower, upper, GGG, in_1, in_2
      REAL(8), INTENT(OUT) :: out_1, out_2, residue

      INTEGER              :: done
      REAL(8), PARAMETER   :: FPP      = 1d-16 ! floating-point precision
      REAL(8), PARAMETER   :: MAX_ITER = 2d2   ! maximum allowed number of iterations
      REAL(8), PARAMETER   :: nearzero = 1d-16
      REAL(8) :: TOL, AA, BB, CC, DD, EE, FA, FB, FC, xm, PP, SS, QQ, RR,&
&                Tol1, Root, Minimum

      TOL    = 1d-12
      done   = 0
      Niter  = 0
!
      AA = lower
      BB = upper
      CALL F_Brent(MODE, FA, out_2, AA, GGG, in_1, in_2)
      IF (abs(FA) .LT. TOL) THEN
       out_1 = AA
       RETURN
      ENDIF
!
      CALL F_Brent(MODE, FB, out_2, BB, GGG, in_1, in_2)
      IF (abs(FB) .LT. TOL) THEN
       out_1 = BB
       RETURN
      ENDIF
!
! Check if root is bracketed
      IF (((FA .GT. 0d0) .AND. (FB .GT. 0d0)) .OR. ((FA .LT. 0d0) .AND.&
     & (FB .LT. 0d0))) THEN
       print*, 'ERROR: Root is not bracketed in MODE', MODE
       STOP
      ELSE
       FC = FB
       DO WHILE ((done .EQ. 0) .AND. (Niter .LT. MAX_ITER))
                Niter =  Niter + 1
!
                 IF (((FC .GT. 0d0) .AND. (FB .GT. 0d0)) .OR.&
     &                        ((FC .LT. 0d0) .AND. (FB .LT. 0d0))) THEN
                 ! rename AA, BB, CC and adjust bounding interval DD
                    CC = AA
                    FC = FA
                    DD = BB - AA
                    EE = DD
                 ENDIF
                 IF (abs(FC) .LT. abs(FB)) THEN
                    AA = BB
                    BB = CC
                    CC = AA
                    FA = FB
                    FB = FC
                    FC = FA
                 ENDIF
                 Tol1 = 2d0 * FPP * abs(BB) + 0.5d0 * TOL  ! convergence check
                 xm   = 0.5d0 * (CC - BB)
                 IF((abs(xm) .LE.Tol1) .OR. (abs(FA) .LT.nearzero)) THEN
                 ! a root has been found
                     Root        = BB
                     done        = 1
                     CALL F_Brent(MODE, residue, out_2,BB,GGG,in_1,in_2)
                 ELSE
                     IF ((abs(EE) .GE. Tol1) .AND.&
     &                                 (abs(FA) .GT. abs(FB))) THEN
                         SS = FB / FA  ! attempt inverse quadratic function
                         IF (abs(AA - CC) .LT. nearzero) THEN
                             PP = 2d0 * xm * SS
                             QQ = 1d0 - SS
                         ELSE
                             QQ = FA / FC
                             RR = FB / FC
                             PP = SS * (2d0 * xm * QQ * (QQ - RR) -&
     &                                           (BB- AA) * (RR - 1d0))
                             QQ = (QQ - 1d0) * (RR - 1d0) * (SS - 1d0)
                         ENDIF
                         IF (PP .GT. nearzero) THEN
                             QQ = -QQ
                         ENDIF
                         PP = abs(PP)
                         ! Evaluation of the minimum
                         IF ((3d0 * xm * QQ - abs(Tol1 * QQ)) .LT.&
     &                                             (abs(EE * QQ))) THEN
                            Minimum = 3d0 * xm * QQ - abs(Tol1 * QQ)
                         ELSE
                            Minimum = abs(EE * QQ)
                         ENDIF
!
                         IF ((2d0 * PP) .LT. Minimum) THEN
                         ! accept interpolation
                             EE = DD
                             DD = PP / QQ
                         ELSE
                         ! interpolation failed, use bisection
                             DD = xm
                             EE = DD
                         ENDIF
                     ELSE
                     ! bounds decreasing too slowly, use bisection
                         DD = xm
                         EE = DD
                     ENDIF
                     AA = BB  ! move last best guess to AA
                     FA = FB
                     IF (abs(DD) .GT. Tol1) THEN
                     ! evaluate new trial root
                         BB = BB + DD
                     ELSE
                         IF (xm .GT. 0d0) THEN
                             BB = BB + abs(Tol1)
                         ELSE
                             BB = BB - abs(Tol1)
                         ENDIF
                     ENDIF
                     CALL F_Brent(MODE, FB, out_2, BB, GGG, in_1, in_2)
                     IF (Niter .GT. MAX_ITER) THEN
                        print*, 'Max N iteration exceeded in Brent'
                        STOP
                     ENDIF
                 ENDIF
       END DO
      ENDIF
!
      out_1 = Root
!
END SUBROUTINE BrentRoots
!
!
!===============================================================================
      SUBROUTINE F_Brent(MODE, F, out_2, XXX, GGG, in_1, in_2)
!===============================================================================
          IMPLICIT NONE
          INTEGER, INTENT(IN)  :: MODE
          REAL(8), INTENT(IN)  :: XXX, GGG, in_1, in_2
          REAL(8), INTENT(OUT) :: F, out_2
          REAL(8)              :: temp, v_v, v_l, e_v, e_l,e,press,v,qual


          IF (MODE .EQ. 1) THEN
!
! In two-phase region, (e,v) --> (psat,x,Tsat)        
        press = XXX
        v = in_1
        e = GGG
        CALL satprop(5, press, temp, v_v, v_l, e_v, e_l)

!        print*, 'iter process',v,v_l, v_v
        out_2 = temp
        qual = (v - v_l)/(v_v-v_l)
        F = (e - qual*e_v - (1_pr - qual)*e_l) * 1e-3_pr
!
!
          ELSE
             STOP 'MODE of BrentRoots unknown'
          ENDIF
!
     END SUBROUTINE F_Brent
!
END MODULE solver_eos
