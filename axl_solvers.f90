MODULE axl_solvers
!
! It contains algorithms and numerical methods used for the inversion of the
! Span-Wagner EoS, necessary for the evaluation of properties in each node of the
! grid.
! 
! It contains a 2D Newton-Raphson algorithm, a 1D Newton-Raphson algorithm, a
! Brent algorithm and the corresponding functions of which we seek the zeros.
!==============================================================================
!============================================================================== 
!
USE def_constants
USE properties
!USE def_variables, ONLY: vL_psat_spline, vV_psat_spline,&
!                        &uL_psat_spline, uV_psat_spline,&
!                        &Tsat_psat_spline, saturP
!USE Interp_table
IMPLICIT NONE
!
!PRIVATE
!PUBLIC  New_Rap2D, New_Rap3D, New_Rap1D, BrentRoots
CONTAINS
!
!
!
!===============================================================================
      SUBROUTINE New_Rap2Daxl(MODE, out1_NewRap, out2_NewRap, &
     &           resnorm, Niter, exitflag, in_1, in_2,in_3,in_4, guess_1, guess_2)
!===============================================================================
!
      IMPLICIT NONE
!
! Global variables
      INTEGER, INTENT(IN)  :: MODE
      INTEGER, INTENT(OUT) :: exitflag, Niter
      REAL(pr),INTENT(IN)  :: guess_1, guess_2, in_1, in_2,in_3,in_4
      REAL(pr),INTENT(OUT) :: out1_NewRap, out2_NewRap, resnorm
!
! Local variables                      
      INTEGER              :: n, nx, nf, INFO, MAX_ITER                 
      INTEGER,POINTER      :: IPIV(:)
      REAL(8) ::   x, TOL_X,  TOL_FUN, ALPHA,  MIN_LAMBDA, MAX_LAMBDA,&
     &          Dx_J, lambda, slope,   fffold, lambda_min, fff,        &
     &          fff2, A_l, lambda2, lambda_OLD, aa, bb, discriminant 
      REAL(8), DIMENSION (2,1)  ::  XXX, F, delta, XXX_forw, XXX_back,&
     &       F_forw, F_back, dF, TYPX, dx_star, dx, XXX_old, C_l, aabb  
      REAL(8), DIMENSION (1,2)  ::  g
      REAL(8), DIMENSION (2,2)  ::  J, J0, Jstar, AAA, B_l
      !
      ! ------------- INITIALIZATION ----------
      XXX(1,1) = guess_1 
      XXX(2,1) = guess_2
      !
      !   Options
      !
      TOL_X      = 1.d-20           ! relative max step tolerance-20
      TOL_FUN    = 1.d-12           ! function tolerance -12
      MAX_ITER   = 500             ! max number of iterations
      ALPHA      = 1.d-4            ! criteria for decrease -4
      MIN_LAMBDA = 0.1d0            ! min lambda 0.1,0.5
      MAX_LAMBDA = 0.5d0 
      TYPX       = abs(XXX)         ! x scaling value, remove zeros
      !
      !
      call F_New_Rap2D (MODE, F, XXX, in_1, in_2,in_3,in_4)
      !
      ! ---- Jacobian estimation
      !
      Dx_J = 0.7d-9        ! finite difference delta
      nx   = 2             ! degrees of freedom
      nf   = 2             ! number of functions
      DO n = 1,nx
         delta(1,1) = 0.d0; delta(2,1)    = 0.d0
         delta(n,1) = delta(n,1) + Dx_J 
         XXX_forw   = XXX + delta 
         XXX_back   = XXX - delta
         call F_New_Rap2D (MODE, F_forw, XXX_forw, in_1, in_2,in_3,in_4)
         call F_New_Rap2D (MODE, F_back, XXX_back, in_1, in_2,in_3,in_4)
         dF      = F_forw - F_back   ! delta F
         J(1, n) = dF(1,1)/Dx_J/2.d0 ! derivatives dF/d_n
         J(2, n) = dF(2,1)/Dx_J/2.d0 ! derivatives dF/d_n
      ENDDO
      J0(1,1) = 1.d0 / TYPX(1,1)  ! Jacobian scaling matrix
      J0(1,2) = 1.d0 / TYPX(2,1)
      J0(2,1) = J0(1,1)
      J0(2,2) = J0(1,2)
      Jstar   = J/J0;             ! scale Jacobian
      !
      resnorm = sqrt(F(1,1)*F(1,1) + F(2,1)*F(2,1))
      exitflag = 1
      !
      !--- SOLVER -------------------
      !
      Niter  = 0     ! start counter
      lambda = 1.d0  ! backtracking
      lambda_OLD = lambda
!
      allocate(IPIV(nf))
!           
 1    DO WHILE (((resnorm .GT. TOL_FUN) .OR. (lambda .LT. 1.d0)) &
     &        .AND.  (exitflag >= 0)   .AND. (Niter <= MAX_ITER)) 
               Niter = Niter + 1
              ! 
               IF (lambda==1.d0) THEN       
              ! Newton-Raphson solver
              ! dx_star  = -Jstar\F      calculate Newton step
                 AAA     = -Jstar
                 dx_star = F
                !call PLU decomposition routine from LAPACK
                 call DGETRF( nf, nf,  AAA,nf, IPIV, INFO)
                 call DGETRS('N', nf,1,AAA,nf, IPIV, dx_star,nf,INFO)   
      	         dx      = dx_star*TYPX
                 g       = MATMUL(TRANSPOSE(F),Jstar)
      	         slope   = g(1,1)*dx_star(1,1) + g(1,2)*dx_star(2,1)
      	         fffold  = SUM(F*F)
                 XXX_old = XXX      
                 lambda_min = TOL_X/MAXVAL(ABS(dx)/ABS(XXX_old))
               ENDIF
!      
!        print*, 'LAMBDA = ', lambda, 'resnorm',resnorm
               IF (lambda < lambda_min) THEN
                  exitflag = -2; ! XXX is too close to XXX_OLD 
               !   print*,'XXX is too close to XXX_OLD'
               !ELSEIF any(isnan(dx)) || any(isinf(dx))
               !    exitflag = -1; ! matrix may be singular
               !    STOP
               ENDIF
! 
! XXX at new step  
!   	
      	 XXX = XXX_old + dx*lambda
!        print*, 'resnorm',resnorm
!        print*,'dx*lambda= ',dx*lambda
!        print*,'xxx', XXX
!        print*, 'lambda_min', lambda_min
         call F_New_Rap2D (MODE, F, XXX, in_1, in_2,in_3,in_4)
         !
         ! ---- Jacobian estimation
         !
         DO n = 1,nx
     	  delta(1,1)    = 0.d0; delta(2,1)    = 0.d0
     	  delta(n,1)    = delta(n,1) + Dx_J 
     	  XXX_forw   = XXX + delta 
     	  XXX_back   = XXX - delta
          call F_New_Rap2D (MODE, F_forw, XXX_forw, in_1, in_2,in_3,in_4)
          call F_New_Rap2D (MODE, F_back, XXX_back, in_1, in_2,in_3,in_4)
     	  dF      = F_forw - F_back    ! delta F
     	  J(1, n) = dF(1,1)/Dx_J/2.d0  ! derivatives dF/d_n
     	  J(2, n) = dF(2,1)/Dx_J/2.d0 ! derivatives dF/d_n
      	ENDDO      
      	Jstar  = J/J0        ! scale Jacobian
      	fff    = SUM(F*F)    ! f to be minimized
!      
! Optimization
!      
          IF (fff > fffold + ALPHA*lambda*slope) THEN
              IF (lambda .eq. 1.d0) THEN
                  lambda = -slope/2.d0/(fff-fffold-slope) ! calculate lambda
              ELSE
                  A_l      =  1.d0/(lambda_OLD - lambda2)
                  B_l(1,1) =  1.d0/lambda_OLD**2.d0
                  B_l(1,2) = -1.d0/lambda2**2.d0
                  B_l(2,1) = -lambda2/lambda_OLD**2.d0
                  B_l(2,2) =  lambda_OLD/lambda2**2.d0
       	          C_l(1,1) =  fff-fffold-lambda_OLD*slope
                  C_l(2,1) =  fff2-fffold-lambda2*slope
                  !
      	          aabb = A_l* MATMUL(B_l,C_l)
                  aa   = aabb(1,1)
                  bb   = aabb(2,1)
                 IF (aa == 0.d0) THEN
                      lambda = -slope/2.d0/bb;
                 ELSE
                      discriminant = bb**2.d0 - 3.d0*aa*slope;
                      IF (discriminant < 0.d0) THEN
                          lambda = MAX_LAMBDA*lambda_OLD;
                      ELSEIF (bb <= 0.d0) THEN
                          lambda = (-bb + sqrt(discriminant))/3.d0/aa
                      ELSE
                          lambda = -slope/(bb + sqrt(discriminant))
                      ENDIF
      	   ENDIF
                 lambda = min(lambda,MAX_LAMBDA*lambda_OLD)! minimum step length
              ENDIF      
          ELSE
              lambda = 1.d0; ! fraction of Newton step
          ENDIF      
          IF (lambda < 1.d0) THEN
              lambda2 = lambda_OLD;
              fff2    = fff;                     ! save 2nd most previous value
              lambda  = max(lambda,MIN_LAMBDA*lambda_OLD) ! minimum step length
              GO TO 1
          ENDIF
      !
      !
          resnorm = sqrt(F(1,1)*F(1,1) + F(2,1)*F(2,1))
      !
      ENDDO
      !
      deallocate(IPIV)
      out1_NewRap = XXX(1,1)
      out2_NewRap = XXX(2,1)     
      END SUBROUTINE New_Rap2Daxl
!
!
!
      SUBROUTINE F_New_Rap2D (MODE, F, XXX, in_1, in_2,in_3,in_4)     
      IMPLICIT NONE 
      INTEGER,  INTENT(IN) :: MODE         
      REAL(pr), INTENT(IN) :: in_1, in_2, in_3,in_4
      REAL(pr), DIMENSION (2,1), INTENT(IN)   :: XXX
      REAL(pr), DIMENSION (2,1), INTENT(OUT)  :: F
      REAL(pr), DIMENSION (2,1)  :: F_out
      REAL(pr) :: F1_MODE, F2_MODE, Helmholtz, T, v, p, e, cv, cp, s, c, &
     &            v_l, p_l, s_l, e_l, g_l, v_v, p_v, s_v, e_v, g_v,&
     &            helmho_v,helmho_l,deriv,dummy1,dummy2,x_out,a_out,entro,enth
      REAL(pr) :: psat, Tsat, vvsat, vlsat, uvsat, ulsat,sv,sl,qual,velo,rhoinv
      REAL(pr) :: p1,p2,p3,duL_dp, duV_dp, dvL_dp, dvV_dp,sound,du_dp_x,dv_dp_x,ratio
!           
      IF     (MODE .EQ. 1) THEN
! Calculation of the saturation specific volumes of vapor 
! and liquid phase.
! The saturation specific volume is obtained by imposing:
!  p_l = p_v         pressure equality
!  g_l = g_v         Gibbs free enthalpy equality
!
! That is the Maxwell criterion for saturation.
!
          v_l = XXX(1,1)
          v_v = XXX(2,1)
          T   = in_1
!          
          CALL pressure(T,v_l,p_l)
          CALL pressure(T,v_v,p_v)
!         
          CALL helmho(T,v_l,helmho_l)
          CALL helmho(T,v_v,helmho_v) 
!          g_l    = Helmholtz(v_l,T) + p_l*v_l*1e-3_pr      ! kJ / kg
!          g_v    = Helmholtz(v_v,T) + p_v*v_v*1e-3_pr      ! kJ / kg
!         
          g_l    = helmho_l + p_l*v_l      ! J / kg
          g_v    = helmho_v + p_v*v_v
          F(1,1) = (p_l - p_v)*1e-6_pr
          F(2,1) = (g_l - g_v)*1e-3_pr
!
      ELSEIF (MODE .EQ. 2) THEN
! Calculation of temperature and specific volume for imposed pressure p
! and specific internal energy, e.
! 
! For this function:
!         in_1     = p
!         in_2     = e
!         XXX(1,1) = T
!         XXX(2,1) = v
!
! 
          T   = XXX(1,1)
          v   = XXX(2,1)
          CALL pressure(T,v,p)
          CALL inter_energy(T,v,e)
          F(1,1) = (in_1 - p)*1e-6_pr         ! Pa
          F(2,1) = (in_2 - e)*1e-3_pr         ! J/kg
!print*, in_1,in_2, p,e
!
!
      ELSEIF (MODE .EQ. 3) THEN
! META liquid: dpdv_T  = 0 and e(t,v) = e_in
! For this function:
!         in_1     = e
!         in_2     = 
!         XXX(1,1) = T
!         XXX(2,1) = v_l
!
! 
          T   = XXX(1,1)
          v_l   = XXX(2,1)
!
        CALL inter_energy(T,v_l,e)
        CALL dpdv_T(T,v_l,deriv)
!
          F(1,1) = deriv * 1e-6_pr
          F(2,1) = (in_1 - e) *1e-3_pr          
!
!
      ELSEIF (MODE .EQ. 4) THEN
! isentropic
          psat      = XXX(1,1)
          rhoinv    = XXX(2,1)
!          s      = in1
!          htot   = in_2 
!          dA     = in_3
!          u(i-1) = in_4
!
          CALL satprop(3, psat, Tsat, vvsat, vlsat, uvsat, ulsat)
          CALL axlpress(Tsat,vvsat,vlsat,p1,p2,p3)
          psat = (psat+p1+p2+p3)/4.0
          CALL satprop(3, psat, Tsat, vvsat, vlsat, uvsat, ulsat)

          CALL entropy(Tsat, vvsat, sv)       
          CALL entropy(Tsat, vlsat, sl)
!     
         qual = abs((rhoinv - vlsat) / (vvsat-vlsat))
!         print*, 'v',rhoinv, 'qual',qual, vlsat
        CALL satderiv(3, psat, duL_dp, duV_dp, dvL_dp, dvV_dp)


!!
        ratio   = ((uvsat - ulsat)/(vvsat - vlsat))   ! (J/kg)/(m3/kg)
        du_dp_x = x_out * duV_dp + (1_pr - x_out) * duL_dp
        dv_dp_x = x_out * dvV_dp + (1_pr - x_out) * dvL_dp

!!
        sound    = SQRT((psat  + ratio)/(du_dp_x - ratio * dv_dp_x)) * vlsat
        velo     = 1.0/((in_4/sound)-1) * in_4*in_3 + in_4 

!        print*, 'velo', velo
!         velo = in_3 * rhoinv

         enth = in_2 - 0.5*velo*velo  
!       print*, 'p,x',psat,qual,'T',Tsat
!       print*, 'Sv',sv, 'Sl',sl
!       print*, ' '
          F(1,1) =   (in_1 - qual*sv + (1.0-qual)*sl) *1e-3_pr
          F(2,1) =  (enth - (uvsat+psat*vvsat)*qual - (ulsat+psat*vlsat)*(1.0-qual))  *1e-3_pr

      ELSE
         print*,'MODE of New_Rap2D unknown'
         STOP
      ENDIF
!
      END SUBROUTINE F_New_Rap2D
!

!
!===============================================================================
      SUBROUTINE New_Rap3Daxl(MODE, out1_NewRap, out2_NewRap, out3_NewRap, &
     &           resnorm, Niter, exitflag, in_1,in_2, guess_1, guess_2, guess_3)
!===============================================================================
!
      IMPLICIT NONE
!
! Global variables
      INTEGER, INTENT(IN)  :: MODE
      INTEGER, INTENT(OUT) :: exitflag, Niter
      REAL(pr),INTENT(IN)  :: guess_1, guess_2, guess_3, in_1, in_2!,in_3
      REAL(pr),INTENT(OUT) :: out1_NewRap, out2_NewRap, out3_NewRap, resnorm
!
! Local variables                      
      INTEGER              :: n, nx, nf, INFO, MAX_ITER                 
      INTEGER,POINTER      :: IPIV(:)
      REAL(8) ::   x, TOL_X,  TOL_FUN, ALPHA,  MIN_LAMBDA, MAX_LAMBDA,&
     &          Dx_J, lambda, slope,   fffold, lambda_min, fff,        &
     &          fff2, A_l, lambda2, lambda_OLD, aa, bb, discriminant 
      REAL(8), DIMENSION (3,1)  ::  XXX, F, delta, XXX_forw, XXX_back,&
     &       F_forw, F_back, dF, TYPX, dx_star, dx, XXX_old, C_l, aabb  
      REAL(8), DIMENSION (1,3)  ::  g
      REAL(8), DIMENSION (3,3)  ::  J, J0, Jstar, AAA, B_l
      !
      ! ------------- INITIALIZATION ----------
      XXX(1,1) = guess_1 
      XXX(2,1) = guess_2
      XXX(3,1) = guess_3
      !
      !   Options
      !
      TOL_X      = 1.d-20           ! relative max step tolerance
      TOL_FUN    = 1.d-12           ! function tolerance
      MAX_ITER   = 500             ! max number of iterations
      ALPHA      = 1.d-4            ! criteria for decrease
      MIN_LAMBDA = 0.1d0            ! min lambda
      MAX_LAMBDA = 0.5d0 
      TYPX       = abs(XXX)         ! x scaling value, remove zeros
      !
      !
      call F_New_Rap3D (MODE, F, XXX, in_1,in_2)
      !
      ! ---- Jacobian estimation
      !
      Dx_J = 1.d-6         ! finite difference delta
      nx   = 3             ! degrees of freedom
      nf   = 3             ! number of functions
      DO n = 1,nx
         delta(1,1) = 0.d0; delta(2,1) = 0.d0; delta(3,1) = 0.d0
         delta(n,1) = delta(n,1) + Dx_J 
         XXX_forw   = XXX + delta 
         XXX_back   = XXX - delta
         call F_New_Rap3D (MODE, F_forw, XXX_forw, in_1,in_2)
         call F_New_Rap3D (MODE, F_back, XXX_back, in_1,in_2)
         dF      = F_forw - F_back   ! delta F
         J(1, n) = dF(1,1)/Dx_J/2.d0 ! derivatives dF/d_n
         J(2, n) = dF(2,1)/Dx_J/2.d0 ! derivatives dF/d_n
         J(3, n) = dF(3,1)/Dx_J/2.d0 ! derivatives dF/d_n
      ENDDO
      J0(1,1) = 1.d0 / TYPX(1,1)  ! Jacobian scaling matrix
      J0(1,2) = 1.d0 / TYPX(2,1)
      J0(1,3) = 1.d0 / TYPX(3,1)
      J0(2,1) = J0(1,1)
      J0(2,2) = J0(1,2)
      J0(2,3) = J0(1,3)
      J0(3,1) = J0(1,1)
      J0(3,2) = J0(1,2)
      J0(3,3) = J0(1,3)
      Jstar   = J/J0;             ! scale Jacobian
      !
      resnorm = sqrt(F(1,1)*F(1,1) + F(2,1)*F(2,1)+ F(3,1)*F(3,1))
      exitflag = 1
      !
      !--- SOLVER -------------------
      !
      Niter  = 0     ! start counter
      lambda = 1.d0  ! backtracking
      lambda_OLD = lambda
!
      allocate(IPIV(nf))
!           
 1    DO WHILE (((resnorm .GT. TOL_FUN) .OR. (lambda .LT. 1.d0)) &
     &        .AND.  (exitflag >= 0)   .AND. (Niter <= MAX_ITER)) 
               Niter = Niter + 1
              ! 
               IF (lambda==1.d0) THEN       
              ! Newton-Raphson solver
              ! dx_star  = -Jstar\F      calculate Newton step
                 AAA     = -Jstar
                 dx_star = F
                !call PLU decomposition routine from LAPACK
                 call DGETRF( nf, nf,  AAA,nf, IPIV, INFO)
                 call DGETRS('N', nf,1,AAA,nf, IPIV, dx_star,nf,INFO)   
      	         dx      = dx_star*TYPX
                 g       = MATMUL(TRANSPOSE(F),Jstar)
      	         slope   = g(1,1)*dx_star(1,1) + g(1,2)*dx_star(2,1) +&
                &          g(1,3)*dx_star(3,1)
      	         fffold  = SUM(F*F)
                 XXX_old = XXX      
                 lambda_min = TOL_X/MAXVAL(ABS(dx)/ABS(XXX_old))
               ENDIF
!      
!        print*, 'LAMBDA = ', lambda, 'resnorm',resnorm
               IF (lambda < lambda_min) THEN
                  exitflag = -2; ! XXX is too close to XXX_OLD 
               !   print*,'XXX is too close to XXX_OLD'
               !ELSEIF any(isnan(dx)) || any(isinf(dx))
               !    exitflag = -1; ! matrix may be singular
               !    STOP
               ENDIF
! 
! XXX at new step  
!   	
      	 XXX = XXX_old + dx*lambda
!        print*, 'resnorm',resnorm
!        print*,'dx*lambda= ',dx*lambda
!        print*,'xxx', XXX
!        print*, 'lambda_min', lambda_min
         call F_New_Rap3D (MODE, F, XXX, in_1,in_2)
         !
         ! ---- Jacobian estimation
         !
         DO n = 1,nx
     	  delta(1,1)    = 0.d0; delta(2,1) = 0.d0; delta(3,1)=0.d0
     	  delta(n,1)    = delta(n,1) + Dx_J 
     	  XXX_forw   = XXX + delta 
     	  XXX_back   = XXX - delta
          call F_New_Rap3D (MODE, F_forw, XXX_forw, in_1,in_2)
          call F_New_Rap3D (MODE, F_back, XXX_back, in_1,in_2)
     	  dF      = F_forw - F_back    ! delta F
     	  J(1, n) = dF(1,1)/Dx_J/2.d0  ! derivatives dF/d_n
     	  J(2, n) = dF(2,1)/Dx_J/2.d0 ! derivatives dF/d_n
     	  J(3, n) = dF(3,1)/Dx_J/2.d0 ! derivatives dF/d_n
      	ENDDO      
      	Jstar  = J/J0        ! scale Jacobian
      	fff    = SUM(F*F)    ! f to be minimized
!      
! Optimization
!      
          IF (fff > fffold + ALPHA*lambda*slope) THEN
              IF (lambda .eq. 1.d0) THEN
                  lambda = -slope/2.d0/(fff-fffold-slope) ! calculate lambda
              ELSE
                  A_l      =  1.d0/(lambda_OLD - lambda2)
                  B_l(1,1) =  1.d0/lambda_OLD**2.d0
                  B_l(1,2) = -1.d0/lambda2**2.d0
                  B_l(2,1) = -lambda2/lambda_OLD**2.d0
                  B_l(2,2) =  lambda_OLD/lambda2**2.d0
       	          C_l(1,1) =  fff-fffold-lambda_OLD*slope
                  C_l(2,1) =  fff2-fffold-lambda2*slope
                  !
      	          aabb = A_l* MATMUL(B_l,C_l)
                  aa   = aabb(1,1)
                  bb   = aabb(2,1)
                 IF (aa == 0.d0) THEN
                      lambda = -slope/2.d0/bb;
                 ELSE
                      discriminant = bb**2.d0 - 3.d0*aa*slope;
                      IF (discriminant < 0.d0) THEN
                          lambda = MAX_LAMBDA*lambda_OLD;
                      ELSEIF (bb <= 0.d0) THEN
                          lambda = (-bb + sqrt(discriminant))/3.d0/aa
                      ELSE
                          lambda = -slope/(bb + sqrt(discriminant))
                      ENDIF
      	   ENDIF
                 lambda = min(lambda,MAX_LAMBDA*lambda_OLD)! minimum step length
              ENDIF      
          ELSE
              lambda = 1.d0; ! fraction of Newton step
          ENDIF      
          IF (lambda < 1.d0) THEN
              lambda2 = lambda_OLD;
              fff2    = fff;                     ! save 2nd most previous value
              lambda  = max(lambda,MIN_LAMBDA*lambda_OLD) ! minimum step length
              GO TO 1
          ENDIF
      !
      !
          resnorm = sqrt(F(1,1)*F(1,1) + F(2,1)*F(2,1) + F(3,1)*F(3,1))
      !
      ENDDO
      !
      deallocate(IPIV)
      out1_NewRap = XXX(1,1)
      out2_NewRap = XXX(2,1)     
      out3_NewRap = XXX(3,1)     
      END SUBROUTINE New_Rap3Daxl
!
!
!
      SUBROUTINE F_New_Rap3D (MODE, F, XXX, in_1,in_2)     
      IMPLICIT NONE 
      INTEGER,  INTENT(IN) :: MODE         
      REAL(pr), INTENT(IN) :: in_1, in_2!, in_3
      REAL(pr), DIMENSION (3,1), INTENT(IN)   :: XXX
      REAL(pr), DIMENSION (3,1), INTENT(OUT)  :: F
      REAL(pr), DIMENSION (3,1)  :: F_out
      REAL(pr) :: F1_MODE, F2_MODE, Helmholtz, T, v, p, e, cv, cp, s, c, &
     &            v_l, p_l, s_l, e_l, g_l, v_v, p_v, s_v, e_v, g_v,&
     &            helmho_v,helmho_l,x,y,z,qual         
!           
      IF (MODE .EQ. 1) THEN
! FOR LL right boundary (liquid saturation line):
!  p_l = p_v         pressure equality
!  g_l = g_v         Gibbs free enthalpy equality
!  el = e            at one interal energy
! That is the Maxwell criterion for saturation.
!
          v_l = XXX(1,1)
          v_v = XXX(2,1)                
          T   = XXX(3,1)
          e = in_1      
!          
          CALL pressure(T,v_l,p_l)
          CALL pressure(T,v_v,p_v)
!         
          CALL helmho(T,v_l,helmho_l)
          CALL helmho(T,v_v,helmho_v)
!
          CALL inter_energy(T,v_l,e_l)         
!         
          g_l    = helmho_l + p_l*v_l      ! J / kg
          g_v    = helmho_v + p_v*v_v
          F(1,1) = (p_l - p_v)*1e-6_pr
          F(2,1) = (g_l - g_v)*1e-3_pr
          F(3,1) = (e_l - e)*1e-3_pr      
!cas test for 3D solver
      ELSEIF (MODE .EQ. 2) THEN
          x=  XXX(1,1)
          y=  XXX(2,1)
          z=  XXX(3,1)
!
          F(1,1) = sin(x)
          F(2,1) = cos(y)
          F(3,1) = z*z*z - 6_pr*z*z + 11_pr*z -6
!
!           
      ELSEIF (MODE .EQ. 3) THEN
!  FOR TP region
!  p=p_l         pressure equality
!  p =p_v        pressure equality
!  g_l = g_v     Gibbs free enthalpy equality
! That is the Maxwell criterion for saturation.
!
          v_l = XXX(1,1)
          v_v = XXX(2,1)
          T   = XXX(3,1)
          p = in_1
!          
          CALL pressure(T,v_l,p_l)
          CALL pressure(T,v_v,p_v)
!         
          CALL helmho(T,v_l,helmho_l)
          CALL helmho(T,v_v,helmho_v)
!
!         
          g_l    = helmho_l + p_l*v_l      ! J / kg
          g_v    = helmho_v + p_v*v_v
!
          F(1,1) = (p_l - p_v)*1e-6_pr
          F(2,1) = (p_v - p)*1e-6_pr
          F(3,1) = (g_l - g_v)*1e-3_pr
      ELSEIF (MODE .EQ. 4) THEN
!  FOR LH right boundary (vapor saturation line)
!  p_l = p_v         pressure equality
!  g_l = g_v         Gibbs free enthalpy equality
!  ev = e            at one interal energy
! That is the Maxwell criterion for saturation.
!
          v_l = XXX(1,1)
          v_v = XXX(2,1)
          T   = XXX(3,1)
          e = in_1
!          
          CALL pressure(T,v_l,p_l)
          CALL pressure(T,v_v,p_v)
!         
          CALL helmho(T,v_l,helmho_l)
          CALL helmho(T,v_v,helmho_v)
!
          CALL inter_energy(T,v_v,e_v)
!         
          g_l    = helmho_l + p_l*v_l      ! J / kg
          g_v    = helmho_v + p_v*v_v
          F(1,1) = (p_l - p_v)*1e-6_pr
          F(2,1) = (g_l - g_v)*1e-3_pr
          F(3,1) = (e_v - e)*1e-3_pr
!
     ELSEIF (MODE .EQ. 5) THEN
! for a given input couple of values of specific internal energy and
! specific volume, the temperature, vapor specific volume, liquid specific
! volume, are evalued at sauration and the corresponding quality is
! computed. (TP region)
!   (v,e)--> (Vv, Vl, T, x)  
!   F1 = pv = pl
!   F2 = gv = gl
!   F3 = e - x*ev - (1-x)*el
!   x = (v-Vl)/(Vv - Vl)
!
          v_l = XXX(1,1)
          v_v = XXX(2,1)
          T   = XXX(3,1)
          e = in_1
          v = in_2
!          
          CALL pressure(T,v_l,p_l)
          CALL pressure(T,v_v,p_v)
!         
          CALL helmho(T,v_l,helmho_l)
          CALL helmho(T,v_v,helmho_v)
!
          CALL inter_energy(T,v_v,e_v)
          CALL inter_energy(T,v_l,e_l)
!         
          qual = (v - v_l)/(v_v-v_l)
          g_l    = helmho_l + p_l*v_l      ! J / kg
          g_v    = helmho_v + p_v*v_v
!
          F(1,1) = (p_l - p_v)*1e-6_pr
          F(2,1) = (g_l - g_v)*1e-3_pr
          F(3,1) = (e - qual*e_v - (1_pr - qual)*e_l)*1e-3_pr        
!
      ELSE
         print*,'MODE of New_Rap3D unknown'
         STOP
      ENDIF
    END SUBROUTINE F_New_Rap3D
!
!===============================================================================
          SUBROUTINE New_Rap1Daxl(MODE, out_1, out_2, resnorm, Niter,&      
     &                             exitflag, GGG, X0, in_1, out3)
!===============================================================================
!         1D         TWO-PHASE or SINGLE-PHASE
!-------------------------------------------------------------------
!
!  Input :
!  -------
!     MODE        =  it will depend on the used function
!     GGG         =  idem
!     in_1        =  idem
!     in_2        =  idem
!     X0          =  idem
!
!
!  Output :
!  -------
!
!     out_1       =  idem
!     out_2       =  idem
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
      
      END SUBROUTINE New_Rap1Daxl
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
          REAL(8) :: temp, e, v
          REAL(8) :: delta_p,press,e_l,e_v,V_v,V_l,qual 
          
          IF (MODE .EQ. 1) THEN
! (v,e) --> T: v=in_1, e=GGG, T=XXX
                v = in_1
                temp = XXX
!                
                CALL inter_energy(temp, v, e)
!                
                F = (e - GGG) * 1e-3_pr
                out_2 = 0.0_pr
                out_3 = 0.0_pr
!
!        
        ELSEIF (MODE .EQ. 2) THEN
! In two-phase region, (e,v) --> (psat,x,Tsat) 
!       
        press = XXX
        v = in_1
        e = GGG 
!        CALL satprop(3, press, temp, v_v, v_l, e_v, e_l)

!        print*, 'iter process',v,v_l, v_v
        out_3 = temp 
        qual = (v - v_l)/(v_v-v_l)
        out_2 = qual
        F = (e - qual*e_v - (1_pr - qual)*e_l) * 1e-3_pr
!
!       
        ELSEIF (MODE .EQ. 3) THEN
        v     = XXX
        temp  = in_1 
        CALL pressure(temp,v,press)
        F = (press - GGG) * 1e-6_pr
              
        ELSE
                print*,'MODE of New_Rap1D unknown'
                STOP
          ENDIF
!
      END SUBROUTINE F_New_Rap1D

!=======================================================================================
!
SUBROUTINE BrentRootsaxl(MODE, out_1, out_2, residue, Niter, GGG, lower, upper, in_1, in_2)
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
END SUBROUTINE BrentRootsaxl
!
!
!===============================================================================
      SUBROUTINE F_Brent(MODE, F, out_2, XXX, GGG, in_1, in_2)
!===============================================================================
          IMPLICIT NONE
          INTEGER, INTENT(IN)  :: MODE
          REAL(8), INTENT(IN)  :: XXX, GGG, in_1, in_2
          REAL(8), INTENT(OUT) :: F, out_2
          REAL(8)              :: temp, v_v, v_l, e_v, e_l,e,press,v,qual,du_dp_x,dv_dp_x,ratio
          REAL(8)              :: enth, htot,velo,sv,sl,entro,duL_dp, duV_dp, dvL_dp, dvV_dp


    IF (MODE .EQ. 1) THEN
!
! In two-phase region, (e,v) --> (psat,x,Tsat)        
        press = XXX
        v = in_1
        e = GGG
!        CALL satprop(3, press, temp, v_v, v_l, e_v, e_l)

!        print*, 'iter process',v,v_l, v_v
        out_2 = temp
        qual = (v - v_l)/(v_v-v_l)
        F = (e - qual*e_v - (1_pr - qual)*e_l) * 1e-3_pr
!
!
    ELSEIF (MODE .EQ. 2) THEN
! isentropic
          press    = XXX
          enth     = GGG
!          velo     = in_1
!
          CALL satprop(3, press, temp, v_v, v_l, e_v, e_l)
! 
!          print*, 'T=',temp,'v_L',v_l, 'P', press
!         velo = in_1 * v_l
!        CALL satderiv(3, press, duL_dp, duV_dp, dvL_dp, dvV_dp)


!!
!        ratio   = ((e_v - e_l)/(v_v - v_l))   ! (J/kg)/(m3/kg)
!        du_dp_x =  duL_dp
!        dv_dp_x =  dvL_dp
!!
!        velo    = SQRT((press  + ratio)/(du_dp_x - ratio * dv_dp_x)) * v_l
         print*, 'sound liquid', velo
          
         F =  (enth - e_l - press*v_l - 0.5*velo*velo )  *1e-3_pr
!         F =  (entro - sl)  *1e-3_pr

    ELSE
         STOP 'MODE of BrentRoots unknown'
    ENDIF
!
     END SUBROUTINE F_Brent
!
!
!
!
!
!
!
END MODULE axl_solvers
