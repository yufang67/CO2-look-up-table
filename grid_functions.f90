!========================================================================
!@ The Module contains the routine to help construct the grid.
!
!@ polyfit: compute the coefficients for the construction of splines.
!  Max_point: compute the point on the sautration line with maximun e
!  corner_A, corner_B, corner_c: compute some points at the intersection
!  between domains.
!=======================================================================
!
!
MODULE grid_functions

     USE def_constants
     USE def_variables
     USE non_linear_solvers
     USE properties
!
      IMPLICIT NONE
!
     PRIVATE
!
     PUBLIC   polyfit, Max_point,&  
&               corner_A, corner_B, corner_C, corner_D, corner_E
!
       CONTAINS
!
! ===================================================================
!
!                              polyfit 
!
! ===================================================================
      SUBROUTINE polyfit(p,x,y,N)
!POLYFIT Fit polynomial to data.
!   P = POLYFIT(X,Y,N) finds the coefficients of a polynomial P(X) of
!   degree N that fits the data Y best in a least-squares sense. P is a
!   row vector of length N+1 containing the polynomial coefficients in
!   descending powers, P(1)*X^N + P(2)*X^(N-1) +...+ P(N)*X + P(N+1).
!
!   [P,S] = POLYFIT(X,Y,N) returns the polynomial coefficients P and a
!   structure S for use with POLYVAL to obtain error estimates for
!   predictions.  S contains fields for the triangular factor (R) from a
!   QR
!   decomposition of the Vandermonde matrix of X, the degrees of freedom
!   (df), and the norm of the residuals (normr).  If the data Y are
!   random,
!   an estimate of the covariance matrix of P is
!   (Rinv*Rinv')*normr^2/df,
!   where Rinv is the inverse of R.
!
!   [P,S,MU] = POLYFIT(X,Y,N) finds the coefficients of a polynomial in
!   XHAT = (X-MU(1))/MU(2) where MU(1) = MEAN(X) and MU(2) = STD(X).
!   This
!   centering and scaling transformation improves the numerical
!   properties
!   of both the polynomial and the fitting algorithm.
!
!   Warning messages result if N is >= length(X), if X has repeated, or
!   nearly repeated, points, or if X might need centering and scaling.
!
        INTEGER, INTENT(IN)  ::  N
        INTEGER           ::  i, j, INFO
        INTEGER,allocatable   ::  IPIV(:)
        REAL(pr), DIMENSION (N+1), INTENT(OUT)     ::  p
        REAL(pr), DIMENSION (N+1), INTENT(IN)      ::  x, y
        REAL(pr), DIMENSION (N+1)     ::  TAU, b
        REAL(pr), DIMENSION (N+1,N+1) ::  Q, WORK, V, QR, R

        IF (size(x) .NE. size(y)) THEN
            print*, 'polyfit: XYSizeMismatch'
        ENDIF

!if nargout > 2
!   mu = [mean(x); std(x)];
!   x = (x - mu(1))/mu(2);
!end

! Construct Vandermonde matrix.
        V         = 0_pr
        V(:, N+1) = 1_pr
        DO j = N,1,-1
           V(:,j) = x*V(:,j+1);
        ENDDO
!
! QR factorization
!
        QR = V
         CALL DGEQRF( N+1, N+1, QR, N+1, TAU, WORK, N+1, INFO)
        IF (INFO .NE. 0) STOP 'ERROR in DGEQRF'
! getting R from QR
        R = 0_pr
        DO I=1,N+1
            R(I,I:N+1) = QR(I,I:N+1)
        ENDDO
! getting Q from QR
        Q = 0_pr
        DO I=1,N+1
        Q(I,I) = 1_pr
        ENDDO
         CALL DORMQR('L','N', N+1, N+1, N+1, QR, N+1, TAU, Q, N+1, WORK, N+1,INFO)
        IF (INFO .NE. 0) STOP 'ERROR in DORMQR'
!
! solving linear system to get coefficients
!
        b = MATMUL(TRANSPOSE(Q),y)
!
! PLU factorization of matrix
!
        allocate(IPIV(N+1))
! maybe not needed because R is triangular
         CALL DGETRF(N+1, N+1, R, N+1, IPIV, INFO)
        IF (INFO .NE. 0) STOP 'ERROR in DGETRF'
!

! Solving the system factorized
!
        p = b
         CALL DGETRS('N' , N+1, 1, R, N+1, IPIV, p, N+1, INFO)
        deallocate(IPIV)
!
      END SUBROUTINE polyfit
!
!
!
!
!
!
!
!===============================================================================
      SUBROUTINE corner_A()
!===============================================================================
!
! It is necessary to evaluate the properties of the point A on the Pmax isobar
! line at the intersection between the LL and the LH domains;
!
!
!
! Declaration of variables
        INTEGER            :: N, flag
        REAL(pr) :: p_corner_A, T_guess, v_guess, u_corner_A, res
!
!
        p_corner_A   = p_max    ! MPa
        u_corner_A   = e_cr  ! kJ/kg
        T_guess      = 319_pr          ! K ???
        v_guess      = 1.9e-3_pr       ! m3/kg
!
!
        call New_Rap2D(2, T_corner_A, v_corner_A, res, N, flag, p_corner_A, u_corner_A, T_guess, v_guess)
        call New_Rap2D(2, T_corner_A, v_corner_A, res, N, flag, p_corner_A, u_corner_A, T_corner_A, v_corner_A)
        call New_Rap2D(2, T_corner_A, v_corner_A, res, N, flag, p_corner_A, u_corner_A, T_corner_A, v_corner_A)
!      
      END SUBROUTINE corner_A
!
!
!
!===============================================================================
      SUBROUTINE corner_B()
!===============================================================================
!
! It is necessary to evaluate the properties of the point B on the Pmax isobar
! line at the intersection between the LH and HT domains;
!
! Declaration of variables
        INTEGER            :: N, flag
        REAL(pr)   :: p_corner_B, T_guess, v_guess, u_corner_B, res
!
!
        p_corner_B   = p_max        ! Pa
        u_corner_B   = e_umax       ! J/kg
        T_guess      = 336.94_pr    ! K
        v_guess      = 3.69e-3_pr   ! m3/kg
!
      CALL New_Rap2D(2, T_corner_B, v_corner_B, res, N, flag, p_corner_B,u_corner_B, T_guess, v_guess)
      CALL New_Rap2D(2, T_corner_B, v_corner_B, res, N, flag, p_corner_B,u_corner_B, T_corner_B, v_corner_B)
      CALL New_Rap2D(2, T_corner_B, v_corner_B, res, N, flag, p_corner_B,u_corner_B, T_corner_B, v_corner_B)
! 
!
     END SUBROUTINE corner_B
!
!
!
!===============================================================================
      SUBROUTINE corner_C()
!===============================================================================
!
! It is necessary to evaluate the properties on the Pmin isobar line at the
! intersection between R and HT domains;
!
! Declaration of variables
        INTEGER   :: N, flag
        REAL(pr)   :: p_corner_C, T_guess, v_guess, u_corner_C, res
!
!
        p_corner_C   = p_tri  ! Pa
        u_corner_C   = e_umax  ! J/kg
        v_guess      = 7.6597e-2_pr  ! K
        T_guess      = 225.553  ! m3/kg
!
      CALL New_Rap2D(2, T_corner_C, v_corner_C, res, N, flag, p_corner_C,u_corner_C, T_guess, v_guess)
      CALL New_Rap2D(2, T_corner_C, v_corner_C, res, N, flag, p_corner_C,u_corner_C, T_corner_C, v_corner_C)
      CALL New_Rap2D(2, T_corner_C, v_corner_C, res, N, flag, p_corner_C,u_corner_C, T_corner_C, v_corner_C)



      END SUBROUTINE corner_C
!
!
!
!
!===============================================================================
      SUBROUTINE corner_D()
!===============================================================================
!
!         USE saturation, ONLY: saturation_curve
         REAL(8)   :: p_corner_D, v_corner_D, T_corner_D, u_corner_D
         REAL(8)   :: press,delta_p, e_l,e_v,v_v,v_l,temp,qual
         integer   :: i,j
!
         p_corner_D = P_tri
         u_corner_D = e_cr 
!
!       CALL saturation_curve()
        press = p_corner_D
        delta_p = saturP(2) - saturP(1)
        i = INT((press - saturP(1))/delta_p) + 1
!
        print*, press,delta_p,i,p_corner_D,u_corner_D 
        e_l = 0_pr
        e_v = 0_pr
        v_l = 0_pr
        v_v = 0_pr
        temp = 0_pr
!
        DO j = 1, ord_spline + 1
        e_l = e_l + uL_psat_spline(ord_spline+2-j,i) * press**(j-1)
        e_v = e_v + uV_psat_spline(ord_spline+2-j,i) * press**(j-1)
        v_l = v_l + vL_psat_spline(ord_spline+2-j,i) * press**(j-1)
        v_v = v_v + vV_psat_spline(ord_spline+2-j,i) * press**(j-1)
        temp = temp + Tsat_psat_spline(ord_spline+2-j,i) * press**(j-1)
!
        ENDDO
        print*, e_v,e_l,v_v,v_l
        qual = (u_corner_D - e_l)/(e_v - e_l)
        v_corner_D = qual*v_v + (1_pr-qual)*v_l
        print*, 'v',v_corner_D, 'x',qual,temp 
!
      END SUBROUTINE corner_D
!
!
!
!===============================================================================
      SUBROUTINE corner_E()
!===============================================================================
!
!         USE saturation, ONLY: saturation_curve
         REAL(8)   :: p_corner_E, v_corner_E, T_corner_E, u_corner_E
         REAL(8)   :: pguess,vlguess,vvguess,Tguess,in_2,res
         REAL(8)   :: Tsat,vv,vl
         integer   :: Niter, flag
         REAL(8)   :: press,delta_p, e_l,e_v,v_v,v_l,temp,qual
         integer   :: i,j
!
         pguess = 6d6
         u_corner_E = e_tri_R
         vlguess    = 0.9d-3
         vvguess    = 1.0475d-2
         Tguess     = 255d0 


         CALL New_Rap3D(4,vl,vv,Tsat,&
              & res, Niter, flag, u_corner_E,in_2,&
              & vlguess,vvguess,Tguess) 
     
      print*,'RES', res
      print*, 'vv',vv,'vl',vl,'T',Tsat
      print*, 'vvguess',vvguess,'vlguess',vlguess,'Tguess',Tguess

!
       CALL pressure(Tsat,vv,press)
!       CALL saturation_curve()
      print*, 'pressure', press
        delta_p = saturP(2) - saturP(1)
        i = INT((press - saturP(1))/delta_p) + 1
!
!        print*, press,delta_p,i,p_corner_D,u_corner_D
        e_l = 0_pr
        e_v = 0_pr
        v_l = 0_pr
        v_v = 0_pr
        temp = 0_pr
!
        DO j = 1, ord_spline + 1
        e_l = e_l + uL_psat_spline(ord_spline+2-j,i) * press**(j-1)
        e_v = e_v + uV_psat_spline(ord_spline+2-j,i) * press**(j-1)
        v_l = v_l + vL_psat_spline(ord_spline+2-j,i) * press**(j-1)
        v_v = v_v + vV_psat_spline(ord_spline+2-j,i) * press**(j-1)
        temp = temp + Tsat_psat_spline(ord_spline+2-j,i) * press**(j-1)
!
        ENDDO
        print*,'spline', e_v,e_l,v_v,v_l,temp
!
      END SUBROUTINE corner_E

!
!===============================================================================
      SUBROUTINE Max_point(u_max, v_umax, T_umax, p_umax, c_umax)
!===============================================================================
! Compute the maximum internal enregy at the saturation line. The solutions
! are very sensible to the initial guesses, so the domaine of the computaion is
! limited from 236K to 268K at saturation line. Once the maximum is found, it is
! considered as a constant point as same as critial point and triple point.
!-------------------------------------------------------------------------------
! 
! Global variables
        REAL(pr), INTENT(OUT)  :: T_umax, u_max, v_umax, p_umax, c_umax
!
! Local variables
        INTEGER  :: N_iter, Niter, exitflag
        REAL(pr)  :: TOL, tau, deltaK, T_right, T_left, &
        &   s, T_c, T_d, p, v, u_c, u_d, u_left, u_right, &
        &   guess_1_c, guess_2_c, v_l_c, v_v_c, resnorm_c, resnorm_d, in_2,&
        &   guess_1_d, guess_2_d, press, sound, energy, v_left,&
        &   v_right, v_l_d, v_v_d
!
!----------------------------------------------------------------------------
! EVALUATION OF THE MAXIMUM POINT USING GSS (Golden Section Search)
! ALGORITHM
!----------------------------------------------------------------------------
!
        tau = 2_pr/(1_pr + 5_pr**(0.5_pr)) ! Golden ratio for the GSS
!        
!        guess_1_c = 1_pr / rho_cr - tau*(1_pr/rho_cr - 1_pr/rho_tri_L)
!        guess_2_c = 1_pr / rho_cr + tau*(1_pr/rho_tri_R - 1_pr/rho_cr)
!        guess_1_d = 1_pr / rho_tri_L + tau*(1_pr/rho_cr - 1_pr/rho_tri_L)
!        guess_1_d = 1_pr / rho_tri_L 
!        guess_2_d = 1_pr / rho_cr + 0.2_pr*(1_pr/rho_tri_R - 1_pr/rho_cr)
!
!        
! 
        in_2 = 0_pr
!
!         
!
        TOL = 1e-12_pr                 ! Tolerance
        N_iter = 0
!
        deltaK  = 1_pr               ! for the definition of the initial extremes
        T_right = 236_pr+ deltaK
        T_left  = 268_pr - deltaK
!
        guess_1_c = 1_pr/1105.12_pr
        guess_2_c = 1_pr/28.935_pr
        guess_1_d = 1_pr/957.04_pr
        guess_2_d = 1_pr/82.965_pr
!
        DO WHILE (abs(T_right - T_left) .GT. TOL)
!         do while (N_iter < 10) 
          N_iter = N_iter + 1
          s      = tau * (T_right - T_left)
! Two additional points are necessary
          T_c = T_left + s
          T_d = T_right - s
! Evaluation of internal energy in c and d
! Iterative routine in the iterative routine, cost!!
!
        CALL New_Rap2D(1, v_l_c, v_v_c, resnorm_c, Niter, exitflag, T_c, in_2, guess_1_c, guess_2_c)
        CALL inter_energy(T_c,v_v_c,energy)
!       
          u_c = energy
!print*, "temp",T_c,"vl_c",v_l_c,"vv_c",v_v_c,"res_c",resnorm_c,"N",Niter                     
!
        CALL New_Rap2D(1, v_l_d, v_v_d, resnorm_d, Niter, exitflag, T_d, in_2, guess_1_d, guess_2_d)
        CALL inter_energy(T_d,v_v_d,energy)
!print*, "temp",T_d,"vl_d",v_l_d,"vv_d",v_v_d,"res_d",resnorm_d,"N",Niter                      
!       
          u_d = energy
!
          IF (u_d .LT. u_c) THEN
              T_left  = T_d
              v_left = v_v_d
              T_right = T_right
              guess_1_d = v_l_c
              guess_2_d = v_left             
!              guess_1_c = v_l_d + tau*(1_pr/rho_tri_L - v_l_d)
!              guess_2_c = v_v_d + tau*(1_pr/rho_tri_R - v_v_d)
          ELSE
              T_left  = T_left
              T_right = T_c
              v_right = v_v_c
                guess_1_c = v_l_d
                guess_2_c = v_right
!              guess_1_d = v_l_c + tau*(1_pr/rho_cr - v_l_c)
!              guess_2_d = v_v_c + tau*(1_pr/rho_cr - v_v_c )
          ENDIF
!print*, N_iter,"res",resnorm_c,resnorm_d, abs(T_right-T_left)/TOL
!print*, "Tleft",T_left,"rho_left",1_pr/v_left
!print*, "Tright",T_right,"rho_right",1_pr/v_right
!print*,"-----------------------------------"
!print*," "
        ENDDO   
! once that the loop ends, the range is sufficiently small; so it is possible to evaluate 
! in which point the maximum is and what is the value of the maximum itself
!
     CALL inter_energy(T_left,v_left,energy)
        u_left = energy
!
     CALL inter_energy(T_right,v_right,energy)
        u_right = energy
!
        IF (u_left .LT. u_right) THEN
            T_umax = T_right      ! K
            u_max  = u_right      ! J/kg
            v_umax = v_right
        ELSE
            T_umax = T_left       ! K
            u_max  = u_left       ! J/kg
            v_umax = v_left
        ENDIF
!
!
        print*, v_l_c
        print*, v_l_d
        CALL pressure(T_umax, v_umax, press)
        p_umax = press   ! Pa
!
        CALL sound_speed(T_umax, v_umax, sound)
        c_umax = sound
!
      END SUBROUTINE Max_point
!
!
END MODULE grid_functions
