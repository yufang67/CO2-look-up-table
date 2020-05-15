! ===================================================================
!
!                      Saturation_Curve 
!
! ===================================================================
! This subroutine create a grid in the 2-phase region by discretizing P_sat
! from P_tri to P_crit to obtain a table of the saturation properties:
!
!                 P_sat, V_v, V_l, e_v, e_l, T_sat 
!
! The spline coefficients are also computed, so for given values of 
! saturation pressure, coefficients of spline are ble to express 
!                   u(p_sat), v(p_sat) , T(p_sat)
!
! ===================================================================
MODULE saturation
!
     USE def_constants, ONLY: pr, ord_spline, P_cr,P_tri, NNN_TP, NNN_sat_TP,&
&                             rho_tri_L, rho_tri_R, T_tri,e_cr,NNN_TP2, NNN_sat_TP2,&
&                             rho_cr, T_cr,e_tri_L,e_tri_R
     USE def_variables, ONLY: saturP, saturP_sat, vL_psat_spline, vV_psat_spline,&
&                             uL_psat_spline, uV_psat_spline, Tsat_psat_spline,  &
&                             saturP2, saturP_sat2, vL_psat_spline2, vV_psat_spline2,&
&                             uL_psat_spline2, uV_psat_spline2, Tsat_psat_spline2
     USE non_linear_solvers, ONLY: New_Rap3D
     USE properties, ONLY: inter_energy
     USE grid_functions, ONLY: polyfit

     PRIVATE
     PUBLIC  saturation_curve, sat_curve
     CONTAINS     
!========================================================================
!
     SUBROUTINE saturation_curve()
!
!========================================================================       
        IMPLICIT NONE        
!
        INTEGER  :: i, j, exitflag,Niter
!
        REAL(pr)  :: delta, delta_sat,v_v,v_l,Tsat,u_l,u_v,&
 &                   resnorm, guess1,guess2,guess3,in_2
!
        REAL(pr), DIMENSION (NNN_sat_TP) :: vL, vV, uL, uV,saturT
!
!
!
        delta      = (P_cr - P_tri) / (NNN_TP-1)
        saturP     = P_tri + (/(i*delta, i=0, NNN_TP-1)/)
        delta_sat  = (P_cr - P_tri) / (NNN_sat_TP-1)
        saturP_sat = P_tri + (/(i*delta_sat, i=0,NNN_sat_TP-1)/)
!
! points on saturation curve 
!
        guess1 = 1_pr/rho_tri_L
        guess2 = 1_pr/rho_tri_R
        guess3 = T_tri
!
        DO i = 1, NNN_sat_TP-1

         CALL New_Rap3D(3,v_l , v_v, Tsat, &
     &   resnorm, Niter, exitflag, saturP_sat(i),in_2, guess1, guess2,guess3)
        
        IF (resnorm > 1e-7_pr) THEN
        print*, "saturation curve", resnorm, i
        STOP
        ENDIF
!
!        print*, v_l,v_v
!        print*, u_l,u_v
!        print*," "
         CALL inter_energy(Tsat,v_v,u_v)
         CALL inter_energy(Tsat,v_l,u_l)
           vL(i)     = v_l
           vV(i)     = v_v
           uL(i)     = u_l
           uV(i)     = u_v
           saturT(i) = Tsat
!
           guess1 = v_l
           guess2 = v_v
           guess3 = Tsat
        ENDDO
!
        vL(NNN_sat_TP)     = 1_pr/rho_cr
        vV(NNN_sat_TP)     = 1_pr/rho_cr
        uL(NNN_sat_TP)     = e_cr
        uV(NNN_sat_TP)     = e_cr
        saturT(NNN_sat_TP) = T_cr
! intervals on saturation curve
!
        DO i = 1, NNN_TP-1
           j = 3 * (i-1) + 1
             CALL polyfit(vL_psat_spline(:,i), saturP_sat (j:j+3), vL(j:j+3),ord_spline)
             CALL polyfit(vV_psat_spline(:,i), saturP_sat (j:j+3), vV(j:j+3),ord_spline)
             CALL polyfit(uL_psat_spline(:,i), saturP_sat (j:j+3), uL(j:j+3),ord_spline)
             CALL polyfit(Tsat_psat_spline(:,i), saturP_sat(j:j+3),saturT(j:j+3), ord_spline)
             CALL polyfit(uV_psat_spline(:,i), saturP_sat (j:j+3), uV(j:j+3),ord_spline)
        ENDDO
!print*,'-------------------------------------------------------------------------'
!print*,'           construction for saturation curve with 3e order spline        '
!print*,'-------------------------------------------------------------------------'
!
!
     END SUBROUTINE saturation_curve
!
!===================================================================================
!
     SUBROUTINE sat_curve()
!
!===================================================================================       
        IMPLICIT NONE
!
        INTEGER  :: i, j, exitflag,Niter
!
        REAL(pr)  :: delta, delta_sat,v_v,v_l,Tsat,u_l,u_v,&
 &                   resnorm, guess1,guess2,guess3,in_2
!
        REAL(pr), DIMENSION (NNN_sat_TP2) :: vL, vV, uL, uV,saturT
!
!
!
        delta      = (P_cr - P_tri) / (NNN_TP2-1)
        saturP2     = P_tri + (/(i*delta, i=0, NNN_TP2-1)/)
        delta_sat  = (P_cr - P_tri) / (NNN_sat_TP2-1)
        saturP_sat2 = P_tri + (/(i*delta_sat, i=0,NNN_sat_TP2-1)/)
!
! points on saturation curve 
!
        guess1 = 1_pr/rho_tri_L
        guess2 = 1_pr/rho_tri_R
        guess3 = T_tri
!
        DO i = 2, NNN_sat_TP2-1

         CALL New_Rap3D(3,v_l , v_v, Tsat, &
     &   resnorm, Niter, exitflag, saturP_sat2(i),in_2, guess1, guess2,guess3)

        IF (resnorm > 1e-12_pr) THEN
        print*, "sat  curve", resnorm, i, Niter
        STOP
        ENDIF
!
!        print*, v_l,v_v
!        print*, u_l,u_v
!        print*," "
         CALL inter_energy(Tsat,v_v,u_v)
         CALL inter_energy(Tsat,v_l,u_l)
           vL(i)     = v_l
           vV(i)     = v_v
           uL(i)     = u_l
           uV(i)     = u_v
           saturT(i) = Tsat
!
           guess1 = v_l
           guess2 = v_v
           guess3 = Tsat
        ENDDO
!
        vL(NNN_sat_TP2)     = 1_pr/rho_cr
        vV(NNN_sat_TP2)     = 1_pr/rho_cr
        uL(NNN_sat_TP2)     = e_cr
        uV(NNN_sat_TP2)     = e_cr
        saturT(NNN_sat_TP2) = T_cr
!
        vL(1)     = 1_pr/rho_tri_L
        vV(1)     = 1_pr/rho_tri_R
        uL(1)     = e_tri_L
        uV(1)     = e_tri_R
        saturT(1) = T_tri

! intervals on saturation curve
!
        DO i = 1, NNN_TP2-1
           j = 5 * (i-1) + 1
             CALL polyfit(vL_psat_spline2(:,i),   saturP_sat2 (j:j+5), vL(j:j+5),     5)
             CALL polyfit(vV_psat_spline2(:,i),   saturP_sat2 (j:j+5), vV(j:j+5),     5)
             CALL polyfit(uL_psat_spline2(:,i),   saturP_sat2 (j:j+5), uL(j:j+5),     5)
             CALL polyfit(Tsat_psat_spline2(:,i), saturP_sat2 (j:j+5), saturT(j:j+5), 5)
             CALL polyfit(uV_psat_spline2(:,i),   saturP_sat2 (j:j+5), uV(j:j+5),     5)
        ENDDO
!print*,'-------------------------------------------------------------------------'
!print*,'construction for saturation curve with 5e order spline                   '
!print*,'-------------------------------------------------------------------------'
!
!
     END SUBROUTINE sat_curve
!
!=================================================================================
!     
!     SUBROUTINE satprop(mode, psat, Tsat, vvsat, vlsat, uvsat, ulsat)
!
!=================================================================================
!      IMPLICIT NONE
!!
!      INTEGER  :: i, j, j_sat
!!
!      REAL(pr), INTENT(in)  :: psat
!      INTEGER , INTENT(in)  :: mode
!!
!      REAL(pr), INTENT(out) :: Tsat, vvsat,vlsat,uvsat,ulsat
!!
!      REAL(pr) :: delta, pp, ratio, c_out, a_out, x_out, v_in, du_dp_x, dv_dp_x, temp
!      REAL(pr) :: duL_dp, duV_dp, dvL_dp, dvV_dp 
!      REAL(pr) :: vL, uL, vV, uV
!
!
!use 3e order spline
!
!      IF (mode == 3) THEN 
!!
!       delta = saturP(2) - saturP(1)
!       j_sat = INT((psat  - saturP(1))/delta) + 1
!!
!!##computing derivatives
!!
!       duL_dp = 0_pr;  duV_dp = 0_pr;  dvL_dp = 0_pr;  dvV_dp = 0_pr
!       DO i = 1, ord_spline
!          pp      = psat**(ord_spline - i)
!          duL_dp  = duL_dp + (ord_spline+1-i) * uL_psat_spline(i,j_sat) *pp
!          duV_dp  = duV_dp + (ord_spline+1-i) * uV_psat_spline(i,j_sat) *pp
!          dvL_dp  = dvL_dp + (ord_spline+1-i) * vL_psat_spline(i,j_sat) *pp
!          dvV_dp  = dvV_dp + (ord_spline+1-i) * vV_psat_spline(i,j_sat) *pp
!       ENDDO
!
!##computing saturation quantities
!
!       vL = 0_pr; uL = 0_pr;  vV = 0_pr;  uV = 0_pr; temp = 0_pr
!       DO i = 1, ord_spline+1
!          pp = psat**(i-1)
!          vL   = vL  + vL_psat_spline  (ord_spline+2-i, j_sat) * pp
!          uL   = uL  + uL_psat_spline  (ord_spline+2-i, j_sat) * pp
!          vV   = vV  + vV_psat_spline  (ord_spline+2-i, j_sat) * pp
!          uV   = uV  + uV_psat_spline  (ord_spline+2-i, j_sat) * pp
!          temp = temp+ Tsat_psat_spline(ord_spline+2-i, j_sat) * pp 
!       ENDDO
!
!       ratio = ((uV - uL)/(vV - vL))   ! (J/kg)/(m3/kg)
!       du_dp_x = x_out * duV_dp + (1_pr - x_out) * duL_dp
!       dv_dp_x = x_out * dvV_dp + (1_pr - x_out) * dvL_dp
!
!       c_out = SQRT((psat  + ratio)/(du_dp_x - ratio * dv_dp_x)) * v_in ! (m/s)
!
!       a_out = x_out*vV/v_in
!
!##for output saturation quantities
!
!       Tsat  = temp 
!       vvsat = vV
!       vlsat = vL
!       uvsat = uV
!       ulsat = uL
!
!use 5e order spline
!
!      ELSEIF (mode == 5) THEN
!
!       delta = saturP2(2) - saturP2(1)
!       j_sat = INT((psat  - saturP2(1))/delta) + 1
!
!##computing derivatives
!
!       duL_dp = 0_pr;  duV_dp = 0_pr;  dvL_dp = 0_pr;  dvV_dp = 0_pr; temp=0_pr
!       DO i = 1, 5
!          pp      = psat**(5 - i)
!          duL_dp  = duL_dp + (5+1-i) * uL_psat_spline2(i,j_sat) *pp
!          duV_dp  = duV_dp + (5+1-i) * uV_psat_spline2(i,j_sat) *pp
!          dvL_dp  = dvL_dp + (5+1-i) * vL_psat_spline2(i,j_sat) *pp
!          dvV_dp  = dvV_dp + (5+1-i) * vV_psat_spline2(i,j_sat) *pp
!       ENDDO
!
!##computing saturation quantities
!
!       vL = 0_pr; uL = 0_pr;  vV = 0_pr;  uV = 0_pr
!       DO i = 1, 6
!          pp = psat**(i-1)
!          vL   = vL  + vL_psat_spline2  (5+2-i, j_sat) * pp
!          uL   = uL  + uL_psat_spline2  (5+2-i, j_sat) * pp
!          vV   = vV  + vV_psat_spline2  (5+2-i, j_sat) * pp
!          uV   = uV  + uV_psat_spline2  (5+2-i, j_sat) * pp
!          temp = temp+ Tsat_psat_spline2(5+2-i, j_sat) * pp
!       ENDDO
!
!       ratio = ((uV - uL)/(vV - vL))   ! (J/kg)/(m3/kg)
!       du_dp_x = x_out * duV_dp + (1_pr - x_out) * duL_dp
!       dv_dp_x = x_out * dvV_dp + (1_pr - x_out) * dvL_dp
!
!       c_out = SQRT((psat  + ratio)/(du_dp_x - ratio * dv_dp_x)) * v_in ! (m/s)
!
!       a_out = x_out*vV/v_in
!
!##for output saturation quantities
!
!       Tsat  = temp
!       vvsat = vV
!       vlsat = vL
!       uvsat = uV
!       ulsat = uL
!
!      ENDIF
!     END SUBROUTINE satprop
!
!
!
END MODULE saturation
