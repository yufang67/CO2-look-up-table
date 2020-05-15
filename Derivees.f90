MODULE derivees
!
!
     USE def_variables     
     USE def_constants
     USE properties, ONLY: dpdu_v, dpdr_u 
!     USE interp_der
!      
     IMPLICIT NONE
!
     PRIVATE
     PUBLIC  CO2DER
!
     CONTAINS
!
!===============================================================================
!
SUBROUTINE CO2DER(dp_dv_u, dp_du_v , u_in,v_in, T_in, p_TP,res2,res3,res4) 
!
!===============================================================================
!
!                       CO2 DERIVATIVES
!
! Input: u_in (the specific internal energy)
!        v_in (the specific volume)
!        T_in (Temperature)
!        P_TP (saturation pressure for two-phase) 
! Output: dp_dv_u, dp_du_v 
!
!===============================================================================
REAL(pr),INTENT(OUT)  :: dp_dv_u, dp_du_v
REAL(pr),INTENT(IN)   :: T_in,p_TP
!
!Local Variables
INTEGER :: i, j, flag_TP, j_sat, i_R, i_L, flag_loca,Niter,exitflag
REAL(pr) :: v_min, v_max, v_sat, v_sat_log, delta,&
&            qual,press,temp,sound,T_guess,out2
REAL(pr) :: duL_dp,  duV_dp,  dvL_dp,  dvV_dp
REAL(pr) :: vL, uL, vV, uV, pp, ratio
REAL(pr) :: du_dp_x, dv_dp_x, du_dr_p
REAL(pr) :: x_check, v_check, f_out
REAL(pr) :: gamma_pg, e_pg, x_out,dT_dv_u, dT_du_vi,deriv
REAL(pr) :: res2,res3,res4,u_in,v_in
!
!
!
flag_TP = 1
flag_loca = 0
f_out   = 0_pr
!
!check input
IF ( (u_in /= u_in) .OR. (v_in /= v_in) .OR. (v_in < 0.0) ) THEN
!
!   write(*,*) 'INPUT ERRORS FOR CO2 TABLE from ', mess2,' CALLED by ',mess1
!   write(*,*) 'u_in=',u_in,' v_in=',v_in, ' at node, x=',xloc,' y=',yloc
    write(*,*) 'INPUT ERRORS FOR CO2 TABLE ','u_in=',u_in,' v_in=',v_in
    STOP
ENDIF
!
IF (u_in .GT. e_umax) THEN
! 
          IF (u_in .GT. u_end) THEN
!             STOP '** Out of range. Too high specific internal energy'
              write(*,*) 'Too high eint for CO2_TABLE ','u_in=',u_in,' v_in=',v_in
              STOP
          ENDIF
!
! To evaluate the interval on the vertical axis, BC of HT
          delta = y_mesh_HT(2) - y_mesh_HT(1)
          i     = INT((u_in    - y_mesh_HT(1))/delta) + 1 !location indice 
!
          v_min = 0_pr
          v_max = 0_pr
          DO j = 1, ord_spline + 1
             v_min = v_min + spline_left_HT (ord_spline+2-j,i) * u_in**(j-1)
             v_max = v_max + spline_right_HT(ord_spline+2-j,i) * u_in**(j-1)
          ENDDO
!
          IF (v_in .LT. v_min) THEN
!             STOP '** Out of range. Pressure higher than 50 MPa in HT'
             write(*,*) 'Pressure higher than 50 MPa in HT','u_in=',u_in,'v_in=',v_in
             STOP
          ELSEIF (v_in .GT. v_max) THEN
! In LP
             flag_loca = 6
          ELSE
! In HT
             flag_loca = 4
          END IF
ELSEIF ((u_in .LE. e_umax) .AND. (v_in .LE. v_umax)) THEN 
! Left part LL or LH or TP
       IF (u_in .LE. e_cr) THEN
! In LL
           IF (u_in .LT. e_tri_L) THEN
!           STOP '** Out of range. Too low specific internal energy value in LL'
             write(*,*) 'Too low eint in LL ','u_in=',u_in,' v_in=',v_in
             STOP
           ENDIF  
!
! To evaluate the interval on the vertical axis, BC of LL
          delta  = y_mesh_LL(2) - y_mesh_LL(1)
          i      = INT((u_in    - y_mesh_LL(1))/delta) + 1
!
          v_min = 0_pr;          v_sat = 0_pr
          DO j = 1, ord_spline + 1
             v_min = v_min + spline_pmax_LL(ord_spline+2-j,i) * u_in**(j-1)
             v_sat = v_sat + spline_Lsat_LL(ord_spline+2-j,i) * u_in**(j-1)
!             v_sat_log = v_sat_log + spline_Lsat_LL(ord_spline+2-j,i) *u_in**(j-1)
!             v_sat     = 10_pr ** v_sat_log
!print*, spline_Lsat_LL(ord_spline+2-j,i)
          ENDDO
!          
          IF (v_in .LT. v_min) THEN
!             STOP '** Out of range. Pressure higher than 50 MPa in LL'
             write(*,*) 'P > 50 MPa in LL','u_in=',u_in,'v_in=',v_in
             STOP
          ELSEIF (v_in .GT. v_sat) THEN
! Region Two-phase (TP)
             x_check = (u_in - e_tri_L)/(e_tri_R - e_tri_L)
             v_check = x_check*v_tri_R + (1_pr-x_check)*v_tri_L
             IF (v_in .GT. v_check) THEN
! In Solid, we put it back to TP
                v_in = v_check - 1e-5_pr
                f_out = 10
             ENDIF
             flag_loca = 5  
          ELSE
! Region Left Low (LL) 
             flag_loca = 1 
          ENDIF
! End in LL and begin in LH
       ELSE
!
! To evaluate the interval on the vertical axis, BC of LH
          delta = y_mesh_LH(2)- y_mesh_LH(1)
          i     = INT((u_in - y_mesh_LH(1))/delta) + 1
!
          v_min = 0_pr;          v_sat = 0_pr
          v_sat_log = 0_pr
          DO j = 1, ord_spline + 1
             v_min     = v_min     + spline_pmax_LH(ord_spline+2-j,i) *u_in**(j-1)
             v_sat     = v_sat     + spline_Lsat_LH(ord_spline+2-j,i) *u_in**(j-1)
!             v_sat_log = v_sat_log + spline_Lsat_LH(ord_spline+2-j,i) *u_in**(j-1)
!             v_sat     = 10_pr ** v_sat_log
          ENDDO
!
          IF (v_in .LT. v_min) THEN
!             STOP '** Out of range. Pressure higher than 10 MPa in LH'
             write(*,*) 'P > 50 MPa in LH ','u_in=',u_in,'v_in=',v_in
             STOP
          ELSEIF (v_in .GT. v_sat) THEN
! Region Two-phase (TPH + TPM)   
             x_check = (u_in - e_tri_L)/(e_tri_R - e_tri_L)
             v_check = x_check*v_tri_R + (1_pr-x_check)*v_tri_L
             IF (v_in .GT. v_check) THEN
! In Solid, we put it back to TP
                v_in = v_check - 1e-5_pr
                f_out = 10
             ENDIF
             flag_loca = 5
          ELSE
! Region Left High (LH)
             flag_loca = 2
          ENDIF
       ENDIF
! Right part R or TP or LP or solid(not considered yet) 
ELSE
! 
       IF (u_in .LT. e_tri_R) THEN
          x_check = (u_in - e_tri_L)/(e_tri_R - e_tri_L)
          v_check = x_check*v_tri_R + (1_pr-x_check)*v_tri_L
          IF (v_in .GT. v_check) THEN
! In Solid, we put it back to TP
             v_in = v_check - 1e-5_pr
             flag_loca = 5
             f_out = 10_pr
          ENDIF
! In TP
          flag_loca = 5
       ELSEIF (u_in .LT. e_tri_L) THEN
          x_check = (v_in - v_tri_L)/(v_tri_R - v_tri_L)
          u_in    = x_check*e_tri_R + (1_pr-x_check)*e_tri_L
          flag_loca = 5
       ELSE      
! BC of R
          delta = y_mesh_R(2) - y_mesh_R(1)
          i     = INT((u_in   - y_mesh_R(1))/delta) + 1
!
          v_sat = 0_pr;       v_max = 0_pr
          DO j = 1, ord_spline + 1
            v_sat = v_sat + spline_Vsat(ord_spline+2-j,i) * u_in**(j-1)
            v_max = v_max + spline_pmin(ord_spline+2-j,i) * u_in**(j-1)
          ENDDO
!
          IF (v_in .GT. v_max) THEN
! In LP
            flag_loca=6
          ELSEIF (v_in .LT. v_sat) THEN
! In TP
            flag_loca = 5
          ELSE
! In Region Right (R)
            flag_loca = 3
          ENDIF
       ENDIF
!
ENDIF
!print*, "location", flag_loca
!
!
!
!SELECT CASE (flag_loca)
!
!CASE( 0 )
!STOP '** Locating the points in the physical domaion failed '
!
!## LL ##
!CASE( 1 )
!
!     delta         = (x_mesh_max - x_mesh_min)/(MMM_LL-1)
!     x_mesh_LL     = x_mesh_min + (/(i*delta, i=0,MMM_LL-1)/)
!
!     CALL Lin_der(dT_dv_u, dT_du_v, dp_dv_u, dp_du_v,&
!&                 NNN_LL,MMM_LL,x_mesh_LL,y_mesh_LL, spline_Lsat_LL, spline_pmax_LL,&
!&                 u_in, v_in, TTT_LL, ppp_LL, vvv_LL)
!
!## LH ##
!CASE( 2 )
!
!     delta         = (x_mesh_max - x_mesh_min)/(MMM_LH-1)
!     x_mesh_LH     = x_mesh_min + (/(i*delta, i=0,MMM_LH-1)/)
!
!     CALL Lin_der_Log10(dT_dv_u, dT_du_v, dp_dv_u, dp_du_v,&
!&                       NNN_LH,MMM_LH,x_mesh_LH,y_mesh_LH, spline_Lsat_LH, spline_pmax_LH,&
!&                       u_in, v_in, TTT_LH, ppp_LH, vvv_LH)
!
!## R ##
!CASE( 3 )
!
!     delta        = (x_mesh_max - x_mesh_min)/(MMM_R-1)
!     x_mesh_R     = x_mesh_min + (/(i*delta, i=0,MMM_R-1)/)
!
!     CALL Lin_der_Log10(dT_dv_u, dT_du_v, dp_dv_u, dp_du_v,&
!&                       NNN_R, MMM_R,x_mesh_R,y_mesh_R, spline_pmin, spline_Vsat,&
!                        u_in, v_in, TTT_R, ppp_R, vvv_R)
!     CALL Lin_der(dT_dv_u, dT_du_v, dp_dv_u, dp_du_v,&
!&                 NNN_R, MMM_R,x_mesh_R,y_mesh_R, spline_pmin, spline_Vsat,&
!                  u_in, v_in, TTT_R, ppp_R, vvv_R)
!
!## HT ## 
!CASE( 4 )
!
!     delta       = (x_mesh_max - x_mesh_min)/(MMM_HT-1)
!     x_mesh_HT   =  x_mesh_min + (/(i*delta, i=0,MMM_HT-1)/)
!
!!     CALL Lin_der(dT_dv_u, dT_du_v, dp_dv_u, dp_du_v,&
!!&                 NNN_HT,MMM_HT,x_mesh_HT,y_mesh_HT, spline_right_HT,spline_left_HT,&
!!&                 u_in, v_in, TTT_HT, ppp_HT, vvv_HT)
!
!     CALL Lin_der_Log10(dT_dv_u, dT_du_v, dp_dv_u, dp_du_v,&
!                        NNN_HT,MMM_HT,x_mesh_HT,y_mesh_HT, spline_right_HT,spline_left_HT,&
!                        u_in, v_in, TTT_HT, ppp_HT, vvv_HT)
!
!## TP ##
IF (flag_loca == 5) THEN
!CASE( 5 )
!! derivatives computed by spline with saturation pressure p_TP
!
       delta = saturP(2) - saturP(1)
       j_sat = INT((p_TP  - saturP(1))/delta) + 1
!
       duL_dp = 0_pr;  duV_dp = 0_pr;  dvL_dp = 0_pr;  dvV_dp = 0_pr
       DO i = 1, ord_spline
          pp      = p_TP**(ord_spline - i)
          duL_dp  = duL_dp + (ord_spline+1-i) * uL_psat_spline(i,j_sat) *pp
          duV_dp  = duV_dp + (ord_spline+1-i) * uV_psat_spline(i,j_sat) *pp
          dvL_dp  = dvL_dp + (ord_spline+1-i) * vL_psat_spline(i,j_sat) *pp
          dvV_dp  = dvV_dp + (ord_spline+1-i) * vV_psat_spline(i,j_sat) *pp
       ENDDO
!
!
       vL = 0_pr; uL = 0_pr;  vV = 0_pr;  uV = 0_pr
       DO i = 1, ord_spline+1
          pp = p_TP**(i-1)
          vL = vL + vL_psat_spline(ord_spline+2-i, j_sat) * pp
          uL = uL + uL_psat_spline(ord_spline+2-i, j_sat) * pp
          vV = vV + vV_psat_spline(ord_spline+2-i, j_sat) * pp
          uV = uV + uV_psat_spline(ord_spline+2-i, j_sat) * pp
       ENDDO
!
       x_out   = (v_in - vL) / (vV - vL)
       ratio   = ((uV - uL)/(vV - vL))   ! (J/kg)/(m3/kg)
       du_dp_x = x_out * duV_dp + (1_pr - x_out) * duL_dp
       dv_dp_x = x_out * dvV_dp + (1_pr - x_out) * dvL_dp
       du_dr_p = -v_in*v_in*ratio
!
!
       dp_du_v = 1_pr / (du_dp_x - ratio*dv_dp_x)
       dp_dv_u = du_dr_p / dp_du_v * 1_pr/(v_in*v_in) 
!
!       dT_du_v = 0.0
!       dT_dv_u = 0.0       
!
!       print*, 'du_dr_p',du_dr_p, 'dp_du_v',dp_du_v
       IF ((dp_du_v /= dp_du_v) .OR. (dp_dv_u /=  dp_dv_u) ) THEN
         STOP '**TP Derivatives problem in CO2DER case (5)'
       ENDIF
!
!
!CASE( 6 )
ELSEIF (flag_loca == 6) THEN
!        
        IF (u_in .LT. -123.74e3_pr) THEN
! In Solid phase, but we push it back to LP
         u_in  = -123.74e3_pr + 1e3_pr
         f_out = 2_pr
         Print*,'** out of range (SOLID) in Interp_table case(6)'
        ENDIF
         gamma_pg = 1.313_pr
         e_pg     = 236.0294e3_pr
!
         dp_du_v  = (gamma_pg-1_pr)/v_in
         dp_dv_u  = -(gamma_pg-1_pr)*(u_in+e_pg)/(v_in*v_in)         
!
!         dT_du_v  = (gamma_pg-1_pr)/R
!         dT_dv_u  = 0.0
!
         IF ((dp_du_v /= dp_du_v) .OR. (dp_dv_u /=  dp_dv_u) ) THEN
           STOP '**PG Derivatives problem in CO2DER case (6)'
         ENDIF
!
!END SELECT
!
ELSE
         CALL dpdu_v(T_in, v_in, deriv)
         dp_du_v = deriv
!
         CALL dpdr_u(T_in, v_in, deriv)
         dp_dv_u = - deriv * 1_pr/(v_in*v_in)
ENDIF
END SUBROUTINE CO2DER
!
END MODULE derivees 
