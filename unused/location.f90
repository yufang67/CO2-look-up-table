MODULE location
!
!
     USE def_variables, ONLY: y_mesh_LL,y_mesh_LH,y_mesh_R,y_mesh_HT,&
&                             spline_Lsat_LL, spline_pmax_LL,spline_meta_LL,&
&                             spline_Lsat_LH, spline_pmax_LH,&
&                             spline_pmin, spline_Vsat,&
&                             spline_left_HT, spline_right_HT
     USE def_constants, ONLY: e_tri_L, e_tri_R,v_tri_R,v_tri_L,&
&                             e_cr,e_umax,v_umax,ord_spline, u_end,pr
!      
     IMPLICIT NONE
!
     PRIVATE
     PUBLIC :: phaseloca
!
!
     CONTAINS
!
!===============================================================================
!
SUBROUTINE phaseloca(x_out, a_out, flag, u_in,v_in)
!
!===============================================================================
!
!                       CO2 BILINEAR LOOK-UP TABLES
!
! Input: u_in (the specific internal energy)
!        v_in (the specific volume v_in) 
! Output:  
!         flag
!         a_out
!         x_out (thermodynamic quality for two-phase regions) 
!
!===============================================================================
REAL(8),INTENT(OUT)  :: x_out, a_out
INTEGER,INTENT(OUT)  :: flag 
REAL(8),INTENT(IN)   :: u_in!, v_in
!
!Local Variables
INTEGER :: i, j, j_sat, i_R, i_L, flag_loca,Niter,exitflag,f_out
REAL(8) :: v_min, v_max, v_sat, v_sat_log, delta,&
&            qual,press,temp,sound,T_guess,out2,res
!REAL(8) :: duL_dp,  duV_dp,  dvL_dp,  dvV_dp
!REAL(8) :: vL, uL, vV, uV, pp, ratio, du_dp_x, dv_dp_x
REAL(8) :: x_check, v_check, v_in
flag_loca = 0
x_out   = 1.0d0
a_out   = 1.0
res     = 0.0
!
IF (u_in .GT. e_umax) THEN
! In HT
          IF (u_in .GT. u_end) THEN
             STOP '** Out of range. Too high specific internal energy'
          ENDIF
!
! To evaluate the interval on the vertical axis
          delta = y_mesh_HT(2) - y_mesh_HT(1)
          i     = INT((u_in    - y_mesh_HT(1))/delta) + 1 !location indice 
!
          v_min = 0.0
          v_max = 0.0
          DO j = 1, ord_spline + 1
             v_min = v_min + spline_left_HT (ord_spline+2-j,i) * u_in**(j-1)
             v_max = v_max + spline_right_HT(ord_spline+2-j,i) * u_in**(j-1)
          ENDDO
!
          IF (v_in .LT. v_min) THEN
             print*, 'u=',u_in,'v=',v_in
             print*,'** Out of range. Pressure higher than 10 MPa in HT'
             STOP
          ELSEIF (v_in .GT. v_max) THEN
             PRINT*,'Small pressure region'
          ! small pressure region
          flag_loca = 6
          ELSE
          ! Then we are in Region High Temperature (HT)
          flag_loca = 4
          END IF
ELSEIF ((u_in .LE. e_umax) .AND. (v_in .LE. v_umax)) THEN
! In LL or LH
        IF (u_in .LE. e_cr) THEN
! In LL
                IF (u_in .LT. e_tri_L) THEN
                print*, 'u=',u_in,'v=',v_in
                STOP '** Out of range. Too low specific internal energy value in LL'
                ENDIF
!
! To evaluate the interval on the vertical axis
          delta  = y_mesh_LL(2) - y_mesh_LL(1)
          i      = INT((u_in    - y_mesh_LL(1))/delta) + 1
!
          v_min = 0.0;          v_sat = 0.0
          DO j = 1, ord_spline + 1
             v_min = v_min + spline_pmax_LL(ord_spline+2-j,i) * u_in**(j-1)
             v_sat = v_sat + spline_Lsat_LL(ord_spline+2-j,i) * u_in**(j-1)
!             v_sat_log = v_sat_log + spline_Lsat_LL(ord_spline+2-j,i) *u_in**(j-1)
!             v_sat     = 10_pr ** v_sat_log
!print*, spline_Lsat_LL(ord_spline+2-j,i)
          ENDDO
!
!print*,'============================='
!print*, "v_min", v_min
!print*, "v_sat", v_sat
                IF (v_in .LT. v_min) THEN
                print*, 'u=',u_in,'v=',v_in
                STOP '** Out of range. Pressure higher than 10 MPa in LL'
                ELSEIF (v_in .GT. v_sat) THEN
                 ! Region Two-phase (TP)
                        flag_loca = 5
                ELSE
                 ! Region Left Low (LL) 
                        flag_loca = 1
                        x_out   = 0.0
                        a_out   = 0.0
                ENDIF
! End in LL and begin in LH
        ELSE
!
!  To evaluate the interval on the vertical axis
          delta = y_mesh_LH(2)- y_mesh_LH(1)
          i     = INT((u_in - y_mesh_LH(1))/delta) + 1
!
          v_min = 0.0;          v_sat = 0.0
          v_sat_log = 0.0
          DO j = 1, ord_spline + 1
             v_min     = v_min     + spline_pmax_LH(ord_spline+2-j,i) *u_in**(j-1)
             v_sat_log = v_sat_log + spline_Lsat_LH(ord_spline+2-j,i) *u_in**(j-1)
             v_sat     = 10.0 ** v_sat_log
          ENDDO
!
          IF (v_in .LT. v_min) THEN
             print*, 'u=',u_in,'v=',v_in
             STOP '** Out of range. Pressure higher than 10 MPa in LH'
          ELSEIF (v_in .GT. v_sat) THEN
          ! Region Two-phase (HT)   
                flag_loca = 5
          ELSE
          ! Region Left High (LH)
                flag_loca = 2
          ENDIF
!End in LH
        ENDIF
! Begin in R or p<5bar
ELSE
! It occurs that:  ((u_in .LE. u_max) .AND. (v_in .GT. v_umax))
!
! To evaluate the interval on the vertical axis
     IF (u_in .LT. e_tri_R) THEN
        x_check = (u_in - e_tri_L)/(e_tri_R - e_tri_L)
        v_check = x_check*v_tri_R + (1_pr-x_check)*v_tri_L
        IF (v_in .GT. v_check) THEN
           IF (v_in .GT. v_tri_R) THEN
!              print*, '==>S-V or small pressure region <===,u_in=',u_in,'v_in=',v_in
              flag_loca = 6
!              v_in = v_check
              f_out= 1_pr
           ELSE
!              print*, '====>ATTENDION!!!, S-L-V region<====,u_in=',u_in,'v_in=',v_in
              v_in = v_check-1d-4
              flag_loca = 5
              f_out = 2_pr
           ENDIF
        ELSE
        ! Region two-phase
           flag_loca = 5
        ENDIF
!     IF (u_in .LT. e_tri_R) THEN
!        x_check = (u_in - e_tri_L)/(e_tri_R - e_tri_L)
!        v_check = x_check*v_tri_R + (1.0-x_check)*v_tri_L
!        IF (v_in .GT. v_check) THEN
!       PRINT*, 'Small pressure region'
!        ! Region small pressure
!        flag_loca = 6
!        ENDIF
        ! Region two-phase
!        flag_loca = 5
     ELSE
       delta = y_mesh_R(2) - y_mesh_R(1)
       i     = INT((u_in   - y_mesh_R(1))/delta) + 1
!
       v_sat = 0.0;       v_max = 0.0
       DO j = 1, ord_spline + 1
          v_sat = v_sat + spline_Vsat(ord_spline+2-j,i) * u_in**(j-1)
          v_max = v_max + spline_pmin(ord_spline+2-j,i) * u_in**(j-1)
       ENDDO
!
! To evaluate if it is in the single phase domain or in the two-phase
!
!       IF ((u_in .LT. e_umax) .AND. (u_in .GT. e_tri_R) .AND. &
!&           (v_in .LT.v_sat)) THEN
        ! Region two-phase  
!        flag_loca = 5   
!        
!       ENDIF
!
       IF (v_in .GT. v_max) THEN
!          PRINT*, 'Small pressure region in R'
       ! Small pressure
          flag_loca=6
       ELSEIF (v_in .LT. v_sat) THEN
        ! Region two-phase
          flag_loca = 5
       ELSE
        ! Region Right (R)
          flag_loca = 3
       ENDIF
     ENDIF
!
ENDIF
flag= flag_loca
!print*, "location=", flag_loca
!
END SUBROUTINE phaseloca
!
!
END MODULE location                                                           
