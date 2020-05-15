! ===================================================================
!
!                   Grid_construction_High_Temperature 
!
! ===================================================================
!
! Region High Temperature (HT): 
!
! From the specific internal energy value referred to the maximum point
! on the saturation curve up to T = 500K (internal energy).
! Left boundary: Maximum pressure 50MPa
! Right boundary: Minimum pressure 0.5MPa
! note: discritization log on vvv_HT to minimize errors 
! ===================================================================
SUBROUTINE grid_construction_high_temperature()
!
     USE def_constants
     USE def_variables, ONLY: y_mesh_HT,x_mesh_HT,y_mesh_sat_HT,v_left_HT,v_right_HT,&
&                             spline_left_HT, spline_right_HT,ppp_HT,vvv_HT,ccc_HT,TTT_HT,cpcp_HT
     USE non_linear_solvers
     USE properties
     USE grid_functions
!
INTEGER  ::  i, j, Niter, exitflag
!INTEGER, DIMENSION (NNN_HT,MMM_HT)  ::  
!
REAL(pr)   ::  u_min, res, delta, T_guess, v_guess, p_min,&
&              temp,press,out_2, sound
!
REAL(pr), DIMENSION (NNN_HT)     :: T_pmax, p_HT, T_pmin
REAL(pr), DIMENSION (NNN_sat_HT) :: T_sat_pmax, T_sat_pmin, p_sat_HT
REAL(pr), DIMENSION (NNN_HT,MMM_HT) :: res_pT_HT, vvv_log
!
! data used for bicubic interpolation
!
!REAL(8) :: v1, v2, v3, v4, u1, u2, u3, u4
!REAL(8) :: F1, F2, F3, F4, G1, G2, G3, G4
!REAL(8) :: a1, a2, a3, a4, b1, b2, b3, b4
!REAL(8) :: dF1v, dF2v, dF3v, dF4v, dF1u, dF2u, dF3u, dF4u, dF1vu, dF2vu,
!dF3vu, dF4vu
!REAL(8) :: dF1vv, dF2vv, dF3vv, dF4vv, dF1uu, dF2uu, dF3uu, dF4uu
!REAL(8) :: dF1x, dF2x, dF3x, dF4x, dF1y, dF2y, dF3y, dF4y, dF1xy, dF2xy,
!dF3xy, dF4xy
!REAL(8), DIMENSION (4, 4) :: A_bicub, C_bicub
!
!--------------------------------------------------
! 
!1) DEFINITION OF ARRAY u:
!
!--------------------------------------------------
u_min = e_umax
!
! arrays x_mesh_HT, y_mesh_HT
!
delta       = (u_end - u_min)/(NNN_HT-1)
y_mesh_HT   =  u_min + (/(i*delta, i=0,NNN_HT-1)/)
delta       = (x_mesh_max - x_mesh_min)/(MMM_HT-1)
x_mesh_HT   =  x_mesh_min + (/(i*delta, i=0,MMM_HT-1)/)
!
! array y_mesh_sat_HT
!
delta         = (u_end - u_min)/(NNN_sat_HT-1)
y_mesh_sat_HT = u_min + (/(i*delta, i=0,NNN_sat_HT-1)/)
!
!-----------------------------------------------------------------------
! 
!2) EVALUATION OF v_max(u) AND v_min(u) FOR EACH ELEMENT OF THE ARRAY u:
!
! Concerning the evaluation of v_min(u): P = 10MPa, vmin = v_B
! V_max: P = 0.5MPa, vmax = v_C
!
!-----------------------------------------------------------------------
!print*, 'start left'
!
!Left boundary, vmin
!
        T_guess = T_B
        v_guess = v_B

        ppp_HT(:,1) = p_max
!
        DO i=1,NNN_HT
!print*,i, 'in', ppp_HT(i,1), y_mesh_HT(i),T_guess, v_guess
!        
              CALL New_Rap2D(2, TTT_HT(i,1), vvv_HT(i,1), &
              & res, Niter, exitflag, ppp_HT(i,1), y_mesh_HT(i),&
              & T_guess, v_guess)
!
             CALL sound_speed(TTT_HT(i,1),vvv_HT(i,1), sound)
             ccc_HT(i,1) = sound
!print*,'out',TTT_HT(i,1), vvv_HT(i,1),res, Niter, exitflag,ccc_HT(i,1)
!
              IF (res > res_ref) THEN
                print*, "left_HT", res,i, "iter", Niter,"flag",exitflag
                STOP
              ENDIF
              T_guess = TTT_HT(i,1)
              v_guess = vvv_HT(i,1)              
         END DO
!
! Right boundary of HT domain. 
!
!print*, 'start right'
        T_guess = T_C
        v_guess = v_C
        ppp_HT(:,MMM_HT) = p_tri
!
        DO i=1,NNN_HT
!        
              CALL New_Rap2D(2, TTT_HT(i,MMM_HT), vvv_HT(i,MMM_HT), &
              & res, Niter, exitflag, ppp_HT(i,MMM_HT), y_mesh_HT(i),&
              & T_guess, v_guess)
!
             CALL sound_speed(TTT_HT(i,MMM_HT),vvv_HT(i,MMM_HT), sound)
             ccc_HT(i,MMM_HT) = sound
!
              IF (res > res_ref) THEN
                print*, "right_HT", res,i, "iter", Niter,"flag",exitflag
                STOP
              ENDIF
              T_guess = TTT_HT(i,MMM_HT)
              v_guess = vvv_HT(i,MMM_HT)
!              print*, vvv_HT(i,MMM_HT) 
         ENDDO
!   
! Middle domain of the two boundaries. v_min-->v_max, u_min-->u_max
! vvv_HT in the physical space is built by the X in the transformed
! space with a LINEAR SCALING TRANSFORMATION 
!print*, 'start middle'
!
        T_guess = TTT_HT(1,1)
        DO j=2,MMM_HT-1
           DO i=1, NNN_HT
!
        vvv_log(i,j) = ((x_mesh_HT(j) - x_mesh_min)/(x_mesh_max - x_mesh_min))*&
     &               (LOG10(vvv_HT(i,MMM_HT))-LOG10(vvv_HT(i,1)))+ LOG10(vvv_HT(i,1))
       vvv_HT(i,j) = 10_pr **vvv_log(i,j)
!                vvv_HT(i,j) = ((x_mesh_HT(j) - x_mesh_min) /&
!                & (x_mesh_max - x_mesh_min)) *(vvv_HT(i,MMM_HT)-vvv_HT(i,1))+&
!                & vvv_HT(i,1)
!                
                CALL New_Rap1D(1, temp, out_2, res, Niter,&
               & exitflag, y_mesh_HT(i), T_guess,&
               & vvv_HT(i,j), out_2)
!                
                TTT_HT(i,j) = temp
                T_guess = temp
!                
!
              CALL pressure(TTT_HT(i,j),vvv_HT(i,j),press)
              ppp_HT(i,j) = press
!
             CALL sound_speed(TTT_HT(i,j),vvv_HT(i,j), sound)
             ccc_HT(i,j) = sound
!
                IF (res > res_ref) THEN
                  print*, "middle_HT", res,i, "iter", Niter
                  STOP
                ENDIF
           END DO
        END DO
!        
!
!
!
! grid for Cp
!
!        DO j = 1,MMM_HT
!           DO  i = 1, NNN_HT
!               CALL heat_cap_p(TTT_HT(i,j),vvv_HT(i,j),cpcp_HT(i,j))
!               
!               IF ( (cpcp_HT(i,j) /= cpcp_HT(i,j)) .OR. (cpcp_HT(i,j) <= 0.0) .OR.&
!&                   (cpcp_HT(i,j) >  cp_cr  ) ) THEN
!                  print*,'cp no values in HT kJ/(kgK)', cpcp_HT(i,j)/1000.0,i,j
!                  cpcp_HT(i,j) = cp_cr
!               ENDIF
!
!             ENDDO
!        ENDDO
!
!

!------------------------------------------------------------------
!       
! Spline coeff construction
!
!---------------------------------------------------------------------
!print*, 'start spline left'
        T_guess = T_B
        v_guess = v_B
!
        DO i=1,NNN_sat_HT
!        
              CALL New_Rap2D(2, T_sat_pmax(i), v_left_HT(i), &
              & res, Niter, exitflag, p_max, y_mesh_sat_HT(i),&
              & T_guess, v_guess)
!
              IF (res > res_ref) THEN
                print*, "left_sp_HT", res,i, "iter", Niter,"flag",exitflag
                STOP
              ENDIF
              T_guess = T_sat_pmax(i)
              v_guess = v_left_HT(i)
!                
         END DO
!
!print*, 'start spline right'
        T_guess = T_C
        v_guess = v_C
!
        DO i=1,NNN_sat_HT
!        
              CALL New_Rap2D(2, T_sat_pmin(i), v_right_HT(i), &
              & res, Niter, exitflag, p_tri, y_mesh_sat_HT(i),&
              & T_guess, v_guess)
!
              IF (res > res_ref) THEN
                print*, "right_sp_HT", res,i, "iter", Niter,"flag",exitflag
                STOP
              ENDIF
              T_guess = T_sat_pmin(i)
              v_guess = v_right_HT(i)
!                
         END DO
!
! construction of spline coefficients associated to the boundary curves
!
DO i = 1, NNN_HT-1
   j = 3 * (i-1) + 1
   call polyfit(spline_left_HT(:,i), y_mesh_sat_HT(j:j+3), v_left_HT(j:j+3),ord_spline)
   call polyfit(spline_right_HT(:,i), y_mesh_sat_HT(j:j+3), v_right_HT(j:j+3),ord_spline)
ENDDO
!print*,'------------------------------------------------------'
print*,           'HT GRID FINISH'
!print*,'------------------------------------------------------'
!
!
END SUBROUTINE Grid_construction_High_Temperature                    
