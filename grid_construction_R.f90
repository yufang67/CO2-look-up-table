! ===================================================================
!
!                      Grid_construction_Right 
!
! ===================================================================
!
! Region Right (R): 
!
! From the Triple Point specific internal energy in the vapour region
! up to the specific internal energy value referred to the maximum point
! on the saturation curve.
! Bounded by the saturated vapour curve and the minimum pressure line
!
! note: 23/11/2017 equidistance on vvv_R for derivatives
! ===================================================================
SUBROUTINE grid_construction_right()

     USE def_constants
     USE def_variables
     USE non_linear_solvers
     USE properties
     USE grid_functions
!        
      IMPLICIT NONE
!
INTEGER                           ::  i, j, Niter, exitflag
INTEGER, DIMENSION (NNN_R,MMM_R)  ::  res_R
!
REAL(pr)   ::  u_min, res, delta, T_guess,v_guess,  vliq,vliq_guess,&
&              sound, press, temp, out_2,in_2 
!
REAL(pr), DIMENSION (NNN_sat_R) :: T_sat_sat, T_pmin,  res_sp
REAL(pr), DIMENSION (NNN_R,MMM_R) ::  vvv_log
REAL(pr), DIMENsION (NNN_R) :: v_max, v_min
!
! data used for bicubic interpolation
!
!REAL(pr) :: v1, v2, v3, v4, u1, u2, u3, u4
!REAL(pr) :: F1, F2, F3, F4, G1, G2, G3, G4
!REAL(pr) :: a1, a2, a3, a4, b1, b2, b3, b4
!REAL(pr) :: dF1v, dF2v, dF3v, dF4v, dF1u, dF2u, dF3u, dF4u, dF1vu, dF2vu,dF3vu, dF4vu
!REAL(pr) :: dF1vv, dF2vv, dF3vv, dF4vv, dF1uu, dF2uu, dF3uu, dF4uu
!REAL(pr) :: dF1x, dF2x, dF3x, dF4x, dF1y, dF2y, dF3y, dF4y, dF1xy, dF2xy,dF3xy, dF4xy
!REAL(pr), DIMENSION (4, 4) :: A_bicub, C_bicub
!
!--------------------------------------------------
! 
!1) DEFINITION OF ARRAY u:
!
!--------------------------------------------------
!
u_min = e_tri_R
!
! arrays x_mesh_R, y_mesh_R
!
delta    = (e_umax - u_min)/(NNN_R-1)
y_mesh_R = u_min + (/(i*delta, i=0,NNN_R-1)/)
delta    = (x_mesh_max - x_mesh_min)/(MMM_R-1)
x_mesh_R = x_mesh_min + (/(i*delta, i=0,MMM_R-1)/)
!
! array y_mesh_sat_R
!
delta         = (e_umax - u_min)/(NNN_sat_R-1)
y_mesh_sat_R  = u_min + (/(i*delta, i=0,NNN_sat_R-1)/)
!
!-----------------------------------------------------------------------
! 
!2) EVALUATION OF v_max(u) AND v_min(u) FOR EACH ELEMENT OF THE ARRAY u:
!
! Since the vapour part is now analyzed, v_max(u) is referred to points 
! corresponding to the p_min curve, while v_min(u) is referred to points
! associated to the saturated vapour line.vmax(1,1) = vmin(1,1)
!
!-----------------------------------------------------------------------
!
!-----------------vmax(pmin = p_tri)-----------------------------------
vvv_R(1,:) = v_tri_R
TTT_R(1,:) = T_tri
ppp_R(1,:) = p_tri
ppp_R(:,MMM_R) = p_tri
ccc_R(1,:) = c_tri_R
!
T_guess = T_tri
v_guess = v_tri_R
!
        DO i=2,NNN_R
!        
              CALL New_Rap2D(2, TTT_R(i,MMM_R), vvv_R(i,MMM_R), &
              & res, Niter, exitflag, ppp_R(i,MMM_R), y_mesh_R(i),&
              & T_guess, v_guess)
!
             CALL sound_speed(TTT_R(i,MMM_R),vvv_R(i,MMM_R), sound)
             ccc_R(i,MMM_R) = sound
!
              IF (res > res_ref) THEN
                print*, "right_R", res, i ,"flag",exitflag
                STOP
              ENDIF
              T_guess = TTT_R(i,MMM_R)
              v_guess = vvv_R(i,MMM_R)  
!       print*, "v_left", vvv_LL(i,1)
         END DO
         
!---------------------------------------------------------------------------
!---------------vmin(saturation line from triple point to umax)-----------
!
v_guess = v_tri_R
T_guess = T_tri
vliq_guess = v_tri_L
!
        DO i=2, NNN_R-1
!         
              CALL New_Rap3D(4,vliq,vvv_R(i,1),TTT_R(i,1),&
              & res, Niter, exitflag, y_mesh_R(i),in_2,&
              & vliq_guess, v_guess,T_guess)
!
              CALL pressure(TTT_R(i,1),vvv_R(i,1),press)
              ppp_R(i,1) = press
!
             CALL sound_speed(TTT_R(i,1),vvv_R(i,1), sound)
             ccc_R(i,1) = sound
!
              IF (res > res_ref) THEN
                print*, "left R", res,i, "iter", Niter
                STOP
              ENDIF
!
        v_guess = vvv_R(i,1)
        T_guess = TTT_R(i,1)
        vliq_guess = vliq
        
        END DO
        vvv_R(NNN_R,1) = v_umax
        ppp_R(NNN_R,1) = p_umax
        TTT_R(NNN_R,1) = T_umax
        ccc_R(NNN_R,1) = c_umax
        
        
!--------------------------------------------------------------------------
!------------vmin < v < vmax ---------------------------------------------
!
DO i = 2,NNN_R 
    DO j = 2,MMM_R - 1
!       vvv_log(i,j) = ((x_mesh_R(j) - x_mesh_min)/(x_mesh_max - x_mesh_min))*&
!     &               (LOG10(vvv_R(i,MMM_R))-LOG10(vvv_R(i,1)))+ LOG10(vvv_R(i,1))
!       vvv_R(i,j) = 10_pr **vvv_log(i,j)
       vvv_R(i,j) = ((x_mesh_R(j) - x_mesh_min)/(x_mesh_max - x_mesh_min))*&
  &               (vvv_R(i,MMM_R)-vvv_R(i,1))+ vvv_R(i,1)



       IF (vvv_R(i,j) .GT. vvv_R(i,MMM_R)) THEN
         print*, vvv_R(i,1),vvv_R(i,MMM_R), x_mesh_R(j)
         print*, ''
         print*, vvv_R(i,j), i, j
         print*, LOG10(vvv_R(i,1)), LOG10(vvv_R(i,MMM_R)) 
         stop
       ENDIF

    ENDDO
ENDDO
!
! ppp_R, TTT_R, ccc_R
!
T_guess = T_tri
!
        DO j=2,MMM_R-1
           DO i=2, NNN_R
!
              CALL New_Rap1D(1, temp, out_2, res, Niter,&
               & exitflag, y_mesh_R(i), T_guess,&
               & vvv_R(i,j), out_2)
!                
                TTT_R(i,j) = temp
!                
!
              CALL pressure(TTT_R(i,j),vvv_R(i,j),press)
              ppp_R(i,j) = press
!
              CALL sound_speed(TTT_R(i,j),vvv_R(i,j), sound)
              ccc_R(i,j) = sound
!
                IF (res > res_ref) THEN
                  print*, "middle R", res, i,j, "iter", Niter
                  STOP
                ENDIF
!          
        T_guess = temp
           END DO
        END DO
!
!
!
!
! grid for Cp
!
        DO j = 1,MMM_R
           DO  i = 1, NNN_R
               CALL heat_cap_p(TTT_R(i,j),vvv_R(i,j),cpcp_R(i,j))
!               
               IF ( (cpcp_R(i,j) /= cpcp_R(i,j)) .OR. (cpcp_R(i,j) <= 0.0) .OR.&
&                   (cpcp_R(i,j) >  cp_cr  ) ) THEN
                  print*,'cp no values in R (kJ/kgK)', cpcp_R(i,j)/1000.0,i,j
                  cpcp_R(i,j) = cp_cr
               ENDIF
!
             ENDDO
        ENDDO
!
!
        
!---------------------------------------------------------------- 
!
! MESH CONSTRUCTION for spline coefficients
!------------------------------------------------------------
!
v_Vsat(1) = v_tri_R
v_Vsat(NNN_sat_R) = v_umax
v_Vpmin(1) = v_tri_R
!Pmin line(right boundary)
!
T_guess = T_tri
v_guess = v_tri_R
!
DO i = 2, NNN_sat_R
!
     CALL New_Rap2D(2, T_pmin(i), v_Vpmin(i), &
     & res_sp(i), Niter, exitflag, p_tri, y_mesh_sat_R(i),&
     & T_guess, v_guess)
!
     T_guess = T_pmin(i)
     v_guess = v_Vpmin(i)
!
   IF (res_sp(i) > res_ref) THEN
       print*,'R-Convergence failed calculating v_Vpmin(',i,')'
       print*,'resnorm    = ',res_sp(i),'N_iter_New = ',Niter
   ENDIF
ENDDO
 
!
!v_Vsat(NNN_sat_R) = v_Vsat(NNN_sat_R-1)  ! near the maximum on the saturation curve
!
! saturation line
v_guess = v_tri_R
T_guess = T_tri
vliq_guess = v_tri_L
!
        DO i=2, NNN_sat_R-1
!         
              CALL New_Rap3D(4,vliq,v_Vsat(i),T_sat_sat(i),&
              & res_sp(i), Niter, exitflag, y_mesh_sat_R(i),in_2,&
              & vliq_guess, v_guess,T_guess)
!
        IF (res_sp(i) > res_ref) THEN
                print*, 'R-Convergence failed calculating v_Vsat(',i,')'
                print*,'resnorm    = ',res_sp(i),'N_iter_New = ',Niter
        ENDIF
!
        v_guess = v_Vsat(i)
        T_guess = T_sat_sat(i)
        vliq_guess = vliq
        END DO


 
      
!
! construction of spline coefficients associated to saturation curve nodes and
! pmin curve nodes
!
DO i = 1, NNN_R-1
   j = 3 * (i-1) + 1
   call polyfit(spline_Vsat(:,i), y_mesh_sat_R(j:j+3), v_Vsat (j:j+3),ord_spline)
   call polyfit(spline_pmin(:,i), y_mesh_sat_R(j:j+3), v_Vpmin(j:j+3),ord_spline)
ENDDO
!print*,'------------------------------------------------------'
print*,           'R GRID FINISH'
!print*,'------------------------------------------------------'
!
!
END SUBROUTINE grid_construction_right
