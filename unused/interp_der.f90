MODULE interp_der
!==================================================================
!Evaluate derivatives is not used any more because of the log grid
!==================================================================
     USE def_constants, ONLY: ord_spline, pr
     IMPLICIT NONE
     PRIVATE
     PUBLIC :: Lin_der, Lin_der_Log10
!
!
     CONTAINS
! ===================================================================
!
!
! ===================================================================
SUBROUTINE Lin_der(dT_dv_u, dT_du_v, dp_dv_u, dp_du_v,&
&          NNN, MMM, x_mesh,y_mesh, spline_pmin, spline_Vsat, u_in, v_in,TTT,ppp,vvv)
!
IMPLICIT NONE
!
!
REAL(pr), INTENT(OUT) :: dT_dv_u, dT_du_v, dp_dv_u, dp_du_v
REAL(pr), INTENT(IN)  :: u_in, v_in
INTEGER, INTENT(IN)  :: NNN, MMM
REAL(pr), DIMENSION(NNN), INTENT(IN)  :: y_mesh
REAL(pr), DIMENSION(MMM), INTENT(IN)  :: x_mesh
REAL(pr), DIMENSION(ord_spline+1,NNN-1), INTENT(IN) :: spline_pmin, spline_Vsat
REAL(pr), DIMENSION(NNN,MMM), INTENT(IN) :: TTT, ppp, vvv
!
!
INTEGER :: i, j
REAL(8) :: x_mesh_test, delta_u, delta_v, v_min_li, v_max_li,    &
&           T_down_left, T_down_right, T_up_left, T_up_right,    &   
&           p_down_left, p_down_right, p_up_left, p_up_right,    &      
&           v1, v2, v3, v4, u1, u2, a2, a3, a4, b3,              &
&  dX_dv_u, dX_du_v, dY_du_v, dT_dX_Y, dT_dY_X, dp_dX_Y, dp_dY_X  
!
!
!
delta_u = y_mesh(2) - y_mesh(1)
i       = INT((u_in - y_mesh(1))/delta_u) + 1
!
IF(( i .LE. 0) .OR. (i .GE. NNN)) THEN
  print*,'In Lin_der'
  print*,'i lower than 0 or greater then NNN',i, u_in 
  stop
ENDIF
!      
!
v_min_li = 0_pr
v_max_li = 0_pr
!
DO j = 1, ord_spline + 1
   v_min_li =v_min_li + spline_Vsat(ord_spline+2-j,i)*u_in**(j-1)
   v_max_li =v_max_li + spline_pmin(ord_spline+2-j,i)*u_in**(j-1)
ENDDO
!
! Evaluation of x_mesh_test:
!
delta_v     = x_mesh(2) - x_mesh(1)
x_mesh_test = ((x_mesh(MMM) - x_mesh(1)) / (v_max_li - &
&                v_min_li)) * (v_in - v_min_li) + x_mesh(1)
!
j = INT((x_mesh_test - x_mesh(1))/delta_v) + 1
!
IF ((j .LT. 1) .OR. (j .GE. MMM)) THEN
   print*,'In Lin_der'
   print*, 'Error in index determination j',j
   stop
ENDIF
!      
T_down_left  = TTT(i,j)
T_down_right = TTT(i,j+1)
T_up_left    = TTT(i+1,j)
T_up_right   = TTT(i+1,j+1)
!
p_down_left  = ppp(i,j)
p_down_right = ppp(i,j+1)
p_up_left    = ppp(i+1,j)
p_up_right   = ppp(i+1,j+1)
!
! Derivatives in transformed space:
!       
dT_dX_Y =((T_down_right - T_down_left) +& 
&         (T_up_right   + T_down_left  - T_down_right - T_up_left)*&
&         (u_in - y_mesh(i)) / delta_u )/delta_v 
!     
dT_dY_X =((T_up_left    - T_down_left) + &
&         (T_up_right   + T_down_left  - T_down_right - T_up_left)*&
&         (x_mesh_test - x_mesh(j)) / delta_v)/delta_u
!
dp_dX_Y =((p_down_right - p_down_left) + &
&         (p_up_right   + p_down_left  - p_down_right - p_up_left)*&
&         (u_in - y_mesh(i)) / delta_u )/delta_v
!     
dp_dY_X =((p_up_left    - p_down_left) + &
&         (p_up_right   + p_down_left  - p_down_right - p_up_left)*&
&         (x_mesh_test - x_mesh(j)) / delta_v)/delta_u    
!
! TRANSFORMED SPACE
!
v1 = vvv(i,   j  )
v2 = vvv(i+1, j  )
v3 = vvv(i+1, j+1)
v4 = vvv(i,   j+1)
!
u1 = y_mesh(i)
u2 = y_mesh(i+1)
!
a2 = v4-v1;  a3 = v2-v1;  a4 = v1+v3-v2-v4
b3 = u2-u1
!      
dX_dv_u = delta_v / ( a2 + a4*(u_in - y_mesh(i))/delta_u )
dX_du_v =-dX_dv_u * ( a3 + a4*(x_mesh_test-x_mesh(j))/delta_v ) / b3
dY_du_v = delta_u / b3 
!
! Derivatives in physical space:
!
dT_dv_u = dX_dv_u * dT_dX_Y
dT_du_v = dX_du_v * dT_dX_Y + dY_du_v * dT_dY_X
!
dp_dv_u = dX_dv_u * dp_dX_Y
dp_du_v = dX_du_v * dp_dX_Y + dY_du_v * dp_dY_X
print*, dX_du_v * dp_dX_Y , dY_du_v * dp_dY_X, dp_du_v
!
IF ( (dT_dv_u /= dT_dv_u) .OR.&
&    (dT_du_v /= dT_du_v) .OR.&
&    (dp_dv_u /= dp_dv_u) .OR.&
&    (dp_du_v /= dp_du_v) )THEN
   print*, 'dT_dv_u',dT_dv_u 
   print*, 'dT_du_v',dT_du_v
   print*, 'dp_dv_u',dp_dv_u
   print*, 'dp_du_v',dp_du_v
   STOP
ENDIF
!      
END SUBROUTINE Lin_der
!
!
! ===================================================================
!
!
! ===================================================================
!
SUBROUTINE Lin_der_Log10(dT_dv_u, dT_du_v, dp_dv_u, dp_du_v,&
&          NNN, MMM, x_mesh,y_mesh, spline_pmin, spline_Vsat, u_in, v_in,TTT,ppp,vvv)
!
IMPLICIT NONE
!
!
REAL(pr), INTENT(OUT) :: dT_dv_u, dT_du_v, dp_dv_u, dp_du_v
REAL(pr), INTENT(IN)  :: u_in, v_in
INTEGER, INTENT(IN)  :: NNN, MMM
REAL(pr), DIMENSION(NNN), INTENT(IN)  :: y_mesh
REAL(pr), DIMENSION(MMM), INTENT(IN)  :: x_mesh
REAL(pr), DIMENSION(MMM)  :: vlog
REAL(pr), DIMENSION(ord_spline+1,NNN-1), INTENT(IN) :: spline_pmin, spline_Vsat
REAL(pr), DIMENSION(NNN,MMM), INTENT(IN) :: TTT, ppp, vvv
!
!
INTEGER :: i, j
REAL(8) :: x_mesh_test, delta_u, delta_v, v_min_li, v_max_li,    &
&           T_down_left, T_down_right, T_up_left, T_up_right,    &
&           p_down_left, p_down_right, p_up_left, p_up_right,    &
&           v1, v2, v3, v4, u1, u2, a2, a3, a4, b3,              &
&  dX_dv_u, dX_du_v, dY_du_v, dT_dX_Y, dT_dY_X, dp_dX_Y, dp_dY_X
!
!
!
delta_u = y_mesh(2) - y_mesh(1)
i       = INT((u_in - y_mesh(1))/delta_u) + 1
!
IF(( i .LE. 0) .OR. (i .GE. NNN)) THEN
  print*,'In Lin_der'
  print*,'i lower than 0 or greater then NNN',i, u_in
  stop
ENDIF
!      
!
v_min_li = 0_pr
v_max_li = 0_pr
!
DO j = 1, ord_spline + 1
   v_min_li =v_min_li + spline_Vsat(ord_spline+2-j,i)*u_in**(j-1)
   v_max_li =v_max_li + spline_pmin(ord_spline+2-j,i)*u_in**(j-1)
ENDDO
!
! Evaluation of x_mesh_test:
!
delta_v     = x_mesh(2) - x_mesh(1)
x_mesh_test = ((x_mesh(MMM) - x_mesh(1)) / (LOG10(v_max_li) - LOG10(v_min_li)))&
&* (LOG10(v_in) - LOG10(v_min_li)) +  x_mesh(1)
!



j = INT((x_mesh_test - x_mesh(1))/delta_v) + 1
!
IF ((j .LT. 1) .OR. (j .GE. MMM)) THEN
   print*,'In Lin_der'
   print*, 'Error in index determination j',j
   stop
ENDIF
!      
T_down_left  = TTT(i,j)
T_down_right = TTT(i,j+1)
T_up_left    = TTT(i+1,j)
T_up_right   = TTT(i+1,j+1)
!
p_down_left  = ppp(i,j)
p_down_right = ppp(i,j+1)
p_up_left    = ppp(i+1,j)
p_up_right   = ppp(i+1,j+1)
!
! Derivatives in transformed space:
!       
dT_dX_Y =((T_down_right - T_down_left) +&
&         (T_up_right   + T_down_left  - T_down_right - T_up_left)*&
&         (u_in - y_mesh(i)) / delta_u )/delta_v
!     
dT_dY_X =((T_up_left    - T_down_left) + &
&         (T_up_right   + T_down_left  - T_down_right - T_up_left)*&
&         (x_mesh_test - x_mesh(j)) / delta_v)/delta_u
!
dp_dX_Y =((p_down_right - p_down_left) + &
&         (p_up_right   + p_down_left  - p_down_right - p_up_left)*&
&         (u_in - y_mesh(i)) / delta_u )/delta_v
!     
dp_dY_X =((p_up_left    - p_down_left) + &
&         (p_up_right   + p_down_left  - p_down_right - p_up_left)*&
&         (x_mesh_test - x_mesh(j)) / delta_v)/delta_u
!
! TRANSFORMED SPACE
!
v1 = LOG10(vvv(i,   j  ))
v2 = LOG10(vvv(i+1, j  ))
v3 = LOG10(vvv(i+1, j+1))
v4 = LOG10(vvv(i,   j+1))
!
u1 = y_mesh(i)
u2 = y_mesh(i+1)
!
a2 = v4-v1;  a3 = v2-v1;  a4 = v1+v3-v2-v4
b3 = u2-u1
!      
!dX_dv_u = delta_v / ( a2 + a4*(u_in - y_mesh(i))/delta_u )
dX_dv_u = delta_v / ( a2 + a4*(u_in - y_mesh(i))/delta_u ) / (vvv(i,j)*2.302585) 
dX_du_v =-dX_dv_u * ( a3 + a4*(x_mesh_test-x_mesh(j))/delta_v ) / b3
dY_du_v = delta_u / b3
!       
! Derivatives in physical space:
!
dT_dv_u = dX_dv_u * dT_dX_Y
dT_du_v = dX_du_v * dT_dX_Y + dY_du_v * dT_dY_X
dp_dv_u = dX_dv_u * dp_dX_Y
dp_du_v = dX_du_v * dp_dX_Y + dY_du_v * dp_dY_X
!
print*, dX_du_v * dp_dX_Y , dY_du_v * dp_dY_X, dp_du_v
IF ( (dT_dv_u /= dT_dv_u) .OR.&
&    (dT_du_v /= dT_du_v) .OR.&
&    (dp_dv_u /= dp_dv_u) .OR.&
&    (dp_du_v /= dp_du_v) )THEN
   print*, 'dT_dv_u',dT_dv_u
   print*, 'dT_du_v',dT_du_v
   print*, 'dp_dv_u',dp_dv_u
   print*, 'dp_du_v',dp_du_v
   STOP
ENDIF
!      
END SUBROUTINE Lin_der_Log10
!
END MODULE interp_der
