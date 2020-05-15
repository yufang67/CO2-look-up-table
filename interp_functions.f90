MODULE interp_functions
!
     USE def_variables
     USE def_constants

     IMPLICIT NONE
     PRIVATE
     PUBLIC :: Lin_int_Left_Low ,Lin_int_Left_High,Lin_int_Log10,&
&              Lin_TP_Log10,Lin_TP,Lin_int
!
!
     CONTAINS
!
! ===================================================================
!
!                         Lin_int_Left_Low
!
! ===================================================================
!
! It receives as input a couple of values (u_in, v_in),
! identifies the cell of the grid within the point is, and, finally,
! evaluates properties using the matrices already computed, by using
! a bilinear interpolation between the vertices of the cell.
!
! * Lin_int_Left_Low - linear scaling transformation and spline 
!                      coefficients "linearly" determined
! ===================================================================
SUBROUTINE Lin_int_Left_Low(T_out, p_out, c_out, u_in, v_in)
!
! Global Variables
REAL(pr), INTENT(OUT) :: T_out, p_out, c_out
REAL(pr), INTENT(IN)  :: u_in, v_in
!
! Local Variables
INTEGER :: i, j
REAL(pr) :: x_mesh_test, delta_u, delta_v, v_min_li, v_max_li,v_max_li_log, &
&           T_down_left, T_down_right, T_up_left, T_up_right, &
&           p_down_left, p_down_right, p_up_left, p_up_right, &
&           c_down_left, c_down_right, c_up_left, c_up_right, &
&           DEN, A, B, D, E
!
! to identify the cell within the test-point is:
!
delta_u = y_mesh_LL(2) - y_mesh_LL(1)
i = INT((u_in - y_mesh_LL(1))/delta_u) + 1
!
! once that the interval on the vertical axis has been identified, 
! it is possible to know the spline coefficients associated to it. 
! By simply using them it is possible to evaluate the v_min_li 
! (related to the pmax curve) and the v_max_li 
! (related to the sat curve) for that u_in:
!
v_min_li = 0_pr
v_max_li = 0_pr
v_max_li_log = 0_pr
!
!
DO j = 1, ord_spline + 1
   v_min_li = v_min_li + spline_pmax_LL(ord_spline+2-j,i) * u_in**(j-1)
!Evaluate Vmin e Vmax, on the row i(high already identified)
   v_max_li = v_max_li + spline_Lsat_LL(ord_spline+2-j,i) * u_in**(j-1)
!   v_max_li_log = v_max_li_log + spline_Lsat_LL(ord_spline+2-j,i) * u_in**(j-1)
!   v_max_li     = 10_pr ** v_max_li_log
!Approximated through the spline 
ENDDO
!
! Evaluation of x_mesh_test: Lin_int_Left_Low
! We have the X to identify
x_mesh_test = ((x_mesh_LL(MMM_LL) - x_mesh_LL(1)) / &   
& (v_max_li - v_min_li)) * (v_in - v_min_li) + x_mesh_LL(1)
!x_mesh_test = ((x_mesh_LL(MMM_LL) - x_mesh_LL(1)) / (LOG10(v_max_li) - LOG10(v_min_li)))&
!&* (LOG10(v_in) - LOG10(v_min_li)) +  x_mesh_LL(1)
!
! to identify the interval on the horizontal axis
delta_v = x_mesh_LL(2) - x_mesh_LL(1)
j = INT((x_mesh_test - x_mesh_LL(1))/delta_v) + 1 !X of the cell identified
! 
! once that the cell has been identified, 
! it is possible to evaluate properties at the vertices:
T_down_left  = TTT_LL(i,j)
T_down_right = TTT_LL(i,j+1)
T_up_left    = TTT_LL(i+1,j)
T_up_right   = TTT_LL(i+1,j+1)
!
p_down_left  = ppp_LL(i,j)
p_down_right = ppp_LL(i,j+1)
p_up_left    = ppp_LL(i+1,j)
p_up_right   = ppp_LL(i+1,j+1)
!
c_down_left  = ccc_LL(i,j)
c_down_right = ccc_LL(i,j+1)
c_up_left    = ccc_LL(i+1,j)
c_up_right   = ccc_LL(i+1,j+1)
!
! Bilinear Interpolation:
!
DEN = (x_mesh_LL(j+1)- x_mesh_LL(j))*(y_mesh_LL(i)- y_mesh_LL(i+1))
A   = (x_mesh_LL(j+1)- x_mesh_test) *(y_mesh_LL(i)- u_in)
B   = (x_mesh_test   - x_mesh_LL(j))*(y_mesh_LL(i)- u_in)
D   = (x_mesh_LL(j+1)- x_mesh_test) *(u_in - y_mesh_LL(i+1))
E   = (x_mesh_test   - x_mesh_LL(j))*(u_in - y_mesh_LL(i+1))
!
T_out = (A/DEN) * T_up_left + (B/DEN) * T_up_right + (D/DEN) * T_down_left + &
&    (E/DEN) * T_down_right
!
p_out = (A/DEN) * p_up_left + (B/DEN) * p_up_right + (D/DEN) * p_down_left + &
&    (E/DEN) * p_down_right
!
c_out = (A/DEN) * c_up_left + (B/DEN) * c_up_right + (D/DEN) * c_down_left + &
&    (E/DEN) * c_down_right
!
END SUBROUTINE Lin_int_Left_Low
!
!
! ===================================================================
!
!                              Lin_int_Left_High
!
! ===================================================================
!
! It receives as input a couple of values (u_in, v_in),
! identifies the cell of the grid within the point is, and, finally,
! evaluates properties using the matrices already computed, by using
! a bilinear interpolation between the vertices of the cell.
!
! * Lin_int_Left_High - linear scaling transformation and spline 
!                       coefficients "logaritmically" determined
!                       at the saturation
! ===================================================================
SUBROUTINE Lin_int_Left_High(T_out, p_out, c_out, u_in, v_in)
!
!
! Global Variables
REAL(pr), INTENT(OUT)  :: T_out, p_out, c_out
REAL(pr), INTENT(IN)    :: u_in, v_in
!
! Local Variables
INTEGER  :: i, j
REAL(pr)  :: x_mesh_test, delta_u, delta_v, v_min_li, v_max_li, v_max_li_log, &
&           T_down_left, T_down_right, T_up_left, T_up_right, &
&           p_down_left, p_down_right, p_up_left, p_up_right, &
&           c_down_left, c_down_right, c_up_left, c_up_right, &
&           DEN, A, B, D, E
!
!
! to identify the cell within the test-point is:
!
delta_u = y_mesh_LH(2)     - y_mesh_LH(1)
i       = INT((u_in - y_mesh_LH(1))/delta_u) + 1
!
! Initialization of v_min_li and v_max_li
v_min_li     = 0_pr
v_max_li     = 0_pr
v_max_li_log = 0_pr
!
DO j = 1, ord_spline + 1
   v_min_li     = v_min_li     + spline_pmax_LH(ord_spline+2-j,i) * u_in**(j-1)
   v_max_li     = v_max_li     + spline_Lsat_LH(ord_spline+2-j,i) * u_in**(j-1)
!   v_max_li_log = v_max_li_log + spline_Lsat_LH(ord_spline+2-j,i) * u_in**(j-1)
!   v_max_li     = 10_pr ** v_max_li_log
ENDDO
!
! Evaluation of x_mesh_test:
! In the middle of the LH region v is discretized logarithmic
x_mesh_test = ((x_mesh_LH(MMM_LH) - x_mesh_LH(1)) / (LOG10(v_max_li) - LOG10(v_min_li)))&
&* (LOG10(v_in) - LOG10(v_min_li)) +  x_mesh_LH(1)
!x_mesh_test = ((x_mesh_LH(MMM_LH) - x_mesh_LH(1)) / &
!& (v_max_li - v_min_li)) * (v_in - v_min_li) + x_mesh_LH(1)
!
! to identify the interval on the horizontal axis
delta_v = x_mesh_LH(2)     - x_mesh_LH(1)
j       = INT((x_mesh_test - x_mesh_LH(1))/delta_v) + 1
!
! once that the cell has been identified,it is possible to evaluate properties at the vertices:
T_down_left  = TTT_LH(i,j)
T_down_right = TTT_LH(i,j+1)
T_up_left    = TTT_LH(i+1,j)
T_up_right   = TTT_LH(i+1,j+1)
!
p_down_left  = ppp_LH(i,j)
p_down_right = ppp_LH(i,j+1)
p_up_left    = ppp_LH(i+1,j)
p_up_right   = ppp_LH(i+1,j+1)
!
c_down_left  = ccc_LH(i,j)
c_down_right = ccc_LH(i,j+1)
c_up_left    = ccc_LH(i+1,j)
c_up_right   = ccc_LH(i+1,j+1)
!
! Bilinear Interpolation:
!
DEN = (x_mesh_LH(j+1)- x_mesh_LH(j))* (y_mesh_LH(i)- y_mesh_LH(i+1))
A   = (x_mesh_LH(j+1)- x_mesh_test) * (y_mesh_LH(i)- u_in)
B   = (x_mesh_test   - x_mesh_LH(j))* (y_mesh_LH(i)- u_in)
D   = (x_mesh_LH(j+1)- x_mesh_test) * (u_in - y_mesh_LH(i+1))
E   = (x_mesh_test   - x_mesh_LH(j))* (u_in - y_mesh_LH(i+1))
!
T_out = (A/DEN) * T_up_left + (B/DEN) * T_up_right + (D/DEN) * T_down_left + &
&    (E/DEN) * T_down_right
!
p_out = (A/DEN) * p_up_left + (B/DEN) * p_up_right + (D/DEN) * p_down_left + &
&    (E/DEN) * p_down_right
!
c_out = (A/DEN) * c_up_left + (B/DEN) * c_up_right + (D/DEN) * c_down_left + &
&    (E/DEN) * c_down_right
!
END SUBROUTINE Lin_int_Left_High
!
!
!
! ===================================================================
!
!                              Lin_int_Log10 
!
! ===================================================================
!
! It receives as input a couple of values (u_in, v_in),
! identifies the cell of the grid within the point is, and, finally,
! evaluates properties using the matrices already computed, by using
! a bilinear interpolation between the vertices of the cell.
!
! * Lin_int_Log10 - log scaling transformation
! note : 23/23/2017 not log scaling any more 
! ===================================================================
!
SUBROUTINE Lin_int_Log10(T_out, p_out, c_out, u_in, v_in, NNN, MMM, x_mesh,&
&y_mesh, spline_pmin, spline_Vsat, TTT, ppp, ccc)
!
! Global Variables
REAL(pr), INTENT(OUT) :: T_out, p_out, c_out
REAL(pr), INTENT(IN)  :: u_in, v_in
INTEGER, INTENT(IN)  :: NNN, MMM
REAL(pr), DIMENSION(NNN), INTENT(IN)  :: y_mesh
REAL(pr), DIMENSION(MMM), INTENT(IN)  :: x_mesh
REAL(pr), DIMENSION(ord_spline+1,NNN-1), INTENT(IN) :: spline_pmin, spline_Vsat
REAL(pr), DIMENSION(NNN,MMM), INTENT(IN) :: TTT, ppp, ccc
!
! Local Variables
INTEGER  :: i, j
REAL(pr)  :: x_mesh_test, delta_u, delta_v, v_min_li, v_max_li, &
&            T_down_left, T_down_right, T_up_left, T_up_right, &
&            p_down_left, p_down_right, p_up_left, p_up_right, &
&            c_down_left, c_down_right, c_up_left, c_up_right, &
&            DEN, A, B, D, E
!
! to identify the cell within the test-point is:
!
delta_u = y_mesh(2) - y_mesh(1)
i       = INT((u_in - y_mesh(1))/delta_u) + 1 ! to identify the interval on the vertical axis
!
! Initialization of v_min_li and v_max_li
v_min_li = 0_pr
v_max_li = 0_pr
!
DO j = 1, ord_spline + 1
   v_min_li = v_min_li + spline_Vsat(ord_spline+2-j,i) * u_in**(j-1)
   v_max_li = v_max_li + spline_pmin(ord_spline+2-j,i) * u_in**(j-1)
ENDDO
!
! Evaluation of x_mesh_test:
!
x_mesh_test = ((x_mesh(MMM) - x_mesh(1)) / (LOG10(v_max_li) - LOG10(v_min_li)))&
&* (LOG10(v_in) - LOG10(v_min_li)) +  x_mesh(1)
!
!x_mesh_test = ((x_mesh(MMM) - x_mesh(1)) / (v_max_li - v_min_li))&
!&* (v_in - v_min_li) +  x_mesh(1)
!
delta_v = x_mesh(2) - x_mesh(1)
j = INT((x_mesh_test - x_mesh(1))/delta_v) + 1 ! to identify the interval on the horizontal axis
!
!
! once that the cell has been identified, 
! it is possible to evaluate properties at the vertices:
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
c_down_left  = ccc(i,j)
c_down_right = ccc(i,j+1)
c_up_left    = ccc(i+1,j)
c_up_right   = ccc(i+1,j+1)
!
! Bilinear Interpolation:
!
DEN = (x_mesh(j+1) - x_mesh(j))   * (y_mesh(i) - y_mesh(i+1))
A   = (x_mesh(j+1) - x_mesh_test) * (y_mesh(i) - u_in)
B   = (x_mesh_test - x_mesh(j))   * (y_mesh(i) - u_in)
D   = (x_mesh(j+1) - x_mesh_test) * (u_in - y_mesh(i+1))
E   = (x_mesh_test - x_mesh(j))   * (u_in - y_mesh(i+1))
!
T_out = (A/DEN) * T_up_left + (B/DEN) * T_up_right + (D/DEN) * T_down_left + &
&    (E/DEN) * T_down_right
!
p_out = (A/DEN) * p_up_left + (B/DEN) * p_up_right + (D/DEN) * p_down_left + &
&    (E/DEN) * p_down_right
!
c_out = (A/DEN) * c_up_left + (B/DEN) * c_up_right + (D/DEN) * c_down_left + &
&    (E/DEN) * c_down_right
!
END SUBROUTINE Lin_int_Log10
!
!=============================================================================
!
!Bilinear_interp
!
!============================================================================
SUBROUTINE Lin_int(T_out, p_out, c_out, u_in, v_in, NNN, MMM, x_mesh,&
&y_mesh, spline_pmin, spline_Vsat, TTT, ppp, ccc)
!
! Global Variables
REAL(pr), INTENT(OUT) :: T_out, p_out, c_out
REAL(pr), INTENT(IN)  :: u_in, v_in
INTEGER, INTENT(IN)  :: NNN, MMM
REAL(pr), DIMENSION(NNN), INTENT(IN)  :: y_mesh
REAL(pr), DIMENSION(MMM), INTENT(IN)  :: x_mesh
REAL(pr), DIMENSION(ord_spline+1,NNN-1), INTENT(IN) :: spline_pmin, spline_Vsat
REAL(pr), DIMENSION(NNN,MMM), INTENT(IN) :: TTT, ppp, ccc
!
! Local Variables
INTEGER  :: i, j
REAL(pr)  :: x_mesh_test, delta_u, delta_v, v_min_li, v_max_li, &
&            T_down_left, T_down_right, T_up_left, T_up_right, &
&            p_down_left, p_down_right, p_up_left, p_up_right, &
&            c_down_left, c_down_right, c_up_left, c_up_right, &
&            DEN, A, B, D, E
!
! to identify the cell within the test-point is:
!
delta_u = y_mesh(2) - y_mesh(1)
i       = INT((u_in - y_mesh(1))/delta_u) + 1 ! to identify the interval on the vertical axis
!
! Initialization of v_min_li and v_max_li
v_min_li = 0_pr
v_max_li = 0_pr
!
DO j = 1, ord_spline + 1
   v_min_li = v_min_li + spline_Vsat(ord_spline+2-j,i) * u_in**(j-1)
   v_max_li = v_max_li + spline_pmin(ord_spline+2-j,i) * u_in**(j-1)
ENDDO
!
! Evaluation of x_mesh_test:
!
!
x_mesh_test = ((x_mesh(MMM) - x_mesh(1)) / (v_max_li - v_min_li))&
&* (v_in - v_min_li) +  x_mesh(1)
!
delta_v = x_mesh(2) - x_mesh(1)
j = INT((x_mesh_test - x_mesh(1))/delta_v) + 1 ! to identify the interval on the horizontal axis
!
!
! once that the cell has been identified, 
! it is possible to evaluate properties at the vertices:
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
c_down_left  = ccc(i,j)
c_down_right = ccc(i,j+1)
c_up_left    = ccc(i+1,j)
c_up_right   = ccc(i+1,j+1)
!
! Bilinear Interpolation:
!
DEN = (x_mesh(j+1) - x_mesh(j))   * (y_mesh(i) - y_mesh(i+1))
A   = (x_mesh(j+1) - x_mesh_test) * (y_mesh(i) - u_in)
B   = (x_mesh_test - x_mesh(j))   * (y_mesh(i) - u_in)
D   = (x_mesh(j+1) - x_mesh_test) * (u_in - y_mesh(i+1))
E   = (x_mesh_test - x_mesh(j))   * (u_in - y_mesh(i+1))
!
T_out = (A/DEN) * T_up_left + (B/DEN) * T_up_right + (D/DEN) * T_down_left + &
&    (E/DEN) * T_down_right
!
p_out = (A/DEN) * p_up_left + (B/DEN) * p_up_right + (D/DEN) * p_down_left + &
&    (E/DEN) * p_down_right
!
c_out = (A/DEN) * c_up_left + (B/DEN) * c_up_right + (D/DEN) * c_down_left + &
&    (E/DEN) * c_down_right
!
END SUBROUTINE Lin_int
! ===================================================================
!
!                         Bilinear_interp_Log10_TP 
!
! ===================================================================
!
! It receives as input a couple of values (u_in, v_in),
! identifies the cell of the grid within the point is, and, finally,
! evaluates properties using the matrices already computed, by using
! a bilinear interpolation between the vertices of the cell.
!
! * Lin_int_Log10 - log scaling transformation
! BUG: v_in ~ vsat => j=0 but 1<j<100 
! ===================================================================
!
SUBROUTINE Lin_TP_Log10(T_out, p_out, c_out, x_out, u_in, v_in, NNN, MMM, x_mesh,&
&y_mesh, spline_pmin, spline_Vsat, TTT, ppp, ccc, xxx)
!
! Global Variables
REAL(pr), INTENT(OUT) :: T_out, p_out, c_out, x_out
REAL(pr), INTENT(IN)  :: u_in, v_in
INTEGER, INTENT(IN)  :: NNN, MMM
REAL(pr), DIMENSION(NNN), INTENT(IN)  :: y_mesh
REAL(pr), DIMENSION(MMM), INTENT(IN)  :: x_mesh
REAL(pr), DIMENSION(ord_spline+1,NNN-1), INTENT(IN) :: spline_pmin, spline_Vsat
REAL(pr), DIMENSION(NNN,MMM), INTENT(IN) :: TTT, ppp, ccc, xxx
!
! Local Variables
INTEGER  :: i, j
REAL(pr)  :: x_mesh_test, delta_u, delta_v, v_min_li, v_max_li, &
&            T_down_left, T_down_right, T_up_left, T_up_right, &
&            p_down_left, p_down_right, p_up_left, p_up_right, &
&            c_down_left, c_down_right, c_up_left, c_up_right, &
&            x_down_left, x_down_right, x_up_left, x_up_right, &
&            DEN, A, B, D, E
!
! to identify the cell within the test-point is:
!
delta_u = y_mesh(2) - y_mesh(1)
i       = INT((u_in - y_mesh(1))/delta_u) + 1 ! to identify the interval on the vertical axis
!
! Initialization of v_min_li and v_max_li
v_min_li = 0_pr
v_max_li = 0_pr
!
DO j = 1, ord_spline + 1
   v_min_li = v_min_li + spline_Vsat(ord_spline+2-j,i) * u_in**(j-1)
   v_max_li = v_max_li + spline_pmin(ord_spline+2-j,i) * u_in**(j-1)
ENDDO
!
! Evaluation of x_mesh_test:
!
!print*, x_mesh(MMM), x_mesh(1)
!print*, LOG10(v_max_li),LOG10(v_min_li)
!print*, LOG10(v_in),LOG10(v_min_li)
x_mesh_test = ((x_mesh(MMM) - x_mesh(1)) / (LOG10(v_max_li) - LOG10(v_min_li)))&
&* (LOG10(v_in) - LOG10(v_min_li)) +  x_mesh(1)
!
delta_v = x_mesh(2) - x_mesh(1)
j = INT((x_mesh_test - x_mesh(1))/delta_v) + 1 ! to identify the interval on the horizontal axis
! for v_in very close to the saturation curve
IF (j==0) THEN 
    j=1
ENDIF
!print*, 'x_mesh_test', x_mesh_test, x_mesh(1),delta_v
!
! once that the cell has been identified, 
! it is possible to evaluate properties at the vertices:
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
!print*, p_down_left,p_down_right,p_up_left,p_up_right
!print*, i,j
!
c_down_left  = ccc(i,j)
c_down_right = ccc(i,j+1)
c_up_left    = ccc(i+1,j)
c_up_right   = ccc(i+1,j+1)
!
x_down_left  = xxx(i,j)
x_down_right = xxx(i,j+1)
x_up_left    = xxx(i+1,j)
x_up_right   = xxx(i+1,j+1)
!
! Bilinear Interpolation:
!
DEN = (x_mesh(j+1) - x_mesh(j))   * (y_mesh(i) - y_mesh(i+1))
A   = (x_mesh(j+1) - x_mesh_test) * (y_mesh(i) - u_in)
B   = (x_mesh_test - x_mesh(j))   * (y_mesh(i) - u_in)
D   = (x_mesh(j+1) - x_mesh_test) * (u_in - y_mesh(i+1))
E   = (x_mesh_test - x_mesh(j))   * (u_in - y_mesh(i+1))
!
T_out = (A/DEN) * T_up_left + (B/DEN) * T_up_right + (D/DEN) * T_down_left + &
&    (E/DEN) * T_down_right
!
p_out = (A/DEN) * p_up_left + (B/DEN) * p_up_right + (D/DEN) * p_down_left + &
&    (E/DEN) * p_down_right
!
c_out = (A/DEN) * c_up_left + (B/DEN) * c_up_right + (D/DEN) * c_down_left + &
&    (E/DEN) * c_down_right
!
x_out = (A/DEN) * x_up_left + (B/DEN) * x_up_right + (D/DEN) * x_down_left + &
&    (E/DEN) * x_down_right
!
END SUBROUTINE Lin_TP_Log10
!
!
! ===================================================================
!
!                         Bilinear_interp_TP
!
! ===================================================================
!
! It receives as input a couple of values (u_in, v_in),
! identifies the cell of the grid within the point is, and, finally,
! evaluates properties using the matrices already computed, by using
! a bilinear interpolation between the vertices of the cell.
!
! ===================================================================
!
SUBROUTINE Lin_TP(T_out, p_out, c_out, x_out, u_in, v_in, NNN, MMM, x_mesh,&
&y_mesh, spline_pmin, spline_Vsat, TTT, ppp, ccc, xxx)
!
! Global Variables
REAL(pr), INTENT(OUT) :: T_out, p_out, c_out, x_out
REAL(pr), INTENT(IN)  :: u_in, v_in
INTEGER, INTENT(IN)  :: NNN, MMM
REAL(pr), DIMENSION(NNN), INTENT(IN)  :: y_mesh
REAL(pr), DIMENSION(MMM), INTENT(IN)  :: x_mesh
REAL(pr), DIMENSION(ord_spline+1,NNN-1), INTENT(IN) :: spline_pmin, spline_Vsat
REAL(pr), DIMENSION(NNN,MMM), INTENT(IN) :: TTT, ppp, ccc, xxx
!
! Local Variables
INTEGER  :: i, j
REAL(pr)  :: x_mesh_test, delta_u, delta_v, v_min_li, v_max_li, &
&            T_down_left, T_down_right, T_up_left, T_up_right, &
&            p_down_left, p_down_right, p_up_left, p_up_right, &
&            c_down_left, c_down_right, c_up_left, c_up_right, &
&            x_down_left, x_down_right, x_up_left, x_up_right, &
&            DEN, A, B, D, E
!
! to identify the cell within the test-point is:
!
delta_u = y_mesh(2) - y_mesh(1)
i       = INT((u_in - y_mesh(1))/delta_u) + 1 ! to identify the interval on the vertical axis
!
! Initialization of v_min_li and v_max_li
v_min_li = 0_pr
v_max_li = 0_pr
!
DO j = 1, ord_spline + 1
   v_min_li = v_min_li + spline_Vsat(ord_spline+2-j,i) * u_in**(j-1)
   v_max_li = v_max_li + spline_pmin(ord_spline+2-j,i) * u_in**(j-1)
ENDDO
!
! Evaluation of x_mesh_test:
!
x_mesh_test = ((x_mesh(MMM) - x_mesh(1)) / &
& (v_max_li - v_min_li)) * (v_in - v_min_li) + x_mesh(1)
!
delta_v = x_mesh(2) - x_mesh(1)
j = INT((x_mesh_test - x_mesh(1))/delta_v) + 1 ! to identify the interval on the horizontal axis
!
!
! once that the cell has been identified, 
! it is possible to evaluate properties at the vertices:
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
c_down_left  = ccc(i,j)
c_down_right = ccc(i,j+1)
c_up_left    = ccc(i+1,j)
c_up_right   = ccc(i+1,j+1)
!
x_down_left  = xxx(i,j)
x_down_right = xxx(i,j+1)
x_up_left    = xxx(i+1,j)
x_up_right   = xxx(i+1,j+1)
!
! Bilinear Interpolation:
!
DEN = (x_mesh(j+1) - x_mesh(j))   * (y_mesh(i) - y_mesh(i+1))
A   = (x_mesh(j+1) - x_mesh_test) * (y_mesh(i) - u_in)
B   = (x_mesh_test - x_mesh(j))   * (y_mesh(i) - u_in)
D   = (x_mesh(j+1) - x_mesh_test) * (u_in - y_mesh(i+1))
E   = (x_mesh_test - x_mesh(j))   * (u_in - y_mesh(i+1))
!
T_out = (A/DEN) * T_up_left + (B/DEN) * T_up_right + (D/DEN) * T_down_left + &
&    (E/DEN) * T_down_right
!
p_out = (A/DEN) * p_up_left + (B/DEN) * p_up_right + (D/DEN) * p_down_left + &
&    (E/DEN) * p_down_right
!
c_out = (A/DEN) * c_up_left + (B/DEN) * c_up_right + (D/DEN) * c_down_left + &
&    (E/DEN) * c_down_right
!
x_out = (A/DEN) * x_up_left + (B/DEN) * x_up_right + (D/DEN) * x_down_left + &
&    (E/DEN) * x_down_right
!
END SUBROUTINE Lin_TP

END MODULE interp_functions
