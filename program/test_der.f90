PROGRAM test_der
!
!
USE def_constants
USE def_variables
USE non_linear_solvers
USE Grid
USE properties
USE Interp_table
USE derivees
!
IMPLICIT NONE
!
INTEGER :: j, i, N, flag, k, p,domain
! Number of nodes for the test
!LL
INTEGER, PARAMETER ::  NNN_test_LL = NNN_LL , MMM_test_LL =  MMM_LL
!LH
INTEGER, PARAMETER ::  NNN_test_LH = NNN_LH , MMM_test_LH =  MMM_LH
!R
INTEGER, PARAMETER ::  NNN_test_R  = NNN_R  , MMM_test_R  =  MMM_R
!HT
INTEGER, PARAMETER ::  NNN_test_HT = NNN_HT , MMM_test_HT =  MMM_HT
!
REAL(8) :: dtest_u, dtest_v, res, time1, time2, p_guess, T_guess                
!
REAL(8) :: delta, a_out,out3
REAL(8) :: v_max_log!, v_sat_li, v_sat_left_log, v_sat_left, v_sat_log, &
!&          v_sat_right
!
!LL
REAL(8)                         :: d_y_test_LL, d_x_test_LL
REAL(8), DIMENSION(NNN_test_LL) :: y_test_LL
REAL(8), DIMENSION(MMM_test_LL) :: x_test_LL
!LH
REAL(8)                         :: d_y_test_LH, d_x_test_LH
REAL(8), DIMENSION(NNN_test_LH) :: y_test_LH
REAL(8), DIMENSION(MMM_test_LH) :: x_test_LH
!R
REAL(8)                         :: d_y_test_R, d_x_test_R
REAL(8), DIMENSION(NNN_test_R)  :: y_test_R
REAL(8), DIMENSION(MMM_test_R)  :: x_test_R
!HT
REAL(8)                         :: d_y_test_HT, d_x_test_HT
REAL(8), DIMENSION(NNN_test_HT) :: y_test_HT
REAL(8), DIMENSION(MMM_test_HT) :: x_test_HT
!
! Matrixes declaration for each region
!LL
REAL(8), DIMENSION(NNN_test_LL, MMM_test_LL) :: v_test_LL_log,v_test_LL,&
&  dT_dv_LL, dT_du_LL, dp_dv_LL, dp_du_LL,c_BL_LL,p_BL_LL,T_BL_LL,x_out_LL
!
!LH
REAL(8), DIMENSION(NNN_test_LH, MMM_test_LH) :: v_test_LH_log,v_test_LH, &
&  dT_dv_LH, dT_du_LH, dp_dv_LH, dp_du_LH,c_BL_LH,p_BL_LH,T_BL_LH,x_out_LH
!
!R
REAL(8), DIMENSION(NNN_test_R, MMM_test_R) :: v_test_R_log, v_test_R, &
&  dT_dv_R, dT_du_R, dp_dv_R, dp_du_R,c_BL_R,p_BL_R,T_BL_R,x_out_R
!HT
REAL(8), DIMENSION(NNN_test_HT, MMM_test_HT) :: v_test_HT_log, v_test_HT, &
&  dT_dv_HT, dT_du_HT, dp_dv_HT, dp_du_HT,c_BL_HT,p_BL_HT,T_BL_HT,x_out_HT
   
REAL(8) :: v_min_li, v_max_li, out_2, p_TP, res2,res3,res4,res1

REAL(8) :: press, temp, sound, u_in, v_in,T_ref,v_l,v_v,&
&                 resnorm,in_2, vl_guess,vv_guess,e_v,e_l,x_out,&
&                 mass
INTEGER :: Niter, exitflag,NN
REAL(8), DIMENSION(500) :: T_tab,p_tab
REAL(8), DIMENSION(500) :: x_tab
REAL(8) :: Perr,Terr,xer
REAL(8) :: dT_dv_TP, dT_du_TP, dp_dv_TP, dp_du_TP
!
!
!
!
!====================================================
! CONSTRUCTION OF THE GRID ON THE (u,v) DIAGRAM
!====================================================
!
CALL cpu_time(time1)
CALL MAKE_GRID() 
CALL cpu_time(time2)
print*, '===================================================='
print*, 'TIME INTERVAL FOR GRID CONSTRUCTION = ', time2-time1
print*, '===================================================='
!
! Introducing test values:
!
dtest_u      = 5d1
dtest_v      = 5d-4
domain       = 1

SELECT CASE(domain)
CASE(1)
!
!
! Left Low LL
!
y_test_LL(1)           = y_mesh_LL(1) + dtest_u
y_test_LL(NNN_test_LL) = y_mesh_LL(NNN_LL) - dtest_u
d_y_test_LL            = (y_test_LL(NNN_test_LL) - y_test_LL(1))/(NNN_test_LL - 1)
y_test_LL              = y_test_LL(1) + (/(j*d_y_test_LL, j=0,NNN_test_LL-1)/)

x_test_LL(1)           = x_mesh_LL(1) + dtest_v
x_test_LL(MMM_test_LL) = x_mesh_LL(MMM_LL) - dtest_v
d_x_test_LL            = (x_test_LL(MMM_test_LL) - x_test_LL(1))/(MMM_test_LL - 1)      
x_test_LL              = x_test_LL(1) + (/(j*d_x_test_LL, j=0,MMM_test_LL-1)/)
!
! I need a different mesh to test some points for each cell
!
!
OPEN (UNIT = 41,             FILE = 'dp_du_LL.txt', &
&     FORM = 'formatted',  ACTION = 'write',         &
&   STATUS = 'replace')
!
OPEN (UNIT = 57,             FILE = 'dp_dv_LL.txt', &
&     FORM = 'formatted',  ACTION = 'write',         &
&   STATUS = 'replace')
!
!OPEN (UNIT = 63,             FILE = 'dT_du_LL.txt', &
!&     FORM = 'formatted',  ACTION = 'write',         &
!&   STATUS = 'replace')
!
!OPEN (UNIT = 64,             FILE = 'dT_dv_LL.txt', &
!&     FORM = 'formatted',  ACTION = 'write',         &
!&   STATUS = 'replace')
!
!
DO k = 1, NNN_test_LL-1  !For each cell of the original grid to test
   delta = y_mesh_LL(2)  - y_mesh_LL(1)
   i = INT((y_test_LL(k) - y_mesh_LL(1))/delta) + 1  !I found the i-level of the tested point
   v_min_li = 0d0
   v_max_li = 0d0
   DO j = 1, ord_spline + 1
      v_min_li = v_min_li + spline_pmax_LL(ord_spline+2-j,i) * y_test_LL(k) ** (j - 1)
!      v_max_log  = v_max_log + spline_Lsat_LL(ord_spline+2-j,i) * y_test_LL(k) ** (j - 1)
!      v_max_li   = 10d0 ** v_max_log
      v_max_li = v_max_li + spline_Lsat_LL(ord_spline+2-j,i) * y_test_LL(k) ** (j - 1)  !I found the min and max volumes correspondent to that energy level
   ENDDO
   DO p = 1, MMM_test_LL-1
   delta = x_mesh_LL(2)  - x_mesh_LL(1)
     j = INT((x_test_LL(p) - x_mesh_LL(1))/delta) + 1   !I found the j-level of the tested point
     v_test_LL(k,p) = ((x_test_LL(p) - x_mesh_LL(1))/(x_mesh_LL(MMM_LL) - x_mesh_LL(1))) * &
&                        (v_max_li - v_min_li) + v_min_li
!     v_test_LL_log(k,p) = ((x_test_LL(p) - x_mesh_LL(1))/(x_mesh_LL(MMM_LL) - x_mesh_LL(1))) * &
!&                        (LOG10(v_max_li) - LOG10(v_min_li)) + LOG10(v_min_li)
!     v_test_LL(k,p) = 10d0 ** v_test_LL_log(k,p)

!
!  CO2 BI-LINEAR LOOK-UP TABLES
       CALL CO2BLLT_EQUI(p_BL_LL(k,p), T_BL_LL(k,p), c_BL_LL(k,p), &
&                        x_out_LL(k,p), a_out, res, y_test_LL(k), v_test_LL(k,p),0.0_pr)

      CALL CO2DER(dp_dv_LL(k,p), dp_du_LL(k,p),y_test_LL(k), v_test_LL(k,p), T_BL_LL(k,p),p_BL_LL(k,p),res2,res3,res4)
!
!
         WRITE(41,*) y_test_LL(k), v_test_LL(k,p), dp_du_LL(k,p)
         WRITE(57,*) y_test_LL(k), v_test_LL(k,p), dp_dv_LL(k,p)
!         WRITE(63,*) y_test_LL(k), v_test_LL(k,p), dT_du_LL(k,p)
!         WRITE(64,*) y_test_LL(k), v_test_LL(k,p), dT_dv_LL(k,p) 
!               
!
   ENDDO
ENDDO   !For all the tested values NNN_test x MMM_test in the Region
   CLOSE(UNIT = 41)
   CLOSE(UNIT = 57)
!   CLOSE(UNIT = 63)
!   CLOSE(UNIT = 64)
   
!
!----------------------------------------------------------------------
print*, 'LL Test finished'
!----------------------------------------------------------------------
STOP
CASE(2)
!======================================================
! Left High
!======================================================
!
y_test_LH(1)           = y_mesh_LH(1) + dtest_u
y_test_LH(NNN_test_LH) = y_mesh_LH(NNN_LH) - dtest_u
d_y_test_LH            = (y_test_LH(NNN_test_LH) - y_test_LH(1))/(NNN_test_LH - 1)

y_test_LH              = y_test_LH(1) + (/(j*d_y_test_LH, j=0,NNN_test_LH-1)/)

x_test_LH(1)           = x_mesh_LH(1) + dtest_v
x_test_LH(MMM_test_LH) = x_mesh_LH(MMM_LH) - dtest_v
d_x_test_LH            = (x_test_LH(MMM_test_LH) - x_test_LH(1))/(MMM_test_LH - 1)
      
x_test_LH              = x_test_LH(1) + (/(j*d_x_test_LH, j=0,MMM_test_LH-1)/)
!
! I need a different mesh to test some points for each cell
!
!
OPEN (UNIT = 41,             FILE = 'dp_du_LH.txt', &
&     FORM = 'formatted',  ACTION = 'write',         &
&   STATUS = 'replace')
!
OPEN (UNIT = 57,             FILE = 'dp_dv_LH.txt', &
&     FORM = 'formatted',  ACTION = 'write',         &
&   STATUS = 'replace')
!
!OPEN (UNIT = 64,             FILE = 'dT_dv_LH.txt', &
!&     FORM = 'formatted',  ACTION = 'write',         &
!&   STATUS = 'replace')
!
!OPEN (UNIT = 63,             FILE = 'dT_du_LH.txt', &
!&     FORM = 'formatted',  ACTION = 'write',         &
!&   STATUS = 'replace')
!
!
DO k = 1, NNN_test_LH-1  !For each cell of the original grid to test
   delta = y_mesh_LH(21)  - y_mesh_LH(20) 
   i     = INT((y_test_LH(k) - y_mesh_LH(1))/delta) + 1  !I found the i-level of the tested point
   v_min_li    = 0d0
   v_max_li    = 0d0
   v_max_log   = 0d0

   DO j = 1, ord_spline + 1
      v_min_li   = v_min_li   + spline_pmax_LH(ord_spline+2-j,i) * y_test_LH(k) ** (j - 1)
      v_max_li   = v_max_li   + spline_Lsat_LH(ord_spline+2-j,i) * y_test_LH(k) ** (j - 1)     
!      v_max_log  = v_max_log + spline_Lsat_LH(ord_spline+2-j,i) * y_test_LH(k) ** (j - 1)
!      v_max_li   = 10d0 ** v_max_log
!      print*, spline_pmax_LH(ord_spline+2-j,i)
   ENDDO

!   stop
   DO p = 1, MMM_test_LH-1

   delta = x_mesh_LH(2)  - x_mesh_LH(1)
     j = INT((x_test_LH(p) - x_mesh_LH(1))/delta) + 1   !I found the j-level of the tested point
!     v_test_LH_log(k,p) = ((x_test_LH(p) - x_mesh_LH(1))/(x_mesh_LH(MMM_LH) - x_mesh_LH(1))) * &
!&                        (LOG10(v_max_li) - LOG10(v_min_li)) + LOG10(v_min_li)
!     v_test_LH(k,p)     = 10d0 ** v_test_LH_log(k,p)
     v_test_LH(k,p) = ((x_test_LH(p) - x_mesh_LH(1))/(x_mesh_LH(MMM_LH) - x_mesh_LH(1))) * &
&                        (v_max_li - v_min_li) + v_min_li
     
     
     IF ( v_test_LH(k,p) .GE. v_max_li) THEN
        print*, 'Error in meshtest', k, p, y_test_LH(k), v_max_li, v_test_LH(k,p)!----->Problem in spline spin_LH
        stop
     ENDIF   
!        print*, ''
!   print*, v_min_li,v_test_LH(k,p), v_max_li
!       
!        print*, k, p, v_test_LH(k,p)
!  CO2 BI-LINEAR LOOK-UP TABLES
!     print*, y_test_LH(k),v_test_LH(k,p),p
!
        CALL CO2BLLT_EQUI (p_BL_LH(k,p), T_BL_LH(k,p), c_BL_LH(k,p), &
&                          x_out_LH(k,p), a_out, res, y_test_LH(k), v_test_LH(k,p),0.0_pr)
        CALL CO2DER(dp_dv_LH(k,p), dp_du_LH(k,p),y_test_LH(k), v_test_LH(k,p), T_BL_LH(k,p),p_BL_LH(k,p),res2,res3,res4)
!
!
         WRITE(41,*) y_test_LH(k), v_test_LH(k,p), dp_du_LH(k,p)
         WRITE(57,*) y_test_LH(k), v_test_LH(k,p), dp_dv_LH(k,p)
!         WRITE(63,*) y_test_LH(k), v_test_LH(k,p), dT_du_LH(k,p)
!         WRITE(64,*) y_test_LH(k), v_test_LH(k,p), dT_dv_LH(k,p) 
!               
!
   ENDDO
ENDDO   !For all the tested values NNN_test x MMM_test in the Region
   CLOSE(UNIT = 41)
   CLOSE(UNIT = 57)
!   CLOSE(UNIT = 63)
!   CLOSE(UNIT = 64)
!
!----------------------------------------------------------------------
print*, 'LH Test finished'
!----------------------------------------------------------------------
STOP
CASE(3)
!======================================================
! Right
!==========================================
!
y_test_R(1)           = y_mesh_R(1) + dtest_u
y_test_R(NNN_test_R)  = y_mesh_R(NNN_R) - dtest_u
d_y_test_R            = (y_test_R(NNN_test_R) - y_test_R(1))/(NNN_test_R - 1)
y_test_R              = y_test_R(1) + (/(j*d_y_test_R, j=0,NNN_test_R-1)/)
x_test_R(1)           = x_mesh_R(1) + dtest_v
x_test_R(MMM_test_R)  = x_mesh_R(MMM_R) - dtest_v
d_x_test_R            = (x_test_R(MMM_test_R) - x_test_R(1))/(MMM_test_R - 1)      
x_test_R              = x_test_R(1) + (/(j*d_x_test_R, j=0,MMM_test_R-1)/)
!
! I need a different mesh to test some points for each cell
!
!
OPEN (UNIT = 41,             FILE = 'dp_du_R.txt', &
&     FORM = 'formatted',  ACTION = 'write',         &
&   STATUS = 'replace')
!
OPEN (UNIT = 57,             FILE = 'dp_dv_R.txt', &
&     FORM = 'formatted',  ACTION = 'write',         &
&   STATUS = 'replace')
!
!OPEN (UNIT = 58,             FILE = 'dT_du_R.txt', &
!&     FORM = 'formatted',  ACTION = 'write',         &
!&   STATUS = 'replace')
!
!OPEN (UNIT = 64,             FILE = 'dT_dv_R.txt', &
!&     FORM = 'formatted',  ACTION = 'write',         &
!&   STATUS = 'replace')
!
!
DO k = 1, NNN_test_R-1  !For each cell of the original grid to test
   delta = y_mesh_R(2)  - y_mesh_R(1) 
   i     = INT((y_test_R(k) - y_mesh_R(1))/delta) + 1  !I found the i-level of the tested point
   v_min_li    = 0d0
   v_max_li    = 0d0

   DO j = 1, ord_spline + 1
      v_max_li   = v_max_li + spline_pmin(ord_spline+2-j,i) * y_test_R(k) ** (j - 1)
      v_min_li   = v_min_li + spline_Vsat(ord_spline+2-j,i) * y_test_R(k) ** (j - 1)  !I found the min and max volumes correspondent to that energy level

   ENDDO
!   print*, ''
!   print*, v_min_li, v_max_li
!   stop
   DO p = 1, MMM_test_R-1
   delta = x_mesh_R(2)  - x_mesh_R(1)
     j = INT((x_test_R(p) - x_mesh_R(1))/delta) + 1   !I found the j-level of the tested point
!     v_test_R_log(k,p) = ((x_test_R(p) - x_mesh_R(1))/(x_mesh_R(MMM_R) - x_mesh_R(1))) * &
!&                        (LOG10(v_max_li) - LOG10(v_min_li)) + LOG10(v_min_li)
!     v_test_R(k,p)     = 10d0 ** v_test_R_log(k,p)

     v_test_R(k,p) = ((x_test_R(p) - x_mesh_R(1))/(x_mesh_R(MMM_R) - x_mesh_R(1))) * &
&                    (v_max_li - v_min_li) + v_min_li
     
     
     IF (( v_test_R(k,p) .LT. v_min_li) .OR. ( v_test_R(k,p) .GT. v_max_li))  THEN
        print*, 'Error in meshtest', k, p, y_test_R(k), v_min_li, v_test_R(k,p)!----->Problem in spline spin_R
        stop
     ENDIF       
!
!  CO2 BI-LINEAR LOOK-UP TABLES
        CALL CO2BLLT_EQUI (p_BL_R(k,p), T_BL_R(k,p), c_BL_R(k,p), &
&                          x_out_R(k,p), a_out, res, y_test_R(k), v_test_R(k,p),0.0_pr)

        CALL CO2DER(dp_dv_R(k,p), dp_du_R(k,p),y_test_R(k), v_test_R(k,p),T_BL_R(k,p),p_BL_R(k,p),res2,res3,res4)
!
         WRITE(41,*) y_test_R(k), v_test_R(k,p), dp_du_R(k,p)
         WRITE(57,*) y_test_R(k), v_test_R(k,p), dp_dv_R(k,p)
!         WRITE(58,*) y_test_R(k), v_test_R(k,p), dT_du_R(k,p)
!         WRITE(64,*) y_test_R(k), v_test_R(k,p), dT_dv_R(k,p) 
!               
!
   ENDDO
ENDDO   !For all the tested values NNN_test x MMM_test in the Region
   CLOSE(UNIT = 41)
   CLOSE(UNIT = 57)
!   CLOSE(UNIT = 58)
!   CLOSE(UNIT = 64)
!
!
!----------------------------------------------------------------------
print*, 'R Test finished'
!----------------------------------------------------------------------
STOP
CASE(4)
!======================================================
! High Temperature
!==========================================
!
y_test_HT(1)           = y_mesh_HT(1) + dtest_u
y_test_HT(NNN_test_HT) = y_mesh_HT(NNN_HT) - dtest_u
d_y_test_HT            = (y_test_HT(NNN_test_HT) - y_test_HT(1))/(NNN_test_HT - 1)
y_test_HT              = y_test_HT(1) + (/(j*d_y_test_HT, j=0,NNN_test_HT-1)/)
x_test_HT(1)           = x_mesh_HT(1) + dtest_v
x_test_HT(MMM_test_HT)  = x_mesh_HT(MMM_HT) - dtest_v
d_x_test_HT            = (x_test_HT(MMM_test_HT) - x_test_HT(1))/(MMM_test_HT - 1)      
x_test_HT              = x_test_HT(1) + (/(j*d_x_test_HT, j=0,MMM_test_HT-1)/)
!
! I need a different mesh to test some points for each cell
!
!
OPEN (UNIT = 41,             FILE = 'dp_du_HT.txt', &
&     FORM = 'formatted',  ACTION = 'write',         &
&   STATUS = 'replace')
!
OPEN (UNIT = 57,             FILE = 'dp_dv_HT.txt', &
&     FORM = 'formatted',  ACTION = 'write',         &
&   STATUS = 'replace')
!
!OPEN (UNIT = 58,             FILE = 'dT_du_HT.txt', &
!&     FORM = 'formatted',  ACTION = 'write',         &
!&   STATUS = 'replace')
!
!OPEN (UNIT = 64,             FILE = 'dT_dv_HT.txt', &
!&     FORM = 'formatted',  ACTION = 'write',         &
!&   STATUS = 'replace')
!
!
!
DO k = 1, NNN_test_HT-1  !For each cell of the original grid to test
   delta = y_mesh_HT(2)  - y_mesh_HT(1) 
   i     = INT((y_test_HT(k) - y_mesh_HT(1))/delta) + 1  !I found the i-level of the tested point
   v_min_li    = 0d0
   v_max_li    = 0d0

   DO j = 1, ord_spline + 1
      v_max_li   = v_max_li + spline_right_HT(ord_spline+2-j,i) * y_test_HT(k) ** (j - 1)
      v_min_li   = v_min_li + spline_left_HT(ord_spline+2-j,i) * y_test_HT(k) ** (j - 1)  !I found the min and max volumes correspondent to that energy level

   ENDDO
!   print*, ''
!   print*, v_min_li
!   stop
   DO p = 1, MMM_test_HT-1
   delta = x_mesh_HT(2)  - x_mesh_HT(1)
     j = INT((x_test_HT(p) - x_mesh_HT(1))/delta) + 1   !I found the j-level of the tested point

!     v_test_HT(k,p) = ((x_test_HT(p) - x_mesh_HT(1))/(x_mesh_HT(MMM_HT) - x_mesh_HT(1))) * &
!&                        (v_max_li - v_min_li) + v_min_li



     v_test_HT_log(k,p) = ((x_test_HT(p) - x_mesh_HT(1))/(x_mesh_HT(MMM_HT) - x_mesh_HT(1))) * &
&                        (LOG10(v_max_li) - LOG10(v_min_li)) + LOG10(v_min_li)
     v_test_HT(k,p)     = 10d0 ** v_test_HT_log(k,p)
     
     
     IF (( v_test_HT(k,p) .LT. v_min_li) .OR. ( v_test_HT(k,p) .GT. v_max_li))  THEN
        print*, 'Error in meshtest', k, p, y_test_HT(k), v_min_li, v_test_HT(k,p)!----->Problem in spline spin_HT
        stop
     ENDIF       
!
!  CO2 BI-LINEAR LOOK-UP TABLES
     CALL CO2BLLT_EQUI (p_BL_HT(k,p), T_BL_HT(k,p), c_BL_HT(k,p), &
&                       x_out_HT(k,p), a_out, res, y_test_HT(k), v_test_HT(k,p),0.0_pr)
     CALL CO2DER( dp_dv_HT(k,p), dp_du_HT(k,p),y_test_HT(k), v_test_HT(k,p),T_BL_HT(k,p),p_BL_HT(k,p),res2,res3,res4)

!
         WRITE(41,*) y_test_HT(k), v_test_HT(k,p), dp_du_HT(k,p)
         WRITE(57,*) y_test_HT(k), v_test_HT(k,p), dp_dv_HT(k,p)
!         WRITE(58,*) y_test_HT(k), v_test_HT(k,p), dT_du_HT(k,p)
!         WRITE(64,*) y_test_HT(k), v_test_HT(k,p), dT_dv_HT(k,p) 
!               
!
   ENDDO
ENDDO   !For all the tested values NNN_test x MMM_test in the Region
   CLOSE(UNIT = 41)
   CLOSE(UNIT = 57)
!  CLOSE(UNIT = 58)
!  CLOSE(UNIT = 64)
!
!----------------------------------------------------------------------
print*, 'HT Test finished'
!----------------------------------------------------------------------
STOP
CASE(5)
!
      NN = 100
      delta = (T_cr-0.1_pr - (T_tri+0.01_pr))/ (NN)
      T_tab = T_tri+0.01_pr + (/(i*delta,i=0,NN-1)/)
      delta = 1.0_pr / (NN)
      x_tab = 0.0001_pr + (/(i*delta,i=0,NN-1)/)

!print*,x_tab
!
      press=0_pr
      temp=0_pr
      sound=0_pr
      mass=0_pr
      res1=0_pr
!print*, "T_tab", T_tab
      vl_guess = 1_pr/rho_tri_L + 1e-4
      vv_guess = 1_pr/rho_tri_R - 1e-4
      T_ref = 260_pr
!
OPEN (UNIT = 22,             FILE = 'deriv_TP.txt', &
&   FORM = 'formatted',     ACTION = 'write',   &
&   STATUS = 'replace')
!
        DO i=1,NN
!        print*, i, T_tab(i)
        CALL New_Rap2D(1, v_l,v_v, &
     &           resnorm, Niter, exitflag, T_tab(i), in_2, vl_guess, vv_guess)
!
       IF (resnorm > 1e-10) THEN
                print*, "res", resnorm, "iter", Niter
                STOP
              ENDIF
!
        vl_guess = v_l
        vv_guess = v_v
        CALL inter_energy(T_tab(i),v_v,e_v)
        CALL inter_energy(T_tab(i),v_l,e_l)
        CALL pressure(T_tab(i),v_v,p_tab(i))!
!print*, T_tab(i),p_tab(i)
!print*,i, v_l, v_v
!print*,i, e_v, e_l
!!    print*,"error v", (v_l-1_pr/998.89_pr)/(1_pr/998.89_pr),&
!!&           (v_v-1_pr/64.417_pr)/(1_pr/64.417_pr)
            p_guess = p_tab(i)
          DO j=1,NN
 ! print*,i,j
            v_in = x_tab(j)*v_v + (1_pr - x_tab(j))* v_l
            u_in = x_tab(j)*e_v + (1_pr - x_tab(j))* e_l
!   print*,x_tab(j), T_tab(i)!,v_v,v_l
!print*, i,j,"v_in",v_in,"u_in",u_in, p_guess

        CALL CO2BLLT_EQUI(press, temp, sound, x_out, a_out, res, u_in,v_in,p_guess)
        p_guess = press
        CALL CO2DER(dp_dv_TP, dp_du_TP,u_in, v_in,temp, press,res2,res3,res4)



        WRITE(22,*) T_tab(i), press,dp_dv_TP, dp_du_TP, v_v ,p,u_in,v_in
          ENDDO
         ENDDO
!!
 CLOSE(UNIT = 22)
print*," ==========>test_TP FINISH<=============="


END SELECT
!
END PROGRAM test_der
