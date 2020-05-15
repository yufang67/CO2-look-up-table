
PROGRAM test_3D
!      
      USE def_constants
      USE def_variables  
      USE non_linear_solvers
      USE properties
      USE grid_functions
      USE Grid      
      USE Interp_table
      IMPLICIT NONE
!      
!     INTEGER :: Niter, exitflag      
!      REAL(pr):: x, y, z, out1, out2, out3, resnorm,in_1, in_2
!      REAL(pr):: p,v,u,T,c
!========================================================================
!
!           Test for the 3D nonlinear solver
!
!========================================================================
!        x = -1.0_pr
!        y = 3.14_pr / 1.1_pr
!        z = -1_pr
!        resnorm = 1e-10_pr
!        Niter = 100

!        CALL  New_Rap3D(2, out1, out2, out3, &
!       &       resnorm, Niter, exitflag, 1.1_pr, x, y, z)

!    print*, "---------------------------------------------------"
!    print*, "x = ", out1, "y = ", out2, "z = ", out3
!    print*, "---------------------------------------------------"
!
!
!======================================================================
!
!               Test for the Max_point routine
!
!======================================================================
!
!        Niter = 5
!        x = 2_pr
!print*, x/Niter
!        u = 0_pr
!        v = 0_pr
!        T = 0_pr
!        p = 0_pr
!        c = 0_pr
!          CALL Max_point(u,v,T,p,c)
!
!    print*, "==================================================="
!    print*, "e_max =", u*1e-3_pr, v, 1_pr/v,c
!    print*, "---------------------------------------------------"
!    print*, "T_umax =", T, "P_umax =", p*1e-6_pr
!    print*, "==================================================="
!======================================================================
!
!              Test corner
!
!=====================================================================
!
!  
      CALL  corner_A()
      CALL  corner_B()
      CALL  corner_C()
!       CALL   corner_D()
!
!     STOP
!    print*, "---------------------------------------------------"
    print*, "v=", v_corner_A,v_corner_B,v_corner_C
    print*, "T=", T_corner_A,T_corner_B,T_corner_C
!    print*, "---------------------------------------------------"
!
!======================================================================
!
!              Test location of interpolation
!
!=====================================================================
!
!
!      REAL(pr) :: p, T, c, x, a, res, u_in, v_in,p_ref,T_ref
!
!      p=0_pr
!      T=0_pr
!      c=0_pr
!      x=0_pr
!      a=0_pr
!      res=0_pr
!
!      u_in = -126.27e3_pr
!      v_in = 1_pr/261.29_pr
!      T_ref = 315_pr
!
!      CALL  MAKE_GRID()
!      CALL  CO2BLLT_EQUI(p,T,c,x,a,res,u_in,v_in)
!      CALL  pressure(T_ref,v_in,p_ref)
!
!    print*, "---------------------------------------------------"
!    print*, "p", (p_ref-p)/p_ref
!    print*, "T", (T-T_ref)/T_ref
!    print*, "---------------------------------------------------"
!======================================================================
!
!       Test points in two-phase region
!
!======================================================================
!      REAL(pr) :: p, T, c, a, res, u_in, v_in,T_ref,v_l,v_v,&
!&                 resnorm,in_2, vl_guess,vv_guess,e_v,e_l,p_guess,x_out,&
!&                 delta
!      INTEGER :: Niter, exitflag,i,j
!!      REAL(pr), DIMENSION(500) :: T_tab,p_tab 
!      REAL(pr), DIMENSION(500) :: x_tab
!      REAL(pr) :: Perr,Terr,xerr       
!!
!      delta = (T_cr-0.1_pr - (T_tri+0.01_pr))/ (500_pr)
!      T_tab = T_tri+0.01_pr + (/(i*delta,i=0,500-1)/)
!      delta = 1.0_pr / (500_pr)
!      x_tab = 0.0001_pr + (/(i*delta,i=0,500-1)/)
!
!print*,x_tab
!
!        CALL MAKE_GRID()
!      p=0_pr
!      T=0_pr
!      c=0_pr
!      a=0_pr
!      res=0_pr
!!print*, "T_tab", T_tab
!      vl_guess = 1_pr/rho_tri_L + 1e-4
!      vv_guess = 1_pr/rho_tri_R - 1e-4
!      T_ref = 260_pr
!
!OPEN (UNIT = 22,             FILE = 'TP_T-x-Perr-p-u-v.txt', &
!&   FORM = 'formatted',     ACTION = 'write',   &
!&   STATUS = 'replace')
!!
!OPEN (UNIT = 42,             FILE = 'TP_error-P-T-x.txt', &
!&   FORM = 'formatted',     ACTION = 'write',   &
!&   STATUS = 'replace')
!
!        DO i=1,500
!        CALL New_Rap2D(1, v_l,v_v, &
!     &           resnorm, Niter, exitflag, T_tab(i), in_2, vl_guess, vv_guess)
!!
!       IF (resnorm > 1e-10) THEN
!                print*, "res", resnorm, "iter", Niter
!                STOP
!              ENDIF
!
!        vl_guess = v_l
!        vv_guess = v_v 
!        CALL inter_energy(T_tab(i),v_v,e_v)
!        CALL inter_energy(T_tab(i),v_l,e_l)
!        CALL pressure(T_tab(i),v_v,p_tab(i))
!!print*, T_tab(i),p_tab(i)
!!print*, v_l, v_v
!!print*, e_v, e_l
!!!    print*,"error v", (v_l-1_pr/998.89_pr)/(1_pr/998.89_pr),&
!!!&           (v_v-1_pr/64.417_pr)/(1_pr/64.417_pr)
!!            p_guess = p_tab(i)
!          DO j=1,500
! ! print*,i,j
!            v_in = x_tab(j)*v_v + (1_pr - x_tab(j))* v_l
!            u_in = x_tab(j)*e_v + (1_pr - x_tab(j))* e_l
!!   print*,x_tab(j), T_tab(i)!,v_v,v_l
!!print*, "v_in",v_in,"u_in",u_in, p_guess
!           
!        CALL CO2BLLT_EQUI(p, T, c, x_out, a, res, u_in,v_in,p_guess)
!!print*, 'output',p,T,c,x_out   
!        Perr = abs(p-p_tab(i))/p_tab(i)
!        Terr = abs(T-T_tab(i))/T_tab(i)
!        xerr = abs(x_out-x_tab(j))/x_tab(j)
!        xerr = 0.0
!        p_guess = p
!
!WRITE(22,*) T_tab(i), x_tab(j),v_v ,p,u_in,v_in,c,a 
!WRITE(42,*) Perr, Terr,xerr,u_in,v_in, abs(T-T_tab(i))
!
!
!!   print*, "error for T=", abs(T-T_tab(i))/T_tab(i) 
!!   print*, "error for p=", p,p_tab(i),abs(p-p_tab(i))/p_tab(i)
!!   print*, "error for x=", x_tab(j), x_out,abs(x_tab(j)-x_out)/x_tab(j)
!          ENDDO
!!print*,"---------------------------------------------------"       
!         ENDDO 
!!!
! CLOSE(UNIT = 22)
! CLOSE(UNIT = 42)
!print*," ==========>test_3D FINISH<=============="
!print*, "---------------------------------------------------"
! !   print*, "T_ref=", T_ref,"T_table",T
! !   print*, "p_ref=", p_ref,"p_table",p
! !   print*, "x_ref=", x,"x_table",x_out
! !   print*, "error for T=", abs(T-T_ref)/T_ref 
! !  print*, "error for p=", abs(p-p_ref)/p_ref
! !  print*, "error for x=", abs(x-x_out)/x
! !  print*, "---------------------------------------------------"
END PROGRAM test_3D
