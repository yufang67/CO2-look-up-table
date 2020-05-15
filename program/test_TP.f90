
PROGRAM test_TP
!      
      USE def_constants
      USE def_variables  
      USE non_linear_solvers
      USE properties
      USE grid_functions
      USE Grid      
      USE Interp_table
!      USE derivees
!      USE interp_functions
!      USE saturation
      IMPLICIT NONE
!======================================================================
!
!       Test points in two-phase region
!
!======================================================================
      REAL(pr) :: p, T, c, a, res, u_in, v_in,T_ref,v_l,v_v,&
&                 resnorm,in_2, vl_guess,vv_guess,e_v,e_l,p_guess,x_out,&
&                 delta
      INTEGER :: Niter, exitflag,i,j,NN
      REAL(pr), DIMENSION(500) :: T_tab,p_tab 
      REAL(pr), DIMENSION(500) :: x_tab
      REAL(pr) :: Perr,Terr,xerr       
!
      NN = 500
      delta = (T_cr-0.01_pr - (T_tri+0.01_pr))/ (NN)
      T_tab = T_tri+0.1_pr + (/(i*delta,i=0,NN-1)/)
      delta = 1.0_pr / (NN)
      x_tab = 0.0001_pr + (/(i*delta,i=0,NN-1)/)

!print*,x_tab
!
        CALL MAKE_GRID()
!STOP
      p=0_pr
      T=0_pr
      c=0_pr
      a=0_pr
      res=0_pr
!print*, "T_tab", T_tab
      vl_guess = 1_pr/rho_tri_L + 1e-4
      vv_guess = 1_pr/rho_tri_R - 1e-4
      T_ref = 260_pr

OPEN (UNIT = 22,             FILE = 'TP2_T-x-vV-p-u-v.txt', &
&   FORM = 'formatted',     ACTION = 'write',   &
&   STATUS = 'replace')
!
OPEN (UNIT = 42,             FILE = 'TP2_error-P-T-x.txt', &
&   FORM = 'formatted',     ACTION = 'write',   &
&   STATUS = 'replace')

        DO i=1,NN
!print*, i, T_tab(i)
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
           
        CALL CO2BLLT_EQUI(p, T, c, x_out, a, res, u_in,v_in,p_guess)
!print*,i,j, 'output',p,T,c,x_out   
        Perr = abs(p-p_tab(i))/p_tab(i)
        Terr = abs(T-T_tab(i))/T_tab(i)
        xerr = abs(x_out-x_tab(j))/x_tab(j)
!        xerr = 0.0
        p_guess = p

WRITE(22,*) T_tab(i), x_tab(j),v_v ,p,u_in,v_in,c,a 
WRITE(42,*) Perr, Terr,xerr,u_in,v_in, Terr*T_tab(i) 
    IF ( (Perr> 0.1_pr) .OR. (Terr >0.1) ) THEN
       print*, i,j, 'v', v_in, 'u', u_in
       print*, "error for T=", abs(T-T_tab(i))/T_tab(i), 'T_tb', T,'T_ref',T_tab(i) 
       print*, "error for p=", p,p_tab(i),abs(p-p_tab(i))/p_tab(i)
       print*," "
!       print*, "error for x=", x_tab(j), x_out,abs(x_tab(j)-x_out)/x_tab(j)
    ENDIF
          ENDDO
         ENDDO 
!!
 CLOSE(UNIT = 22)
 CLOSE(UNIT = 42)
print*," ==========>test_TP FINISH<=============="
END PROGRAM test_TP
