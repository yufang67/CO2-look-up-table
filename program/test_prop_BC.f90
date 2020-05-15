PROGRAM test_pressure

! Short test program for the function: 'Helmholtz.f90'.
!
! Pressure corresponds to: -(df / dv)_T.
!
! Let us check if pressure values from tables at the end of the paper
! (from page 1562 on) are obtained using 'Helmholtz.f90'.
!
! The derivative is approximated by a 2nd order finite difference.
!
      USE def_constants
      USE Grid
      USE non_linear_solvers
      USE properties
      USE solver_eos
      IMPLICIT NONE
!
      INTEGER :: Niter, exitflag      
      REAL(pr):: T, v, p, e, cv, cp, s, c, v_l, v_v, v_l_Span, v_v_Span, &
     &           resnorm, guess_1, guess_2, T_span, v_span, e_span,helmho1
      REAL(pr):: p2, e2, cv2, cp2, s2, c2,helmho2,Helmholtz,eguess,out3,out2,press
      REAL(pr):: T1,v1,s_entro1,e1,p1,u1,press1,h1
      REAL(pr):: T2,v2,s_entro2,u2,press2,h2
!
! TEST 1 pressure, e, cv, cp, c, s, helmholtz
!     v = 0.00124937  ! corresponding to 1 MPa (Table 35, page 1568)
!      v = 1_pr/ 801.62_pr
     T = 10+273.15!294.95!309.65!308.95!299.85!310.15
     p = 3.0e6!7.1e6 !9.1e6
!     eguess = -262219.0
     eguess = 1.0/900!524.06
!!!      CALL property (T,v,p,e,cv,cp,s,c)
!!!      helmho1 = Helmholtz(v,T)
     T1= 40+273.15
     v1= 1/485.5
     p1= 9.0e6
     u1= 0 

     T2= 10.0+273.15
     v2= 1/161.8
     p2= 6.1e6
     u2= 0
!
!!!        print*, 'test pressure -- T,v', T,v
!     CALL pressure(T,v,p2)
!      CALL entropy(T,v,s2)
!      CALL heat_cap_v(T,v,cv2)
!      CALL heat_cap_p(T,v,cp2)
!      CALL inter_energy(T,v,e2)
!      CALL sound_speed(T,v,c2)
!      CALL helmho(T,v,helmho2)

!       CALL MAKE_GRID()
!       CALL eos_1d(3, e, e2, resnorm, Niter,&
!     &             exitflag, p, eguess, v, p2) 

       CALL New_Rap1D(3, v, out2, resnorm, Niter,&
     &                             exitflag, p, eguess, T, out3)     
       CALL inter_energy(T,v,e2)
       CALL pressure(T,v,press)
      print*, 'v',v,'rho', 1.0/v
      print*, 'e ',e2
      print*, '---------------------------------------'
      print*, 'p in MPa', press*1e-6


       CALL entropy(T1,v1,s_entro1)
       CALL inter_energy(T1,v1,e1)
       CALL pressure(T1,v1,press1)
       h1 = e1+p1*v1


       CALL entropy(T2,v2,s_entro2)
       CALL inter_energy(T2,v2,e2)
       CALL pressure(T2,v2,press2)
       h2 = e2+p2*v2


      print*,'------------------Primary inlet-------------------------------------'
      print*, 's1= ', s_entro1,' h1= ', h1, "p1= ", press1, " ",p1
      print*,'------------------Secondary inlet----------------------------------------'
      print*, 's2= ', s_entro2,' h2= ', h2, "p2= ", press2, " ",p2
      print*,'----------------------------------------------------------'
      print*, 'exergy primary inlet= ', (h1-h2)+0.5*(u1*u1-u2*u2)-T2*(s_entro1-s_entro2)
      print*,'----------------------------------------------------------'

!      print*, ' cv ', cv2*1e-3_pr, ' error ', cv2*1e-3_pr - 0.68217_pr 
      print*,'----------------------------------------------------------'
!      print*, ' cp  ' , cp2*1e-3_pr,' error ', cp2*1e-3_pr - 0.92089_pr
      print*,'----------------------------------------------------------'
!      print*, ' s  ' , s2*1e-3_pr,' error ', s2*1e-3_pr + 0.44964
      print*,'----------------------------------------------------------'
!      print*, ' c  ' , c2,' error ', c2 - 262.43_pr
      print*,'----------------------------------------------------------'
!!      print*, ' helmholtz  ' , helmho2,' error ', helmho1*1e3_pr-helmho2
!!      print*,'----------------------------------------------------------'
!             
! TEST 3 - MARCO   saturation curve
!
!     T = 300_pr                                  ! K
! from table 34, page 1561. 
! At 300 K the saturated specific volumes are:
!     v_l_Span = 1_pr / 679.24_pr 
!     v_v_Span = 1_pr / 268.58_pr
!     
!     guess_1 =  v_l_Span*0.8_pr 
!     guess_2 =  v_v_Span*1.2_pr 
!
!
!      CALL New_Rap2D(1, v_l, v_v, &
!     &           resnorm, Niter, exitflag, T, 0d0, guess_1, guess_2)
!     
!      CALL pressure(T,v_l,p)
!
!      
      print*, '----------------TEST 3-----------------------'
!      print*, 'resnorm, Niter, exitflag', resnorm, Niter, exitflag
      print*,'-------------------------------------------------------------------------------'
!      print*, 'GUESS:        rho_l = ', 1_pr/guess_1,  ' rho_v= ',1_pr/guess_2
!      print*, 'From article: rho_l = ', 1_pr/v_l_Span, ' rho_v= ',1_pr/v_v_Span
!      print*, 'Newton:       rho_l = ', 1_pr/v_l,      ' rho_v= ',1_pr/v_v
      print*,'--------------------------------------------------------------------------------'
!      print*, 'pressure in MPa, test 3', p*1e-6_pr
!             
! TEST 4 - MARCO    calculation of p and e corresponding to assigned T and V
!
! from table 34, page 1567. 
! At 300 K, 1 MPa the internal energy and the specific volumes are:
!     e_Span  = -61.765e3_pr 
!     v_Span  = 1_pr / 18.579_pr
!     T_Span       = 300_pr
!     
!     guess_1 = T*0.8_pr 
!     guess_2 = v_Span*0.8_pr 
!
!
!      CALL New_Rap2D(2, T, v, &
!     &           resnorm, Niter, exitflag, 1e6_pr, e_Span, guess_1, guess_2)
!     
!      CALL pressure(T,v,p)
!
!      
!      print*, '----------------TEST 4-----------------------'
!      print*, 'resnorm, Niter, exitflag', resnorm, Niter, exitflag
!     print*,'-------------------------------------------------------------------------------'
!      print*, 'GUESS:        T = ', guess_1,  ' v= ', guess_2
!      print*, 'From article: T = ', T_Span,   ' v= ', v_Span
!      print*, 'Newton:       T = ', T,        ' v= ', v
!      print*,'--------------------------------------------------------------------------------'
!      print*, 'pressure in MPa, test 4', p*1e-6_pr
!
!
END PROGRAM test_pressure
