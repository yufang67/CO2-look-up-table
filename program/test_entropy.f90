PROGRAM test_entropy
!
      USE def_constants
      USE def_variables 
      USE Grid
!      USE Interp_table 
!      USE non_linear_solvers
      USE properties
      IMPLICIT NONE
!
!      INTEGER :: flag!, Niter, exitflag,i,j
!      REAL(pr):: T, v, p, e, cv, cp, s, c, v_l, v_v, v_l_Span,v_v_Span,&
!&           resnorm, guess_1, guess_2, T_span, v_span, e_span,helmho1
!      REAL(pr):: p2, e2, cv2, cp2, s2, c2, helmho2,Helmholtz
      REAL(pr):: v_in1,v_in2,T_in1,T_in2,v_in3,v_in4,T_in3,T_in4
      REAL(pr):: s_out1,s_out2,s_out3,s_out4,p_out2,p_out4
      REAL(pr):: p_in1, T1, vv1, vl1, uv1, ul1
      REAL(pr):: p_in2, T2, vv2, vl2, uv2, ul2

!
!
!
!
        CALL MAKE_GRID()

        p_in1 = 4.3339e6
        p_in2 = 5.785e6

        CALL satprop(3, p_in1, T1, vv1, vl1, uv1, ul1)
        CALL satprop(3, p_in2, T2, vv2, vl2, uv2, ul2)

        print*, 'p1= ',p_in1, 'T1= ',T1, 'rho1= ',1.0/vv1, 'rhol1= ',1.0/vl1
        print*, 'p2= ',p_in2, 'T2= ',T2, 'rho2= ',1.0/vv2, 'rhol2= ',1.0/vl2
!STOP  
        v_in1=1.0/128.42686 
        T_in1=281.62 
        
        v_in2= 1.0/872.23144
        T_in2= 281.62


        v_in3= 1.0/197.49158
        T_in3= 293.56

        v_in4= 1.0/768.93522
        T_in4= 293.56

        CALL entropy(T_in1,v_in1,s_out1)
        CALL entropy(T_in2,v_in2,s_out2)
        CALL entropy(T_in3,v_in3,s_out3)
        CALL entropy(T_in4,v_in4,s_out4)

        
        CALL pressure(T_in2,v_in2,p_out2)
        CALL pressure(T_in4,v_in4,p_out4)

        print*, 'v1', v_in1, 'T1=',T_in1, 's1=',s_out1
        print*, 'v2', v_in2, 'T2=',T_in2, 's2=',s_out2,p_out2
        print*, 'v3', v_in3, 'T3=',T_in3, 's3=',s_out3
        print*, 'v4', v_in4, 'T4=',T_in4, 's4=',s_out4,p_out4




END PROGRAM test_entropy
