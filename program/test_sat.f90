
PROGRAM test_TP
!      
      USE def_constants, ONLY: pr, T_cr,T_tri,P_cr,P_tri
!      USE def_variables  
!      USE non_linear_solvers
      USE properties, ONLY: satprop
!      USE grid_functions
!      USE Grid      
!      USE Interp_table
      USE saturation, ONLY: saturation_curve, sat_curve
      IMPLICIT NONE
!
      INTEGER :: i,j,NN
      REAL(pr), DIMENSION(600) :: p_tab
      REAL(pr) :: Terr,xerr       
      REAL(pr) :: delta, Tsat, vvsat, vlsat, uvsat, ulsat
      REAL(pr) :: Tsat5, vvsat5, vlsat5, uvsat5, ulsat5
!
        CALL saturation_curve()
        CALL sat_curve()
!
      NN = 600
      delta = (P_cr - P_tri)/ (NN)
      p_tab = p_tri + (/(i*delta,i=0,NN-1)/)

  
      DO i = 1,600
           CALL satprop(3, p_tab(i), Tsat, vvsat, vlsat, uvsat, ulsat)
           CALL satprop(5, p_tab(i), Tsat5, vvsat5, vlsat5, uvsat5, ulsat5)
          IF ( (abs(Tsat-Tsat5)  > 1e-12_pr).OR.&
&              (abs(vvsat-vvsat5)> 1e-12_pr).OR.&
&              (abs(vlsat-vlsat5)> 1e-12_pr).OR.&
&              (abs(uvsat-uvsat5)> 1e-12_pr).OR.&
&              (abs(ulsat-ulsat5)> 1e-12_pr) ) THEN
           print*, 'i',i,'p', p_tab(i)
           print*, 'T', Tsat,Tsat5, abs(Tsat-Tsat5)
           print*, 'v-vapor', vvsat,vvsat5, abs(vvsat-vvsat5)
           print*, 'v-liquid', vlsat,vlsat5, abs(vlsat-vlsat5)
           print*, 'u-vapor', uvsat,uvsat5, abs(uvsat-uvsat5)
           print*, 'u-liquid', ulsat,ulsat5, abs(ulsat-ulsat5)
          ENDIF
!       
      ENDDO
END PROGRAM test_TP
