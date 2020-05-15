!---------------------------------------------------------------------------------------------------
! All the formulations  have been taken from the original
! article, THERMODYNAMICS PROPERTIES INVOLVING DERIVATIVES using the peng-robinson equation of state
!---------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------
! MODULE   peng_robinson
! @routine 
! @authors Yu Fang
! @date    07-2017
!--------------------------------------------------------------------------------------------------

MODULE peng_robinson
!
!
       IMPLICIT NONE
!
! CONSTANTS for CO2

       REAL,PARAMETER,PRIVATE :: R_u   =  8.31445d0
       REAL,PARAMETER,PRIVATE :: T_c   =  304.1282d0
       REAL,PARAMETER,PRIVATE :: P_c   =  7377300.0d0
       REAL,PARAMETER,PRIVATE :: omega =  0.228d0
       REAL,PARAMETER,PRIVATE :: m_m   =  0.044d0
       REAL,PARAMETER,PRIVATE :: gam   =  1.3d0       !1.4? 
       REAL,PRIVATE :: m     =  0.37464d0 + 1.54226d0*omega - 0.26992d0*omega*omega
       REAL,PRIVATE :: a_c   =  0.45723553d0*R_u*R_u*T_c*T_c/P_c
       REAL,PRIVATE :: b     =  0.077796074d0*R_u*T_c/P_c
  
!      
!       
       PUBLIC :: pressure_pr, interenergy_pr, soundspeed_pr
       
!
       CONTAINS
       

!
!
!===============================================================================================
       SUBROUTINE pressure_pr(T,v,p)  !Pa
!===============================================================================================
        IMPLICIT NONE
!
        REAL(8) :: T, v, a
        REAL(8) :: p
!LOCAL
!       
        REAL(8) :: v_mol
!
!       
        v_mol = m_m*v
        a = a_c * ( 1.d0 + m * (1d0 - sqrt(T/T_c)) ) 
!        
        p = R_u*T/(v_mol-b) - a/( v_mol*(v_mol+b) + b*(v_mol-b) )
!
!
       END SUBROUTINE pressure_pr

!===============================================================================================
       SUBROUTINE interenergy_pr(T,v,e)  !J/kg
!===============================================================================================
        IMPLICIT NONE
!
!IN/OUT
        REAL(8) :: T,v,a
        REAL(8) :: e
!LOCAL
        REAL(8) :: v_mol,u0,ui,a_div1, z, BB,p
!
!
       u0 = R_u/( m_m *(gam-1.d0) ) * T
!
       v_mol = m_m*v
       a = a_c * ( 1.d0 + m * (1d0 - sqrt(T/T_c)) ) 

       CALL pressure_pr(T,v,p)

       z = p*v_mol/(R_u*T)
       BB = b*p/(R_u*T)
       a_div1 = -a_c*m*0.5d0*sqrt(1.d0/(T_c*T))
!
       ui = (T*a_div1 - a)/(b*sqrt(8.0d0)) * log(( z + BB*(1.d0+sqrt(2.d0)) )/ &
&                                           ( z + BB*(1.d0-sqrt(2.d0)) ))    
!unit from J/mol --> J/kg          
       ui = ui/m_m

       e = u0 + ui
       END SUBROUTINE interenergy_pr

!===============================================================================================
       SUBROUTINE heatcap_v_pr(T,v,cv)  !J/(kg.K)
!===============================================================================================
        IMPLICIT NONE
!
!IN/OUT
        REAL(8) :: T,v,a
        REAL(8) :: cv
!LOCAL
        REAL(8) :: v_mol,cv0,cvi,a_div2, z, BB,p
!
!
       cv0 = R_u/( m_m *(gam-1.d0) )
!
       v_mol = m_m*v
       a = a_c * ( 1.d0 + m * (1d0 - sqrt(T/T_c)) )

       CALL pressure_pr(T,v,p)

       z = p*v_mol/(R_u*T)
       BB = b*p/(R_u*T)
       a_div2 = a_c*m*0.25d0*sqrt(1.d0/(T_c*T*T*T))
!
       cvi = (T*a_div2)/(b*sqrt(8.0d0)) * log(( z + BB*(1.d0+sqrt(2.d0)) )/ &
&                                           ( z + BB*(1.d0-sqrt(2.d0)) ))
!unit from J/mol --> J/kg          
       cvi = cvi/m_m

       cv = cv0 + cvi
       END SUBROUTINE heatcap_v_pr

!===============================================================================================
       SUBROUTINE heatcap_p_pr(T,v,cp)  !J/(kg.K)
!===============================================================================================
        IMPLICIT NONE
!
!IN/OUT
        REAL(8) :: T,v,a
        REAL(8) :: cp
!LOCAL
        REAL(8) :: v_mol,cv0,cp0,cvi,cpi,a_div1,a_div2, z, BB,p,AA
!
        REAL(8) :: dp_dT, dv_dT, dz_dT, dAA_dT, dBB_dT
!
       cv0 = R_u/( m_m *(gam-1.d0) ) 
       cp0 = cv0*gam
!
       v_mol = m_m*v
       a = a_c * ( 1.d0 + m * (1d0 - sqrt(T/T_c)) )

       CALL pressure_pr(T,v,p)

       z = p*v_mol/(R_u*T)
       AA = a*p / (R_u*T*T*R_u)
       BB = b*p/(R_u*T)
       a_div1 = -a_c*m*0.5d0*sqrt(1.d0/(T_c*T))
       a_div2 = a_c*m*0.25d0*sqrt(1.d0/(T_c*T*T*T))
!
       cvi = (T*a_div2)/(b*sqrt(8.0d0)) * log(( z + BB*(1.d0+sqrt(2.d0)) )/ &
&                                           ( z + BB*(1.d0-sqrt(2.d0)) ))
!       
       dAA_dT = p/(R_u*R_u*T*T) * (a_div1 - 2.0d0*a/T)
       dBB_dT = -b*p / (R_u*T*T)
!      
       dz_dT = dAA_dT*(BB-z) + dBB_dT*(6.0*BB*z+2.0*z-3.0*BB*BB-2.0*BB+AA-z*z)/&
&              (3.0*z*z+2.0*(BB-1.0)*z+(AA-2.0*BB-3.0*BB*BB)) 


       dp_dT = R_u/(v_mol-b) - a_div1/(v_mol*(v_mol+b) + b*(v_mol-b))
       dv_dT = R_u/p * ( T*dz_dT + z ) 
!
       cpi = cvi + T*dp_dT*dv_dT - R_u
       cpi = cpi/ m_m
!       
       cp = cp0 + cpi
!
!
       END SUBROUTINE heatcap_p_pr


!===============================================================================================
       SUBROUTINE soundspeed_pr(T,v,c)  !m/s
!===============================================================================================
        IMPLICIT NONE
!
!IN/OUT
        REAL(8) :: T,v,a
        REAL(8) :: c
!LOCAL
         REAL(8) :: v_mol,cp,cv,dp_dv,c2
!
       v_mol = m_m*v
       a = a_c * ( 1.d0 + m * (1d0 - sqrt(T/T_c)) )
!
       dp_dv = -R_u*T/((v_mol-b)*(v_mol-b)) + 2.0*a*(v_mol+b)/&
&             ( (v_mol*(v_mol+b)+b*(v_mol-b)) * (v_mol*(v_mol+b)+b*(v_mol-b)) )
!       
      CALL heatcap_v_pr(T,v,cv)
      CALL heatcap_p_pr(T,v,cp)
!
      c2 = -v_mol*v_mol*cp/cv*dp_dv
!convert unit
      c = sqrt(c2/m_m)
!
       END SUBROUTINE soundspeed_pr


END MODULE peng_robinson  
