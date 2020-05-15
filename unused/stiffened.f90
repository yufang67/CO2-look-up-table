!--------------------------------------------------------------------------------------------------
! MODULE   EOS stiffened gas (pure phase)
! @note     input flag, e, v 
! @authors Yu Fang
! @date    07-2017
!--------------------------------------------------------------------------------------------------

MODULE stiffened
!
!
       IMPLICIT NONE

       REAL,PARAMETER,PRIVATE :: e_refv=290.0d3 
       REAL,PARAMETER,PRIVATE :: e_refl=1040.0d3   !880.0d3
!!------------------------------------------------
       REAL,PARAMETER,PRIVATE :: p0=  2.0d6
       REAL,PARAMETER,PRIVATE :: T0=  225.0  !360.0
       REAL,PARAMETER,PRIVATE :: c0l= 923.07 !799.17d0
       REAL,PARAMETER,PRIVATE :: c0v= 276.67
       REAL,PARAMETER,PRIVATE :: h0vv= 27.237d3
       REAL,PARAMETER,PRIVATE :: h0ll= -409.89d3 !-369.55d3
       REAL,PARAMETER,PRIVATE :: rho0l=1150.73 !1080.22
       REAL,PARAMETER,PRIVATE :: rho0v=66.092
!!------------------------------------------------
       CONTAINS


!===============================================================================================
       SUBROUTINE coeff_st(gamv,gaml,cvv,cvl,piv,pil)  !sk,Qk? mass transfer
!===============================================================================================
        IMPLICIT NONE

!Input: choose one saturation state
!        REAL(8) :: p0,T0,c0v,c0l,h0v,h0l,rho0l
!Output: 
        REAL(8) :: gamv,gaml,cvv,cvl,piv,pil
!LOCAL
!       
        REAL(8) :: e0l,h0v,h0l,e0v
!
!       
        h0v = h0vv+e_refv
        h0l = h0ll+e_refl
!
        e0l = h0l-p0/rho0l
        e0v = h0v-p0/rho0v
!
        gamv = 1.0d0 + c0v*c0v/h0v
        gaml = 1.0d0 + c0l*c0l/h0l
!
!        cvv = h0v / (T0*gamv)
        cvv = e0v / T0
        cvl = h0l / (T0*gaml)
!   
        pil = (rho0l * e0l * (gaml-1.0d0) - p0) / gaml
!        pil = (rho0l * cvl * T0 * (gaml-1.0d0) - p0) / gaml
        piv = 0.0d0
!        cvl = (e0l+pil/rho0l) / T0
!
!        PRINT*, 'gas',gamv,cvv,piv
!        PRINT*, 'liq',gaml,cvl,pil
!
       END SUBROUTINE coeff_st  
!===============================================================================================
       SUBROUTINE pressure_st(flag_phase, v,e, p)  !Pa
!===============================================================================================
        IMPLICIT NONE
!Input:
        INTEGER,INTENT(IN) :: flag_phase 
        REAL(8),INTENT(IN) :: v,e
!Output:     
        REAL(8),INTENT(OUT) :: p
!LOCAL
!        REAL(8) :: p0,T0,c0v,c0l,h0v,h0l,rho0l
        REAL(8) :: gamv,gaml,cvv,cvl,piv,pil,energy
   
       CALL coeff_st(gamv,gaml,cvv,cvl,piv,pil)


       
! flag_phase=1 liquid, 2 gas 

!        energy = e_ref + e
        IF (flag_phase==1) THEN
!          
          energy = e_refl + e
          p = (energy)/v * (gaml-1.0d0) -  pil*gaml
!       
!          print*, e/v*(gaml-1d0),  pil*gaml
        ELSEIF (flag_phase==5) THEN         
!
          Print*, 'two phases exist for stiffened gas EOS'
!          STOP 
!
        ELSE
!
          energy = e_refv + e
         
          p = (energy)/v * (gamv-1.0d0) - piv*gamv       
!
        ENDIF
!
       END SUBROUTINE pressure_st
!===============================================================================================
       SUBROUTINE temperature_st(flag_phase,v,e, T)  
!===============================================================================================
        IMPLICIT NONE
!Input:
        INTEGER,INTENT(IN) :: flag_phase
        REAL(8),INTENT(IN) :: v,e
!Output:       
        REAL(8),INTENT(OUT) :: T
!LOCAL
!       
        REAL(8) :: p,energy
!        REAL(8) :: p0,T0,c0v,c0l,h0v,h0l,rho0l
        REAL(8) :: gamv,gaml,cvv,cvl,piv,pil

       CALL coeff_st(gamv,gaml,cvv,cvl,piv,pil)
       CALL pressure_st(flag_phase,v,e, p)
!       print*, 'pressure', p,e,v
!
! flag_phase=1 liquid, 2,3,4 gas
!        energy=e+e_ref
        IF (flag_phase==1) THEN
!         
          energy=e+e_refl
!          p = 3.0d6
          
!          T = (energy  -  pil*v) / cvl
          T = (p  +  pil)*v / (cvl*(gaml-1.0d0))
!    
    
        ELSEIF (flag_phase==5) THEN
!
          Print*, 'two phases exist for stiffened gas EOS'   
!          STOP
!
        ELSE 
!
!          T = (energy) / cvv
          energy = e_refv + e
          T = (p  +  piv)*v / (cvv*(gamv-1.0d0))
!         
!       print*, flag_phase,p,T
        ENDIF
!
       END SUBROUTINE temperature_st

!===============================================================================================
       SUBROUTINE sound_speed_st(flag_phase,v,e_in, c)
!===============================================================================================
        IMPLICIT NONE
!Input:
        INTEGER,INTENT(IN) :: flag_phase
        REAL(8),INTENT(IN) :: v,e_in
!Output:       
        REAL(8),INTENT(OUT) :: c
!LOCAL
        REAL(8) :: T, s_tt, s_te, s_ee,p,e
!       
!        REAL(8) :: p0,T0,c0v,c0l,h0v,h0l,rho0l
        REAL(8) :: gamv,gaml,cvv,cvl,piv,pil

        CALL coeff_st(gamv,gaml,cvv,cvl,piv,pil)
        CALL pressure_st(flag_phase,v,e_in, p)
!        CALL temperature_st(flag_phase,v,e, T)
  

!        e = e_in+e_ref
! flag_phase=1 liquid, 2,3,4 gas 
        IF (flag_phase==1) THEN
!       
!          s_tt = cvl*(e*(gaml-1.0)*(gaml-2.0)*v**(gaml-3.0) - gaml*pil*(gaml-1.0)*v**(gaml-2.0))/&
!&                     (e*v**(gaml-1.0) - pil*v**gaml) - &
!&                cvl*(e*(gaml-1.0)*v**(gaml-2.0) - gaml*pil*v**(gaml-1.0))**(2.0) /&
!&                     ((e*v**(gaml-1.0) - pil*v**gaml)**(2.0))
!          s_te = cvl*(gaml-1.0)*v**(gaml-2.0) / (e*v**(gaml-1.0)-pil*v**gaml) - &
!&                cvl*(e*(gaml-1.0)*v**(gaml-2.0) - gaml*pil*v**(gaml-1.0))*v**(gaml-1.0)/&
!&                    (e*v**(gaml-1.0)-pil*v**gaml)**(2.0)
!          s_ee = -cvl*v**(2.0*(gaml-1.0)) / ((e*v**(gaml-1.0)-pil*v**gaml)**(2.0))
!                
!          c = T*v*v* (-s_tt + 2.0*p*s_te - p*p*s_ee)
          c = gaml*(p+pil)*v
!          c = (gaml-1.0)*cvl*gaml*T
  
        ELSEIF (flag_phase==5) THEN
!
          Print*, 'two phases exist for stiffened gas EOS'   
!          STOP
!
        ELSE 
!           s_tt = cvv*(e*(gamv-1.0)*(gamv-2.0)*v**(gamv-3.0) - gamv*piv*(gamv-1.0)*v**(gamv-2.0))/&
!&                     (e*v**(gamv-1.0) - piv*v**gamv) - &
!&                cvv*(e*(gamv-1.0)*v**(gamv-2.0) - gamv*piv*v**(gamv-1.0))**(2.0) /&
!&                     ((e*v**(gamv-1.0) - piv*v**gamv)**(2.0))
!          s_te = cvv*(gamv-1.0)*v**(gamv-2.0) / (e*v**(gamv-1.0)-piv*v**gamv) - &
!&                cvv*(e*(gamv-1.0)*v**(gamv-2.0) - gamv*piv*v**(gamv-1.0))*v**(gamv-1.0)/&
!&                    (e*v**(gamv-1.0)-piv*v**gamv)**(2.0)
!          s_ee = -cvv*v**(2.0*(gamv-1.0)) / ((e*v**(gamv-1.0)-piv*v**gamv)**(2.0))
!        
!          c = T*v*v* (-s_tt + 2.0*p*s_te - p*p*s_ee)
          c = gamv*(p+piv)*v
!          c = (gamv-1.0)*cvv*gamv*T

!        
!         
        ENDIF
          c = sqrt(c)
!
       END SUBROUTINE sound_speed_st
!
!
!
END MODULE stiffened
