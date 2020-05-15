MODULE var_const
!  USE mod_param_defs

  IMPLICIT NONE
  
!-----------------------------------------------------------------
!
!-----------------------------------------------------------------
!    real(8),parameter :: e_const = 1.0d6 
!    integer,parameter :: Ngost   = 2
!    integer,parameter :: Nxmax   = 100
!    integer,parameter :: Nymax   = 40
!    integer,parameter :: N_1d    = 1000




!-----------------------------------------------------------------
!
!-----------------------------------------------------------------

!    real(8),dimension(1-Ngost:Nxmax+Ngost,1-Ngost:Nymax+Ngost) :: Tguess
!    real(8),dimension(1-Ngost:Nxmax+Ngost,1-Ngost:Nymax+Ngost) :: guessP
!    real(8),dimension(1-Ngost:Nxmax+Ngost,1-Ngost:Nymax+Ngost) :: c_BC
!    real(8),dimension(1-Ngost:N_1d+Ngost) :: guessp_1d
!    real(8),dimension(1-Ngost:Nymax+Ngost) :: guessp_BC
!    real(8),dimension(1-Ngost:Nymax+Ngost) :: vguess_r
!    real(8),dimension(1-Ngost:Nymax+Ngost) :: guessE
    integer :: i_bc
    integer :: flag_diagno
!    integer :: flag_loca        
    integer :: f_out
    REAL(8) :: xloc, yloc
    CHARACTER(LEN=240) :: mess1, mess2


END MODULE var_const
