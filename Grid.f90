!===========================================================================
!@ Module to create the grid in (e,v) domaine for LL, LH, R, HT and TP 
!@ The values of the properties (e,v,T,P,c) are computed by EoS Span-Wagner
!@ and available at each node of the grid.
!==========================================================================
MODULE Grid
!
!        USE def_constants
!        USE def_variables
!        USE non_linear_solvers
!      
        IMPLICIT NONE
!
        PRIVATE
!
        PUBLIC MAKE_GRID
CONTAINS
!       
!
!
!==========================================================================
        SUBROUTINE MAKE_GRID()
!==========================================================================
!
         USE saturation, ONLY: sat_curve, saturation_curve
!         USE mod_param_defs
!         USE mod_error, ONLY: print_error_and_quit, print_message
         IMPLICIT NONE

!         CHARACTER(LEN=strl) :: message



                CALL grid_construction_left_low
                CALL grid_construction_left_high
                CALL grid_construction_right
                CALL grid_construction_high_temperature
                CALL saturation_curve
!                CALL sat_curve
                CALL grid_construction_TPL
                CALL grid_construction_TPM
                CALL grid_construction_TPH
!
         print*,  " ----> Look-up table CO2 MESH CONSTRUCTION ----"
!         CALL print_message ( message )
!
        END SUBROUTINE MAKE_GRID
!
!
END MODULE Grid
