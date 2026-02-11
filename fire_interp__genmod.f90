        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:53 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FIRE_INTERP__genmod
          INTERFACE 
            SUBROUTINE FIRE_INTERP(INDX,LON,LAT,SIBG)
              USE MODULE_SIB, ONLY :                                    &
     &          GRIDCELL_TYPE
              INTEGER(KIND=4), INTENT(IN) :: INDX
              REAL(KIND=4), INTENT(IN) :: LON
              REAL(KIND=4), INTENT(IN) :: LAT
              TYPE (GRIDCELL_TYPE), INTENT(INOUT) :: SIBG
            END SUBROUTINE FIRE_INTERP
          END INTERFACE 
        END MODULE FIRE_INTERP__genmod
