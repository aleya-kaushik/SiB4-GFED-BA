        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:30 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RADABS__genmod
          INTERFACE 
            SUBROUTINE RADABS(DLWBOT,GDIAGT,RADT)
              USE MODULE_SIB, ONLY :                                    &
     &          GDIAG_TYPE,                                             &
     &          RAD_TYPE
              REAL(KIND=8), INTENT(IN) :: DLWBOT
              TYPE (GDIAG_TYPE), INTENT(IN) :: GDIAGT
              TYPE (RAD_TYPE), INTENT(INOUT) :: RADT
            END SUBROUTINE RADABS
          END INTERFACE 
        END MODULE RADABS__genmod
