        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:46:52 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INIT_SSCOLT__genmod
          INTERFACE 
            SUBROUTINE INIT_SSCOLT(NSNOW,NSOIL,SSCOLT)
              USE MODULE_SIB, ONLY :                                    &
     &          SSCOL_TYPE
              INTEGER(KIND=4), INTENT(IN) :: NSNOW
              INTEGER(KIND=4), INTENT(IN) :: NSOIL
              TYPE (SSCOL_TYPE), INTENT(INOUT) :: SSCOLT
            END SUBROUTINE INIT_SSCOLT
          END INTERFACE 
        END MODULE INIT_SSCOLT__genmod
