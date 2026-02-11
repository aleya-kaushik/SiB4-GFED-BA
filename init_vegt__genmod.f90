        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:46:52 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INIT_VEGT__genmod
          INTERFACE 
            SUBROUTINE INIT_VEGT(NSOIL,VEGT)
              USE MODULE_SIB, ONLY :                                    &
     &          VEG_TYPE
              INTEGER(KIND=4), INTENT(IN) :: NSOIL
              TYPE (VEG_TYPE), INTENT(INOUT) :: VEGT
            END SUBROUTINE INIT_VEGT
          END INTERFACE 
        END MODULE INIT_VEGT__genmod
