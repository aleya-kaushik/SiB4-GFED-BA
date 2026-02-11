        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:46:52 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INIT_EQUIBLT__genmod
          INTERFACE 
            SUBROUTINE INIT_EQUIBLT(NPOOLPFT,EQUIBLT)
              USE MODULE_SIB, ONLY :                                    &
     &          EQUIBL_TYPE
              INTEGER(KIND=4), INTENT(IN) :: NPOOLPFT
              TYPE (EQUIBL_TYPE), INTENT(INOUT) :: EQUIBLT
            END SUBROUTINE INIT_EQUIBLT
          END INTERFACE 
        END MODULE INIT_EQUIBLT__genmod
