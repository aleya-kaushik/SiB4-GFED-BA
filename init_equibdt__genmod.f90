        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:46:52 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INIT_EQUIBDT__genmod
          INTERFACE 
            SUBROUTINE INIT_EQUIBDT(NPOOLLU,EQUIBDT)
              USE MODULE_SIB, ONLY :                                    &
     &          EQUIBD_TYPE
              INTEGER(KIND=4), INTENT(IN) :: NPOOLLU
              TYPE (EQUIBD_TYPE), INTENT(INOUT) :: EQUIBDT
            END SUBROUTINE INIT_EQUIBDT
          END INTERFACE 
        END MODULE INIT_EQUIBDT__genmod
