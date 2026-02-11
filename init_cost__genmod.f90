        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:46:52 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INIT_COST__genmod
          INTERFACE 
            SUBROUTINE INIT_COST(COST,NSOIL)
              USE MODULE_SIB, ONLY :                                    &
     &          COS_TYPE
              TYPE (COS_TYPE), INTENT(INOUT) :: COST
              INTEGER(KIND=4), INTENT(IN) :: NSOIL
            END SUBROUTINE INIT_COST
          END INTERFACE 
        END MODULE INIT_COST__genmod
