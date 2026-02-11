        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:46:52 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INIT_HYDROVT__genmod
          INTERFACE 
            SUBROUTINE INIT_HYDROVT(NSOIL,HYDROVT)
              USE MODULE_SIB, ONLY :                                    &
     &          HYDROV_TYPE
              INTEGER(KIND=4), INTENT(IN) :: NSOIL
              TYPE (HYDROV_TYPE), INTENT(INOUT) :: HYDROVT
            END SUBROUTINE INIT_HYDROVT
          END INTERFACE 
        END MODULE INIT_HYDROVT__genmod
