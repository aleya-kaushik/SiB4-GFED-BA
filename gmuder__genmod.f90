        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:23 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GMUDER__genmod
          INTERFACE 
            SUBROUTINE GMUDER(LAT,DOY,CHIL,GMUDMU)
              REAL(KIND=4), INTENT(IN) :: LAT
              INTEGER(KIND=4), INTENT(IN) :: DOY
              REAL(KIND=4), INTENT(IN) :: CHIL
              REAL(KIND=8), INTENT(OUT) :: GMUDMU
            END SUBROUTINE GMUDER
          END INTERFACE 
        END MODULE GMUDER__genmod
