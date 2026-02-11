        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:33 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MAX_DAY_LENGTH__genmod
          INTERFACE 
            SUBROUTINE MAX_DAY_LENGTH(SUBCOUNT,LAT,MAX_DLENGTH)
              INTEGER(KIND=4), INTENT(IN) :: SUBCOUNT
              REAL(KIND=4), INTENT(IN) :: LAT(SUBCOUNT)
              REAL(KIND=4), INTENT(INOUT) :: MAX_DLENGTH(SUBCOUNT)
            END SUBROUTINE MAX_DAY_LENGTH
          END INTERFACE 
        END MODULE MAX_DAY_LENGTH__genmod
