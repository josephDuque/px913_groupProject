MODULE functions
USE ISO_FORTRAN_ENV

CONTAINS

!-------------------------------------------------------------
! Subroutine to convert a particles's real position into
! an integer value corresponding to the cell the particle is 
! in.
SUBROUTINE convertPos(particlePos, cellPos, dx, dy)
INTEGER, DIMENSION(2), INTENT(OUT) :: cellPos
REAL(REAL64), INTENT(IN) :: dx, dy
REAL(REAL64), DIMENSION(2), INTENT(IN) :: particlePos

cellPos(1) = FLOOR((particlePos(1) - 1.0_REAL64)/dx) + 1
cellPos(2) = FLOOR((particlePos(2) - 1.0_REAL64)/dy) + 1

END SUBROUTINE convertPos
!-------------------------------------------------------------


!-------------------------------------------------------------
! Verlet algorithm subroutine:
SUBROUTINE verletAlg(particlePos, particleVel, particleAcc, Ex, Ey, q, dx, dy, dt)
INTEGER, DIMENSION(2) :: indexPos
REAL(REAL64), INTENT(IN) :: q, dx, dy, dt
REAL(REAL64), DIMENSION(2) :: tempPos, tempVel, tempAcc
REAL(REAL64), DIMENSION(2), INTENT(INOUT) :: particlePos, particleVel, particleAcc
REAL(REAL64), DIMENSION(:, :), INTENT(IN), ALLOCATABLE :: Ex, Ey

!-------------------------------
! Update positions into a temporary array:
tempPos(1) = particlePos(1) + particleVel(1)*dt + 0.5_REAL64*particleAcc(1)*(dt**2)
tempPos(2) = particlePos(2) + particleVel(2)*dt + 0.5_REAL64*particleAcc(2)*(dt**2)
!-------------------------------

!-------------------------------
! Convert positions into indexes:
CALL convertPos(particlePos, indexPos, dx, dy)
!-------------------------------

!-------------------------------
! Update acceleration into a temporary array:
tempAcc(1) = q * Ex(indexPos(1), indexPos(2)) 
tempAcc(2) = q * Ey(indexPos(1), indexPos(2))
!-------------------------------

!-------------------------------
! Update velocity into a temporary array:
tempVel(1) = particleVel(1) + 0.5_REAL64*dt*(tempAcc(1) + particleAcc(1))
tempVel(2) = particleVel(2) + 0.5_REAL64*dt*(tempAcc(2) + particleAcc(2))
!-------------------------------

!-------------------------------
! Save the updated values:
particlePos(1) = tempPos(1)
particlePos(2) = tempPos(2)

particleAcc(1) = tempAcc(1)
particleAcc(2) = tempAcc(2)

particleVel(1) = tempVel(1)
particleVel(2) = tempVel(2)
!-------------------------------

END SUBROUTINE verletAlg
!-------------------------------------------------------------




END MODULE functions

PROGRAM velocityVerlet
USE functions
USE command_line
USE domain_tools
IMPLICIT NONE
! Initial particle position and velocity:
INTEGER :: Nx, Ny, i, j
INTEGER, PARAMETER :: iterations = 1000
INTEGER, DIMENSION(2) :: pPosIndex
LOGICAL :: successNx, successNy, successPx, successPy, successVx, successVy
REAL(REAL64) :: dx, dy
REAL(REAL64), PARAMETER :: dt = (1.0_REAL64)*1.0e-2, eCharge = -1.0_REAL64 !eCharge = (1.602176634_REAL64)*1.0e-19 
! Particle position, velocity and acceleration:
REAL(REAL64), DIMENSION(2) :: pPos, pVel, pAcc, axisSize,
REAL(REAL64), DIMENSION(:), ALLOCATABLE :: xAxis, yAxis
REAL(REAL64), DIMENSION(:, :), ALLOCATABLE :: Ex, Ey
CHARACTER(LEN=*), PARAMETER :: ExFilename = 'xE_field.txt' 
CHARACTER(LEN=*), PARAMETER :: EyFilename = 'yE_field.txt'
CHARACTER(LEN=40) :: formatName

CALL parse_args
axisSize(1) = -1.0_REAL64
axisSize(2) = 1.0_REAL64

successNx = get_arg('Nx', Nx)
successNy = get_arg('Ny', Ny)

IF (successNx) THEN
  IF (Nx <= 0 .AND. KIND(Nx) /= KIND(INT32)) THEN
    PRINT *, 'Please input a positive integer for the value of Nx.'
    STOP
  ENDIF
ELSE
  PRINT *, 'Command line argument for Nx was not parsed. Please make sure Nx is a positive integer.'
  STOP
END IF
IF (successNy) THEN
  IF (Ny <= 0 .AND. KIND(Ny) /= KIND(INT32)) THEN
    PRINT *, 'Please input a positive integer for the value of Ny.'
    STOP
  ENDIF
ELSE
  PRINT *, 'Command line argument for Ny was not parsed. Please make sure Ny is a positive integer.'
  STOP
END IF

CALL create_axis(xAxis, Nx, axisSize, delta=dx)
CALL create_axis(yAxis, Ny, axisSize, delta=dy)

! Initialise particle position and velocity:
successPx = get_arg('Px', pPos(1))
successPy = get_arg('Py', pPos(2))
successVx = get_arg('Vx', pVel(1))
successVy = get_arg('Vy', pVel(2))
IF (successPx) THEN
  IF (pPos(1) < xAxis(1) .OR. pPos(1) > xAxis(UBOUND(xAxis, 1)) .OR. KIND(pPos(1)) /= 2*KIND(REAL64)) THEN
    PRINT *, "Please input a real value of the particle's x-position within the domain: ", xAxis(1), 'to ', &
      & xAxis(UBOUND(xAxis,1))
    STOP
  ENDIF
ELSE
  PRINT *, "Command line argument for the particle's x-position was not parsed. Please make sure it is a real value within the &
  & domain."
  STOP
END IF
IF (successPy) THEN
  IF (pPos(2) < yAxis(1) .OR. pPos(2) > yAxis(UBOUND(yAxis, 1)) .OR. KIND(pPos(2)) /= 2*KIND(REAL64)) THEN
    PRINT *, "Please input a real value of the particle's y-position within the domain: ", yAxis(1), 'to ', &
      & yAxis(UBOUND(yAxis,1))
    STOP
  END IF
ELSE
  PRINT *, "Command line argument for the particle's y-position was not parsed. Please make sure it is a real value within the &
  & domain."
  STOP
END IF
IF (successVx) THEN
  IF (KIND(pVel(1)) /= 2*KIND(REAL64)) THEN
    PRINT *, "Please input a real number for the value of the particle's x-velocity. This can be positive or negative depending on &
    & the direction, or it can have no component to the x-velocity."
    STOP
  END IF
ELSE
  PRINT *, "Command line argument for the particle's x-velocity not parsed. Please make sure it is the correct data type and kind."
  STOP
END IF
IF (successVy) THEN
  IF (KIND(pVel(2)) /= 2*KIND(REAL64)) THEN
    PRINT *, "Please input a real number for the value of the particle's y-velocity. This can be positive or negative depending on &
    & the direction, or it can have no component to the y-velocity."
    STOP
  END IF
ELSE
  PRINT *, "Command line argument for the particle's y-velocity not parsed. Please make sure it is the correct data type and kind."
  STOP
END IF

ALLOCATE(Ex(Ny, Nx))
ALLOCATE(Ey(Ny, Nx))

WRITE (formatName, '(A1,I4,A13)') "(",Nx,"(ES23.12,2X))"
OPEN (9, file=ExFilename, form='formatted')
OPEN (10, file=EyFilename, form='formatted')

DO i = 1, Ny
  READ(9, formatName) (Ex(i, j), j = 1, Nx)
END DO

DO i = 1, Ny
  READ(10, formatName) (Ey(i, j), j = 1, Nx)
END DO

CLOSE(9)
CLOSE(10)

!DO i = 1, Ny
!  PRINT *, (Ex(i, j), j = 1, Nx)
!END DO

CALL convertPos(pPos, pPosIndex, dx, dy)

!---------------------------------------------------
! Initial acceleration:
pAcc(1) = eCharge * Ex(pPosIndex(1), pPosIndex(2))
pAcc(2) = eCharge * Ey(pPosIndex(1), pPosIndex(2))
!---------------------------------------------------




DEALLOCATE(Ex); DEALLOCATE(Ey)

END PROGRAM velocityVerlet
