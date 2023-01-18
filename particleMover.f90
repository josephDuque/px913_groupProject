MODULE subroutines
    USE iso_fortran_env
    USE command_line
    IMPLICIT NONE
    CONTAINS
    
    SUBROUTINE params(Ny, Nx, rho)
        LOGICAL :: success, exists
        INTEGER, INTENT(OUT) :: Nx,Ny
        CHARACTER(len=25) :: rho

        CALL parse_args !CALL parse_args function to 'grab and store' arguments specified in command line.
        success = get_arg("Nx", Nx,exists=exists)
        IF(success) THEN
            IF(Nx <= 0_INT32) THEN    !The number of grid cells in x-direction must be a positive integer. Thus, if a either a negative number or zero are specified, an error relating to this is printed, and the program is stopped.
                PRINT*,'Error: Nx must be set to a positive value'
                STOP
            END IF
        ELSE
            IF(exists) THEN
                PRINT*,"Value for 'Nx' of incorrect type/kind passed, ensure command line input for 'Nx' is suitable for type and&
                 & kind specified by INT32."
                STOP
            ELSE
                PRINT*,"Value for 'Nx' not parsed, check command line input."
                STOP
            END IF
        END IF

        success = get_arg("Ny", Ny,exists=exists)
        IF(success) THEN
            IF(Ny <= 0_INT32) THEN    !The number of grid cells in y-direction must be a positive integer. Thus, if a either a negative number or zero are specified, an error relating to this is printed, and the program is stopped.
                PRINT*,'Error: Ny must be set to a positive value'
                STOP
            END IF
        ELSE
            IF(exists) THEN
                PRINT*,"Value for 'Ny' of incorrect type/kind passed, ensure command line input for 'Ny' is suitable for type and&
                 & kind specified by INT32."
                STOP
            ELSE
                PRINT*,"Value for 'Ny' not parsed, check command line input."
                STOP
            END IF
        END IF

        success = get_arg("rho", rho, exists=exists)
        IF(success) THEN
        ELSE
            IF(exists) THEN
                PRINT*,"Value for 'rho' of incorrect type passed, ensure command line input for 'init' is one of &
                &'null','single' or 'double'"
                STOP
            ELSE
                PRINT*,"Value for 'rho' not parsed, check command line input"
                STOP
            END IF
        END IF

    END SUBROUTINE params

    
    SUBROUTINE null(lower, rho_grid, initPos, initVel)
        INTEGER :: lower
        REAL(REAL64), DIMENSION(2), INTENT(OUT) :: initPos, initVel
        REAL(REAL64), DIMENSION(lower:,lower:), INTENT(INOUT) :: rho_grid
                
        rho_grid = 0.0_REAL64

        ! Setting initial position and velocity:
        initPos = 0.0_REAL64
        initVel = 0.1_REAl64

    END SUBROUTINE null


    SUBROUTINE single(lower, rho_grid,x_axis,y_axis,Nx, Ny, initPos, initVel)
        INTEGER, INTENT(IN) :: lower, Nx, Ny
        REAL(REAL64), DIMENSION(2), INTENT(OUT) :: initPos, initVel
        REAL(REAL64), DIMENSION(lower:,lower:), INTENT(INOUT) :: rho_grid
        REAL(REAL64), DIMENSION(:), INTENT(IN) :: x_axis, y_axis
        REAL(REAL64) :: x_div, y_div, z_1
        INTEGER :: i,j


        z_1 = 0.1_REAL64


            
        DO i=1,Nx
            x_div = ((x_axis(i))/z_1)        
            DO j=1,Ny
                y_div = ((y_axis(j))/z_1)
                rho_grid(j,i) = EXP(-(x_div**2)-(y_div**2))
            
            END DO
        
        END DO
        

        ! Setting initial position and velocity:
        initPos(1) = 0.1_REAL64; initPos(2) = 0.0_REAL64
        initVel = 0.0_REAl64
          

    END SUBROUTINE single


    SUBROUTINE double(lower, rho_grid,x_axis,y_axis,Nx, Ny, initPos, initVel)
        INTEGER, INTENT(IN) :: lower, Nx, Ny
        REAL(REAL64), DIMENSION(2), INTENT(OUT) :: initPos, initVel
        REAL(REAL64), DIMENSION(lower:,lower:), INTENT(INOUT) :: rho_grid
        REAL(REAL64), DIMENSION(:), INTENT(IN) :: x_axis, y_axis
        REAL(REAL64) :: x_div, y_div, extra_divx, extra_divy, z_1, z_2
        INTEGER :: i, j
        
        z_1 = 0.1_REAL64
        z_2 = 0.2_REAL64


        DO i=1,Nx
            x_div = ((x_axis(i)+0.25_REAL64)/z_1)

            DO j=1,Ny
                y_div = ((y_axis(j)+0.25_REAL64)/z_1)
                extra_divy = ((y_axis(j)-0.75_REAL64)/z_2)
                extra_divx = ((x_axis(i)-0.75_REAL64)/z_2)
                rho_grid(j,i) = EXP(-(x_div**2)-(y_div**2)) + EXP(-(extra_divx**2)-(extra_divy**2))
            
            END DO

        END DO
        
        ! Setting initial position and velocity:
        initPos(1) = 0.0_REAL64; initPos(2) = 0.5_REAL64
        initVel = 0.0_REAl64

    END SUBROUTINE double


    SUBROUTINE phi(lower, phi_grid, Nx, Ny, rho_grid, dx, dy)
        INTEGER, INTENT(IN) :: lower, Nx, Ny
        REAL(REAL64), INTENT(IN) :: dx, dy
        REAL(REAL64), DIMENSION(lower:,lower:), INTENT(INOUT) :: rho_grid, phi_grid
        REAL(REAL64) :: dens, div1, div2, two, denom, N
        REAL(REAL64) :: etot, conv1, conv2, drms_sum, drms,drms1
        INTEGER :: i, j, zed

        two = 2.0_REAL64

        etot = 10.0_REAL64
        drms = 1.0_REAL64

        zed=0_INT32
        DO i=1,Nx
            DO j=1,Ny                
                dens=rho_grid(j,i)                                      
                div1=((phi_grid(j,i+1) + phi_grid(j,i-1))/(dx**2))
                div2=((phi_grid(j+1,i) + phi_grid(j-1,i))/(dy**2))
                denom=((two/(dx**2))+(two/(dy**2)))
                phi_grid(j,i)=(-(ABS(dens-div1-div2)/denom))            !Calculate one initial set of elements for phi grid.
                                                                        !ABS() in final line of DO loop initializes grid to give           
            END DO                                                      !positive values of 'dmrs', allowing convergence criteria to work.

        END DO


        DO WHILE(etot/drms > 0.00001_REAL64)
            
            DO i=1,Nx
                DO j=1,Ny                
                    dens=rho_grid(j,i)                                      !Calculates phi values for whole grid once
                    div1=((phi_grid(j,i+1) + phi_grid(j,i-1))/(dx**2))
                    div2=((phi_grid(j+1,i) + phi_grid(j-1,i))/(dy**2))
                    denom=((two/(dx**2))+(two/(dy**2)))
                    phi_grid(j,i)=-((dens-div1-div2)/denom)

                END DO

            END DO

            N = REAL(Nx, kind=REAL64)*REAL(Ny, kind=REAL64)    !Make REAL to allow use of N in generation of REALs below.

            etot = 0.0_REAL64
            drms_sum = 0.0_REAL64       !(Re-)Initialise values
            drms = 0.0_REAL64

        
            DO i=1,Nx
                DO j=1,Ny
                    
                    dens=rho_grid(j,i)
                    conv1 = ((phi_grid(j,i-1) - (two*(phi_grid(j,i))) + phi_grid(j,i+1))/(dx**2))
                    conv2 = ((phi_grid(j-1,i) - (two*(phi_grid(j,i))) + phi_grid(j+1,i))/(dy**2))
                    
                    etot = etot + ABS(conv1 + conv2 - dens)

                    drms_sum = drms_sum + conv1 + conv2

                    !PRINT*,ABS(drms_sum)
                    !PRINT*,1/N
                    !PRINT*,(1/N)*(ABS(drms_sum))
                
                END DO
            
            END DO

            drms = SQRT((1/N)*((drms_sum)))
            drms1 = SQRT((1/N)*(ABS(drms_sum)))        
            zed=zed+1
        
        END DO
        !PRINT*,'------------------------------------------------------------------------------'
        !PRINT*,'------------------------------------------------------------------------------'
        !PRINT*,'------------------------------------------------------------------------------'
        !PRINT*,'zed=',zed
        !PRINT*,etot/drms
        !PRINT*, etot/drms1
        !PRINT*,'------------------------------------------------------------------------------'
        !PRINT*,'------------------------------------------------------------------------------'
        !PRINT*,etot
        !PRINT*,'------------------------------------------------------------------------------'
        !PRINT*,drms_sum
        !PRINT*,'------------------------------------------------------------------------------'
        !PRINT*,drms
   


    END SUBROUTINE phi


    SUBROUTINE E_fields(E_field_x, E_field_y, phi_grid, dx, dy, Nx, Ny,lower)
        INTEGER, INTENT(IN) :: Nx, Ny, lower
        REAL(REAL64), DIMENSION(:,:), INTENT(INOUT) :: E_field_x, E_field_y
        REAL(REAL64), DIMENSION(lower:,lower:), INTENT(IN) :: phi_grid
        REAL(REAL64), INTENT(IN) :: dx, dy
        REAL(REAL64) :: two
        INTEGER :: i, j

        two = 2.0_REAL64
        DO i=1,Nx
            Do j=1,Ny
                E_field_x(j,i) = ((phi_grid(j,i+1) - phi_grid(j,i-1))/(two*dx))
                E_field_y(j,i) = ((phi_grid(j+1,i) - phi_grid(j-1,i))/(two*dy))

            END DO
        
        END DO
    END SUBROUTINE

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
    tempAcc(1) = q * Ex(indexPos(2), indexPos(1))
    tempAcc(2) = q * Ey(indexPos(2), indexPos(1))
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


END MODULE subroutines


PROGRAM main                !REMEMBER KING, I HAVE INDEXED IN STRANGE WAYS (J,I)
    USE iso_fortran_env
    USE subroutines
    USE domain_tools
    USE write_netcdf
    IMPLICIT NONE
    INTEGER :: Ny, Nx, i, lower
    INTEGER, PARAMETER :: iterations = 1000
    INTEGER, DIMENSION(2) :: pPosIndex
    REAL(REAL64) :: dx, dy
    REAL(REAL64), PARAMETER :: dt = (1.0_REAL64)*1.0e-2, eCharge = -1.0_REAL64 !eCharge = (1.602176634_REAL64)*1.0e-19
    REAL(REAL64), DIMENSION(2) :: pPos, pVel, pAcc, axis_range
    REAL(REAL64), DIMENSION(2, 0:iterations) :: positions, velocities, accels
    REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: grid, rho_grid, phi_grid, E_field_x, E_field_y, &
      & netcdf_rho, netcdf_phi, netcdf_Ex, netcdf_Ey
    REAL(REAL64), DIMENSION(:), ALLOCATABLE :: x_axis, y_axis
    CHARACTER(len=25) :: rho
    CHARACTER(LEN=*), PARAMETER :: netcdf_filename = 'particle-move.nc'
    CHARACTER(LEN=*), PARAMETER :: code_filename = 'particleMover.f90'
    INTEGER :: neterr
    
    CALL params(Ny, Nx, rho) !This subroutine retrieves values inputted into command line for Nx and Ny parameters.

    ALLOCATE(grid(0:Ny+1,0:Nx+1))
!-------------------------------------------------Initialise All Cells, Then Set Ghost Cells-----------------------------------------------------
    grid=1.0_REAL64
    
    DO i=0,Nx+1
        grid(0,i) = 0.0_REAL64
        grid((Ny+1),i) = 0.0_REAL64
    END DO

    DO i=0,Ny+1
        grid(i,0) = 0.0_REAL64
        grid(i,(Nx+1)) = 0.0_REAL64
    END DO

!-------------------------------------------------Set axes to be indexed later-------------------------------------------

    axis_range(1)=-1.0_REAL64
    axis_range(2)=1.0_REAL64
    CALL create_axis(x_axis, Nx, axis_range, delta=dx)
    CALL create_axis(y_axis, Ny, axis_range, delta=dy)

!--------------------------------------------------Generate rho_grid-----------------------------------------------------
    
    lower=0_INT32
    ALLOCATE(rho_grid(0:Ny+1,0:Nx+1))
    rho_grid = grid
    IF(rho == 'null' .OR. rho == 'NULL' .OR. rho == 'Null') THEN !It is plausible that users will use upper and lower case values of init interchangeably. 
      CALL null(lower,rho_grid, pPos, pVel)                    !Thus, the code is written to allow this. 
        
        
    ELSE IF(rho == 'single' .OR. rho == 'SINGLE' .OR. rho == 'Single' ) THEN
      CALL single(lower, rho_grid,x_axis,y_axis, Nx, Ny, pPos, pVel)
            
    
    ELSE IF(rho == 'double' .OR. rho == 'DOUBLE' .OR. rho=='Double') THEN
      CALL double(lower, rho_grid,x_axis,y_axis,Nx, Ny, pPos, pVel)
    
            
    ELSE
      PRINT*,"Error: Invalid state initialisation requested, choose one of 'null', 'single' or 'double'. &
        &Write in lowercase only."           !If no valid state of initialisation is requested, this program will not work as 
      STOP                                   !intended. Thus, in response to an invalid input, a related error is printed, 
                                             !and the program is stopped.To minimise confusion, it is suggested to write in 
                                             !lowercase, even though other cases would work.
    END IF
    
    !---------------------------------------Generate phi_grid------------------------------------------------------------
    ALLOCATE(phi_grid(0:Ny+1,0:Nx+1))
    phi_grid = grid
    
    ALLOCATE(E_field_x(Ny,Nx))
    ALLOCATE(E_field_y(Ny,Nx))
    
    
    CALL phi(lower, phi_grid, Nx, Ny, rho_grid, dx, dy)
    CALL E_fields(E_field_x, E_field_y, phi_grid, dx, dy, Nx, Ny,lower)


    !-----------------------------------------------------------------
    ! Using the electric fields to set the inital acceleration
    CALL convertPos(pPos, pPosIndex, dx, dy)

    !---------------------------------------------------
    ! Initial acceleration:
    pAcc(1) = eCharge * E_field_x(pPosIndex(2), pPosIndex(1))
    pAcc(2) = eCharge * E_field_y(pPosIndex(2), pPosIndex(1))
    !---------------------------------------------------
    !-----------------------------------------------------------------
    
    
    !-----------------------------------------------------------------
    ! Saving the inital position, velocity and acceleration into
    ! the arrays that are outputted to netcdf.
    positions(:, 0) = pPos; velocities(:, 0) = pVel; accels(:, 0) = pAcc
    !-----------------------------------------------------------------


    !----------------------------------------------------------------- 
    ! Perform the verlet algorithm for 1000 iterations, and at each step
    ! save the positions, velocities and accelerations into the netcdf
    ! arrays:
    DO i = 1, iterations
      CALL verletAlg(pPos, pVel, pAcc, E_field_x, E_field_y, eCharge, dx, dy, dt)
      positions(:, i) = pPos; velocities(:, i) = pVel; accels(:, i) = pAcc
    END DO
    !----------------------------------------------------------------- 

    ALLOCATE(netcdf_rho(Nx, Ny))
    ALLOCATE(netcdf_phi(Nx, Ny))
    ALLOCATE(netcdf_Ex(Nx, Ny))
    ALLOCATE(netcdf_Ey(Nx, Ny))


    ! Transposing the arrays for netcdf output as per the specification:
    netcdf_rho = TRANSPOSE(rho_grid(1:Ny, 1:Nx))
    netcdf_phi = TRANSPOSE(phi_grid(1:Ny, 1:Nx))
    netcdf_Ex = TRANSPOSE(E_field_x(1:Ny, 1:Nx))
    netcdf_Ey = TRANSPOSE(E_field_y(1:Ny, 1:Nx))



    !----------------------------------------------------------------- 
    ! Call the netcdf writer subroutine:
    CALL write_sub(positions, velocities, accels, netcdf_rho, netcdf_phi, netcdf_Ex, netcdf_Ey, Nx, Ny, &
     & iterations, netcdf_filename, neterr, code_filename, rho)
    !----------------------------------------------------------------- 
    
    
    DEALLOCATE(netcdf_Ey)
    DEALLOCATE(netcdf_Ex)
    DEALLOCATE(netcdf_phi)
    DEALLOCATE(netcdf_rho)
    DEALLOCATE(rho_grid)
    DEALLOCATE(grid)
    DEALLOCATE(phi_grid)
    DEALLOCATE(E_field_x)
    DEALLOCATE(E_field_y)




END PROGRAM main



!-------------------Graveyard-----------------------------

    !print*,lbound(rho_grid)
    !print*,ubound(rho_grid)
    !PRINT*, rho_grid


    !PRINT*,rho_grid
    !PRINT*,'---------------------------------------------------------------------------'
    !PRINT*,phi_grid
    !PRINT*,'---------------------------------------------------------------------------'
    !Print*,E_field_x
    !PRINT*,E_field_y
!   ! Variable format depending on the number of columns.
!    write (formatName,'(A1,I4,A13)' ) "(",Nx,"(ES23.12,2X))"
!   ! print*,formatName
!
!    open (10 , file = "xE_field.txt" , form = 'formatted')
!   
!   
!    do i = 1, Ny
!        write (10,formatName) (E_field_x(i, j), j = 1, Nx)
!    end do 
!!  !  close(10)
!
!    open (11 , file = "yE_field.txt" , form = 'formatted')
!
!    !23 FORMAT (3 ( ES23 .12 E3 ) )
!    do i = 1, Ny
!        write (11,formatName) (E_field_y(i, j), j = 1, Nx)
!    end do
!    
!    close(10)
!    close(11)
