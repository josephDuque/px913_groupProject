MODULE import_params
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

END MODULE import_params


MODULE initialization
    USE iso_fortran_env
    USE command_line
    IMPLICIT NONE
    CONTAINS
    
    SUBROUTINE null(lower, rho_grid)
        INTEGER :: lower

        REAL(REAL64), DIMENSION(lower:,lower:), INTENT(INOUT) :: rho_grid
                

        rho_grid = 0_REAL64

        !Need to set initial velocity and position

    END SUBROUTINE null


    SUBROUTINE single(lower, rho_grid,x_axis,y_axis,Nx, Ny)
        INTEGER, INTENT(IN) :: lower, Nx, Ny
        REAL(REAL64), DIMENSION(lower:,lower:), INTENT(INOUT) :: rho_grid
        REAL(REAL64), DIMENSION(:), INTENT(IN) :: x_axis, y_axis
        REAL(REAL64) :: x_div, y_div, z_1
        INTEGER :: i,j


        z_1 = 0.1_REAL64


            
        DO i=1,Nx
            x_div = ((x_axis(i))/z_1)        
            rho_grid(j,i) = EXP(-(x_div**2)-(y_div**2))
            DO j=1,Ny
                y_div = ((y_axis(j))/z_1)
            
            END DO
        
        END DO
        !Need to set initial velocity and position        
          
        

    END SUBROUTINE single


    SUBROUTINE double(lower, rho_grid,x_axis,y_axis,Nx, Ny)
        INTEGER, INTENT(IN) :: lower, Nx, Ny
        REAL(REAL64), DIMENSION(lower:,lower:), INTENT(INOUT) :: rho_grid
        REAL(REAL64), DIMENSION(:), INTENT(IN) :: x_axis, y_axis
        REAL(REAL64) :: x_div, y_div, extra_divx, extra_divy, z_1, z_2
        INTEGER :: i, j
        
        z_1 = 0.1_REAL64
        z_2 = 0.2_REAL64


        DO i=1,Nx
            x_div = ((x_axis(i)+0.25)/z_1)

            DO j=1,Ny
                y_div = ((y_axis(j)+0.25)/z_1)
                extra_divy = ((y_axis(j)-0.75)/z_2)
                extra_divx = ((x_axis(i)-0.75)/z_2)
                rho_grid(j,i) = EXP(-(x_div**2)-(y_div**2)) + EXP(-(extra_divx**2)-(extra_divy**2))
            
            END DO

        END DO
        !Need to set initial velocity and position

    END SUBROUTINE double

END MODULE initialization

MODULE phi_calc
    USE iso_fortran_env
    IMPLICIT NONE
    CONTAINS

    SUBROUTINE phi(lower, phi_grid, Nx, Ny, rho_grid, dx, dy)
        INTEGER, INTENT(IN) :: lower, Nx, Ny
        REAL(REAL64), INTENT(IN) :: dx, dy
        REAL(REAL64), DIMENSION(lower:,lower:), INTENT(INOUT) :: rho_grid, phi_grid
        REAL(REAL64) :: dens, div1, div2, two, denom
        REAL(REAL64) :: etot, conv1, conv2, drms_sum, drms
        INTEGER :: i, j, N, zed

        two = 2.0_REAL64

        etot = 10.0_REAL64
        drms = 1.0_REAL64

        zed=0

        DO WHILE(etot/drms > 0.00001)
            
            DO i=1,Nx
                DO j=1,Ny                
                    dens=rho_grid(j,i)                                      !Calculates phi values for whole grid once
                    div1=((phi_grid(j,i+1) + phi_grid(j,i-1))/(dx**2))
                    div2=((phi_grid(j+1,i) + phi_grid(j-1,i))/(dy**2))
                    denom=((two/(dx**two))+(two/(dy**two)))
                    phi_grid(i,j)=-((dens-div1-div2)/denom)

                END DO

            END DO

            N = Nx*Ny                                                        !Total no. of elements

            etot = 0.0_REAL64
            drms_sum = 0.0_REAL64       !Initialise
            drms = 0.0_REAL64
            
            DO i=1,Nx
                DO j=1,Ny
                                    dens=rho_grid(j,i)
                    conv1 = ((phi_grid(j,i-1) - (two*(phi_grid(j,i))) + phi_grid(j,i+1))/(dx**2))
                    conv2 = ((phi_grid(j-1,i) - (two*(phi_grid(j,i))) + phi_grid(j+1,i))/(dy**2))
                    
                    etot = etot + ABS(conv1 + conv2 - dens)

                    drms_sum = drms_sum + conv1 + conv2
                
                END DO
            
            END DO

            drms = SQRT((1/N)*(drms_sum))
            zed=zed+1
        
        END DO
        PRINT*,'zed=', zed
        PRINT*,etot
        PRINT*,'------------------------------------------------------------------------------'
        PRINT*,drms
        PRINT*,'------------------------------------------------------------------------------'
   


    END SUBROUTINE phi

END MODULE phi_calc


PROGRAM main
    USE iso_fortran_env
    USE import_params
    USE initialization
    !USE ascii_display
    USE domain_tools
    USE phi_calc
    IMPLICIT NONE
    INTEGER :: Ny, Nx, i, lower
    REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: grid, rho_grid, phi_grid
    REAL(REAL64), DIMENSION(:), ALLOCATABLE :: x_axis, y_axis
    REAL(REAL64), DIMENSION(2) :: axis_range
    REAL(REAL64) :: dx, dy
    CHARACTER(len=25) :: rho
    
    CALL params(Ny, Nx, rho) !This subroutine retrieves values inputted into command line for Nx and Ny parameters.

    ALLOCATE(grid(0:Ny+1,0:Nx+1))
!-------------------------------------------------Initialise All Cells, Then Set Ghost Cells-----------------------------------------------------
    grid=1
    
    DO i=0,Nx+1
        grid(0,i) = 0.0
        grid((Ny+1),i) = 0.0
    END DO

    DO i=0,Ny+1
        grid(i,0) = 0.0
        grid(i,(Nx+1)) = 0.0
    END DO

!-------------------------------------------------Set axes to be indexed later-------------------------------------------

    axis_range(1)=-1
    axis_range(2)=1
    CALL create_axis(x_axis, Nx, axis_range)
    !PRINT*,x_axis
    CALL create_axis(y_axis, Ny, axis_range)
    !PRINT*,y_axis
    dx=x_axis(2)-x_axis(1)
    dy=y_axis(2)-y_axis(1)
    !print*,dx
    !print*,dy

!--------------------------------------------------Generate rho_grid-----------------------------------------------------
    
    lower=0
    ALLOCATE(rho_grid(0:Ny+1,0:Nx+1))
    rho_grid = grid
    IF(rho == 'null' .OR. rho == 'NULL' .OR. rho == 'Null') THEN !It is plausible that users will use upper and lower case values of init interchangeably. Thus, the code is written to allow this. 
        CALL null(lower,rho_grid)
        
        
    ELSE
        IF(rho == 'single' .OR. rho == 'SINGLE' .OR. rho == 'Single' ) THEN
            CALL single(lower, rho_grid,x_axis,y_axis, Nx, Ny)
            
    
        ELSE 
            IF(rho == 'double' .OR. rho == 'DOUBLE' .OR. rho=='Double') THEN
                CALL double(lower, rho_grid,x_axis,y_axis,Nx, Ny)
    
            
            ELSE
                PRINT*,"Error: Invalid state initialisation requested, choose one of 'null', 'single' or 'double'. &
                &Write in lowercase only."   !If no valid state of initialisation is requested, this program will not work as 
                STOP                         !intended. Thus, in response to an invalid input, a related error is printed, 
                                             !and the program is stopped.To minimise confusion, it is suggested to write in 
            END IF                           !lowercase, even though other cases would work.
        
        END IF
    
    END IF
    
    !---------------------------------------Generate phi_grid------------------------------------------------------------
    ALLOCATE(phi_grid(0:Ny+1,0:Nx+1))
    phi_grid = grid
    CALL phi(lower, phi_grid, Nx, Ny, rho_grid, dx, dy)
    PRINT*,rho_grid
    PRINT*,'---------------------------------------------------------------------------'
    PRINT*,phi_grid
    PRINT*,'---------------------------------------------------------------------------'
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    !print*,lbound(rho_grid)
    !print*,ubound(rho_grid)
    !PRINT*, rho_grid
    DEALLOCATE(rho_grid)
    DEALLOCATE(grid)





END PROGRAM main