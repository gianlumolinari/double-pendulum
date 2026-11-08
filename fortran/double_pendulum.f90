program double_pendulum
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    integer :: i, nsteps
    real(dp) :: t, dt, tf
    real(dp) :: yA(4), yB(4)
    real(dp) :: x1A, y1A, x2A, y2A
    real(dp) :: x1B, y1B, x2B, y2B
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    real(dp), parameter :: L1 = 1.0_dp, L2 = 1.0_dp

    tf = 20.0_dp
    dt = 0.01_dp
    nsteps = int(tf / dt)

    ! State vector:
    ! y(1) = theta1
    ! y(2) = omega1
    ! y(3) = theta2
    ! y(4) = omega2

    yA = [120.0_dp*pi/180.0_dp, 0.0_dp, -10.0_dp*pi/180.0_dp, 0.0_dp]

    ! Slightly different initial condition
    yB = [120.001_dp*pi/180.0_dp, 0.0_dp, -10.0_dp*pi/180.0_dp, 0.0_dp]

    open(unit=10, file="double_pendulum_output.csv", status="replace")

    write(10,*) "t,theta1A,theta2A,theta1B,theta2B,x1A,y1A,x2A,y2A,x1B,y1B,x2B,y2B"

    do i = 0, nsteps
        t = i * dt

        call positions(yA, L1, L2, x1A, y1A, x2A, y2A)
        call positions(yB, L1, L2, x1B, y1B, x2B, y2B)

        write(10,'(13(ES20.10,:,","))') t, yA(1), yA(3), yB(1), yB(3), &
                                      x1A, y1A, x2A, y2A, &
                                      x1B, y1B, x2B, y2B

        call rk4_step(yA, dt)
        call rk4_step(yB, dt)
    end do

    close(10)

    print *, "Simulation complete."
    print *, "Output written to double_pendulum_output.csv"

contains

    subroutine derivatives(y, dydt)
        implicit none
        real(dp), intent(in) :: y(4)
        real(dp), intent(out) :: dydt(4)

        real(dp) :: theta1, theta2, omega1, omega2
        real(dp) :: delta
        real(dp) :: den1, den2
        real(dp), parameter :: g = 9.81_dp
        real(dp), parameter :: m1 = 1.0_dp, m2 = 1.0_dp
        real(dp), parameter :: L1 = 1.0_dp, L2 = 1.0_dp

        theta1 = y(1)
        omega1 = y(2)
        theta2 = y(3)
        omega2 = y(4)

        delta = theta1 - theta2

        dydt(1) = omega1
        dydt(3) = omega2

        den1 = L1 * (2.0_dp*m1 + m2 - m2*cos(2.0_dp*delta))
        den2 = L2 * (2.0_dp*m1 + m2 - m2*cos(2.0_dp*delta))

        dydt(2) = ( &
            -g*(2.0_dp*m1 + m2)*sin(theta1) &
            -m2*g*sin(theta1 - 2.0_dp*theta2) &
            -2.0_dp*sin(delta)*m2*(omega2**2*L2 + omega1**2*L1*cos(delta)) &
            ) / den1

        dydt(4) = ( &
            2.0_dp*sin(delta)*( &
            omega1**2*L1*(m1 + m2) &
            + g*(m1 + m2)*cos(theta1) &
            + omega2**2*L2*m2*cos(delta) &
            ) ) / den2

    end subroutine derivatives


    subroutine rk4_step(y, dt)
        implicit none
        real(dp), intent(inout) :: y(4)
        real(dp), intent(in) :: dt
        real(dp) :: k1(4), k2(4), k3(4), k4(4)
        real(dp) :: yt(4)

        call derivatives(y, k1)

        yt = y + 0.5_dp * dt * k1
        call derivatives(yt, k2)

        yt = y + 0.5_dp * dt * k2
        call derivatives(yt, k3)

        yt = y + dt * k3
        call derivatives(yt, k4)

        y = y + dt * (k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4) / 6.0_dp

    end subroutine rk4_step


    subroutine positions(y, L1, L2, x1, yy1, x2, yy2)
        implicit none
        real(dp), intent(in) :: y(4)
        real(dp), intent(in) :: L1, L2
        real(dp), intent(out) :: x1, yy1, x2, yy2

        x1 = L1 * sin(y(1))
        yy1 = -L1 * cos(y(1))

        x2 = x1 + L2 * sin(y(3))
        yy2 = yy1 - L2 * cos(y(3))

    end subroutine positions

end program double_pendulum