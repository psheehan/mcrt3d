program polarmc
    use omp_lib
    use grid_general
    use grid_geometry
    use photon_packet
    use lucy_iteration
    use bw_iteration
    use thermal_mc
    use dust
    use sources
    use misc

    implicit none

    type(grid) :: g
    integer :: nphot
    logical :: bw
    real :: t1,t2

    call init_grid(g)

    nphot = 1000000
    bw = .true.

    call cpu_time(t1)
    call do_thermal_mc(g,nphot,bw)
    call cpu_time(t2)

    print*, t2-t1

    call clean_grid(g)
    
end program polarmc
