module bw_iteration
    use grid_general
    use grid_geometry
    use photon_packet
    use dust
    use sources
    use misc

    private
    public :: do_bw_iteration

contains

    subroutine do_bw_iteration(g,nphot,pcount)
        implicit none
        type(grid) :: g
        type(photon) :: p
        real(8) :: tau
        integer :: nphot, i, j, maxniter
        real(8), dimension(:,:,:) :: pcount
        logical, parameter :: bw = .true.
        real(8) :: t1, t2, t3

        maxniter = 1.0e8
        t3 = 0.0
        do i=1, nphot
            if (mod(i,nphot/10).eq.0) then
                print*, i
            end if

            call emit(p,g)

            j = 0
            do while (in_grid(p,g).and.(j.le.maxniter))
                tau = -log(1-random())

                call propagate_photon(p,g,pcount,tau,.false.)

                if (in_grid(p,g)) then
                    if (random().lt.&
                            &p%current_albedo(g%dust(p%l1,p%l2,p%l3))) then
                        call dust_isoscatt(p)
                    else
                        pcount(p%l1,p%l2,p%l3) = pcount(p%l1,p%l2,p%l3) + 1
                        call update_grid(g,nphot,pcount,p%l1,p%l2,p%l3)
                        call dust_absorb(g%dust_species(g%dust(p%l1,p%l2,&
                            &p%l3)),p,g%temp(p%l1,p%l2,p%l3),bw,g%dust_species)
                    end if
                end if

                j = j+1

                if (j.eq.maxniter) then
                    print*, "Killing Photon"
                end if
            end do
            call clean_photon(p)
        end do
        print*, t3

    end subroutine do_bw_iteration

end module bw_iteration
