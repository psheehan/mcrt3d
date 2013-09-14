module lucy_iteration
    use grid_general
    use grid_geometry
    use photon_packet
    use dust
    use sources
    use misc

    private
    public :: do_lucy_iteration

contains

    subroutine do_lucy_iteration(g,nphot,pcount)
        implicit none
        type(grid) :: g
        type(photon) :: p
        real(8) :: tau
        integer :: nphot, i, j, k, maxniter
        real(8), dimension(:,:,:) :: pcount
        logical, parameter :: bw = .false.

        maxniter = 1.0e8
        do i=1, nphot
            if (mod(i,nphot/10).eq.0) then
                print*, i
            end if

            call emit(p,g)

            j = 0
            do while (in_grid(p,g).and.(j.le.maxniter))
                tau = -log(1-random())

                call propagate_photon(p,g,pcount,tau,.true.)

                if (in_grid(p,g)) then
                    if (random().lt.p%current_albedo(g%dust(p%l1,p%l2,p%l3))) then
                        call dust_isoscatt(p)
                    else
                        call dust_absorb(g%dust_species(g%dust(p%l1,p%l2,&
                            &p%l3)),p,g%temp(p%l1,p%l2,p%l3),bw,g%dust_species)
                    end if
                endif

                j = j+1

                if (j.eq.maxniter) then
                    print*, "Killing photon"
                end if
            end do
            call clean_photon(p)
        end do

        do i=1, size(g%w1)-1
            do j=1, size(g%w2)-1
                do k=1, size(g%w3)-1
                    call update_grid(g,nphot,pcount,i,j,k)
                end do
            end do
        end do

        pcount(:,:,:) = 0.0

    end subroutine do_lucy_iteration

end module lucy_iteration
