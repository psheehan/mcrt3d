module grid_general
    use photon_packet
    use grid_geometry
    use dust
    use sources
    use misc

    implicit none
    private
    public :: init_grid, emit, propagate_photon, cell_lum
    public :: update_grid, write_temperature, clean_grid

contains

    subroutine init_grid(g)
        type(grid) ::  g
        real(8), allocatable, dimension(:) :: g1, g2, g3
        real(8), allocatable, dimension(:,:,:) :: dens
        real(8), allocatable, dimension(:) :: Bnu
        integer :: i,j,k,n1,n2,n3,nspecies

        open(unit=2,file="grid.dat")
        read(2,*) n1, g%mirrorx
        read(2,*) n2, g%mirrory
        read(2,*) n3, g%mirrorz

        g%nw1 = n1
        g%nw2 = n2
        g%nw3 = n3

        allocate (g1(n1))
        allocate (g2(n2))
        allocate (g3(n3))
        allocate (g%w1(n1))
        allocate (g%w2(n2))
        allocate (g%w3(n3))

        allocate (dens(n1-1,n2-1,n3-1))
        allocate (g%dens(n1-1,n2-1,n3-1))
        allocate (g%temp(n1-1,n2-1,n3-1))
        allocate (g%mass(n1-1,n2-1,n3-1))
        allocate (g%dust(n1-1,n2-1,n3-1))

        do i=1, n1
            read(2,*) g1(i)
        enddo
        do i=1, n2
            read(2,*) g2(i)
        enddo
        do i=1, n3
            read(2,*) g3(i)
        enddo
        close(2)

        g%w1 = g1
        g%w2 = g2
        g%w3 = g3

        open(unit=2,file="dens.dat")
        do i=1, n1-1
            do j=1, n2-1
                do k=1, n3-1
                    read(2,*) dens(i,j,k)
                enddo
            enddo
        enddo
        close(2)

        g%dens = dens

        g%temp(:,:,:) = 1.

        do i=1, n1-1
            do j=1, n2-1
                do k=1, n3-1
                    g%mass(i,j,k) = cell_mass(g,i,j,k)
                enddo
            enddo
        enddo

        deallocate (g1)
        deallocate (g2)
        deallocate (g3)
        deallocate (dens)

        g%dust(:,:,:) = 1

        open(unit=2,file="dust.dat")
        read(2,*) nspecies
        close(2)

        g%nspecies = nspecies

        allocate (g%dust_species(nspecies))

        do i=1, nspecies
            call init_dust(g%dust_species(i),i)
        enddo

        allocate (g%sources(1))

        do i=1, size(g%sources)
            call init_source(g%sources(i),g%dust_species(1)%nu,M_sun,R_sun,&
                &4000.d0,0d0,0d0,0d0)
        end do

    end subroutine init_grid

    subroutine emit(p,g)
        implicit none
        type(photon) :: p
        type(grid) :: g

        call source_emit(g%sources(1),p,g%nspecies,g%dust_species)

        p%l1 = 0
        p%l2 = 0
        p%l3 = 0
        call update_phot_loc(p,g)

    end subroutine emit

    subroutine propagate_photon(p,g,pcount,tau,lucy)
        implicit none
        type(photon) :: p
        type(grid) :: g
        real(8) :: tau, s, s1, s2
        real(8), dimension(:,:,:) :: pcount
        logical :: lucy

        do while ((tau.gt.1.0d-6).and.in_grid(p,g))

            s1 = next_wall_distance(g,p)

            s2 = tau/(p%current_kext(g%dust(p%l1,p%l2,p%l3))*&
                &g%dens(p%l1,p%l2,p%l3))

            s = min(s1,s2)

            call move(p,s)

            if (lucy) then
                pcount(p%l1,p%l2,p%l3) = pcount(p%l1,p%l2,p%l3) + &
                    &s*p%current_kext(g%dust(p%l1,p%l2,p%l3))*&
                    &g%dens(p%l1,p%l2,p%l3)
            end if

            tau = tau - s*p%current_kext(g%dust(p%l1,p%l2,p%l3))*&
                &g%dens(p%l1,p%l2,p%l3)

            if (s1.lt.s2) then
                call update_phot_loc(p,g)
            end if
        end do

    end subroutine propagate_photon

    real(8) function cell_lum(g,l1,l2,l3)
        implicit none
        type(grid) :: g
        integer :: l1, l2, l3

        cell_lum = 4.0d0*pi*g%mass(l1,l2,l3)*&
            &planck_mean_opacity(g%dust_species(g%dust(l1,l2,l3)),&
            g%temp(l1,l2,l3))*sigma*g%temp(l1,l2,l3)**4/pi

    end function cell_lum

    subroutine update_grid(g,nphot,pcount,l1,l2,l3)
        implicit none
        type(grid) :: g
        integer :: l1, l2, l3, nphot
        real(8) :: T_old
        real(8), dimension(:,:,:) :: pcount
        logical :: not_converged

        not_converged = .true.

        do while (not_converged)
            T_old = g%temp(l1,l2,l3)

            g%temp(l1,l2,l3) = (pcount(l1,l2,l3)*g%sources(1)%luminosity/&
                &(4.0d0*sigma*nphot*&
                &planck_mean_opacity(g%dust_species(g%dust(l1,l2,l3)),T_old)*&
                &g%mass(l1,l2,l3)))**(0.25)

            if ((abs(T_old-g%temp(l1,l2,l3))/T_old.lt.1.0d-2)) then
                not_converged = .false.
            end if
        end do

    end subroutine update_grid

    subroutine write_temperature(g)
        implicit none
        type(grid) :: g
        integer :: i, j, k

        open(unit=2,file="temp.dat")
        do i=1, size(g%w1)-1
            do j=1, size(g%w2)-1
                do k=1, size(g%w3)-1
                    write(2,*) g%temp(i,j,k)
                enddo
            enddo
        enddo
        close(2)

    end subroutine write_temperature

    subroutine clean_grid(g)
        implicit none
        type(grid) :: g
        integer :: i

        do i=1, size(g%dust_species)
            call clean_dust(g%dust_species(i))
        end do

        do i=1, size(g%sources)
            call clean_source(g%sources(i))
        end do

        deallocate (g%w1)
        deallocate (g%w2)
        deallocate (g%w3)
        deallocate (g%dens)
        deallocate (g%temp)
        deallocate (g%mass)
        deallocate (g%dust)
        deallocate (g%dust_species)
        deallocate (g%sources)

    end subroutine clean_grid

end module grid_general
