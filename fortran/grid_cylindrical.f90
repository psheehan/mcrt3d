module grid_geometry
    use photon_packet
    use dust
    use sources
    use misc

    implicit none
    private
    public :: grid, in_grid, update_phot_loc, cell_mass
    public :: next_wall_distance

    type grid
        real(8), allocatable, dimension(:) :: w1, w2, w3
        integer :: nw1, nw2, nw3
        integer :: mirrorx, mirrory, mirrorz
        real(8), allocatable, dimension(:,:,:) :: dens, temp, mass
        integer, allocatable, dimension(:,:,:) :: dust
        type(dust_grain), allocatable, dimension(:) :: dust_species
        integer :: nspecies
        type(spherical_source), allocatable, dimension(:) :: sources
        real(8) :: total_lum
    end type grid

contains

    real(8) function next_wall_distance(g,p)
        implicit none
        type(grid) :: g
        type(photon) :: p
        real(8) :: r, a, b, c, d, sr, sp, sz, gnx, gny
        real(8), dimension(2) :: spa
        integer :: lr, lp, lz, i, l, u
        real(8), dimension(5) :: carr, darr, sr1, sr2, sp1, sz1

        r = sqrt(p%x**2+p%y**2)

        ! Calculate distance to the intersection with the next radial wall.

        a = p%nx**2+p%ny**2
        b = p%x*p%nx+p%y*p%ny

        carr = r**2 - g%w1(p%l1-2:p%l1+2)**2
        darr = b**2 - a*carr

        where(darr.ge.0) sr1 = (-b + sqrt(darr))/a
        where(darr.ge.0) sr2 = (-b - sqrt(darr))/a
        where(darr.lt.0) sr1 = huge(b)
        where(darr.lt.0) sr2 = huge(b)

        where(sr1.lt.1.0e-6*au) sr1 = huge(b)
        where(sr2.lt.1.0e-6*au) sr2 = huge(b)

        sr = min(minval(sr1),minval(sr2))

        ! Calculate distance to intersection with the nearest theta wall.

        if (g%nw2.eq.2) then
            sp = huge(b)
        else
            carr = p%x*sin(g%w2(p%l2-2:p%l2+2))-p%y*cos(g%w2(p%l2-2:p%l2+2))
            darr = p%nx*sin(g%w2(p%l2-2:p%l2+2))-p%ny*cos(g%w2(p%l2-2:p%l2+2))

            sp1 = carr/darr

            where(sp1.lt.1.0e-6*au) sp1 = huge(b)

            sp = minval(sp1)
        end if

        ! Calculate distance to nearest z wall.

        sz1 = (g%w3(p%l3-2:p%l3+2)-p%z)/p%nz

        where(sz1.lt.1.0e-6*au) sz1 = huge(sz)

        sz = minval(sz1)

        next_wall_distance = min(sr,sp,sz)
    end function next_wall_distance

    subroutine update_phot_loc(p,g)
        implicit none
        type(photon) :: p
        type(grid) :: g
        real(8) :: r, phi, gnx, gny

        !call get_cell(g,p%x,p%y,p%z,lr,lp,lz)

        r = sqrt(p%x**2+p%y**2)
        phi = mod(atan2(p%y,p%x)+2*pi,2*pi)

        ! Find the location in the radial grid.

        gnx = cos(phi)
        gny = sin(phi)

        if (r.ge.g%w1(g%nw1)) then
            p%l1 = g%nw1
        else if (r.lt.g%w1(1)) then
            p%l1 = 1
        else
            !if ((p%l1.eq.0).and.(p%l2.eq.0).and.(p%l3.eq.0)) then
                p%l1 = find_in_arr(r,g%w1)
            !else
            !    if (p%l1.eq.1) then
            !        p%l1 = p%l1 + find_in_arr(r,g%w1(p%l1:p%l1+3))-1
            !    else if (p%l1.eq.g%nw1-1) then
            !        p%l1 = p%l1 + find_in_arr(r,g%w1(p%l1-2:p%l1+1))-3
            !    else
            !        p%l1 = p%l1 + find_in_arr(r,g%w1(p%l1-2:p%l1+3))-3
            !    end if
            !end if
        end if

        if ((equal(r,g%w1(p%l1),1.0d-6)).and.((p%nx*gnx+p%ny*gny).le.0)) then
            p%l1 = p%l1-1
        else if ((equal(r,g%w1(p%l1+1),1.0d-6)).and.&
            &((p%nx*gnx+p%ny*gny).ge.0)) then
            p%l1 = p%l1+1
        end if

        ! Find the location in the phi grid.

        if (g%nw2.eq.2) then
            p%l2 = 1
        else
            if (phi.ge.g%w2(g%nw2)) then
                p%l2 = 1
            else if (phi.lt.g%w2(1)) then
                p%l2 = g%nw2-1
            else
                p%l2 = find_in_arr(phi,g%w2)
            end if

            gnx = -sin(phi)
            gny = cos(phi)

            if (equal(phi,g%w2(p%l2),1.0d-3).and.&
                &((p%nx*gnx+p%ny*gny).lt.0).and.(p%l2.ne.1)) then
                p%l2 = p%l2-1
            else if (equal(phi,g%w2(p%l2+1),1.0d-3).and.&
                &((p%nx*gnx+p%ny*gny).ge.0).and.&
                &(p%l2.ne.g%nw2-1)) then
                p%l2 = p%l2+1
            end if
        end if

        ! Find the location in the z grid.
        
        if (p%z.ge.g%w3(g%nw3)) then
            p%l3 = g%nw3
        else if (p%z.le.g%w3(1)) then
            p%l3 = 1
        else
            p%l3 = find_in_arr(p%z,g%w3)
        end if

        if ((equal(p%z,g%w3(p%l3),1.0d-6)).and.(p%nz.le.0)) then
            p%l3 = p%l3-1
        else if ((equal(p%z,g%w3(p%l3+1),1.0d-6)).and.(p%nz.ge.0)) then
            p%l3 = p%l3+1
        end if

    end subroutine update_phot_loc

    logical function in_grid(p,g)
        implicit none
        type(photon) :: p
        type(grid) :: g
        real(8) :: r

        r = sqrt(p%x**2+p%y**2)

        if ((r.ge.g%w1(g%nw1)).or.(r.le.g%w1(1)).or.&
            &(equal(r,g%w1(g%nw1),1.0d-6)).or.&
            &(equal(r,g%w1(1),1.0d-6)).or.(p%z.ge.g%w3(g%nw3)).or.&
            &(p%z.le.g%w3(1)).or.(equal(p%z,g%w3(1),1.0d-6)).or.&
            &(equal(p%z,g%w3(g%nw3),1.0d-6))) then
            in_grid = .false.
        else
            in_grid = .true.
        end if
    end function in_grid

    real(8) function cell_mass(g,lr,lp,lz)
        implicit none
        type(grid) :: g
        integer :: lr, lp, lz

        cell_mass = 0.5*g%dens(lr,lp,lz)*(g%w1(lr+1)**2-g%w1(lr)**2)*&
            &(g%w2(lp+1)-g%w2(lp))*(g%w3(lz+1)-g%w3(lz))
    end function cell_mass

end module grid_geometry
