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
        real(8) :: r, a, b, c, d, sr, st, sp, gnx, gny, gnz
        real(8), dimension(4) :: sta
        real(8), dimension(2) :: spa
        integer :: lr, lt, lp, i
        real(8), dimension(5) :: aarr, barr, carr, darr, sr1, sr2, st1, st2, sp1

        r = sqrt(p%x**2+p%y**2+p%z**2)

        ! Calculate distance to the intersection with the next radial wall.

        b = p%x*p%nx+p%y*p%ny+p%z*p%nz

        carr = r**2 - g%w1(p%l1-2:p%l1+2)**2
        darr = b**2 - carr

        where(darr.ge.0) sr1 = -b + sqrt(darr)
        where(darr.ge.0) sr2 = -b - sqrt(darr)
        where(darr.lt.0) sr1 = huge(b)
        where(darr.lt.0) sr2 = huge(b)
        
        where(sr1.lt.1.0e-6*au) sr1 = huge(b)
        where(sr2.lt.1.0e-6*au) sr2 = huge(b)

        sr = min(minval(sr1),minval(sr2))

        ! Calculate distance to intersection with the nearest theta wall.

        if (g%nw2.eq.2) then
            st = huge(b)
        else
            aarr = p%nx**2+p%ny**2+p%nz**2*tan(g%w2(p%l2-2:p%l2+2))**2
            barr = 2*(p%x*p%nx+p%y*p%ny+p%z*p%nz*tan(g%w2(p%l2-2:p%l2+2))**2)
            carr = p%x**2+p%y**2+p%z**2*tan(g%w2(p%l2-2:p%l2+2))**2
            darr = b**2-4*aarr*carr

            where(darr.ge.0) st1 = (-b + sqrt(darr))/(2*a)
            where(darr.ge.0) st2 = (-b - sqrt(darr))/(2*a)
            where(darr.lt.0) st1 = huge(b)
            where(darr.lt.0) st2 = huge(b)

            where(st1.lt.1.0e-6*au) st1 = huge(b)
            where(st2.lt.1.0e-6*au) st2 = huge(b)

            st = min(minval(st1),minval(st2))
        end if

        ! Calculate distance to intersection with the nearest phi wall.

        if (g%nw3.eq.2) then
            sp = huge(b)
        else
            aarr = p%x*sin(g%w3(p%l3-2:p%l3+2))-p%y*cos(g%w3(p%l3-2:p%l3+2))
            barr = p%nx*sin(g%w3(p%l3-2:p%l3+2))-p%ny*cos(g%w3(p%l3-2:p%l3+2))

            sp1 = aarr/barr

            where(sp1.lt.1.0e-6*au) sp1 = huge(b)

            sp = minval(sp1)
        end if

        next_wall_distance = min(sr,st,sp)
    end function next_wall_distance

    subroutine update_phot_loc(p,g)
        implicit none
        type(photon) :: p
        type(grid) :: g
        real(8) :: r, theta, phi, gnx, gny, gnz

        r = sqrt(p%x**2+p%y**2+p%z**2)
        theta = acos(p%z/r)
        phi = mod(atan2(p%y,p%x)+2*pi,2*pi)

        ! Find the location in the radial grid.

        gnx = sin(theta)*cos(phi)
        gny = sin(theta)*sin(phi)
        gnz = cos(theta)

        if (r.ge.g%w1(g%nw1)) then
            p%l1 = g%nw1
        else if (r.lt.g%w1(1)) then
            p%l1 = 1
        else
            if ((p%l1.eq.0).and.(p%l2.eq.0).and.(p%l3.eq.0)) then
                p%l1 = find_in_arr(r,g%w1)
            else
                if (p%l1.eq.1) then
                    p%l1 = p%l1 + find_in_arr(r,g%w1(p%l1:p%l1+3))-1
                else if (p%l1.eq.g%nw1-1) then
                    p%l1 = p%l1 + find_in_arr(r,g%w1(p%l1-2:p%l1+1))-3
                else
                    p%l1 = p%l1 + find_in_arr(r,g%w1(p%l1-2:p%l1+3))-3
                end if
            end if
        end if

        if ((equal(r,g%w1(p%l1),1.0d-6)).and.&
            &((p%nx*gnx+p%ny*gny+p%nz*gnz).lt.0)) then
            p%l1 = p%l1-1
        else if ((equal(r,g%w1(p%l1+1),1.0d-6)).and.&
            &((p%nx*gnx+p%ny*gny+p%nz*gnz).ge.0)) then
            p%l1 = p%l1+1
        end if

        ! Find the location in the theta grid.

        if (g%nw2.eq.2) then
            p%l2 = 1
        else
            if (theta.ge.g%w2(g%nw2)) then
                p%l2 = g%nw2-1
            else if (theta.lt.g%w2(1)) then
                p%l2 = 1
            else
                p%l2 = find_in_arr(theta,g%w2)
            end if

            gnx = cos(theta)*cos(phi)
            gny = cos(theta)*sin(phi)
            gnz = -sin(theta)

            if (equal(theta,g%w2(p%l2),1.0d-3).and.&
                &((p%nx*gnx+p%ny*gny+p%nz*gnz).lt.0).and.(p%l2.ne.1)) then
                p%l2 = p%l2-1
            else if (equal(theta,g%w2(p%l2+1),1.0d-3).and.&
                &((p%nx*gnx+p%ny*gny+p%nz*gnz).ge.0).and.&
                &(p%l2.ne.g%nw2-1)) then
                p%l2 = p%l2+1
            end if
        end if

        ! Find location in phi grid.

        if (g%nw3.eq.2) then
            p%l3 = 1
        else
            if (phi.ge.g%w3(g%nw3)) then
                p%l3 = 1
            else if (phi.lt.g%w3(1)) then
                p%l3 = g%nw3-1
            else
                p%l3 = find_in_arr(phi,g%w3)
            end if

            gnx = -sin(phi)
            gny = cos(phi)
            gnz = 0.0

            if (equal(phi,g%w3(p%l3),1.0d-3).and.&
                &((p%nx*gnx+p%ny*gny+p%nz*gnz).lt.0).and.(p%l3.ne.1)) then
                p%l3 = p%l3-1
            else if (equal(phi,g%w3(p%l3+1),1.0d-3).and.&
                &((p%nx*gnx+p%ny*gny+p%nz*gnz).ge.0).and.&
                &(p%l3.ne.g%nw3-1)) then
                p%l3 = p%l3+1
            end if
        end if
    end subroutine update_phot_loc

    logical function in_grid(p,g)
        implicit none
        type(photon) :: p
        type(grid) :: g
        real(8) :: r

        r = sqrt(p%x**2+p%y**2+p%z**2)

        if ((r.ge.g%w1(g%nw1)).or.(r.le.g%w1(1)).or.&
            &(equal(r,g%w1(g%nw1),1.0d-6)).or.(equal(r,g%w1(1),1.0d-6)))&
            then
            in_grid = .false.
        else
            in_grid = .true.
        end if
    end function in_grid

    real(8) function cell_mass(g,lr,lt,lp)
        implicit none
        type(grid) :: g
        integer :: lr, lt, lp

        cell_mass = 1./3.*g%dens(lr,lt,lp)*(g%w1(lr+1)**3-g%w1(lr)**3)*&
            &(g%w3(lp+1)-g%w3(lp))*(cos(g%w2(lt))-cos(g%w2(lt+1)))

        if (g%mirrorz.eq.2) then
            cell_mass = 2*cell_mass
        end if
    end function cell_mass

end module grid_geometry
