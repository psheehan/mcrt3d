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
        sequence
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
        real(8) :: sx, sy, sz
        real(8), dimension(4) :: sx1, sy1, sz1

        ! Calculate distance to the intersection with the next x wall.

        sx1 = (g%w1(p%l1-1:p%l1+2)-p%x)*p%invnx

        where(sx1*invau.lt.1.0d-6) sx1 = huge(sx)

        sx = minval(sx1)

        ! Calculate distance to the intersection with the next y wall.

        sy1 = (g%w2(p%l2-1:p%l2+2)-p%y)*p%invny

        where(sy1*invau.lt.1.0e-6) sy1 = huge(sy)

        sy = minval(sy1)

        ! Calculate distance to the intersection with the next z wall.

        sz1 = (g%w3(p%l3-1:p%l3+2)-p%z)*p%invnz

        where(sz1*invau.lt.1.0e-6) sz1 = huge(sz)

        sz = minval(sz1)

        next_wall_distance = min(sx,sy,sz)

    end function next_wall_distance

    subroutine update_phot_loc(p,g)
        implicit none
        type(photon) :: p
        type(grid) :: g
        integer :: l, u

        ! Find photon location in the x grid.

        if (p%x.ge.g%w1(g%nw1)) then
            p%l1 = g%nw1
        else if (p%x.le.g%w1(1)) then
            p%l1 = 1
        else
            if (p%l1.eq.0) then
                p%l1 = find_in_arr(p%x,g%w1)
            else
                l = max(1,p%l1-1)
                u = min(g%nw1,p%l1+2)

                p%l1 = find_in_arr(p%x,g%w1(l:u))+l-1
            endif
        end if

        if ((equal(p%x,g%w1(p%l1),1.0d-6)).and.(p%nx.le.0)) then
            p%l1 = p%l1-1
        else if ((equal(p%x,g%w1(p%l1+1),1.0d-6)).and.(p%nx.ge.0)) then
            p%l1 = p%l1+1
        end if

        ! Find photon location in the y grid.

        if (p%y.ge.g%w2(g%nw2)) then
            p%l2 = g%nw2
        else if (p%y.le.g%w2(1)) then
            p%l2 = 1
        else
            if (p%l2.eq.0) then
                p%l2 = find_in_arr(p%y,g%w2)
            else
                l = max(1,p%l2-1)
                u = min(g%nw2,p%l2+2)

                p%l2 = find_in_arr(p%y,g%w2(l:u))+l-1
            endif
        end if

        if ((equal(p%y,g%w2(p%l2),1.0d-6)).and.(p%ny.le.0)) then
            p%l2 = p%l2-1
        else if ((equal(p%y,g%w2(p%l2+1),1.0d-6)).and.(p%ny.ge.0)) then
            p%l2 = p%l2+1
        end if

        ! Find photon location in the z grid.

        if (p%z.ge.g%w3(g%nw3)) then
            p%l3 = g%nw3
        else if (p%z.le.g%w3(1)) then
            p%l3 = 1
        else
            if (p%l3.eq.0) then
                p%l3 = find_in_arr(p%z,g%w3)
            else
                l = max(1,p%l3-2)
                u = min(g%nw3,p%l3+2)

                p%l3 = find_in_arr(p%z,g%w3(l:u))+l-1
            endif
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

        if ((p%x.ge.g%w1(g%nw1)).or.(p%y.ge.g%w2(g%nw2)).or.&
            &(p%z.ge.g%w3(g%nw3)).or.(equal(p%x,g%w1(g%nw1),1.0d-6)).or.&
            &(equal(p%y,g%w2(g%nw2),1.0d-6)).or.&
            &(equal(p%z,g%w3(g%nw3),1.0d-6))) then
            in_grid = .false.
        else if (((p%x.le.g%w1(1)).or.(equal(p%x,g%w1(1),1.0d-6)))&
            &.and.(g%mirrorx.eq.1)) then
            in_grid = .false.
        else if (((p%y.le.g%w2(1)).or.(equal(p%y,g%w2(1),1.0d-6)))&
            &.and.(g%mirrory.eq.1)) then
            in_grid = .false.
        else if (((p%z.le.g%w3(1)).or.(equal(p%z,g%w3(1),1.0d-6)))&
            &.and.(g%mirrorz.eq.1)) then
            in_grid = .false.
        else
            in_grid = .true.
        end if

        return

    end function in_grid

    real(8) function cell_mass(g,lx,ly,lz)
        implicit none
        type(grid) :: g
        integer :: lx, ly, lz

        cell_mass = g%dens(lx,ly,lz)*(g%w1(lx+1)-g%w1(lx))*&
            &(g%w2(ly+1)-g%w2(ly))*(g%w3(lz+1)-g%w3(lz))

        if (g%mirrorx.eq.2) then
            cell_mass = 2*cell_mass
        end if
        if (g%mirrory.eq.2) then
            cell_mass = 2*cell_mass
        end if
        if (g%mirrorz.eq.2) then
            cell_mass = 2*cell_mass
        end if
    end function cell_mass

end module grid_geometry
