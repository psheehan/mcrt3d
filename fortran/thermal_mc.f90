module thermal_mc
    use lucy_iteration
    use bw_iteration
    use grid_general
    use grid_geometry
    use dust
    use sources

    private
    public :: do_thermal_mc
contains

    subroutine do_thermal_mc(g,nphot,bw)
        implicit none
        type(grid) :: g, gold, greallyold
        integer :: nphot, i, maxniter
        logical :: bw
        real(8), allocatable, dimension(:,:,:) :: pcount

        allocate (pcount(g%nw1-1,g%nw2-1,g%nw3-1))

        pcount(:,:,:) = 0.0

        if (bw) then
            call do_bw_iteration(g,nphot,pcount)
        else
            gold = g
            maxniter = 20
            i = 1
            do while (i.le.maxniter)
                print*, "Starting iteration #", i
                print*, ""

                greallyold = gold
                gold = g

                call do_lucy_iteration(g,nphot,pcount)

                if (i.gt.2) then
                    if (converged(g,gold,greallyold)) then
                        i = maxniter
                    end if
                end if

                i = i + 1
                print*, ""
            end do
        end if
        call write_temperature(g)

        deallocate (pcount)

    end subroutine do_thermal_mc

    logical function converged(g,gold,greallyold)
        implicit none
        type(grid) :: g, gold,greallyold
        real(8), allocatable, dimension(:,:,:) :: R, Rold
        real(8) :: Qthresh, Delthresh, p, Q, Qold, Del
        integer :: i, j, k

        Qthresh = 2.0
        Delthresh = 1.1
        p = 0.99

        allocate(R(size(g%w1)-1,size(g%w2)-1,size(g%w3)-1))
        allocate(Rold(size(g%w1)-1,size(g%w2)-1,size(g%w3)-1))

        do i=1, size(g%w1)-1
            do j=1, size(g%w2)-1
                do k=1, size(g%w3)-1
                    R(i,j,k) = delta(gold%temp(i,j,k),g%temp(i,j,k))
                    Rold(i,j,k) = delta(greallyold%temp(i,j,k),gold%temp(i,j,k))
                end do
            end do
        end do

        Q = quantile(R,p)
        Qold = quantile(Rold,p)
        print*, Q, Qold

        Del = delta(Qold,Q)
        print*, Del

        converged = (Q.lt.Qthresh).and.(Del.lt.Delthresh)

        deallocate (R)
        deallocate (Rold)
        
    end function converged

    real(8) function delta(x1,x2)
        implicit none
        real(8) :: x1, x2

        delta = max(x1/x2,x2/x1)
    end function delta

    real(8) function quantile(R,p)
        implicit none
        real(8), dimension(:,:,:) :: R
        real(8), allocatable, dimension(:) :: Rline
        real(8) :: p
        integer :: i, j, k, nx, ny, nz

        allocate (Rline(size(R)))

        nx = size(R(:,1,1))
        ny = size(R(1,:,1))
        nz = size(R(1,1,:))

        do i=1, nx
            do j=1, ny
                do k=1, nz
                    Rline((i-1)*ny*nz+(j-1)*nz+k) = R(i,j,k)
                end do
            end do
        end do

        call sort(Rline)

        quantile = Rline(int(p*size(Rline)))

        deallocate (Rline)

    end function quantile

    subroutine sort(arr)
        real(8), dimension(:) :: arr
        real(8) :: a
        
        do j=2, size(arr)
            a=arr(j)
            do i=j-1,1,-1
                if (arr(i)<=a) goto 10
                arr(i+1)=arr(i)
            end do
            i=0
10          arr(i+1)=a
        end do
    end subroutine sort

end module thermal_mc
