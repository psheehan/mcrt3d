module misc
    implicit none
    
    real(8), parameter :: pi = 3.1415927
    real(8), parameter :: AU = 1.496d13
    real(8), parameter :: pc = 3.086d18
    real(8), parameter :: R_sun = 6.955d10
    real(8), parameter :: M_sun = 1.98892d33
    real(8), parameter :: L_sun = 3.84e33
    real(8), parameter :: c_l = 2.99792458d10
    real(8), parameter :: h = 6.6260755d-27
    real(8), parameter :: k = 1.380658d-16
    real(8), parameter :: sigma = 5.67051d-5

    real(8), parameter :: hk = h/k
    real(8), parameter :: hh = 2*h
    real(8), parameter :: c_lc_l = c_l**2
    real(8), parameter :: hhc_lc_l = hh/c_lc_l
    real(8), parameter :: hhhh = h**2
    real(8), parameter :: c_lc_lk = c_lc_l*k
    real(8), parameter :: hhhhc_lc_lk = hhhh/(c_lc_lk)
    real(8), parameter :: invau = 1.0d0/au
contains
    real function random()
        random = rand()
    end function random

    real(8) function B_nu(T,nu)
        real(8) :: T, nu
        
        B_nu = 2.*h*nu**3/c_l**2*1./(exp(h*nu/(k*T))-1.)
    end function B_nu

    subroutine B_nu_arr(T,nu,n,Bnu)
        real(8) :: T
        integer :: n
        real(8), dimension(n) :: nu, Bnu
        real(8) :: hkT

        hkT = hk/T

        Bnu = hhc_lc_l*nu**3*1.0/(exp(hkT*nu)-1.)
    end subroutine B_nu_arr

    real(8) function dB_nu(T,nu)
        real(8) :: T, nu

        dB_nu = -2.*h**2*nu**4/(c_l**2*k*T**2)*1./(exp(h*nu/(k*T))-1.)/&
            &(1.-exp(-h*nu/(k*T)))
    end function dB_nu

    subroutine dB_nu_arr(T,nu,n,dBnu)
        real(8) :: T
        real(8), dimension(n) :: nu, dBnu
        real(8), dimension(n) :: expon, invexpon
        real(8) :: hkT, hhhhc_lc_lkTT
        integer :: n

        hkT = hk/T
        expon = exp(hkT*nu)
        invexpon = 1.0/(expon-1.)
        hhhhc_lc_lkTT = -2.0*hhhhc_lc_lk/(T**2)

        dBnu = -hhhhc_lc_lkTT*nu*nu*nu*nu*invexpon/(1.-exp(-hkT*nu))
    end subroutine dB_nu_arr

    real(8) function integrate(y,x)
        implicit none
        real(8), dimension(:) :: y, x
        integer :: n, n1

        n = size(y)
        n1 = n-1
        integrate = sum(0.5d0*(y(2:n)+y(1:n1))*(x(2:n)-x(1:n1)))
    end function integrate

    logical function equal(x,y,tol)
        implicit none
        real(8) :: x, y, tol

        if (abs(x-y).lt.abs(y)*tol) then
            equal = .true.
        else
            equal = .false.
        end if
    end function equal

    integer function find_in_arr(val,arr)
        implicit none
        real(8), dimension(:) :: arr
        real(8) :: val
        integer :: lmin, lmax, ltest
        logical :: not_found

        lmin = 1
        lmax = size(arr)
        not_found = .true.
        do while (not_found)
            ltest = (lmax-lmin)/2+lmin

            if ((val.ge.arr(ltest)).and.(val.lt.arr(ltest+1))) then
                find_in_arr = ltest
                not_found = .false.
            else
                if (val.lt.arr(ltest)) then
                    lmax = ltest
                else
                    lmin = ltest
                end if
            end if
        end do
    end function find_in_arr
    
end module misc
