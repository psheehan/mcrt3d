module dust
    use misc
    use photon_packet

    implicit none
    private
    public :: dust_grain, init_dust, clean_dust, planck_mean_opacity
    public :: dust_emit, dust_absorb, dust_isoscatt, opacity, random_nu
    public :: albedo, int_dBnu_knu

    type dust_grain
        sequence
        real(8), allocatable, dimension(:) :: nu, lam, kabs, ksca, kext, albedo
        real(8), allocatable, dimension(:) :: temp, planck_opacity
        real(8), allocatable, dimension(:) :: int_dBnu_knu
        real(8), allocatable, dimension(:) :: dplanck_opacity_dT, &
            &dint_dBnu_knu_DT
        real(8), allocatable, dimension(:,:) :: Bnu, dBnu, dBnudT, ddBnudT
        integer :: nlam
    end type dust_grain

contains
    
    subroutine init_dust(d,i)
        implicit none
        type(dust_grain) :: d
        real(8), allocatable, dimension(:) :: lam, kabs, ksca, dkextdnu, &
            &dalbedodnu
        real(8), allocatable, dimension(:) :: temp, planck_opacity, int_dBnu_knu
        real(8), allocatable, dimension(:) :: dplanck_opacity_dT, &
            &dint_dBnu_knu_DT
        real(8), allocatable, dimension(:,:) :: Bnu, dBnu, dBnudT, ddBnudT
        integer :: i, j, k, dummy, nlam, nspecies, ntemp
        character*512 :: dust_file

        open(unit=2,file="dust.dat")
        read(2,*) nspecies
        do j=1, i
            read(2,*) dust_file
        end do
        close(2)

        open(unit=2,file=trim(dust_file))
        read(2,*) dummy
        read(2,*) nlam

        d%nlam = nlam

        allocate (lam(nlam))
        allocate (kabs(nlam))
        allocate (ksca(nlam))
        allocate (dkextdnu(nlam))
        allocate (dalbedodnu(nlam))
        allocate (d%lam(nlam))
        allocate (d%nu(nlam))
        allocate (d%kabs(nlam))
        allocate (d%ksca(nlam))
        allocate (d%kext(nlam))
        allocate (d%albedo(nlam))

        do j=1, nlam
            read(2,*) lam(j), kabs(j), ksca(j)
        end do

        close(2)

        d%lam = lam*1.0e-4
        d%nu = c_l/(lam*1.0e-4)
        d%kabs = kabs
        d%ksca = ksca
        d%kext = kabs+ksca
        d%albedo = ksca/(kabs+ksca)

        open(unit=2,file="lookup_table3.txt")
        do j=1, nlam
            read(2,*) dkextdnu(j), dalbedodnu(j)
        end do
        close(2)

        deallocate (lam)
        deallocate (kabs)
        deallocate (ksca)
        deallocate (dkextdnu)
        deallocate (dalbedodnu)

        open(unit=2,file="lookup_table1.txt")
        read(2,*) ntemp

        allocate (d%temp(ntemp))
        allocate (d%planck_opacity(ntemp))
        allocate (d%int_dBnu_knu(ntemp))
        allocate (d%dplanck_opacity_dT(ntemp))
        allocate (d%dint_dBnu_knu_dT(ntemp))
        allocate (temp(ntemp))
        allocate (planck_opacity(ntemp))
        allocate (int_dBnu_knu(ntemp))
        allocate (dplanck_opacity_dT(ntemp))
        allocate (dint_dBnu_knu_dT(ntemp))

        do j=1, ntemp
            read(2,*) temp(j), planck_opacity(j), int_dBnu_knu(j), &
                &dplanck_opacity_dT(j), dint_dBnu_knu_dT(j)
        end do
        d%temp = temp
        d%planck_opacity = planck_opacity
        d%int_dBnu_knu = int_dBnu_knu
        d%dplanck_opacity_dT = dplanck_opacity_dT
        d%dint_dBnu_knu_dT = dint_dBnu_knu_dT

        close(2)

        allocate (d%Bnu(ntemp,nlam))
        allocate (d%dBnu(ntemp,nlam))
        allocate (d%dBnudT(ntemp,nlam))
        allocate (d%ddBnudT(ntemp,nlam))
        allocate (Bnu(ntemp,nlam))
        allocate (dBnu(ntemp,nlam))
        allocate (dBnudT(ntemp,nlam))
        allocate (ddBnudT(ntemp,nlam))

        open(unit=2,file="lookup_table2.txt")

        do j=1, ntemp
            do k=1, nlam
                read(2,*), Bnu(j,k), dBnu(j,k), dBnudT(j,k), ddBnudT(j,k)
            enddo
        enddo

        close(2)

        d%Bnu = Bnu
        d%dBnu = dBnu
        d%dBnudT = dBnudT
        d%ddBnudT = ddBnudT

        deallocate (temp, planck_opacity, int_dBnu_knu, dplanck_opacity_dT,&
            &dint_dBnu_knu_dT)
        deallocate(Bnu)
        deallocate(dBnu)
        deallocate(dBnudT)
        deallocate(ddBnudT)

    end subroutine init_dust

    subroutine dust_emit(d,p,x,y,z,T,bw,species)
        implicit none
        type(photon) :: p
        type(dust_grain) :: d
        real(8) :: x, y, z, T, cost, sint, phi, sinp, cosp
        logical :: bw
        type(dust_grain), allocatable, dimension(:) :: species
        integer :: i, nspecies

        p%x = x
        p%y = y
        p%z = z

        cost = 2*random()-1
        sint = sqrt(1.-cost**2)
        phi = 2.*3.1415927*random()
        sinp = sin(phi)
        cosp = cos(phi)

        p%nx = sint*cosp
        p%ny = sint*sinp
        p%nz = cost

        p%invnx = 1.0d0/p%nx
        p%invny = 1.0d0/p%ny
        p%invnz = 1.0d0/p%nz

        call random_nu(d,T,p,bw)

        nspecies = size(species)

        allocate (p%current_kext(nspecies))
        allocate (p%current_albedo(nspecies))

        do i=1, nspecies
            p%current_kext(i) = opacity(species(i),p%nu)
            p%current_albedo(i) = albedo(species(i),p%nu)
        end do
    end subroutine dust_emit

    subroutine dust_absorb(d,p,T,bw,species)
        implicit none
        type(photon) :: p
        type(dust_grain) :: d
        real(8) :: T, cost, sint, phi, cosp, sinp
        logical :: bw
        type(dust_grain), dimension(:) :: species
        integer :: i

        cost = 2*random()-1.
        sint = sqrt(1.-cost**2)
        phi = 2.*3.1415927*random()
        cosp = cos(phi)
        sinp = sin(phi)

        p%nx = sint*cosp
        p%ny = sint*sinp
        p%nz = cost

        p%invnx = 1.0d0/p%nx
        p%invny = 1.0d0/p%ny
        p%invnz = 1.0d0/p%nz

        call random_nu(d,T,p,bw)

        do i=1, size(species)
            p%current_kext(i) = opacity(d,p%nu)
            p%current_albedo(i) = albedo(d,p%nu)
        end do
    end subroutine dust_absorb

    subroutine dust_isoscatt(p)
        implicit none
        type(photon) :: p
        real(8) :: cost, sint, phi, cosp, sinp

        cost = 2*random()-1.
        sint = sqrt(1.-cost**2)
        phi = 2.*3.1415927*random()
        cosp = cos(phi)
        sinp = sin(phi)

        p%nx = sint*cosp
        p%ny = sint*sinp
        p%nz = cost

        p%invnx = 1.0d0/p%nx
        p%invny = 1.0d0/p%ny
        p%invnz = 1.0d0/p%nz
    end subroutine dust_isoscatt

    subroutine random_nu(d,T,p,bw)
        implicit none
        type(photon) :: p
        type(dust_grain) :: d
        real(8) :: nu, Prob, T, norm, ksi, tot
        real(8), allocatable, dimension(:) :: F
        integer :: i, n
        logical :: bw

        n = size(d%nu)

        allocate (F(n))
        if (bw) then
            call dBnu_arr(d,T,n,F)
            norm = 1.0d0/int_dBnu_knu(d,T)
        else
            !call B_nu_arr(T,d%nu,n,F)
            call Bnu_arr(d,T,n,F)
            norm = -pi/(planck_mean_opacity(d,T)*sigma*T**4)
        end if

        ksi = random()

        tot = 0.0d0
        do i=1, n-1
            tot = tot + 0.5d0*(d%nu(n-i+1)-d%nu(n-i))*&
                &(d%kext(n-i)*F(n-i)+d%kext(n-i+1)*F(n-i+1))
            Prob = tot*norm

            if (Prob.gt.ksi) then
                p%nu = random()*(d%nu(n-i+1)-d%nu(n-i))+&
                    &d%nu(n-i)
                return
            end if
        end do
        deallocate (F)

    end subroutine random_nu

    real(8) function opacity(d,nu)
        implicit none
        type(dust_grain) :: d
        real(8) :: nu
        integer :: l

        call find_freq_bin(d,nu,l)

        opacity = (d%kext(l+1)-d%kext(l))/(d%nu(l+1)-d%nu(l))*(nu-d%nu(l))+&
            &d%kext(l)

    end function opacity

    real(8) function albedo(d,nu)
        implicit none
        type(dust_grain) :: d
        real(8) :: nu
        integer :: l

        call find_freq_bin(d,nu,l)

        albedo = (d%albedo(l+1)-d%albedo(l))/(d%nu(l+1)-d%nu(l))*(nu-d%nu(l))+&
            &d%albedo(l)

    end function albedo

    real(8) function planck_mean_opacity(d,T)
        implicit none
        type(dust_grain) :: d
        real(8) :: T
        integer :: n

        n = find_in_arr(T,d%temp)

        planck_mean_opacity = d%dplanck_opacity_dT(n)*(T-d%temp(n))+&
            &d%planck_opacity(n)

    end function planck_mean_opacity

    real(8) function int_dBnu_knu(d,T)
        implicit none
        type(dust_grain) :: d
        real(8) :: T
        integer :: n

        n = find_in_arr(T,d%temp)

        int_dBnu_knu = d%dint_dBnu_knu_dT(n)*(T-d%temp(n))+d%int_dBnu_knu(n)

    end function int_dBnu_knu

    subroutine Bnu_arr(d,T,n,Bnu)
        implicit none
        type(dust_grain) :: d
        real(8) :: T
        integer :: n, i
        real(8), dimension(n) :: Bnu

        i = find_in_arr(T,d%temp)

        Bnu = d%dBnudT(i,:)*(T-d%temp(i))+d%Bnu(i,:)

    end subroutine Bnu_arr

    subroutine dBnu_arr(d,T,n,dBnu)
        implicit none
        type(dust_grain) :: d
        real(8) :: T
        integer :: n, i
        real(8), dimension(n) :: dBnu

        i = find_in_arr(T,d%temp)

        dBnu = d%ddBnudT(i,:)*(T-d%temp(i))+d%dBnu(i,:)

    end subroutine dBnu_arr

    subroutine find_freq_bin(d,nu,l)
        implicit none
        type(dust_grain) :: d
        real(8) :: nu
        integer :: l, i, lmin, lmax, ltest
        logical :: not_found

        lmin = 1
        lmax = size(d%nu)
        not_found = .true.
        i = 1
        do while (not_found.and.(i.le.10))
            ltest = (lmax-lmin)/2+lmin

            if (in_freq_bin(d,nu,ltest)) then
                l = ltest
                not_found = .false.
            else
                if (nu.lt.d%nu(ltest)) then
                    lmin = ltest
                else
                    lmax = ltest
                endif
            endif
            i = i+1
        end do
    end subroutine find_freq_bin

    logical function in_freq_bin(d,nu,l)
        implicit none
        type(dust_grain) :: d
        real(8) :: nu
        integer :: l

        if ((nu.le.d%nu(l)).and.(nu.gt.d%nu(l+1))) then
            in_freq_bin = .true.
        else
            in_freq_bin = .false.
        end if
    end function in_freq_bin

    subroutine clean_dust(d)
        implicit none
        type(dust_grain) :: d

        deallocate(d%lam)
        deallocate(d%nu)
        deallocate(d%kabs)
        deallocate(d%ksca)
        deallocate(d%kext)
        deallocate(d%albedo)
        deallocate(d%temp)
        deallocate(d%planck_opacity)
        deallocate(d%int_dBnu_knu)
        deallocate(d%dplanck_opacity_dT)
        deallocate(d%dint_dBnu_knu_dT)
        deallocate(d%Bnu)
        deallocate(d%dBnu)
        deallocate(d%dBnudT)
        deallocate(d%ddBnudT)

    end subroutine clean_dust

end module dust
