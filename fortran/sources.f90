module sources
    use misc
    use photon_packet
    use dust

    implicit none
    private
    public :: spherical_source, init_source, clean_source
    public :: source_emit

    type spherical_source
        sequence
        real(8) :: mass, radius, temperature, luminosity, x, y, z
        real(8), allocatable, dimension(:) :: nu, Bnu
        integer :: nnu
    end type spherical_source

contains

    subroutine init_source(s,nu,mass,radius,temperature,x,y,z)
        implicit none
        real(8) :: mass, radius, temperature, x, y, z
        real(8), allocatable, dimension(:) :: nu
        type(spherical_source) :: s
        integer :: i

        s%x = x
        s%y = y
        s%z = z
        s%mass = mass
        s%radius = radius
        s%temperature = temperature
        s%luminosity = 4*pi*s%radius**2*sigma*s%temperature**4

        s%nnu = size(nu)
        allocate (s%nu(s%nnu))
        allocate (s%Bnu(s%nnu))

        s%nu = nu
        call B_nu_arr(s%temperature,s%nu,s%nnu,s%Bnu)

    end subroutine init_source

    subroutine source_emit(s,p,nspecies,species)
        implicit none
        type(photon) :: p
        type(spherical_source) :: s
        real(8) :: theta, phi, cost, sint, cosp, sinp
        integer :: i, ev, nspecies
        type(dust_grain), dimension(nspecies) :: species

        theta = pi*random()
        phi = 2*pi*random()

        p%x = s%radius*sin(theta)*cos(phi)
        p%y = s%radius*sin(theta)*sin(phi)
        p%z = s%radius*cos(theta)

        cost = 2*random()-1
        sint = sqrt(1.-cost**2)
        phi = 2*pi*random()
        sinp = sin(phi)
        cosp = cos(phi)

        p%nx = sint*cosp
        p%ny = sint*sinp
        p%nz = cost
        p%invnx = 1.0d0/p%nx
        p%invny = 1.0d0/p%ny
        p%invnz = 1.0d0/p%nz

        call get_random_nu(s,p)

        allocate(p%current_kext(nspecies))
        allocate(p%current_albedo(nspecies))

        do i=1, nspecies
            p%current_kext(i) = opacity(species(i),p%nu)
            p%current_albedo(i) = albedo(species(i),p%nu)
        end do

    end subroutine source_emit

    subroutine get_random_nu(s,p)
        implicit none
        type(spherical_source) :: s
        type(photon) :: p
        real(8) :: tot, Prob, norm, ksi
        integer :: i, n

        ksi = random()
        norm = -pi/(sigma*s%temperature**4)
        n = s%nnu
        tot = 0.0
        do i=1, n-1
            tot = tot + 0.5*(s%nu(n-i+1)-s%nu(n-i))*(s%Bnu(n-i)+s%Bnu(n-i+1))
            Prob = tot*norm

            if (Prob.gt.ksi) then
                p%nu = random()*(s%nu(n-i+1)-s%nu(n-i))+s%nu(n-i)
                return
            end if
        end do

    end subroutine get_random_nu

    subroutine clean_source(s)
        implicit none
        type(spherical_source) :: s

        deallocate (s%nu)
        deallocate (s%Bnu)

    end subroutine clean_source

end module sources
