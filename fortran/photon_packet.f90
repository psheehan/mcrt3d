module photon_packet
    use misc

    implicit none
    private
    public :: photon,move,clean_photon
    
    type photon
        real(8) :: nu,x,y,z,nx,ny,nz,invnx,invny,invnz
        real(8), allocatable, dimension(:) :: current_kext, current_albedo
        integer :: l1, l2, l3
    end type photon
contains
    subroutine move(p,s)
        type(photon) :: p
        real(8), intent(in) :: s
        
        p%x = p%x + s*p%nx
        p%y = p%y + s*p%ny
        p%z = p%z + s*p%nz
    end subroutine move

    subroutine clean_photon(p)
        type(photon) :: p

        deallocate (p%current_kext)
        deallocate (p%current_albedo)

    end subroutine clean_photon
    
end module photon_packet
