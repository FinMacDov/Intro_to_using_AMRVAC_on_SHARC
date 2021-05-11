module mod_usr
  use mod_hd

  implicit none

contains

  subroutine usr_init()

    usr_init_one_grid => rm1d_init_one_grid
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output 

    call set_coordinate_system("Cartesian")
    call hd_activate()

  end subroutine usr_init

  ! Initialize one grid
  subroutine rm1d_init_one_grid(ixG^L,ix^L,w,x)
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    ! iprob==1 rarefaction wave & shock
    if (iprob==1) then
        where (abs(x(ix^S,1))<0.5d0)
           w(ix^S,rho_)   = 1.5d0
           w(ix^S,mom(1)) = 0.0d0
           w(ix^S,e_)     = 1.5d0
        elsewhere
           w(ix^S,rho_)   = 0.5d0
           w(ix^S,mom(1)) = 0.0d0
           w(ix^S,e_)     = 0.5d0
        end where
    ! iprob==2  shock & shock
    else if (iprob==2) then
        where (abs(x(ix^S,1))<0.5d0)
           w(ix^S,rho_)   = 1.0d0
           w(ix^S,mom(1)) = 2.0d0
           w(ix^S,e_)     = 1.5d0
        elsewhere
           w(ix^S,rho_)   = 1.0d0
           w(ix^S,mom(1)) = -2.0d0
           w(ix^S,e_)     = 1.5d0
        end where
    ! iprob==3  rarefaction wave
    else if (iprob==3) then
        where (abs(x(ix^S,1))<0.5d0)
           w(ix^S,rho_)   = 1.0d0
           w(ix^S,mom(1)) = -0.5d0
           w(ix^S,e_)     = 1.0d0
        elsewhere
           w(ix^S,rho_)   = 1.0d0
           w(ix^S,mom(1)) = 0.5d0
           w(ix^S,e_)     = 1.0d0
        end where
    else
        call mpistop("iprob not available!")
    end if

    call hd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine rm1d_init_one_grid

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: pth(ixI^S), wlocal(ixI^S,1:nw)
    integer :: idirmin,idir,ix^D

    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    ! output temperature
    call hd_get_pthermal(wlocal,x,ixI^L,ixO^L,pth)
    w(ixO^S,nw+1)=dsqrt(hd_gamma*pth(ixO^S)/w(ixO^S,rho_))
  
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames
    varnames='cs'

  end subroutine specialvarnames_output


end module mod_usr
