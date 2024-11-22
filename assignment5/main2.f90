!  This file is an outline of a code to perform a calculation of a
!    charged particle scattered by a spherically-symmetric short-range local
!    potential.
!  There are a number of comments throughout the file beginning with !>>> which
!    indicate sections of code you will need to complete.

! Last modified: May 25 2024
!                Liam Scarlett

program main

  use constants
  implicit none

  real*8, allocatable :: &
    Vmat(:,:),       & !V-matrix elements Vmat(kf,ki)
    rgrid(:),        & !radial grid
    rweights(:),     & !radial integration weights
    kgrid(:),        & !momentum-space grid
    kweights(:),     & !momentum-space integration weights (and Green's func)
    wf1s(:)          !wavefunction for 1s
    
  real*8 :: &
    rmax,   & !max value of radial grid
    dr,     & !radial grid step size
    energy, & !projectile energy
    k,x,theta,    & !projectile momentum
    kg_A, kg_B, kg_P !some parameters for setting up kgrid

  complex*16, allocatable :: Ton(:) !on-shell T-matrix element per l
  character(len=20) :: filename

  integer :: &
    nrmax,      & !number of rgrid points
    nkmax,      & !number of kgrid points
    zproj,      & !projectile charge 
    l,          & !partial-wave angular momentum
    lmin, lmax, & !min and max values of l
    iounit,     & !a unit number for input/ouput
    kg_Na, kg_Nb, kg_Np, & !some more parameters for setting up kgrid
    n, i
  !set kgrid parameters - leave this as is
    kg_Na = 4; kg_Nb = 6; kg_Np = 2
    kg_a = 0.9; kg_b = 2.5; kg_p = 4.0
    nkmax=kg_Na+kg_Nb+kg_Np+1

  !>>> open data.in file and read input parameters
  !    note: energy should be read in electron volts
  !      and grid parameters in atomic units
    open(unit=1, file='data.in', status='old')
    read(1,*) energy
    read(1,*) nrmax, dr
    read(1,*) zproj, lmin, lmax, theta
    close(1)
  !>>> do any input validation you think is necessary here
    print*, 'Eev=', energy, 'nr=', nrmax, 'dr=',dr,'z=',zproj,'lmaxlmin=',lmax,lmin
  !>>> convert the energy to atomic units and calculate the
  !      projectile momentum
    energy = energy / eV  
    k = sqrt(2.0 * energy)    
    print*, 'k-momentum=', k

  !>>> determine number of rgrid points nrmax
  !    note: nrmax should be even for simpson's integration
  !          to take into account that the r=0 point has been omitted
    if (mod(nrmax, 2) /= 0) then
        nrmax = nrmax+1
    end if
    print*, 'new nrmax = ', nrmax

  !allocate memory
    allocate(rgrid(nrmax),rweights(nrmax))
    allocate(kgrid(nkmax),kweights(nkmax))
    allocate(Ton(lmin:lmax))
    allocate(Vmat(nkmax,nkmax))

  !setup grids
    call setup_rgrid(nrmax, dr, rgrid, rweights)
    call setup_kgrid(k, nkmax, kg_Na, kg_a, kg_Nb, kg_b, kg_Np, kg_p, kgrid, kweights)
  ! read wavefunction from file
    allocate(wf1s(nrmax))    
    open(unit=10, file='wf.txt', status='old', action='read')
    do i = 1, nrmax
         read(10, *, iostat=iounit) x, wf1s(i) !x is a dummy variable
        if (iounit /= 0) then
              print *, 'Error reading file at line', i
              stop
        end if
    end do
    

  !begin loop over angular momenta
  do l=lmin, lmax
    !calculate V matrix
      call calculate_Vmatrix(nkmax,kgrid,nrmax,rgrid,rweights,wf1s,energy-0.5,theta,Vmat)

    !solve the Lippman-Schwinger equation for the on-shell T-matrix
      call tmatrix_solver(nkmax,kgrid,kweights,Vmat,Ton(l),theta)
  enddo

   deallocate(rgrid,rweights,kgrid,kweights,Ton,Vmat, wf1s)
end program main


subroutine setup_rgrid(nrmax, dr, rgrid, rweights)
  implicit none
  integer, intent(in) :: nrmax
  real*8, intent(in) :: dr
  real*8, intent(out) :: rgrid(nrmax), rweights(nrmax)
  integer :: ir !index to iterate over r

   do ir=1, nrmax
        rgrid(ir) = ir * dr
        if (mod(ir, 2) == 0) then
                rweights(ir) = 2.0d0 * dr / 3.0d0
        else 
                rweights(ir) = 4.0d0 * dr / 3.0d0
        end if
   end do
   

end subroutine setup_rgrid


subroutine calculate_Vmatrix(nkmax,kgrid,nrmax,rgrid,rweights,wf,energy,theta,Vmat)
  use constants
  implicit none
  integer, intent(in) :: nkmax, nrmax
  real*8, intent(in) :: kgrid(nkmax), rgrid(nrmax), rweights(nrmax), wf(nrmax), energy, theta
  real*8, intent(out) :: Vmat(nkmax,nkmax)
  real*8 :: Vex(nkmax,nkmax), Vdir(nkmax,nkmax),E
  integer :: i,j

  print*, "total energy=", energy,'Ha'
  E = energy*(1.0d0 - 2.0d0*theta)
  Vmat = 0.0d0
  Vdir = 0.0d0
  Vex=0.0d0
  call calc_direct_V_mat_elements(nrmax, nkmax, kgrid, rgrid, wf, rweights, Vdir)      
  call calc_exchange_V_mat_elements(nrmax, nkmax, kgrid, rgrid, wf, rweights, Vex, E)           
  
  Vmat = Vdir+Vex
end subroutine calculate_Vmatrix

subroutine calc_direct_V_mat_elements(nrmax, nkmax, kgrid, r, wf, rweights, V)
        use constants
        integer, intent(in) :: nrmax, nkmax
        real*8, intent(in) :: kgrid(nkmax), r(nrmax), wf(nrmax), rweights(nrmax)
        real*8, intent(out) :: V(nkmax,nkmax)
        integer :: i, j, ki, kj
        real*8 :: f(nrmax), g(nrmax)
        real*8 :: A1(nrmax), A2(nrmax)

        ! Check if there is a transition

        V = 0.0d0
        ! Calculate g(r2)
        g(:) = wf(:) * wf(:)

        ! Calculate A1 array
        A1(1) = rweights(1) * g(1)
        do i = 2, nrmax
                A1(i) = A1(i-1) + rweights(i) * g(i)
        end do
        ! Calculate A2 array
        A2(nrmax) = rweights(nrmax) * g(nrmax) / r(nrmax)
        do i = nrmax-1, 1, -1
                A2(i) = A2(i+1) + rweights(i) * g(i) / r(i)
        end do
        ! Loop over k' and k
        do ki = 1, nkmax
                do kj = 1, nkmax
                        ! Calculate f(r1)
                        f(:) = sin(kgrid(ki) * r) * sin(kgrid(kj) * r)

                        ! Loop over r1
                        V(ki, kj) = 2.0d0 / pi * sum(f(:) * (A1(:) / r(:) + A2(:) - 1.0d0 / r(:)) * rweights(:))

                end do
        end do
end subroutine calc_direct_V_mat_elements

subroutine calc_exchange_V_mat_elements(nrmax, nkmax, kgrid, r, wf, rweights, V, E)
        use constants                
        integer, intent(in) :: nrmax, nkmax
        real*8, intent(in) :: kgrid(nkmax), r(nrmax), wf(nrmax), rweights(nrmax), E
        real*8, intent(out) :: V(nkmax,nkmax)
        integer :: i, j, ki, kj
        real*8 :: f(nkmax,nrmax), g(nkmax,nrmax)
        real*8 :: A1(nkmax,nrmax), A2(nkmax,nrmax)
        real*8 :: Oi(nkmax), Of(nkmax), V1(nkmax), V2(nkmax)

        ! Check if there is a transition

        V = 0.0d0
        ! Calculate g(r2), f(r1), A1, A2 and other matrix element terms
        do i = 1, nkmax
                f(i,:) = sin(kgrid(i) * r(:)) * wf(:)
                Of(i) = sum(sin(kgrid(i) * r(:)) * wf(:) * rweights(:))
                V1(i) = sum(-sin(kgrid(i) * r(:)) *wf(:) * rweights(:) / r(:))
                Oi(i) = sum(sin(kgrid(i) * r(:)) * wf(:) * rweights(:))
                V2(i) = sum(-sin(kgrid(i) * r(:)) * wf(:) * rweights(:) / r(:))
                g(i,:) = wf(:) * sin(kgrid(i) * r)
                ! Calculate A1 array
                A1(i,1) = rweights(1) * g(i,1)
                do j = 2, nrmax
                        A1(i,j) = A1(i,j-1) + rweights(j) * g(i,j)
                end do
                ! Calculate A2 array
                A2(i,nrmax) = rweights(nrmax) * g(i,nrmax) / r(nrmax)
                do j = nrmax-1, 1, -1
                        A2(i,j) = A2(i,j+1) + rweights(j) * g(i,j) / r(j)
                end do
        end do

        ! Loop over k' and k
        do ki = 1, nkmax
                do kj = 1, nkmax
                        ! Calculate V12 term
                        V(ki, kj) = sum(f(ki,:) * (A1(kj,:) / r(:) + A2(kj,:)) * rweights(:))
                        ! Calculate V(k',k)
                        V(ki,kj) = (E - kgrid(ki)**2 / 2.0d0 - kgrid(kj)**2 / 2.0d0) * Oi(kj) * Of(ki) &
                               -V1(ki) * Oi(kj) - Of(ki) * V2(kj) - V(ki,kj)
                end do
        end do
        V = 2.0d0 / pi * V
end subroutine calc_exchange_V_mat_elements
 
    
subroutine tmatrix_solver(nkmax,kgrid,kweights,Vmat,Ton, theta)
  use constants
  implicit none
  integer, intent(in) :: nkmax
  real*8, intent(in) :: kgrid(nkmax), kweights(nkmax), Vmat(nkmax,nkmax), theta
  complex*16, intent(out) :: Ton !on-shell T-matrix element
  real*8 :: &
    Koff(nkmax-1), & !half-off-shell K-matrix elements
    Kon,           & !on-shell K-matrix element
    Von,           & !on-shell V-matrix element
    A(nkmax-1,nkmax-1)!Coefficient matrix for the linear system Ax=b
  integer :: j, ipiv(nkmax-1), info,n
  character(len=20) :: filename

  !>>> store the on-shell V-matrix element in Von
  Von = Vmat(1,1)
  !>>> populate the matrix A according to Eq (142) in the slides
  A = 0.0d0
  do j = 1, nkmax-1
        do n = 1, nkmax-1
            if (j == n) then
                A(j, n) = 1.0d0 - kweights(n+1) * Vmat(n+1, n+1)
            else
                A(j, n) = -kweights(n+1) * Vmat(j+1, n+1)
            end if
        end do
  end do 

  !>>> populate the vector Koff with the half-on-shell V-matrix elements (RHS of Eq (141)) 
  do n=1, nkmax-1
        Koff(n) = Vmat(n+1,1)
  end do
  !Here is the call to DGESV
  call dgesv( nkmax-1, 1, A, nkmax-1, ipiv, Koff, nkmax-1, info )
  if(info /= 0) then
    print*, 'ERROR in dgesv: info = ', info
  endif

  !>>> Now use the half-on-shell K matrix which has been stored in Koff to get the on-shell K-matrix element Kon
  Kon = Von + sum(kweights(2:nkmax) * Vmat(1, 2:nkmax) * Koff)
  !>>> And then use Kon to get the on-shell T-matrix element Ton
  Ton = Kon / (1.0d0 + (0.0d0, 1.0d0) * pi / kgrid(1) * Kon)
  print*, 'Kon',Kon
  print*, 'Ton',Ton
   
  !write to file for plotting
  write(filename, '(A, I0, A)') "k1V_thet", int(theta), ".txt"
  call write_to_file(nkmax, kgrid, Vmat(2:nkmax,1), filename, Von)
  write(filename, '(A, I0, A)') "k1K_thet", int(theta), ".txt"
  call write_to_file(nkmax, kgrid, Koff, filename, Kon)

end subroutine tmatrix_solver

subroutine setup_kgrid(k,nkmax,Na,a,Nb,b,Np,p,kgrid,kweights)
  implicit none
  real*8, intent(in) :: k
  integer, intent(in) :: Na, Nb, Np
  integer, intent(in) :: nkmax
  real*8, intent(out) :: kgrid(nkmax), kweights(nkmax)
  integer :: nk
  real*8, intent(in) ::  a, b, p
  real*8 :: grid1(nkmax-1), weight1(nkmax-1)

  call kgrid_igor(0.0,k,a,Na,b,Nb,p,Np,nkmax-1,grid1,weight1)
  
  kgrid(1) = k
  kgrid(2:nkmax) = grid1
  kweights(1) = 0.0d0
  kweights(2:nkmax) = weight1
  do nk=2, nkmax
    kweights(nk) = 2.0d0* kweights(nk) / (k**2 - kgrid(nk)**2)
  enddo
end subroutine setup_kgrid

subroutine write_to_file(nkmax, kgrid, M, filename, Mon)
        integer, intent(in) :: nkmax
        real*8, intent(in) :: kgrid(nkmax), M(nkmax-1), Mon
        character(len=20), intent(in) :: filename
        integer :: k
        open(unit=3, file=filename, status='replace')
        do k=1, 8
                write(3,'(F10.4, F10.4)') kgrid(k+1), M(k)
        enddo
        write(3, '(F10.4,F10.4)') kgrid(1), Mon
        do k=9,nkmax -1   
                write(3,'(F10.4, F10.4)') kgrid(k+1), M(k)
        end do
        close(3)
end subroutine write_to_file
