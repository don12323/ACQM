program main
        use constants
        implicit none
        real*8 :: alpha, dr, dk, energy
        integer :: N,l,i,j, nrmax, nkmax, fstate, ierr
        real*8, allocatable :: r(:), basis_fns(:,:), wf(:,:), w(:), z(:,:), B(:,:), H(:,:), rweights(:), &
                kgrid(:), Vdir(:,:), Vex(:,:), Von(:), kprime(:)

        character(len=20) :: filename
        !read input file
        open(unit=1, file='input.txt', status='old', action='read')
        read(1,*) alpha, l, N, nkmax, dk, fstate
        close(1)
        
        ! Read data.in file
        open(unit=2, file='data.in', status='old', action='read')
        read(2,*) energy
        read(2,*) nrmax, dr
        close(2)
        
        !convert to Ha
        energy = energy / eV
        print*, 'incident E', energy, 'kmax', nkmax, 'dk', dk, 'maxr', nrmax*dr, 'maxk', nkmax*dk
        !allocate memory for r and phikl
        allocate(r(nrmax), basis_fns(nrmax, N), rweights(nrmax), kgrid(nkmax))

        !create radial coordinates and setup rweights
        call setup_rgrids(nrmax, dr, r, rweights)

        !calculate radial Laguerre basis functions and store them in basis_fns matrix
        call calc_phikl (alpha, l, N, r, basis_fns)
        
        !-----------------------------------------------------------
        !allocate memory
        allocate(B(N,N),H(N,N),w(N),wf(nrmax,N),z(N,N))

        B = 0.0d0
        
        
        do i=1, N-1
                B(i,i)=1.0d0
                B(i,i+1) = -0.5*SQRT(1.0-l*(l+1.0)/((i+l)*(i+l+1.0)))
                B(i+1,i) = B(i,i+1)
        end do
        B(N,N) = 1.0d0
        
        H = -0.5*alpha**2 * B 

        do i=1, N
                H(i,i) = -alpha/(i+l) + alpha**2 + H(i,i)
        end do

        !solve eigenvalue problem
        call rsg(N,N,H,B,w,1,z,ierr)

        !calculating hydrogen radial wavefunctions for a given l
        wf = 0.0d0
        do j=1, N
                do i=1, N
                        wf(:,j) = wf(:,j) + z(i,j)*basis_fns(:,i)
                end do 
        end do
   
        print*, 'Normalisation 1s =', sum(rweights(:)*wf(:,1)**2)
        print*, 'Normalisation final =', sum(rweights(:)*wf(:,fstate)**2)
        !-------------------------------------------------------------- 
       !setup kgrid
        do i = 1, nkmax
            kgrid(i) = i * dk
        end do      

       !calculate kprime values 
       allocate(kprime(nkmax))
       print*, 'printing kprime'
        do i = 1, nkmax
            kprime(i) = sqrt(2.0d0 * (kgrid(i)**2 / 2.0d0 + w(1) - w(fstate)))   
            print*, kprime(i)
        end do 
        print*, '3sE', w(fstate) 
        print*, '1sE', w(1)

       !calculate direct V-mat elements
        allocate(Vdir(nkmax,nkmax), Von(nkmax))
        call calc_direct_V_mat_elements(nrmax, nkmax, fstate, dr, kgrid, r, wf, rweights, Vdir,Von, kprime)
         
       !write direct V-mat to file
        write(filename, '(A,I0,A)') 'Vdir_1s',fstate,'s.txt'
        call write_V_file(nkmax, Vdir, filename, kgrid)
       
       !write on shell values to file
        write(filename, '(A,I0,A)') 'OnSh_dir1s',fstate,'s.txt'
        call write_onshell_tofile(nkmax, kgrid, Von, filename)
        
        !calculate exchange V-mat elements
        allocate(Vex(nkmax,nkmax))
        call calc_exchange_V_mat_elements(nrmax, nkmax, fstate, dr, kgrid, r, wf, rweights, Vex, Von, kprime, energy-0.5)
        write(filename, '(A,I0,A)') 'Vex_1s',fstate,'s.txt'
        call write_V_file(nkmax, Vex, filename, kgrid)
        !write on shell values to file
        write(filename, '(A,I0,A)') 'OnSh_ex1s',fstate,'s.txt'
        call write_onshell_tofile(nkmax, kgrid, Von, filename)
        !write 1s wavefunction to file for q2
        write(filename, '(A)') 'wf.txt'
        call write_onshell_tofile(nrmax, r, wf(:,1), filename)

        deallocate(r,basis_fns,H,B,wf,w,z, rweights, kgrid, Vdir,Vex,Von, kprime)
        contains

                subroutine calc_phikl(alpha, l, N, r,basis_fns)
                        integer, intent(in) :: N,l
                        real,intent(in) :: alpha, r(:)
                        real*8, dimension(:,:), intent(out) :: basis_fns

                        integer :: k,i
                        real*8 :: fact
                        
                        !calculate phi_1,l phi_2,l
                        basis_fns(:,1) = (2*alpha*r)**(l+1) * exp(-alpha*r)
                        basis_fns(:,2) = 2*(l+1-alpha*r)*(2*alpha*r)**(l+1) * exp(-alpha*r) 
                        
                        !loop to make remaining basis functions
                        do k = 3, N
                                basis_fns(:,k) = (2*(k-1+l-alpha*r)*basis_fns(:,k-1) - (k+2*l-1)*basis_fns(:,k-2))/(k-1) 
                        end do
                        !calculate normalisation constants and find non_orth basis functions
                        do k = 1, N
                                fact = k+2*l
                                do i = 1, 2*l
                                        fact = fact*(k+2*l-i) 
                                end do
                                basis_fns(:,k) = SQRT(alpha/((k+l)*fact))*basis_fns(:,k)
                        end do
                end subroutine calc_phikl

                subroutine setup_rgrids(nrmax, dr, r, rweights)
                        integer, intent(in) :: nrmax
                        real*8, intent(in) :: dr
                        real*8, intent(out) :: rweights(nrmax), r(nrmax)
                        integer :: ir

                        do ir = 1, nrmax
                            r(ir) = dr*ir
                            if (mod(ir, 2) == 0) then
                                rweights(ir) = 2.0d0 * dr / 3.0d0
                            else
                                rweights(ir) = 4.0d0 * dr / 3.0d0
                            end if
                        end do
                end subroutine setup_rgrids
                
                subroutine calc_direct_V_mat_elements(nrmax, nkmax, fstate, dr, kgrid, r, wf, rweights, V, Von, kprime)
                        integer, intent(in) :: nrmax, nkmax, fstate
                        real*8, intent(in) :: kgrid(:), r(:), wf(:, :), rweights(:), dr, kprime(:)
                        real*8, intent(out) :: V(:,:), Von(:)
                        integer :: i, j, ki, kj
                        real*8 :: f(nrmax), g(nrmax), krond
                        real*8 :: A1(nrmax), A2(nrmax)
                        
                        ! Check if there is a transition
                        krond = 1.0d0
                        if (fstate /= 1) krond = 0.0d0

                        ! Calculate g(r2)
                        g(:) = wf(:,fstate) * wf(:, 1)
                        
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
                        print*, 'krond', krond
                        ! Loop over k' and k
                        do ki = 1, nkmax
                                do kj = 1, nkmax
                                        ! Calculate f(r1) 
                                        f(:) = sin(kgrid(ki) * r) * sin(kgrid(kj) * r)
                                        
                                        ! Loop over r1
                                        V(ki, kj) = 2.0d0 / pi * sum(f(:) * (A1(:) / r(:) + A2(:) - krond / r(:)) * rweights(:))
                                        
                                end do
                        end do

                        !loop to calculate onshell Values
                        do ki = 1, nkmax
                                        f(:) = sin(kgrid(ki) * r) * sin(kprime(ki) * r)
                                        Von(ki) = 2.0d0 / pi * sum(f(:) * (A1(:) / r(:) + A2(:) - krond / r(:)) * rweights(:))
                                        print*, Von(ki)
                        end do
                end subroutine calc_direct_V_mat_elements
                
                subroutine calc_exchange_V_mat_elements(nrmax, nkmax, fstate, dr, kgrid, r, wf, rweights, V, Von, kprime, E)
                        integer, intent(in) :: nrmax, nkmax, fstate
                        real*8, intent(in) :: kgrid(:), r(:), wf(:, :), rweights(:), dr, kprime(:), E
                        real*8, intent(out) :: V(:,:), Von(:)
                        integer :: i, j, ki, kj
                        real*8 :: f(nkmax,nrmax), g(nkmax,nrmax), V11, Oii 
                        real*8 :: A1(nkmax,nrmax), A2(nkmax,nrmax)
                        real*8 :: Oi(nkmax), Of(nkmax), V1(nkmax), V2(nkmax)
                        
                        ! Calculate g(r2), f(r1), A1, A2 and other matrix element terms
                        do i = 1, nkmax
                                f(i,:) = sin(kgrid(i) * r(:)) * wf(:,1)
                                Of(i) = sum(sin(kgrid(i) * r(:)) * wf(:,1) * rweights(:))
                                V1(i) = sum(-sin(kgrid(i) * r(:)) *wf(:,1) * rweights(:) / r(:))
                                Oi(i) = sum(sin(kgrid(i) * r(:)) * wf(:,fstate) * rweights(:))
                                V2(i) = sum(-sin(kgrid(i) * r(:)) * wf(:,fstate) * rweights(:) / r(:))
                                g(i,:) = wf(:,fstate) * sin(kgrid(i) * r)
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
                        
                        !loop to calculate onshell Values
                        do i = 1, nkmax
                                        !we want to recalculate the k' terms from before
                                        f(i,:) = sin(kprime(i) * r(:)) * wf(:,1)
                                        Of(i) = sum(sin(kprime(i) * r(:)) * wf(:,1) * rweights(:))
                                        V1(i) = sum(-sin(kprime(i) * r(:)) *wf(:,1) * rweights(:) / r(:))

                                        Von(i) = sum(f(i,:) * (A1(i,:) / r(:) + A2(i,:)) * rweights(:))
                                        Von(i) = (kgrid(i)**2 / 2.0d0 -0.5 - kprime(i)**2 / 2.0d0 - kgrid(i)**2 / 2.0d0) * Oi(i) * Of(i) &
                                                -V1(i) * Oi(i) - Of(i) * V2(i) - Von(i)
                        end do
                        Von = 2.0/pi * Von
                end subroutine calc_exchange_V_mat_elements
                
                subroutine write_V_file(nkmax, V, filename, kgrid)
                        integer, intent(in) :: nkmax
                        real*8, intent(in) :: V(:,:), kgrid(:)
                        character(len=20), intent(in) :: filename
                        integer :: ki, kj
                        open(unit=3, file=filename, status='replace')

                        do ki = 1, nkmax
                            do kj = 1, nkmax
                                write(3,*) kgrid(ki), kgrid(kj), V(ki, kj)
                            end do
                            write(3,*) ''
                        end do

                        close(3)
                end subroutine write_V_file
                
                subroutine write_onshell_tofile(nkmax, kgrid, Von, filename)
                        integer, intent(in) :: nkmax
                        real*8, intent(in) :: kgrid(:), Von(:)
                        character(len=20), intent(in) :: filename
                        integer :: k
                        open(unit=3, file=filename, status='replace')
                        do k=1, nkmax
                                write(3,*) kgrid(k), Von(k)
                        enddo
                        close(3)
                end subroutine write_onshell_tofile


end program main

