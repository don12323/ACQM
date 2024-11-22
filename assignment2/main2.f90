program main2
        implicit none
        real*8 :: alpha, rmax, dr, mu
        integer :: N,l,i,j, rn, ierr, num_points
        real*8, allocatable :: r(:), integrand(:), basis_fns(:,:), wf(:,:), w(:), z(:,:), B(:,:), H(:,:), V(:,:), Vin(:), Vout(:), Rin(:), weights(:)
        character(len=20) :: file_name, line, mu_str

        !read input file
        open(unit=1, file='inputq2.txt', status='old', action='read')
        read(1,*) alpha, l, N, dr, rmax, mu
        close(1)

        !allocate memory for r and phikl
        rn = int(rmax/dr) +1
        allocate(r(rn))
        allocate(basis_fns(rn, N))

        !create radial coordinates
        do i = 1, rn
                r(i) = dr*real(i-1)
        end do
        !read data and interpolate data
        file_name = 'PEC.1ssg'
        num_points=0
        open(unit=10, file=file_name, status='old', action='read', iostat=ierr)
        if (ierr /= 0) then
                print *, 'Error opening file ', file_name
                stop
        endif
        do
                read(10, '(A)', iostat=ierr) line
                if (ierr /= 0) exit
                !count number of data
                num_points = num_points + 1
        end do
        close(10)
        
        allocate(Vin(num_points), Rin(num_points), Vout(rn),weights(rn))
        i=0
        open(unit=10, file=file_name, status='old', action='read', iostat=ierr)
        if (ierr /= 0) then
                print *, 'Error opening file ', file_name
                stop
        endif
        do      
                read(10, '(A)', iostat=ierr) line
                if (ierr /= 0) exit
                i = i + 1
                read(line, *) Rin(i), Vin(i)
        end do
        close(10)
        
        !interpolate data
        Vout=0.0d0
        call intrpl(num_points, Rin, Vin, rn, r, Vout)
        deallocate(Vin, Rin)
        
        print*, 'intrp V'
        do i=1, rn
                print*, Vout(i)
        end do
        
        !create weights for simpsons rule
        weights(1) = 1.0d0
        do i=2, rn-1
                weights(i) = 2.0d0 +2.0d0*mod(i+1,2)
        end do 
        weights(rn) = 1.0d0
        weights(:) = weights(:)*dr/3.0d0
        !calculate eigen functions
        call calc_phikl (alpha, l, N, r, basis_fns)
        
        !-----------------------------------------------------------
        !allocate memory
        allocate(B(N,N),H(N,N),V(N,N),w(N),wf(rn,N),z(N,N))

        B = 0.0d0
        
        
        do i=1, N-1
                B(i,i)=1.0d0
                B(i,i+1) = -0.5*SQRT(1.0-l*(l+1.0)/((i+l)*(i+l+1.0)))
                B(i+1,i) = B(i,i+1)
        end do
        B(N,N) = 1.0d0

        !calculate V matrix
        V = 0.0d0
        allocate(integrand(rn))
        integrand = 0.0d0
        do i=1, N
                do j=1, N
                        integrand = basis_fns(:,i) * Vout(:) * basis_fns(:,j)
                        V(i,j)= sum(integrand(:)*weights(:))
                end do
        end do
        !print Vmatrix
        print*, 'V'
        do i=1, N
                do j=1, N
                        write(*, '(F10.4, 1X)', advance='no') V(i,j)
                end do
                print*
        end do
         
        H = -0.5*alpha**2 * B / mu

        do i=1, N
                H(i,i) = alpha**2 / mu + H(i,i)
        end do

        H = H+V

        !solve eigenvalue problem
        call rsg(N,N,H,B,w,1,z,ierr)

        !write energies to file for given N, l and alpha
        write(file_name, '(A,I0,A,I0,A,F4.2,A)') 'EvibN', N,'l',l,'a',alpha,'.txt'
        open(unit=3, file= file_name, status='replace')
        do i=1, N
                write(3,'(I2, F10.6)') N, w(i)
        end do
        close(3)
        
        !calculating radial wavefunctions
        wf = 0.0d0
        do j=1, N
                do i=1, N
                        wf(:,j) = wf(:,j) + z(i,j)*basis_fns(:,i)
                end do 
        end do
        !write v=0 wfs to file
        write(mu_str, '(F12.5)') mu
        !write(file_name, '(A,F12.5,A)') 'mu_',mu,'.txt'
        file_name = 'mu_' // trim(adjustl(mu_str)) // '.txt'
        open(unit=1, file=file_name, status='replace')
        do i=1, rn
                write(1,*) r(i), wf(i,1)
        end do
        close(1)
        !scale and square wfs for poltting
        do i=1,N
                wf(:,i) = 0.005*wf(:,i)**2 + w(i)
        end do
        !write to file
        write(file_name, '(A,I1,A)') 'vibfns_l',l,'.txt' 
        open(unit=4, file=file_name, status='replace')
        do i=1, rn
                write(4,*) r(i),Vout(i),wf(i,:)
        end do
        close(4)
        




        deallocate(r,basis_fns,H,B,V,wf,w,z,Vout,weights,integrand)

        contains
                subroutine calc_phikl(alpha, l, N, r,basis_fns)
                        implicit none
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
end program main2

