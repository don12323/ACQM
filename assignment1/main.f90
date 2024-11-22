program main
        implicit none
        real :: alpha, rmax, dr
        integer :: N,l,i,j, rn, ierr
        real*8, allocatable :: r(:), basis_fns(:,:), wf(:,:), w(:), z(:,:), B(:,:), H(:,:)

        character(len=20) :: file_name

        !read input file
        open(unit=1, file='input.txt', status='old', action='read')
        read(1,*) alpha, l, N, dr, rmax
        close(1)

        !allocate memory for r and phikl
        rn = int(rmax/dr) +1
        allocate(r(rn))
        allocate(basis_fns(rn, N))

        !create radial coordinates
        do i = 1, rn
                r(i) = dr*real(i-1)
        end do

        !calculate eigen functions
        call calc_phikl (alpha, l, N, r, basis_fns)
        
        !write eigenfucntions to output file
        open(unit=2, file='eigenfns.txt', status='replace')
        do i=1, rn
               
                write(2,*) r(i), basis_fns(i,:) 
        end do
        close(unit=2) 
        !-----------------------------------------------------------
        !allocate memory
        allocate(B(N,N),H(N,N),w(N),wf(rn,N),z(N,N))

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

        !write energies to file for given N, l and alpha
        write(file_name, '(A,I0,A,I0,A,F4.2,A)') 'EN', N,'l',l,'a',alpha,'.txt'
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

        !write to file
        write(file_name, '(A,I1,A)') 'radfns_l',l,'.txt' 
        open(unit=4, file=file_name, status='replace')
        do i=1, rn
                write(4,*) r(i), wf(i,:) 
        end do
        close(4)
        




        deallocate(r,basis_fns,H,B,wf,w,z)
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
end program main

