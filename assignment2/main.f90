program main
        implicit none
        real*8 :: alpha, rmax, dr,Rn, Y,integral, temp, E1,E2
        real*8, external :: Yint
        integer :: N,lmax,l,i,j,li,lj, k,nr, num_func,m_max,par,m,lamd,nstates,nstate,ierr
        real*8, allocatable :: r(:), basis_fns(:,:), wf(:,:), w(:), z(:,:), B(:,:), H(:,:), V(:,:), weights(:),f(:), energies(:)
        integer, allocatable :: k_list(:), l_list(:)
        character(len=20) :: file_name
        logical :: sorted


        !read input file
        open(unit=1, file='input.txt', status='old', action='read')
        read(1,*) alpha, lmax, N, dr, rmax, Rn
        close(1)
        
        !allocate memory for r
        nr = int(rmax/dr) +1
        if (mod(nr,2) == 0) nr = nr + 1
        print*, "nr=", nr
        allocate(r(nr), weights(nr),f(nr))

        !create radial coordinates
        do i = 1, nr
                r(i) = dr*real(i-1)
        end do

        !weights 
        weights(1) = 1.0d0
        do i=2, nr-1
                weights(i) = 2.0d0 + 2.0d0*mod(i+1,2)
        end do
        weights(nr) = 1.0d0
        weights(:) = weights(:)*dr/3.0d0

        m_max = lmax
        nstates = 2*m_max +1
        allocate(energies(nstates))
        nstate = 0


        !loop over symmetry pairs
        do par=-1, 1, 2
                do m=0, m_max
                        write(*,'(A,I2,A,I2)') "par =", par, ", m =", m                        
                       !calculate number of basis functions for each symmetry pair
                        num_func = 0
                        do l=0, lmax
                                if((-1)**l /=par .or. l<m) cycle
                                num_func = num_func + N 
                                write(*,'(A,I1)') "l =",l 
                        end do 
                        
                        write(*,'(A,I2)') "num_func=",num_func
                        
                        if (num_func == 0) then
                                write(*,'(A, I1)') "skipping m=", m
                                cycle        !for cases when num_func == 0 cycle to next m
                        end if
                        
                        !keep track of indexes in basis function matrix
                        allocate(basis_fns(nr,num_func), k_list(num_func), l_list(num_func))
                        i=0
                        do l=0, lmax
                                if((-1)**l /= par .or. l < m) cycle
                                do k=1, N
                                        i = i + 1
                                        k_list(i) = k
                                        l_list(i) = l
                                enddo
                        enddo
                        !calculate basis functions
                        do i=0, num_func/N-1
                                call calc_phikl (alpha, l_list(N*i+1), N, r, basis_fns(:,N*i+1:N*i+N))         
                        end do
                        
                        print*, "Basis" 
                        do i=1, nr
                                do j=1,num_func
                                        write(*, '(F10.4, 1X)', advance='no') basis_fns(i,j)
                               end do
                               print*
                        end do
                        
                        !allocate H, V, B matrices
                        allocate(H(num_func,num_func), V(num_func,num_func), B(num_func, num_func))
                        
                        B = 0.0d0
                        do i=1, num_func-1
                                B(i,i) = 1.0d0
                                l = l_list(i)
                                k = k_list(i)
                                if(l_list(i+1) /= l) cycle
                                B(i,i+1) = -0.5*SQRT(1.0-l*(l+1.0)/((i+l)*(i+l+1.0)))
                                B(i+1,i) = B(i,i+1)
                        end do
                        B(num_func,num_func) = 1.0d0       

                        print*, "B" 
                        do i=1, num_func
                                do j=1,num_func
                                        write(*, '(F10.4, 1X)', advance='no') B(i,j)
                               end do
                               print*
                               print*
                        end do
                        
                       !V matrix calculation
                       V = 0.0d0
                       do i=1, num_func
                                do j=1, num_func
                                        li = l_list(i)          
                                        lj = l_list(j)    

                                        do lamd=0, 2*lmax, 2
                                                        
                                                f = basis_fns(:,j) * min(r(:), Rn/2.0d0)**lamd/max(r(:),Rn/2.0d0)**(lamd+1) * basis_fns(:,i)
                                                integral = sum(f(:)*weights(:))                                         
                                                V(i,j) = -2.0d0*integral*yint(dble(li),dble(m),dble(lamd),0.0d0,dble(lj),dble(m)) +V(i,j)
                                        end do

                                end do
                        end do 


                        print*, "V"
                        do i=1, num_func
                                do j=1,num_func
                                        write(*, '(F10.4, 1X)', advance='no') V(i,j)
                               end do
                               print*
                        end do

                        
                        !initialise H matrix
                        H = -0.5*alpha**2 * B
                        
                        
                        do i=1, num_func
                               H(i,i) =  alpha**2 + H(i,i)
                        end do

                        H = H+V
                        
                         
                        print*, "H" 
                        do i=1, num_func
                                do j=1,num_func
                                        write(*, '(F10.4, 1X)', advance='no') H(i,j)
                               end do
                               print*
                        end do

                        !calculate energies for each symmetry
                        allocate(z(num_func,num_func),w(num_func))

                        !solve eigenvalue problem

                        call rsg(num_func,num_func,H,B,w,1,z,ierr)
                        nstate = nstate +1
                        energies(nstate) = w(1)
                        
                        
                        deallocate(k_list, l_list, basis_fns, H, B, V,z,w)
                        print*, "-----------------------------------------------------------------------------"

                end do
        end do
     
 
               
        

        sorted = .false.
        do while(.not.sorted)
                sorted = .true.
                do i=1, nstates-1
                        E1 = energies(i)
                        E2 = energies(i+1)
                        if(E2 < E1) then
                            sorted = .false.
                            temp = energies(i)
                            energies(i) = energies(i+1)
                            energies(i+1) = temp
                        end if
                end do
        end do 
                
        
        
       file_name = "Energies.txt"

      open(unit=11, file= file_name, position='APPEND')
      !write energy states to file
      write(11,'(F10.4)',advance='no') Rn
      do i=1, nstates
                
                write(11, '(1X, F10.7)',advance='no')  energies(i)+1/Rn  
                
      end do
      close(11)      

        

      deallocate(r,weights,f,energies)

      contains
              subroutine calc_phikl(alpha, l, N, r,basis_fns)
                        implicit none
                        integer, intent(in) :: N,l
                        real*8,intent(in) :: alpha, r(:)
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

