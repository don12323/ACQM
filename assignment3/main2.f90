program main2
        implicit none
        real*8, parameter :: pi = 3.14159265358979323846
        real*8 :: alpha, rmax, dr, mu, D, E,Ekmin,Ekmax,s,amp 
        integer :: N,l,i,j,k, startidx,rn, ierr, num_Ek, vi
        real*8, allocatable :: r(:), integrand(:), basis_fns(:,:), wf(:,:), psi2psu(:,:), Ek(:),z(:,:), B(:,:), H(:,:), V(:,:),Vout(:), w(:), weights(:), KED(:,:)
        character(len=20) :: file_name

        !read input file
        open(unit=1, file='input.txt', status='old', action='read')
        read(1,*) alpha, l, N, dr, rmax, mu, Ekmin, Ekmax
        close(1)

        !allocate memory for r and phikl
        rn = int(rmax/dr) +1
        allocate(r(rn))
        allocate(basis_fns(rn, N), weights(rn), Vout(rn))

        !create radial coordinates
        do i = 1, rn
                r(i) = dr*real(i-1)
        end do
        !read data and interpolate data
        file_name = 'PEC.1ssg'
        call read_and_interpolate(file_name, rn, r, Vout)

        print*, 'interpolated PEC.1ssg'
        
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
         
        H = -0.5*alpha**2 * B / mu

        do i=1, N
                H(i,i) = alpha**2 / mu + H(i,i)
        end do

        H = H+V

        !solve eigenvalue problem
        call rsg(N,N,H,B,w,1,z,ierr)

        !write energies to file for given N, l and alpha
        !write(file_name, '(A,I0,A,I0,A,F4.2,A)') 'EvibN', N,'l',l,'a',alpha,'.txt'
        !open(unit=3, file= file_name, status='replace')
        !do i=1, N
        !        write(3,'(I2, F10.6)') N, w(i)
        !end do
        !close(3)
        
        !calculating radial wavefunctions
        wf = 0.0d0
        do j=1, N
                do i=1, N
                        wf(:,j) = wf(:,j) + z(i,j)*basis_fns(:,i)
                end do 
        end do

        !scale and square bound wfs for poltting
        !write to file 
        write(file_name, '(A,I1,A)') 'vibfns_l',l,'.txt' 
        open(unit=4, file=file_name, status='replace')
        do i=1, rn
                write(4,*) r(i),Vout(i), 0.005*wf(i,:) + w(:)
        end do
        close(4)
        !write v=0 wfs to file
       ! write(mu_str, '(F12.5)') mu
        !write(file_name, '(A,F12.5,A)') 'mu_',mu,'.txt'
       ! file_name = 'mu_' // trim(adjustl(mu_str)) // '.txt'
       ! open(unit=1, file=file_name, status='replace')
       ! do i=1, rn
       !         write(1,*) r(i), wf(i,1)
       ! end do
       ! close(1)
        
        !---------------------------------------------------------------------------------------------- 
        ! read and interpolate data for PEC.2psu
        file_name = 'PEC.2psu'
        Vout=0.0
        call read_and_interpolate(file_name, rn, r, Vout)
        print*, 'interpolated data for PEC.2psu'
        D = -0.5d0
        
        !convert eV to Ha
        Ekmax = Ekmax / 27.21136
        Ekmin = Ekmin / 27.21136
        num_Ek = 200 
        s=1.0e-5
        

        allocate(psi2psu(rn, num_Ek), Ek(num_Ek))
        psi2psu = 0.0
        !calculate Ek grid

        do i=1, num_Ek
                Ek(i) = Ekmin + (Ekmax - Ekmin) *(i-1) / (num_Ek - 1)
        end do

        !loop over Energies
        do i=1, num_Ek
                E = Ek(i) + D
                !print*, E
                call numerov_forward(psi2psu(:,i), Vout, E, rn, s , 1, dr, mu)
                !calculate asymptotic amplitude of the wf
                startidx = rn - rn / 10 
                amp = psi2psu(startidx,i)
                do j=startidx, rn
                        if (abs(psi2psu(j,i)) > amp) amp = abs(psi2psu(j,i))                  
                end do

                print*, 'amp', amp

                k = sqrt(2.0*mu*Ek(i))
                psi2psu(:,i) = psi2psu(:,i)*sqrt(2*mu / (k*pi)) / amp
                 

        end do         
        
        !write continuum wavefunctions to file
        file_name = 'contwfs.txt'
        open(unit=4, file=file_name, status='replace')
        do i=1, rn
                write(4,*) r(i), 0.002*psi2psu(i,:) + D + Ek(:)
        end do
        close(4)
        
        !franck condon approximation
        allocate(KED(rn, 4))
        vi=1
        do j= 1,10,3
                do i=1, num_Ek
                        integrand = psi2psu(:,i) * wf(:,j) 
                        KED(i,vi) = (abs(sum(integrand(:)*weights(:))))**2 
                                                
                end do
                vi=vi+1
        end do
        file_name = 'KED.txt'
        open(unit=4, file=file_name, status='replace')
        do i=1, num_Ek
                write(4,*) 27.21136*Ek(i), KED(i,:)
        end do
        close(4)
        


        deallocate(r,basis_fns,H,B,V,wf,w,z,Vout,weights,integrand,Ek, psi2psu, KED)

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
                
                subroutine numerov_forward(psi, V, E, N, s, nodes,dx,mu)
                       integer, intent(in) :: N, nodes
                       real*8, intent(in) :: V(N), E, s, dx, mu
                       real*8, intent(out) :: psi(N)
                       real*8, allocatable :: g(:)
                       integer :: i
                       allocate(g(N))
                       psi(1) = 0.0
                       psi(2) = s * (-1.0)**nodes
                       do i=1, N
                       g(i) = 2.0 * mu * (V(i) - E)
                       end do
                       
                       do i = 3, N
                       psi(i) = ((2.0 + 5.0 * dx**2 / 6.0 * g(i-1)) * psi(i-1) - &
                               (1.0 - dx**2 / 12.0 * g(i-2)) * psi(i-2)) / &
                               (1.0 - dx**2 / 12.0 * g(i))
                       end do
                       deallocate(g)
                end subroutine numerov_forward
                
                subroutine read_and_interpolate(file_name, rn, r, Vout)
                        character(len=*), intent(in) :: file_name
                        integer, intent(in) :: rn                    
                        real*8, intent(in) :: r(rn)
                        real*8, intent(out) :: Vout(rn)
                        integer :: num_points, i, ierr
                        real*8, allocatable :: Rin(:), Vin(:)
                        character(len=100) :: line
                                             
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
                        
                        allocate(Vin(num_points), Rin(num_points))
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
                                    
                end subroutine read_and_interpolate


end program main2

