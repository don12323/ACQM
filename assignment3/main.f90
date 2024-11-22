program main
        implicit none

        real*8, parameter :: xmin = -5.0d0, xmax = 5.0d0, s=1.0e-6, eps=1.0e-6 
        real*8, allocatable :: x(:), V(:), weights(:), psiR(:),psiL(:),psi(:)
        real*8 :: dx, E, Emax, Emin, deltaE, integral
        integer :: i, nodes,N ,iter, max_iter, node_count
        logical converged

        N =100
        max_iter=10000

        if (mod(N,2) == 0) N = N + 1

        print*, 'N=', N

        dx = (xmax-xmin)/(N-1)
        allocate(x(N), V(N), weights(N), psiL(N), psiR(N), psi(N))
        print*, 'dx =' , dx 
        !calculate weights
        weights(1) = 1.0d0
        do i=2, N-1
                weights(i) = 2.0d0 + 2.0d0*mod(i+1,2)
        end do
        weights(N) = 1.0d0
        weights(:) = weights(:)*dx/3.0d0
        print*, '------------------------------------------'
         
        !initialise grid and potential
        do i = 1, N
            x(i) = xmin + (i - 1) * dx
            V(i) = 0.5 * x(i)**2
        end do
        
        !loop over eignestates
        do nodes=0, 3
              !estimate energies
              Emin = 0.0d0
              Emax = 10.0d0
              converged = .false.
              iter = 0

              do while (.not. converged .and. iter < max_iter) 
                  !calculate the new energy
                  E = (Emin + Emax) / 2.0d0
                  print*, 'E', E
                  !shoot from left boundary to right boundary to generate wf
                  call numerov_forward(psiL, V, E, N, s, nodes,dx)
                  call numerov_backward(psiR, V, E, N, s, dx)
                  call match_wavefunctions(psiL, psiR, psi, N)  
                  !count nodes
                  node_count = count_nodes(psi, N)
                  print*, 'nodes = ', node_count 
                  !adjust energies to get correct nodes
                  if (node_count < nodes) then
                     Emin = E
                  else if (node_count > nodes) then
                     Emax = E
                  else 
                       print*, '**calculating cooleys energies' 
                       do while (.not. converged)
                          call numerov_forward(psiL, V, E, N, s, nodes,dx)
                          call numerov_backward(psiR, V, E, N, s, dx)
                          call match_wavefunctions(psiL, psiR, psi, N)
                          call cooley_correction(psiL, psiR, V, E, psi, deltaE, N, dx)              
                          !check for convergence
                          if (abs(deltaE) > eps) then
                                  E = E + deltaE
                                  print*, 'Energy correction:', deltaE, 'New E:', E
                          else if (abs(deltaE) <= eps) then
                                  converged = .true.
                                  print*, '**converged**'
                                  print*, 'Energy correction:', deltaE, 'New E:', E
                          end if
                          
                          iter = iter+1
                          if (iter >= max_iter) exit
                       end do
                  end if     
                  iter=iter+1
              end do
              !normalise wf
              write(*,'(A, I0)') '#Iterations till convergence: ', iter
              integral = sum((psi(:)**2) * weights(:))
              psi(:) = psi(:) / sqrt(integral)
              integral=0.0d0
              do i=1, N-1
                        integral = integral + 0.5*(abs(psi(i)**2) + abs(psi(i+1)**2))*dx
              end do
              print*, 'Area under the curve:', integral
              !write wf to file
              call write_wave_function(nodes, x, psi, E)                        
              print*, '----------------------------------' 
        end do
         

        !write to file
        open(unit=1, file='output.txt', status='replace', action='write')
        do i = 1, N
                write(1, '(F10.4, 1X, F10.4)') x(i), V(i)
        end do


        deallocate(x,V, weights,psiL, psiR, psi)
        contains

                subroutine numerov_forward(psi, V, E, N, s, nodes,dx)
                    integer, intent(in) :: N, nodes
                    real*8, intent(in) :: V(N), E, s, dx
                    real*8, intent(out) :: psi(N)
                    real*8, allocatable :: g(:)
                    integer :: i
                    allocate(g(N))
                    psi(1) = 0.0
                    psi(2) = s * (-1.0)**nodes
                    
                    do i=1, N
                      g(i) = 2.0 * (V(i) - E)
                    end do

                    do i = 3, N
                      psi(i) = ((2.0 + 5.0 * dx**2 / 6.0 * g(i-1)) * psi(i-1) - &
                                (1.0 - dx**2 / 12.0 * g(i-2)) * psi(i-2)) / &
                               (1.0 - dx**2 / 12.0 * g(i))
                    end do
                    deallocate(g)
                end subroutine numerov_forward

                subroutine numerov_backward(psi,V, E, N, s, dx)
                    integer, intent(in) :: N
                    real*8, intent(in) :: V(N), E, s, dx
                    real*8, intent(out) :: psi(N)
                    real*8, allocatable :: g(:)
                    integer :: i
                    allocate(g(N))

                    psi(N) = 0.0
                    psi(N-1) = s
                    
                    do i=1, N
                        g(i) = 2.0 * (V(i) - E)
                    end do

                    do i = N-2, 1, -1
                      psi(i) = ((2.0 + 5.0 * dx**2 / 6.0 * g(i+1)) * psi(i+1) - &
                                (1.0 - dx**2 / 12.0 * g(i+2)) * psi(i+2)) / &
                               (1.0 - dx**2 / 12.0 * g(i))
                    end do
                    deallocate(g)
                end subroutine numerov_backward
                  
                subroutine match_wavefunctions(psiL, psiR, psi, N)
                    integer, intent(in) :: N
                    real*8, intent(in) :: psiL(N), psiR(N)
                    real*8, intent(out) :: psi(N)
                    integer :: i, m

                    m = N / 2
                    !ensure psiL(m) and psiR(m) are non-zero and change m
                    do while((psiL(m) == 0.0 .or. psiR(m) == 0.0) .and. m<N)
                        m=m+1
                    end do
                    do i = 1, m
                      psi(i) = psiL(i) / psiL(m)
                    end do
                    do i = m+1, N
                      psi(i) = psiR(i) / psiR(m)
                    end do
                end subroutine match_wavefunctions 


                subroutine cooley_correction(psiL, psiR, V, E, psi, deltaE, N, dx)
                    integer, intent(in) :: N
                    real*8, intent(in) :: psiL(N), psiR(N), V(N), E, dx
                    real*8, intent(out) :: deltaE
                    real*8, intent(inout) :: psi(N)
                    real*8 :: numerator, denominator, Ym, Ym_minus1, Ym_plus1
                    integer :: i, m

                    !Find the matching point where both psiL(m) and psiR(m) are non-zero
                    m = N / 2
                    do while ((psiL(m) == 0.0 .or. psiR(m) == 0.0) .and. m < N)
                      m = m + 1
                    end do

                    !Calculate psi at the matching point

                    !Calculate Yi terms
                    Ym = (1.0 - dx**2 / 12.0 * 2.0 * (V(m) - E)) * psi(m)
                    Ym_minus1 = (1.0 - dx**2 / 12.0 * 2.0 * (V(m-1) - E)) * psi(m-1)
                    Ym_plus1 = (1.0 - dx**2 / 12.0 * 2.0 * (V(m+1) - E)) * psi(m+1)

                    !Calculate numerator and denominator for deltaE
                    numerator = psi(m) * (-0.5 * (Ym_plus1 - 2.0 * Ym + Ym_minus1) / dx**2 + (V(m) - E) * psi(m))
                    denominator = 0.0

                    !Compute the normalization factor
                    do i = 1, N
                        denominator = denominator + psi(i)**2 
                    end do

                    deltaE = numerator / denominator
                end subroutine cooley_correction

                 subroutine write_wave_function(n, x, psi, E)
                    implicit none
                    integer, intent(in) :: n
                    real*8, intent(in) :: x(:), psi(:), E
                    integer :: i
                    character(len=20) :: filename
                    write(filename, '(A,I0,A)') 'wf_', n,'.txt'

                    open(unit=10, file=filename, status='replace')
                    do i = 1, size(x)
                      write(10,*) x(i), 0.75*psi(i)+E
                    end do
                    close(10)
                  end subroutine write_wave_function





                integer function count_nodes(psi, N)
                    real*8, intent(in) :: psi(N)
                    integer, intent(in) :: N
                    integer :: i

                    count_nodes = 0
                    do i = 2, N
                      if (psi(i-1) * psi(i) < 0.0) then
                        count_nodes = count_nodes + 1
                      end if
                    end do
                end function count_nodes

end program main
       
