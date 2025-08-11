! Ising Model

! 1D case 
!   Finite Size Chains - 2^N particles for the chain (Large N's possible)
!   Varying Temperatures

! E = -J Sum(Si * Si+1) - h Sum (Si)
! PBC : S(n+1) = S1

program Ising1D
    implicit none
    integer :: N, i
    integer, dimension(:, :), allocatable :: ising_lattice
    
    real :: temp = 0, k = 10, J, h
    real :: Q, Z, A, F, E, U
    real :: temp_frac    

    print *, "Enter the number of particles ="
    read *, N

    print *, "Enter the interaction energy ="
    read *, J

    print *, "Enter the magnetic field ="
    read *, h

    !print *, "Enter the temperature ="
    !read *, temp

    allocate(ising_lattice(2**N,N+1))

    call array(N, ising_lattice)
    call hamiltonian(N, ising_lattice, J, h)

    ! Graphing A vs T
    OPEN(unit=30, file="AvsT.dat", status="unknown")

    do i = 0, 100
        temp = temp + 1
        temp_frac = (k*temp) / J
        Z = Q(N, temp, k)
        F = A(Z, temp, k)
        E = U(F, temp, N, k)
        print *, "Iteration ", i, temp, Z!/(float(N)*J)
        write (30, *) temp, Z!/(float(N)*J)
        
        rewind 19
    end do

    deallocate(ising_lattice)

    !rewind 100

end program Ising1D


subroutine array(N, ising_lattice)
    implicit none
    integer :: i, j, temp
    integer, intent(in) :: N
    integer, dimension(2**N, N+1), intent(out) :: ising_lattice
    integer, dimension(2) :: spin_array

    spin_array = (/ -1, 1 /)

    !OPEN(unit=100, file="isingArray.dat", status="unknown")

    do i = 0, (2**N)-1
        temp = i
        do j = 1, N
            if (MOD(temp,2) == 1) THEN
                ising_lattice(i+1,j) = spin_array(2)
            else
                ising_lattice(i+1,j) = spin_array(1)
            endif
            
            temp = temp / 2
        end do

        ising_lattice(i+1, N+1) = ising_lattice(i+1,1)
        
    end do

    !do i = 1, N
    !    print *, ising_lattice(6, i)
    !end do
    !write(100, *) ising_lattice


end 

subroutine hamiltonian(N, ising_lattice, E, B)
    implicit none
    integer, intent(in) :: N
    integer, dimension(2**N, N+1):: ising_lattice

    integer :: k, l, nnsum = 0, sum = 0 
    real :: E, H, B
    
    ! Simple Hamiltonian
    open(unit=19, file="Hdata.dat", status="unknown")

    do k = 1, 2**N
         
        if (N > 2) THEN 
            do l = 1, N
                nnsum = nnsum + (ising_lattice(k, l)*ising_lattice(k,l+1))
                sum = sum + ising_lattice(k, l)
            end do
        
        else 
            nnsum = (ising_lattice(k, 1)*ising_lattice(k,2))
            sum = ising_lattice(k, 1) + ising_lattice(k, 2)
        endif

        H = (-1*E*nnsum) - (B*sum)
        write (19, *) k, H
        
        nnsum = 0
        sum = 0
    end do

    rewind 19

end

real function Q(N, T, k)
    implicit none
    real :: h, T, k
    integer :: i, iter, N

    Q = 0.0

    ! Partition function
    do i = 1, 2**N
        read (19, *) iter, h
        Q = Q + EXP(-h/(k*T))

        print *, "Parition = ", Q 
    end do

end

real function A(part, T, k)
    implicit none
    real :: part, T, k

    A = -k*T*LOG(part)
end

real function U(F, T, N, k)
    implicit none
    real :: F, T, S, k
    integer :: N

    S = k * LOG(float(N))

    ! Internal Energy
    U = F + (T*S)

end



