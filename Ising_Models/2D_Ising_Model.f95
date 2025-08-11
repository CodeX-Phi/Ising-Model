! Ising Model

! 2D case 
!   Finite Size - 2^(N^2) particles for length size N
!   Varying Temps

program Ising2D
    implicit none
    integer :: N, i
    integer, dimension(:, :, :), allocatable :: ising_lattice

    real :: temp = 1, Tc, k = 10, J, h
    real :: Q, Z, A, F, E
    real :: temp_frac   

    print *, "Enter the number of particles ="
    read *, N

    print *, "Enter the interaction energy ="
    read *, J

    Tc = (4*J)/k
    print *, "Critical Temperature given by = ", Tc

    allocate(ising_lattice(2**(N*N), N+1, N+1))

    call array(N, ising_lattice)
    call hamiltonian(N, ising_lattice, J)

    ! call distro(N, J)
    
    ! call U(N, k, Tc, J)

    ! Graphing A vs T
    ! OPEN(unit=31, file="AvsT2.dat", status="unknown")

    OPEN(unit=41, file="QvsT2.dat", status="unknown")

    do i = 0, 100
        temp = temp + 1
        temp_frac = (k*temp) / J
        Z = Q(N, temp, k)
        F = A(Z, temp, k)
        print *, "Iteration ", i, temp, Z ! (float(N)*J)
        write (41, *) temp, Z ! (float(N)*J)
        
        rewind 9
    end do

    deallocate(ising_lattice)

end program Ising2D

subroutine array(N, ising_lattice) 
    implicit none
    integer :: i, j, k, temp
    integer, intent(in) :: N  ! side of square lattice
    integer, dimension(2**(N*N), N+1, N+1), intent(out) :: ising_lattice
    integer, dimension(2) :: spin_array

    spin_array = (/ -1, 1 /)

    !OPEN(unit=100, file="isingArray.dat", status="unknown")

    do i = 0, (2**(N*N))-1
        temp = i
        do j = 1, N

            do k = 1, N
                if (MOD(temp,2) == 1) THEN
                    ising_lattice(i+1, j, k) = spin_array(2)
                else
                    ising_lattice(i+1, j, k) = spin_array(1)
                endif
                
                temp = temp / 2
            end do
            
            ising_lattice(i+1, j, N + 1) = ising_lattice(i+1, j, 1)

            ! print *, ising_lattice(i+1, j, 1), ising_lattice(i+1, j, 2), ising_lattice(i+1, j, 3)
        end do

        do k = 1, N
            ising_lattice(i+1, N + 1, k) = ising_lattice(i+1, 1, k)
        end do

        !print *, ising_lattice(i+1, N+1, 1), ising_lattice(i+1, N+1, 2), ising_lattice(i+1, N+1, 3)
    end do

    print *, i

end

subroutine hamiltonian(N, ising_lattice, E)
    implicit none
    integer, intent(in) :: N
    integer, dimension(2**(N*N), N+1, N+1):: ising_lattice

    integer :: k, l, m, nnsum = 0, sum = 0 
    real :: E, H, B
    
    ! Simple Hamiltonian
    open(unit=9, file="H2data.dat", status="unknown")

    do k = 1, 2**(N*N)
         
        if (N > 2) THEN

            do l = 1, N
                
                do m = 1, N
                    nnsum = nnsum + (ising_lattice(k, l , m) * ising_lattice(k, l, m + 1))
                    nnsum = nnsum + (ising_lattice(k, l , m) * ising_lattice(k, l + 1, m))
                    ! sum = sum + ising_lattice(k, l)
                end do

            end do
        
        else 
            nnsum = nnsum + ((ising_lattice(k, 1, 1) * (ising_lattice(k, 1, 2)+ising_lattice(k, 2, 1))))
            nnsum = nnsum + ((ising_lattice(k, 2, 2) * (ising_lattice(k, 2, 3)+ising_lattice(k, 3, 2))))
            ! sum = ising_lattice(k, 1) + ising_lattice(k, 2)
        endif

        H = (-1*E*nnsum)
        write (9, *) H
        
        nnsum = 0
        sum = 0
    end do

    rewind 9

end

subroutine distro(N, J)
    implicit none
    integer :: i, N

    real :: h, J
    real :: tab(N*N*4 + 1,2)

    OPEN(unit=29, file="distro.dat", status="unknown")

    do i = 1, N*N*4 + 1
        tab (i, 1) = -(N*N*2) + i - 1
        tab(i, 2) = 0
    end do

    print *, tab

    do i = 1, 2**(N*N)
        read(9, *) h
        tab(int(h/J) + (N*N*2) + 1, 2) = tab(int(h/J) + (N*N*2) + 1, 2) + 1
    end do

    print *, tab

    do i = 1, N*N*4 + 1
        if (tab(i, 2) /= 0) THEN
            write(29, *) tab(i, 1), tab(i, 2)
        endif
    end do

end

real function Q(N, T, k)
    implicit none
    real :: h, T, k
    integer :: i, N

    Q = 0.0

    ! Partition function
    do i = 1, 2**(N*N)
        read (9, *) h
        Q = Q + EXP(-h/(k*T))
        
    end do

    rewind 9

end

real function A(Z, T, k)
    implicit none
    real :: Z, T, k

    A = -k*T*LOG(Z)
end

subroutine U(N, k, Tc, J)
    implicit none
    integer :: i
    real :: temp
    
    real :: Tc, k, J
    real :: Z, T, h
    real :: Q, E
    integer :: N

    T = Tc/4.0
    h = Tc/100.0

    OPEN(unit=11, file="2EvsT.dat", status="unknown")

    temp = LOG(Q(N, T + (0.5*h), k))

    Z = (temp - LOG(Q(N, T - (0.5*h), k)))/h
    E = k * (T**2) * Z

    write(11, *) T*k/(4*J), E/float(N)
    print *, "Iteration ", T*k/(4*J), E/float(N)

    do i = 2, 1000
        T = T + h

        Z = (LOG(Q(N, T + (0.5*h), k)) - temp)/ h
        E = k * (T**2) * Z

        write(11, *) ((T + h)*k)/(4*J), E/float(N)
        print *, "Iteration ", ((T + h)*k)/(4*J), E/float(N)

        temp = Z*h + temp

        rewind 9

    end do
    
end