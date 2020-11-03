!Clase Serendipity
module MSerendipity
    implicit none

    type Serendipity
        real :: X(8)
        real :: Y(8)
        real :: t
        real :: l
        real :: l0
        real :: PUNTOS(3)
        real :: PESOS(3)
        integer, allocatable :: nolocales(:)
    contains
        procedure, pass :: psi
        procedure, pass :: dzpsi
        procedure, pass :: dnpsi
        procedure, pass :: TX
        procedure, pass :: TY
        procedure, pass :: JAC
        procedure, pass :: ML
        procedure, pass :: MNL

    end type Serendipity

contains
    real function dzpsi(this,i,z,n)
        class(Serendipity), intent(inout) :: this
        integer :: i
        real :: z
        real :: n
        dzpsi = -1.0
        if (i==1) then
            dzpsi = -1.0/4.0*(n-1.0)*(2.0*z+n)
        else if (i==2) then
            dzpsi = -1.0/4.0*(n-1.0)*(2.0*z-n)
        else if (i==3) then
            dzpsi = 1.0/4.0*(n+1.0)*(2.0*z+n)
        else if (i==4) then
            dzpsi = 1.0/4.0*(n+1.0)*(2.0*z-n)
        else if (i==5) then
            dzpsi = (n-1.0)*z
        else if (i==6) then
            dzpsi = -1.0/2.0*(n*n-1.0)
        else if (i==7) then
            dzpsi = -(n+1.0)*z
        else
            dzpsi = 1.0/2.0*(n*n-1.0)
        endif
    end function dzpsi

    real function dnpsi(this,i,z,n)
        class(Serendipity), intent(inout) :: this
        integer :: i
        real :: z
        real :: n
        dnpsi = -1.0
        if (i==1) then
            dnpsi = -1.0/4.0*(z-1.0)*(2.0*n+z)
        else if (i==2) then
            dnpsi = 1.0/4.0*(z+1.0)*(2.0*n-z)
        else if (i==3) then
            dnpsi = 1.0/4.0*(z+1.0)*(2.0*n+z)
        else if (i==4) then
            dnpsi = -1.0/4.0*(z-1.0)*(2.0*n-z)
        else if (i==5) then
            dnpsi = 1.0/2.0*(z*z-1.0)
        else if (i==6) then
            dnpsi = -n*(z+1.0)
        else if (i==7) then
            dnpsi = -1.0/2.0*(z*z-1.0)
        else
            dnpsi = n*(z-1.0)
        endif
    end function dnpsi

    real function psi(this,i,z,n)
        class(Serendipity), intent(inout) :: this
        integer :: i
        real :: z
        real :: n
        psi = -1.0
        if (i==1) then
            psi = 1.0/4.0*(1.0-z)*(1.0-n)*(-1.0-z-n)
        else if (i==2) then
            psi = 1.0/4.0*(1.0+z)*(1.0-n)*(-1.0+z-n)
        else if (i==3) then
            psi = 1.0/4.0*(1.0+z)*(1.0+n)*(-1.0+z+n)
        else if (i==4) then
            psi = 1.0/4.0*(1.0-z)*(1.0+n)*(-1.0-z+n)
        else if (i==5) then
            psi = 1.0/2.0*(1.0-z*z)*(1.0-n)
        else if (i==6) then
            psi = 1.0/2.0*(1.0+z)*(1.0-n*n)
        else if (i==7) then
            psi = 1.0/2.0*(1.0-z*z)*(1.0+n)
        else
            psi = 1.0/2.0*(1.0-z)*(1.0-n*n)
        endif
    end function psi

    real function TX(this,z,n)
        class(Serendipity), intent(inout) :: this
        integer :: i
        real :: z
        real :: n
        TX = 0.0
        do i=1,8
            TX = TX + this%X(i)*this%psi(i,z,n)
        enddo
    end function TX

    real function TY(this,z,n)
        class(Serendipity), intent(inout) :: this
        integer :: i
        real :: z
        real :: n
        TY = 0.0
        do i=1,8
            TY = TY + this%Y(i)*this%psi(i,z,n)
        enddo
    end function TY


    subroutine JAC(this,matriz,z,n)
        implicit none
        real :: matriz(2,2)
        real :: z
        real :: n
        integer :: i
        integer :: j
        class(Serendipity), intent(inout) :: this
        do i=1,2
            do j=1,2
                matriz(i,j)=0
            enddo
        enddo
        do i=1,8
            matriz(1,1)=matriz(1,1)+this%X(i)*this%dzpsi(i,z,n)
            matriz(1,2)=matriz(1,2)+this%Y(i)*this%dzpsi(i,z,n)
            matriz(2,1)=matriz(2,1)+this%X(i)*this%dnpsi(i,z,n)
            matriz(2,2)=matriz(2,2)+this%Y(i)*this%dnpsi(i,z,n)
        enddo
    end
    subroutine ML(e,matriz,C11,C12,C66)
        implicit none
        real :: matriz(16,16)
        integer :: gdli
        integer :: gdlj
        integer :: i
        integer :: j
        integer :: inl
        integer :: jnl

        real :: z
        real :: n
        real :: x
        real :: y
        real :: detjac
        real :: jacobiano(2,2)
        real :: jacobiano_(2,2)
        real :: dz_i
        real :: dn_i
        real :: dfdx_i
        real :: dfdy_i
        real :: dz_j
        real :: dn_j
        real :: dfdx_j
        real :: dfdy_j
        real :: t
        real :: C11
        real :: C12
        real :: C66
        real :: an
        real :: PESOS(3)

        class(Serendipity), intent(inout) :: e

        PESOS = e%PESOS
        t = e%t
        do i=1,16
            do j=1,16
                matriz(i,j)=0
            enddo
        enddo

        do gdli=1,8
            do gdlj=1,8
                do i=1,3
                    z = e%PUNTOS(i)
                    do j=1,3
                        n = e%PUNTOS(j)
                        
                        x = e%TX(z,n)
                        y = e%TY(z,n)

                        call e%JAC(jacobiano,z,n)
                        detjac = jacobiano(1,1)*jacobiano(2,2)-jacobiano(1,2)*jacobiano(2,1)

                        jacobiano_(1,1)=jacobiano(2,2)/detjac
                        jacobiano_(1,2)=-jacobiano(1,2)/detjac
                        jacobiano_(2,1)=-jacobiano(2,1)/detjac
                        jacobiano_(2,2)=jacobiano(1,1)/detjac

                        dz_i = e%dzpsi(gdli,z,n)
                        dn_i = e%dnpsi(gdli,z,n)

                        dfdx_i = dz_i * jacobiano_(1,1) + dn_i* jacobiano_(1,2)
                        dfdy_i = dz_i * jacobiano_(2,1) + dn_i* jacobiano_(2,2)

                        dz_j = e%dzpsi(gdlj,z,n)
                        dn_j = e%dnpsi(gdlj,z,n)

                        dfdx_J = dz_j * jacobiano_(1,1) + dn_j* jacobiano_(1,2)
                        dfdy_J = dz_j * jacobiano_(2,1) + dn_j* jacobiano_(1,1)
                        
                        an = detjac*PESOS(j)*PESOS(i)

                        matriz(gdli,gdlj)=matriz(gdli,gdlj)+t*(C11*dfdx_i*dfdx_j+C66*dfdy_i*dfdy_j)*an
                        matriz(gdli,gdlj+8)=matriz(gdli,gdlj+8)+t*(C12*dfdx_i*dfdy_j+C66*dfdy_i*dfdx_j)*an
                        matriz(gdli+8,gdlj)=matriz(gdli+8,gdlj)+t*(C12*dfdy_i*dfdx_j+C66*dfdx_i*dfdy_j)*an
                        matriz(gdli+8,gdlj+8)=matriz(gdli+8,gdlj+8)+t*(C11*dfdy_i*dfdy_j+C66*dfdx_i*dfdx_j)*an
                    enddo
                enddo
            enddo
        enddo
    end
    subroutine MNL(e,enl,matriz,C11,C12,C66)
        implicit none
        real :: matriz(16,16)
        integer :: gdli
        integer :: gdlj
        integer :: i
        integer :: j
        integer :: inl
        integer :: jnl

        real :: z
        real :: n
        real :: x
        real :: y
        real :: detjac
        real :: jacobiano(2,2)
        real :: jacobiano_(2,2)
        real :: dz_i
        real :: dn_i
        real :: dfdx_i
        real :: dfdy_i

        real :: znl
        real :: nnl
        real :: xnl
        real :: ynl
        real :: detjacnl
        real :: jacobianonl(2,2)
        real :: jacobiano_nl(2,2)
        real :: dznl_j
        real :: dnnl_j
        real :: dfdxnl_j
        real :: dfdynl_j
        real :: AZN
        real :: t
        real :: distancia
        real :: C11
        real :: C12
        real :: C66
        real :: an
        real :: PESOS(3)


        class(Serendipity), intent(inout) :: e
        class(Serendipity), intent(inout) :: enl

        PESOS = e%PESOS
        t = e%t
        do i=1,16
            do j=1,16
                matriz(i,j)=0
            enddo
        enddo

        do gdli=1,8
            do gdlj=1,8
                do i=1,3
                    z = e%PUNTOS(i)
                    do j=1,3
                        n = e%PUNTOS(j)
                        
                        x = e%TX(z,n)
                        y = e%TY(z,n)

                        call e%JAC(jacobiano,z,n)
                        detjac = jacobiano(1,1)*jacobiano(2,2)-jacobiano(1,2)*jacobiano(2,1)

                        jacobiano_(1,1)=jacobiano(2,2)/detjac
                        jacobiano_(1,2)=-jacobiano(1,2)/detjac
                        jacobiano_(2,1)=-jacobiano(2,1)/detjac
                        jacobiano_(2,2)=jacobiano(1,1)/detjac

                        dz_i = e%dzpsi(gdli,z,n)
                        dn_i = e%dnpsi(gdli,z,n)

                        dfdx_i = dz_i * jacobiano_(1,1) + dn_i* jacobiano_(1,2)
                        dfdy_i = dz_i * jacobiano_(2,1) + dn_i* jacobiano_(2,2)

                        do inl=1,3
                            znl = enl%PUNTOS(inl)
                            do jnl=1,3
                                nnl = enl%PUNTOS(jnl)

                                xnl = enl%TX(zNL,nNL)
                                ynl = enl%TY(zNL,nNL)

                                call enl%JAC(jacobianonl,znl,nnl)
                                detjacnl = jacobianonl(1,1)*jacobianonl(2,2)-jacobianonl(1,2)*jacobianonl(2,1)

                                jacobiano_nl(1,1)=jacobianonl(2,2)/detjacnl
                                jacobiano_nl(1,2)=-jacobianonl(1,2)/detjacnl
                                jacobiano_nl(2,1)=-jacobianonl(2,1)/detjacnl
                                jacobiano_nl(2,2)=jacobianonl(1,1)/detjacnl

                                dznl_j = enl%dzpsi(gdlj,znl,nnl)
                                dnnl_j = enl%dnpsi(gdlj,znl,nnl)

                                dfdxnl_J = dznl_j * jacobiano_nl(1,1) + dnnl_j* jacobiano_nl(1,2)
                                dfdynl_J = dznl_j * jacobiano_nl(2,1) + dnnl_j* jacobiano_nl(1,1)
                                distancia = sqrt((x-xnl)**2.0+(y-ynl)**2.0)
                                AZN = e%L0*exp(-distancia/e%l)
                                ! print *, AZN,distancia,x,xnl,y,ynl,"=",z,n,znl,nnl
                                ! print *, "========================="
                                an = detjac*detjacnl*PESOS(jnl)*PESOS(inl)*PESOS(j)*PESOS(i)
                                matriz(gdli,gdlj)=matriz(gdli,gdlj)+t*t*AZN*(C11*dfdx_i*dfdxnl_j+C66*dfdy_i*dfdynl_j)*an
                                matriz(gdli,gdlj+8)=matriz(gdli,gdlj+8)+t*t*AZN*(C12*dfdx_i*dfdynl_j+C66*dfdy_i*dfdxnl_j)*an
                                matriz(gdli+8,gdlj)=matriz(gdli+8,gdlj)+t*t*AZN*(C12*dfdy_i*dfdxnl_j+C66*dfdx_i*dfdynl_j)*an
                                matriz(gdli+8,gdlj+8)=matriz(gdli+8,gdlj+8)+t*t*AZN*(C11*dfdy_i*dfdynl_j+C66*dfdx_i*dfdxnl_j)*an
                            enddo
                        enddo
                    enddo
                enddo

            enddo
        enddo

    end

end module MSerendipity

program NLFEM
    use MSerendipity
    implicit none
!   Variables globales
    real :: PI
    real,allocatable :: GDLS(:,:)
    real,allocatable :: ELEMENTOS_NO_LOCALES(:)
    real :: E, V, L, T, L0, C11, C12, C66
    real :: PUNTOS(3)
    real :: PESOS(3)
    integer :: NGAUSS
    integer :: NGDL
    integer :: NUMERO_GDL
    integer :: NUMERO_ELEMENTOS
    integer :: NUMERO_NODOS_ELEMENTO
    real :: EQUIS(8)
    real :: YE(8)
    integer :: gdl_e(8)
    integer :: longitudLinea
    real :: start, finish

    integer :: i
    integer :: j
    integer :: k
    integer :: m
    integer :: kkkkk
    character(len=255) :: X1
    character(len=255) :: X2

    character(len=255) format
    real :: MATRIZ_NKLNM(16,16)
    real :: MATRIZ_L(16,16)

    type(Serendipity), allocatable :: ELEMENTOS(:)

    E = 2.1*10.0**6
    V = 0.2
    L = 0.1
    T = 0.5
    C11 = E/(1.0-V**2)
    C12 = V*E/(1.0-V**2)
    C66 = E/2/(1.0+V)
    PI = 3.141592653589793238463
    L0 = 1.0/2.0/PI/L**2/T


    PUNTOS(1) = -sqrt(3.0/5.0)
    PUNTOS(2) = 0
    PUNTOS(3) = sqrt(3.0/5.0)

    NGAUSS = size(PUNTOS)
    NGDL = 16

    PESOS(1) = 5.0/9.0
    PESOS(2) = 8.0/9.0
    PESOS(3) = 5.0/9.0

    open(unit=10, file='input.txt')
    read(10,*) NUMERO_GDL, NUMERO_ELEMENTOS

    allocate(ELEMENTOS(NUMERO_ELEMENTOS))
    allocate(GDLS(NUMERO_GDL,2))

    do i=1,NUMERO_GDL
        read(10,*) GDLS(i,1),GDLS(i,2)
    enddo

    do i=1,NUMERO_ELEMENTOS
        read(10,*) gdl_e
        do j=1,8
            YE(j)=GDLS(gdl_e(j),2)
            EQUIS(j)=GDLS(gdl_e(j),1)
        enddo
        ELEMENTOS(i)%X = EQUIS
        ELEMENTOS(i)%Y = YE
        ELEMENTOS(i)%t = T
        ELEMENTOS(i)%l = L
        ELEMENTOS(i)%l0 = L0

        ELEMENTOS(i)%PUNTOS = PUNTOS
        ELEMENTOS(i)%PESOS = PESOS

    enddo
    do i=1,NUMERO_ELEMENTOS
        read(10, '(A)') format  
        read(format,*) longitudLinea
        allocate(ELEMENTOS(i)%nolocales(longitudLinea+1))
        read(format,*) ELEMENTOS(i)%nolocales
    enddo
    close(10)
    call execute_command_line ('mkdir MATRICES')
    kkkkk = chdir('MATRICES')
    do i=1,NUMERO_ELEMENTOS
        call cpu_time(start)
        write(X1,*) i
        call execute_command_line ('mkdir Elemento'// trim(adjustl(X1)))
        call ELEMENTOS(i)%ML(MATRIZ_L,C11,C12,C66)
        open(10,file='Elemento'//trim(adjustl(X1))//'/KL_'//trim(adjustl(X1))//'.csv')
        write (10, "(*(G0,:,','))") MATRIZ_L
        close(10)
        open(11,file='Elemento'//trim(adjustl(X1))//'/KNLS'//'.csv')
        do j=2,ELEMENTOS(i)%nolocales(1)+1
            call ELEMENTOS(i)%MNL(ELEMENTOS(ELEMENTOS(i)%nolocales(j)),MATRIZ_NKLNM,C11,C12,C66)
            write (11, "(*(G0,:,','))") MATRIZ_NKLNM
        enddo
        close(11)
        print *, "------------------------------------------------------------"
        call cpu_time(finish)
        print *, "Elemento ",i, " - Tiempo: ",(finish-start)*1000, " ms"
    enddo
end program NLFEM
