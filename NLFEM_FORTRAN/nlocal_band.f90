!==============================================================================
!                             PROGRAM NLFEM2D.f
!==============================================================================
!   - A 2D-FEM code for Nonlocal Elasticty
!   - 4-node and 8-node elements can be use for the approximations.
!   - Integration with 3x3 GQP
!   - FERNANDO RAMIREZ SUMMER 2017
!==============================================================================
!   Initial reading file to dimension variables
    implicit real*16(a-h,o-z)
    integer nnd,nel,nne,nelm,ndf
    integer nbce,nbcn
    open(unit=10,file='input.txt',status='unknown')
    read(10,*) nnd,nel,nne,nelm,ndf,nebc,nnbc,nla,wtloc,wtnloc,t,rl,IB,nem
    close(10)
    nem=nem+1;
    write(*,*) nnd,nel,nne,nelm,ndf,nebc,nnbc,nla,wtloc,wtnloc,t,rl,IB,nem
    call maincode(nnd,nel,nne,nelm,ndf,nebc,nnbc,nla,wtloc,wtnloc,t,rl,IB,nem)
    end
!==============================================================================
    subroutine maincode(nnd,nel,nne,nelm,ndf,nebc,nnbc,nla,wtloc,wtnloc,t,rl,IB,nem)
    implicit real*16(a-h,o-z)
    dimension c11(nelm),c22(nelm),c12(nelm),c66(nelm)
    dimension x(nnd),y(nnd),F(2*nnd),U(2*nnd),GK(2*nnd,IB)
    dimension nod(nel,9),nle(nel,nem)
    dimension ebc(nebc,3),bcn(nnbc,3)
    dimension gs9(3),wt9(3)
    dimension sf(nne,9),dsfxi(nne,9),dsfeta(nne,9)
    dimension gj(2,2),gji(2,2),elmk(2*nne,2*nne),elx(nne),ely(nne),wg(9)
    dimension dx(nne),dy(nne),stra(nel,5),stre(nel,5)
    dimension elxn(nne),elyn(nne),elmknl(2*nne,2*nne)
    open(unit=10,file='input.txt',status='unknown')
    read(10,*) nnd,nel,nne,nelm,ndf,nebc,nnbc,nla,wtloc,wtnloc,t,rl,IB,nem
    write(*,*) nnd,nel,nne,nelm,ndf,nebc,nnbc,nla,wtloc,wtnloc,t,rl,IB,nem
    do i=1,nnd
        read(10,*) x(i),y(i)
    enddo
    do i=1,nel
        read(10,*) (nod(i,j),j=1,nne+1)
    enddo
    do i=1,nel
        read(10,*) nle(i,1),(nle(i,j),j=2,(nle(i,1)+1))
        write(*,*) i,nle(i,1)
    enddo
    do i=1,nelm
        read(10,*) c11(i),c22(i),c12(i),c66(i)
    enddo
    do i=1,nebc
        read(10,*) (ebc(i,j),j=1,3)
    enddo
    do i=1,nnbc
        read(10,*) (bcn(i,j),j=1,3)
    enddo
    close(10)
!==========================================================================
    do i=1,2*nnd
        do j=1,IB
            GK(i,j)=0.D0
        enddo
        F(i)=0.D0
    enddo
!==========================================================================
!   Gauss Points and Weights 3x3
    ntpg=9
    gs9(1)=-0.774596669241483
    gs9(2)=0.D0
    gs9(3)=0.774596669241483
    wt9(1)=5.D0/9.D0
    wt9(2)=8.D0/9.D0
    wt9(3)=5.D0/9.D0
!==========================================================================
! Evaluation of interpolation functions and their derivatives at Gauss-Points
    call shape(nne,gs9,wt9,sf,dsfxi,dsfeta,wg)
!==========================================================================
!   Local contribution to stiffness matrix - BANDED
!==========================================================================
    write(*,*) 'Klocal begins'
    do i=1,nel
        write(*,*) i
        do j=1,nne
            ndi=nod(i,j)
            elx(j)=x(ndi)
            ely(j)=y(ndi)
        enddo
        ntm=nod(i,nne+1)
        c11e=c11(ntm)
        c22e=c22(ntm)
        c12e=c12(ntm)
        c66e=c66(ntm)
        call locelmat(nne,ntpg,sf,dsfxi,dsfeta,elx,ely,c11e,c22e,c12e,c66e,elmk,wtloc,wg,t)
        do j=1,2*nne
            if (j.le.nne) then
                ndi=nod(i,j)
                nr=(ndi-1)*ndf+1
            else
                ndi=nod(i,j-nne)
                nr=(ndi-1)*ndf+2
            endif
            do k=1,2*nne
                if (k.le.nne) then
                    ndi=nod(i,k)
                    nc=(ndi-1)*ndf+1
                else
                    ndi=nod(i,k-nne)
                    nc=(ndi-1)*ndf+2
                endif
                if (nc.ge.nr) then
                    ncc=nc-nr+1
                    if (ncc.le.IB) then
                        GK(nr,ncc)=GK(nr,ncc)+elmk(j,k)
                    endif
                endif
             enddo
        enddo
    enddo
!==========================================================================
!   NONLOCAL contribution to stiffness matrix - BANDED
!=========================================================================
    if (nla.eq.1) then
    write(*,*) 'Knonlocal begins'
    do i=1,nel
        write(*,*) i
        do k=1,nne
            ndi=nod(i,k)
            elx(k)=x(ndi)
            ely(k)=y(ndi)
        enddo
        ntm=nod(i,nne+1)
        c11e=c11(ntm)
        c22e=c22(ntm)
        c12e=c12(ntm)
        c66e=c66(ntm)
        nnle=nle(i,1)
        do j=1,nnle
            nli=nle(i,j+1)
            do k=1,nne
                ndi=nod(nli,k)
                elxn(k)=x(ndi)
                elyn(k)=y(ndi)
            enddo
            ntm=nod(nli,nne+1)
            c11en=c11(ntm)
            c22en=c22(ntm)
            c12en=c12(ntm)
            c66en=c66(ntm)
            call nlocelmat(nne,ntpg,sf,dsfxi,dsfeta,wg,elx,ely,c11e,c22e,c12e,c66e,wtloc, &
            elxn,elyn,c11en,c22en,c12en,c66en,wtnloc,elmknl,t,rl)
            do k=1,2*nne
                if (k.le.nne) then
                    ndi=nod(i,k)
                    nr=(ndi-1)*ndf+1
                else
                    ndi=nod(i,k-nne)
                    nr=(ndi-1)*ndf+2
                endif
                do m=1,2*nne
                    if (m.le.nne) then
                        ndi=nod(nli,m)
                        nc=(ndi-1)*ndf+1
                    else
                        ndi=nod(nli,m-nne)
                        nc=(ndi-1)*ndf+2
                    endif
                    if (nc.ge.nr) then
                        ncc=nc-nr+1
                        if (ncc.le.IB) then
                            GK(nr,ncc)=GK(nr,ncc)+elmknl(k,m)
                        endif
                    endif
                enddo
            enddo
        enddo
    enddo
    endif
!==========================================================================
! Considering GK banded
!==========================================================================
    do i=1,nnbc
        nid=(bcn(i,1)-1)*ndf+bcn(i,2)
        F(nid)=bcn(i,3)
    enddo
    do i=1,nebc
        nid=(ebc(i,1)-1)*ndf+ebc(i,2)
        do j=1,2*nnd
            if (j.ge.nid) then
                k=j-nid+1
                if (k.le.IB) then
                    F(j)=F(j)-GK(nid,k)*ebc(i,3)
                    GK(nid,k)=0.D0
                endif
            else
                k=nid-j+1
                if (k.le.IB) then
                    F(j)=F(j)-GK(j,k)*ebc(i,3)
                    GK(j,k)=0.D0
                endif
            endif
         enddo
    enddo
    do i=1,nebc
        nid=(ebc(i,1)-1)*ndf+ebc(i,2)
        GK(nid,1)=1.0D0
        F(nid)=ebc(i,3)
    enddo
!==========================================================================
    write(*,*) 'Solving System of Equations'
!   call solgauss(nnd,GK,F,U)
    numeq=2*nnd
    LU=1
    call GELT(GK,F,LU,numeq,IB,U)
    open(unit=10,file='disp.txt',status='unknown')
!==========================================================================
! Saves displacements into a file
!==========================================================================
    do i=1,nnd
        nid=(i-1)*2
        write(10,*) U(nid+1),U(nid+2)
    enddo
    close(10)
!==========================================================================
! Strain and Stress Calculation - Element Central Point
!==========================================================================
    do i=1,nel
        do j=1,nne
            ndi=nod(i,j)
            elx(j)=x(ndi)
            ely(j)=y(ndi)
            ndi=(ndi-1)*ndf
            dx(j)=U(ndi+1)
            dy(j)=U(ndi+2)
        enddo
        xp=(elx(2)+elx(1))/2.0
        yp=(ely(2)+ely(3))/2.0
        stre(i,4)=xp
        stre(i,5)=yp
        stra(i,4)=xp
        stra(i,5)=yp
        rlx=elx(2)-elx(1)
        rly=ely(3)-ely(2)
        ntm=nod(i,nne+1)
        c11e=c11(ntm)
        c22e=c22(ntm)
        c12e=c12(ntm)
        c66e=c66(ntm)
        ne=i
        call stress(nel,ne,nne,elx,ely,dx,dy,c11e,c22e,c12e,c66e,stra,stre,xp,yp,rlx,rly)
    enddo
    open(unit=10,file='stress.txt',status='unknown')
    do i=1,nel
        write(10,*) (stra(i,j),j=1,5)
    enddo
    do i=1,nel
        write(10,*) (stre(i,j),j=1,5)
    enddo
    close(10)
!==========================================================================
    end
!==========================================================================
! Evaluation of local part element stiffnes matrix
!==========================================================================
    subroutine locelmat(nne,ntpg,sf,dsfxi,dsfeta,elx,ely,c11e,c22e,c12e,c66e,elmk,wtloc,wg,t)
    implicit real*16(a-h,o-z)
    dimension sf(nne,9),dsfxi(nne,9),dsfeta(nne,9),elx(nne),ely(nne)
    dimension elmk(2*nne,2*nne),gj(2,2),gji(2,2),wg(9)
    c66e=c66e
    do i=1,2*nne
        do j=1,2*nne
            elmk(i,j)=0.D0
        enddo
    enddo
    do i=1,ntpg
        do m=1,2
            do n=1,2
                gj(m,n)=0.D0
            enddo
        enddo
        do j=1,nne
            gj(1,1)=gj(1,1)+elx(j)*dsfxi(j,i)
            gj(1,2)=gj(1,2)+ely(j)*dsfxi(j,i)
            gj(2,1)=gj(2,1)+elx(j)*dsfeta(j,i)
            gj(2,2)=gj(2,2)+ely(j)*dsfeta(j,i)
        enddo
        detJ=gj(1,1)*gj(2,2)-gj(1,2)*gj(2,1)
        gji(1,1)=gj(2,2)/detJ
        gji(2,2)=gj(1,1)/detJ
        gji(1,2)=-gj(1,2)/detJ
        gji(2,1)=-gj(2,1)/detJ
        do j=1,nne
            dsxj=gji(1,1)*dsfxi(j,i)+gji(1,2)*dsfeta(j,i)
            dsyj=gji(2,1)*dsfxi(j,i)+gji(2,2)*dsfeta(j,i)
            do k=1,nne
                dsxk=gji(1,1)*dsfxi(k,i)+gji(1,2)*dsfeta(k,i)
                dsyk=gji(2,1)*dsfxi(k,i)+gji(2,2)*dsfeta(k,i)
                elmk(j,k)=elmk(j,k)+wtloc*(c11e*dsxj*dsxk+c66e*dsyj*dsyk)*detJ*wg(i)*t
                elmk(nne+j,nne+k)=elmk(nne+j,nne+k)+wtloc*(c22e*dsyj*dsyk+c66e*dsxj*dsxk)*detJ*wg(i)*t
                elmk(j,k+nne)=elmk(j,k+nne)+wtloc*(c12e*dsxj*dsyk+c66e*dsyj*dsxk)*detJ*wg(i)*t
                elmk(j+nne,k)=elmk(j+nne,k)+wtloc*(c12e*dsyj*dsxk+c66e*dsxj*dsyk)*detJ*wg(i)*t
            enddo
        enddo
    enddo
    end
!==========================================================================
    subroutine shape(nne,gs9,wt9,sf,dsfxi,dsfeta,wg)
    implicit real*16(a-h,o-z)
    dimension gs9(3),wt9(3),wg(9)
    dimension sf(nne,9),dsfxi(nne,9),dsfeta(nne,9)
    if (nne.eq.8) then
        npg=0
        do i=1,3
            xi=gs9(i)
            do j=1,3
                eta=gs9(j)
                npg=npg+1.0
                sf(1,npg)=0.25*(1.0-xi)*(1.0-eta)*(-xi-eta-1.0)
                dsfxi(1,npg)=-0.25*(1.0-eta)*(-2.0*xi-eta)
                dsfeta(1,npg)=-0.25*(1.0-xi)*(-2.0*eta-xi)
                sf(2,npg)=0.25*(1.0+xi)*(1.0-eta)*(xi-eta-1.0)
                dsfxi(2,npg)=0.25*(1.0-eta)*(2.0*xi-eta)
                dsfeta(2,npg)=-0.25*(1.0+xi)*(-2.0*eta+xi)
                sf(3,npg)=0.25*(1.0+xi)*(1.0+eta)*(xi+eta-1.0)
                dsfxi(3,npg)=0.25*(1.0+eta)*(2.0*xi+eta)
                dsfeta(3,npg)=0.25*(1.0+xi)*(2.0*eta+xi)
                sf(4,npg)=0.25*(1.0-xi)*(1.0+eta)*(-xi+eta-1.0)
                dsfxi(4,npg)=-0.25*(1.0+eta)*(-2.0*xi+eta)
                dsfeta(4,npg)=0.25*(1.0-xi)*(2.0*eta-xi)
                sf(5,npg)=0.5*(1.0-xi*xi)*(1.0-eta)
                dsfxi(5,npg)=-1.0*(1.0-eta)*xi
                dsfeta(5,npg)=-0.5*(1.0-xi*xi)
                sf(6,npg)=0.5*(1.0-eta*eta)*(1.0+xi)
                dsfxi(6,npg)=0.5*(1.0-eta*eta)
                dsfeta(6,npg)=-1.0*(1.0+xi)*eta
                sf(7,npg)=0.5*(1.0-xi*xi)*(1.0+eta)
                dsfxi(7,npg)=-1.0*xi*(1.0+eta)
                dsfeta(7,npg)=0.5*(1.0-xi*xi)
                sf(8,npg)=0.5*(1.0-xi)*(1.0-eta*eta)
                dsfxi(8,npg)=-0.5*(1.0-eta*eta)
                dsfeta(8,npg)=-1.0*eta*(1.0-xi)
                wg(npg)=wt9(i)*wt9(j)
            enddo
        enddo
    elseif (nne.eq.4) then
        npg=0
        do i=1,3
            xi=gs9(i)
            do j=1,3
                eta=gs9(j)
                npg=npg+1.0
                sf(1,npg) =0.25*(xi-1.0)*(eta-1.0)
                dsfxi(1,npg)=0.25*(eta-1.0)
                dsfeta(1,npg)=0.25*(xi-1.0)
                sf(2,npg)=-0.25*(eta-1.0)*(xi+1.0)
                dsfxi(2,npg)=-0.25*(eta-1.0)
                dsfeta(2,npg)=-0.25*(xi+1.0)
                sf(3,npg)=0.25*(eta+1.0)*(xi+1.0)
                dsfxi(3,npg)=0.25*(eta+1.0)
                dsfeta(3,npg)=0.25*(xi+1.0)
                sf(4,npg)=-0.25*(eta+1.0)*(xi-1.0)
                dsfxi(4,npg)=-0.25*(eta+1.0)
                dsfeta(4,npg)=-0.25*(xi-1.0)
                wg(npg)=wt9(i)*wt9(j)
            enddo
        enddo
    endif
    end
!==========================================================================
    subroutine solgauss(nnd,GK,F,U)
    implicit real*16(a-h,o-z)
    dimension GK(2*nnd,2*nnd),F(2*nnd),U(2*nnd)
    do i=1,2*nnd-1
        do j=i+1,2*nnd
            fac=GK(j,i)/GK(i,i)
            do k=1,2*nnd
                GK(j,k)=GK(j,k)-fac*GK(i,k)
            enddo
            F(j)=F(j)-fac*F(i)
        enddo
    enddo
    do i=1,2*nnd
        ii=2*nnd-i+1
        do j=ii+1,2*nnd
            F(ii)=F(ii)-GK(ii,j)*F(j)
        enddo
        F(ii)=F(ii)/GK(ii,ii)
        U(ii)=F(ii)
    enddo
    end
!==========================================================================
    subroutine stress(nel,ne,nne,elx,ely,dx,dy,c11e,c22e,c12e,c66e,stra,stre,xp,yp,rlx,rly)
    implicit real*16(a-h,o-z)
    dimension ssf(nne),sfxi(nne),sfeta(nne),elx(nne),ely(nne)
    dimension stra(nel,5),stre(nel,5)
    dimension gj(2,2),gji(2,2),dx(nne),dy(nne)
    do i=1,2
        do j=1,2
            gj(i,j)=0.D0
        enddo
    enddo
    do i=1,3
        stra(ne,i)=0.D0
        stre(ne,i)=0.D0
    enddo
    xi=2.D0*(xp-elx(1))/rlx-1
    eta=2.D0*(yp-ely(1))/rly-1
    call shape2(nne,ssf,sfxi,sfeta,xi,eta)
    do j=1,nne
        gj(1,1)=gj(1,1)+elx(j)*sfxi(j)
        gj(1,2)=gj(1,2)+ely(j)*sfxi(j)
        gj(2,1)=gj(2,1)+elx(j)*sfeta(j)
        gj(2,2)=gj(2,2)+ely(j)*sfeta(j)
    enddo
    detJ=gj(1,1)*gj(2,2)-gj(1,2)*gj(2,1)
    gji(1,1)=gj(2,2)/detJ
    gji(2,2)=gj(1,1)/detJ
    gji(1,2)=-gj(1,2)/detJ
    gji(2,1)=-gj(2,1)/detJ
    do j=1,nne
        dsx=gji(1,1)*sfxi(j)+gji(1,2)*sfeta(j)
        dsy=gji(2,1)*sfxi(j)+gji(2,2)*sfeta(j)
        stra(ne,1)=stra(ne,1)+dx(j)*dsx
        stra(ne,2)=stra(ne,2)+dy(j)*dsy
        stra(ne,3)=stra(ne,3)+0.5D0*(dx(j)*dsy+dy(j)*dsx)
    enddo
    stre(ne,1)=c11e*stra(ne,1)+c12e*stra(ne,2)
    stre(ne,2)=c12e*stra(ne,1)+c22e*stra(ne,2)
    stre(ne,3)=2*c66e*stra(ne,3)
    end
!==========================================================================
!   Evaluate shape functions at specific (xi,eta) point
!==========================================================================
    subroutine shape2(nne,ssf,sfxi,sfeta,xi,eta)
    implicit real*16(a-h,o-z)
    dimension ssf(nne),sfxi(nne),sfeta(nne)
    if (nne.eq.8) then
        ssf(1)=0.25*(1.0-xi)*(1.0-eta)*(-xi-eta-1.0)
        sfxi(1)=-0.25*(1.0-eta)*(-2.0*xi-eta)
        sfeta(1)=-0.25*(1.0-xi)*(-2.0*eta-xi)
        ssf(2)=0.25*(1.0+xi)*(1.0-eta)*(xi-eta-1.0)
        sfxi(2)=0.25*(1.0-eta)*(2.0*xi-eta)
        sfeta(2)=-0.25*(1.0+xi)*(-2.0*eta+xi)
        ssf(3)=0.25*(1.0+xi)*(1.0+eta)*(xi+eta-1.0)
        sfxi(3)=0.25*(1.0+eta)*(2.0*xi+eta)
        sfeta(3)=0.25*(1.0+xi)*(2.0*eta+xi)
        ssf(4)=0.25*(1.0-xi)*(1.0+eta)*(-xi+eta-1.0)
        sfxi(4)=-0.25*(1.0+eta)*(-2.0*xi+eta)
        sfeta(4)=0.25*(1.0-xi)*(2.0*eta-xi)
        ssf(5)=0.5*(1.0-xi*xi)*(1.0-eta)
        sfxi(5)=-1.0*(1.0-eta)*xi
        sfeta(5)=-0.5*(1.0-xi*xi)
        ssf(6)=0.5*(1.0-eta*eta)*(1.0+xi)
        sfxi(6)=0.5*(1.0-eta*eta)
        sfeta(6)=-1.0*(1.0+xi)*eta
        ssf(7)=0.5*(1.0-xi*xi)*(1.0+eta)
        sfxi(7)=-1.0*xi*(1.0+eta)
        sfeta(7)=0.5*(1.0-xi*xi)
        ssf(8)=0.5*(1.0-xi)*(1.0-eta*eta)
        sfxi(8)=-0.5*(1.0-eta*eta)
        sfeta(8)=-1.0*eta*(1.0-xi)
    elseif (nne.eq.4) then
        ssf(1) =0.25*(xi-1.0)*(eta-1.0)
        sfxi(1)=0.25*(eta-1.0)
        sfeta(1)=0.25*(xi-1.0)
        ssf(2)=-0.25*(eta-1.0)*(xi+1.0)
        sfxi(2)=-0.25*(eta-1.0)
        sfeta(2)=-0.25*(xi+1.0)
        ssf(3)=0.25*(eta+1.0)*(xi+1.0)
        sfxi(3)=0.25*(eta+1.0)
        sfeta(3)=0.25*(xi+1.0)
        ssf(4)=-0.25*(eta+1.0)*(xi-1.0)
        sfxi(4)=-0.25*(eta+1.0)
        sfeta(4)=-0.25*(xi-1.0)
    endif
    end
!==========================================================================
    subroutine nlocelmat(nne,ntpg,sf,dsfxi,dsfeta,wg,elx,ely,c11e,c22e,c12e,c66e,wtloc, &
    elxn,elyn,c11en,c22en,c12en,c66en,wtnloc,elmknl,t,rl)
    implicit real*16(a-h,o-z)
    dimension sf(nne,9),dsfxi(nne,9),dsfeta(nne,9),elx(nne),ely(nne)
    dimension elxn(nne),elyn(nne),gj(2,2),gji(2,2),wg(9)
    dimension elmknl(2*nne,2*nne),gjn(2,2),gjin(2,2)
    do i=1,2*nne
        do j=1,2*nne
            elmknl(i,j)=0.D0
        enddo
    enddo
    do il=1,ntpg
        do i=1,2
            do j=1,2
                gj(i,j)=0.D0
            enddo
        enddo
        xl=0.D0
        yl=0.D0
        do j=1,nne
            gj(1,1)=gj(1,1)+elx(j)*dsfxi(j,il)
            gj(1,2)=gj(1,2)+ely(j)*dsfxi(j,il)
            gj(2,1)=gj(2,1)+elx(j)*dsfeta(j,il)
            gj(2,2)=gj(2,2)+ely(j)*dsfeta(j,il)
            xl=xl+elx(j)*sf(j,il)
            yl=yl+ely(j)*sf(j,il)
        enddo
        detJ=gj(1,1)*gj(2,2)-gj(1,2)*gj(2,1)
        gji(1,1)=gj(2,2)/detJ
        gji(2,2)=gj(1,1)/detJ
        gji(1,2)=-gj(1,2)/detJ
        gji(2,1)=-gj(2,1)/detJ
        do inl=1,ntpg
            do i=1,2
                do j=1,2
                    gjn(i,j)=0.D0
                enddo
            enddo
            xnl=0.D0
            ynl=0.D0
            do j=1,nne
                gjn(1,1)=gjn(1,1)+elxn(j)*dsfxi(j,inl)
                gjn(1,2)=gjn(1,2)+elyn(j)*dsfxi(j,inl)
                gjn(2,1)=gjn(2,1)+elxn(j)*dsfeta(j,inl)
                gjn(2,2)=gjn(2,2)+elyn(j)*dsfeta(j,inl)
                xnl=xnl+elxn(j)*sf(j,inl)
                ynl=ynl+elyn(j)*sf(j,inl)
            enddo
            detJn=gjn(1,1)*gjn(2,2)-gjn(1,2)*gjn(2,1)
            gjin(1,1)=gjn(2,2)/detJn
            gjin(2,2)=gjn(1,1)/detJn
            gjin(1,2)=-gjn(1,2)/detJn
            gjin(2,1)=-gjn(2,1)/detJn
            dnd=sqrt((xnl-xl)**2+(ynl-yl)**2)
            rrl=rl/6.D0
            atf=(1/(2*3.14159*rrl*rrl*t))*exp(-dnd/rrl)
            do jl=1,nne
                dsxjl=gji(1,1)*dsfxi(jl,il)+gji(1,2)*dsfeta(jl,il)
                dsyjl=gji(2,1)*dsfxi(jl,il)+gji(2,2)*dsfeta(jl,il)
                do jnl=1,nne
                    dsxjnl=gjin(1,1)*dsfxi(jnl,inl)+gjin(1,2)*dsfeta(jnl,inl)
                    dsyjnl=gjin(2,1)*dsfxi(jnl,inl)+gjin(2,2)*dsfeta(jnl,inl)
                    elmknl(jl,jnl)=elmknl(jl,jnl)+wtnloc*atf*(c11e*dsxjl*dsxjnl+c66e*dsyjl*dsyjnl) &
                    *detJ*detJn*wg(il)*wg(inl)*t*t
                    elmknl(nne+jl,nne+jnl)=elmknl(nne+jl,nne+jnl)+wtnloc*atf*(c22e*dsyjl*dsyjnl+c66e*dsxjl*dsxjnl) &
                    *detJ*detJn*wg(il)*wg(inl)*t*t
                    elmknl(jl,jnl+nne)=elmknl(jl,jnl+nne)+wtnloc*atf*(c12e*dsxjl*dsyjnl+c66e*dsyjl*dsxjnl) &
                    *detJ*detJn*wg(il)*wg(inl)*t*t
                    elmknl(jl+nne,jnl)=elmknl(jl+nne,jnl)+wtnloc*atf*(c12e*dsyjl*dsxjnl+c66e*dsxjl*dsyjnl) &
                    *detJ*detJn*wg(il)*wg(inl)*t*t
                enddo
            enddo
        enddo
    enddo
    end
!==========================================================================
!   gauss elimination by Dr. Erik G. Thompson - Banded
!==========================================================================
    subroutine GELT(A,Y,LU,NUMEQ,IB,X)
    implicit real*16(a-h,o-z)
    integer NEM1,NUMEQ,IB,JEND,I,J,K1,J1,I1
    dimension A(NUMEQ,IB),Y(NUMEQ),X(NUMEQ)
    NEM1=NUMEQ-1
!==========================================================================
!     BEGIN FORWARD ELIMINATION
!==========================================================================
    write(*,*) 'Forward'
    do I=1,NEM1
         JEND=NUMEQ-I+1
         if (JEND.gt.IB)  then
           JEND=IB;
         endif
         do J=2,JEND
            J1=I+J-1
            if (LU.eq.1)     then
               FAC=A(I,J)/A(I,1)
               K1=0
               do K=J,JEND
                 K1=K1+1;
                 A(J1,K1)=A(J1,K1)-A(I,K)*FAC
               enddo
               A(I,J)=FAC;
            endif
            Y(J1)=Y(J1)-Y(I)*A(I,J)
         enddo
      enddo
!==========================================================================
!     BEGIN BACK SUBSTITUTION
!==========================================================================
    write(*,*) 'Back 2'
    Y(NUMEQ)=Y(NUMEQ)/A(NUMEQ,1)
    do I=1,NEM1
         I1=NUMEQ-I
         Y(I1)=Y(I1)/A(I1,1)
         JEND=NUMEQ-I1+1
         if (JEND.gt.IB) then
            JEND=IB;
         endif
        do J=2,JEND
           J1=I1+J-1
            Y(I1)=Y(I1)-A(I1,J)*Y(J1)
         enddo
      enddo
    do I=1,NUMEQ
        if (ABS(Y(I)).lt.1.D-50) then
            X(I)=0.D0
        else
            X(I)=Y(I)
        endif
    enddo
    return
    end
