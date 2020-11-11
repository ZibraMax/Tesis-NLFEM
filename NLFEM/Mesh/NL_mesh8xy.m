%%% Meshing of a rectangular domain using 8-node elements identifying
%%% neighbor elements within a distance lr
clear all
lr=6*0.1;
%lr=0.75;
%%%%%%
lx=5;
ly=5;
t=0.1;
nex=50;
ney=50;
nne=8;
hx=lx/nex;
hy=ly/ney;
wtloc=0.5;
wtnloc=0.5;
nla=1;
nnd=(ney+1)*(2*nex+1)+(ney)*(nex+1);
coor=zeros(nnd,2);
x=zeros(nnd,1);
y=zeros(nnd,1);
nel=nex*ney;
nle=zeros(nel,400);
elm=zeros(nel,9);
ntel=1;
input=[nnd,nel,nne,ntel];
%Coordinate Generation
nd=0;
for i=1:ney
    cy=(i-1)*hy;
    for j=1:2*nex+1
        nd=nd+1;
        y(nd)=cy;
        x(nd)=(j-1)*hx/2;
    end
    cy=(i-1)*hy+hy/2;
    for j=1:nex+1
        nd=nd+1;
        y(nd)=cy;
        x(nd)=(j-1)*hx;
    end
end
cy=ly;
for j=1:2*nex+1
    nd=nd+1;
    y(nd)=cy;
    x(nd)=(j-1)*hx/2;
end
%plot(x,y,'*')
%% Element Node Connectivty    
ne=0;
for i=1:ney
    ne=ne+1;
    elm(ne,1)=(i-1)*(3*nex+2)+1;
    elm(ne,2)=elm(ne,1)+2;
    elm(ne,4)=elm(ne,1)+3*nex+2;
    elm(ne,3)=elm(ne,4)+2;
    elm(ne,5)=elm(ne,1)+1;    
    elm(ne,8)=elm(ne,1)+2*nex+1;
    elm(ne,7)=elm(ne,4)+1;
    elm(ne,6)=elm(ne,8)+1;
    elm(ne,9)=1;
    for j=1:(nex-1)
        ne=ne+1;
        elm(ne,1)=elm(ne-1,2);
        elm(ne,2)=elm(ne,1)+2;
        elm(ne,4)=elm(ne-1,3);
        elm(ne,3)=elm(ne,4)+2;
        elm(ne,5)=elm(ne,1)+1;
        elm(ne,8)=elm(ne-1,6);
        elm(ne,7)=elm(ne,4)+1;
        elm(ne,6)=elm(ne,8)+1;
        elm(ne,9)=1;
    end
end
%% Identify neighbors within influence distance
nem=0;
for i=1:nel
    ne=0;
    nii=elm(i,1);
    nfi=elm(i,3);
    xci=(x(nii)+x(nfi))/2;
    yci=(y(nii)+y(nfi))/2;  
    ne=ne+1;
    nle(i,ne+1)=i;    
    for j=1:nel
        if j~=i
            nij=elm(j,1);
            nfj=elm(j,3);
            xcj=(x(nij)+x(nfj))/2;
            ycj=(y(nij)+y(nfj))/2;
            dist=sqrt((xci-xcj)^2+(yci-ycj)^2);
            if dist<(lr+0.5*sqrt(hx^2+hy^2))
                ne=ne+1;
                nle(i,ne+1)=j;
            end
        end
    end
    nle(i,1)=ne;
    if ne>nem
        nem=ne;
    end
end
%%% Coordinates for saving a file
coor(:,1)=x(:);
coor(:,2)=y(:);
%%% Material Properties
E=21000000.0;
nu=0.20;
nelm=1;
for i=1:nelm
    c11(i)=E/(1-nu^2);
    c22(i)=c11(i);
    c12(i)=nu*E/(1-nu^2);
    c66(i)=E/2/(1+nu);
end
ndf=2;
nebc=0;
nnbc=0;
for i=1:nnd
    if x(i)==0
        nebc=nebc+1;
        ebc(nebc,1)=i;
        ebc(nebc,2)=1;
        ebc(nebc,3)=0.00;
%         nebc=nebc+1;
%         ebc(nebc,1)=i;
%         ebc(nebc,2)=2;
%         ebc(nebc,3)=0.00;
    elseif x(i)==lx
        nebc=nebc+1;
        ebc(nebc,1)=i;
        ebc(nebc,2)=1;
        ebc(nebc,3)=0.5;
%         if y(i)<=ly/2
%            ebc(nebc,3)=0.005/(ly/2)*y(i);
%         else
%             ebc(nebc,3)=0.005-0.005/(ly/2)*(y(i)-ly/2);
%         end
    end
end
%% Bandwidth Calculation
IB=0;
ndf=2;
for i=1:nel
    for j=1:8
        ni=elm(i,j);
        for k=1:8
            nf=elm(i,k);
            nIB=abs(nf-ni+1)*ndf;
            if nIB>IB
                IB=nIB
            end
        end
    end
end
            
        
fid=fopen('input.txt','wt');
fprintf(fid, '%i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %12.8f\t %12.8f\t %12.8f\t %12.8f\t %i\t %i\n', nnd,nel,nne,nelm,ndf,nebc,nnbc,nla,wtloc,wtnloc,t,lr,IB,nem);
for i=1:nnd
    fprintf(fid, '%12.8f\t %12.8f\n',x(i),y(i));
end
for i=1:nel
    fprintf(fid, '%i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\n',elm(i,1),elm(i,2),elm(i,3),elm(i,4),elm(i,5),elm(i,6),elm(i,7),elm(i,8));
end
for i=1:nel
    fprintf(fid, '%i\t',nle(i,1));
    for j=1:nle(i,1)-1
        fprintf(fid, '%i\t',nle(i,j+1));
    end
    fprintf(fid, '%i\n',nle(i,nle(i,1)+1));
end
for i=1:nelm
    fprintf(fid, '%12.8f\t %12.8f\t %12.8f\t %12.8f\n',c11(i), c22(i), c12(i), c66(i));
end
for i=1:nebc
    fprintf(fid, '%i\t %i\t %12.8f\n',ebc(i,1),ebc(i,2),ebc(i,3));
end
for i=1:nnbc
    fprintf(fid, '%i\t %i\t %12.8f\n',nbc(i,1),nbc(i,2),nbc(i,3));
end
fclose(fid);
for ii=1:nel
plot(x,y,'*')
hold on    
i=ii;
    xpl(1)=x(elm(i,1));
    xpl(2)=x(elm(i,2));
    xpl(3)=x(elm(i,3));
    xpl(4)=x(elm(i,4));
    ypl(1)=y(elm(i,1));
    ypl(2)=y(elm(i,2));
    ypl(3)=y(elm(i,3));
    ypl(4)=y(elm(i,4));
    fill(xpl,ypl,'b')
    hold on
    for j=1:nle(i,1)-1
        k=nle(i,j+2)
        xpl(1)=x(elm(k,1));
        xpl(2)=x(elm(k,2));
        xpl(3)=x(elm(k,3));
        xpl(4)=x(elm(k,4));
        ypl(1)=y(elm(k,1));
        ypl(2)=y(elm(k,2));
        ypl(3)=y(elm(k,3));
        ypl(4)=y(elm(k,4));
        fill(xpl,ypl,'r')
        hold on
    end
    axis equal
pause
hold off
close ALL
end

        