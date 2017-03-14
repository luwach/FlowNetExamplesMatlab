clc;
clear all;


%Parameters of wall

Lx=49;
Ly=17;

%Number of nodes

nx=Lx+1;
ny=Ly+1;
    
h=ones(ny,nx);

for i = 1:ny
    for j = 1:nx
        
        h(i,j)= -Inf;
        
    end    
end

%Right nodes

for i = 14:ny
    for j = 1:nx
        
        h(i,j)=NaN;
        
    end    
end

%BOUNDARY CONDITIONS

%Pressure top

h(13,1:25)=linspace(3,3,25);

%Pressure bottom

h(13,26:nx)=linspace(0,0,25);

%Left side of the wall

h(8:13,25)=linspace(1.5,3,6);

%Right side of the wall

h(8:13,26)=linspace(1.5,0,6);

%Left boundary of the model

h(1:13,1)=linspace(4.5,3,ny-5);

%Bottom boundary of the model

h(1,1:nx)=linspace(4.5,4.5,nx);

%Right boundary of the model

h(1:13,nx)=linspace(4.5,0,ny-5);

u_to_w = zeros(ny, nx);
w_to_u = zeros(2, nx*ny);
m = 0;

for i = 1:ny
    for j = 1:nx
        if h(i,j) == -Inf
            m=m+1;
            u_to_w(i,j) = m;
            w_to_u(:, m) = [i, j]';
        end
    end
end

%Matrix building M

M0=zeros(m,m);
M1=zeros(m,m);

for m1=1:m
    
    c = w_to_u(:,m1);
    
    if h(c(1)-1,c(2)) < 10000 && h(c(1)-1,c(2)) > -10000
    
        M1(m1,m1)=-1;
        
    elseif h(c(1)-1,c(2)) == -Inf
        
        M1(m1,m1)=-1;
        M0(m1, u_to_w(c(1)-1, c(2)))=1;
        
    end
    
    if h(c(1)+1,c(2)) < 10000 && h(c(1)+1,c(2)) > -10000
    
        M1(m1,m1)=M1(m1,m1)-1;
      
    elseif h(c(1)+1,c(2)) == -Inf
        
        M1(m1,m1)=M1(m1,m1)-1;
        M0(m1, u_to_w(c(1)+1, c(2)))=1;
        
    end
    
    if h(c(1),c(2)-1) < 10000 && h(c(1),c(2)-1) > -10000
    
        M1(m1,m1)=M1(m1,m1)-1;
       
    elseif h(c(1),c(2)-1) == -Inf
        
        M1(m1,m1)=M1(m1,m1)-1;
        M0(m1, u_to_w(c(1), c(2)-1))=1;
        
    end
    
    if h(c(1),c(2)+1) < 10000 && h(c(1),c(2)+1) > -10000
    
        M1(m1,m1)=M1(m1,m1)-1;
        
    elseif h(c(1),c(2)+1) == -Inf
        
        M1(m1,m1)=M1(m1,m1)-1;
        M0(m1, u_to_w(c(1), c(2)+1))=1;
        
    end
    
end

M=M0+M1;

%Calculation of nodes - h

%Replace 0 in corner with NaN

for i = 14:ny
    for j = 1:nx
        
        h(i,j)=0;
        
    end    
end

tol=1d-6;
err=1;
t = waitbar(0,'h is counting...');

while err > tol
    
    waitbar(tol/err)
    
    for m1=1:m
    
        c = w_to_u(:,m1);
   
        h(c(1), c(2))= -Inf;
    
    end

    b=zeros(m,1);
   
    for m1=1:m
    
        c = w_to_u(:,m1);
    
        if h(c(1)-1,c(2)) < 10000 && h(c(1)-1,c(2)) > -10000
    
            b(m1,1)=-h(c(1)-1,c(2));
        
        end
    
        if h(c(1)+1,c(2)) < 10000 && h(c(1)+1,c(2)) > -10000
    
            b(m1,1)=b(m1,1)-h(c(1)+1,c(2));
        
        end
    
        if h(c(1),c(2)-1) < 10000 && h(c(1),c(2)-1) > -10000
    
            b(m1,1)=b(m1,1)-h(c(1),c(2)-1);
        
        end
    
        if h(c(1),c(2)+1) < 10000 && h(c(1),c(2)+1) > -10000

            b(m1,1)=b(m1,1)-h(c(1),c(2)+1);
        
        end
    
    end
   
    u=M\b;

    for m1=1:m
    
        c = w_to_u(:,m1);
   
        h(c(1), c(2))=u(m1,1);
    
    end

    hkp1=h;
   
    %Boundary conditions

    %Left boundary of the model
    
    for i = 2:12
        
        hkp1(i,1) = 0.25*(h(i-1,1)+2*h(i,2)+h(i+1,1));
        
    end
    
    %Left corner of the model
    
    hkp1(1,1) = 0.5*(h(2,1)+h(1,2));
    
    %Bottom boundary of the model
  
    for j = 2:nx-1
        
        hkp1(1,j) = 0.25*(h(1,j-1)+h(1,j+1)+2*h(2,j));
        
    end
    
    %Right corner of the model

    hkp1(1,nx) = 0.5*(h(1,nx-1)+h(2,nx));
    
    %Right boundary of the model

    for i = 2:12
        
        hkp1(i,nx) = 0.25*(h(i+1,nx)+2*h(i,nx-1)+h(i-1,nx));
        
    end
    
    %Left side of the wall
    
    for i = 8:12
        
        hkp1(i,25) = 0.25*(h(i+1,25)+h(i-1,25)+2*h(i,24));
        
    end
    
    %Right side of the wall
    
    for i = 8:12
        
        hkp1(i,26) = 0.25*(h(i+1,26)+h(i-1,26)+2*h(i,27));
        
    end

    err = sqrt(sum(sum((hkp1-h).^2)));
    
    h = hkp1;
    
end

close(t)

for i = 14:ny
    for j = 1:nx
        
        h(i,j)=NaN;
        
    end    
end

h;

q1=0;

for j=2:24
    
    q1 = q1 + (h(12,j)-h(10,j));
    
end

q_per_unit = (1/(2*1))*(0.5*(h(12,1)-h(10,1))+q1+0.5*(h(12,25)-h(10,25)));

s=ones(ny,nx);

for i = 1:ny
    for j = 1:nx
        
        s(i,j)= -Inf;
        
    end    
end

%Right nodes

for i = 14:ny
    for j = 1:nx
        
        s(i,j)=NaN;
        
    end    
end

%Boundary conditions

%Pressure top

s(13,1:25)=linspace(q_per_unit,0,25);

%Bottom pressure

s(13,26:nx)=linspace(0,q_per_unit,25);

%Left side of the wall

s(7:13,25)=linspace(0,0,7);

%Right side of the wall

s(7:13,26)=linspace(0,0,7);

%Left boundary of the model

s(1:13,1)=linspace(q_per_unit,q_per_unit,ny-5);

%Bottom boundary of the model

s(1,1:nx)=linspace(q_per_unit,q_per_unit,nx);

%Right boundary of the model

s(1:13,nx)=linspace(q_per_unit,q_per_unit,ny-5);

v_to_z = zeros(ny, nx);
z_to_v = zeros(2, nx*ny);
n = 0;

for i = 1:ny
    for j = 1:nx
        if s(i,j) == -Inf
            n=n+1;
            v_to_z(i,j) = n;
            z_to_v(:, n) = [i, j]';
        end
    end
end

%Calculating matrix N

N0=zeros(n,n);
N1=zeros(n,n);

for n1=1:n
    
    d = z_to_v(:,n1);
    
    if s(d(1)-1,d(2)) < 10000 && s(d(1)-1,d(2)) > -10000
    
        N1(n1,n1)=-1;
        
    elseif s(d(1)-1,d(2)) == -Inf
        
        N1(n1,n1)=-1;
        N0(n1, v_to_z(d(1)-1, d(2)))=1;
        
    end
    
    if s(d(1)+1,d(2)) < 10000 && s(d(1)+1,d(2)) > -10000
    
        N1(n1,n1)=N1(n1,n1)-1;
      
    elseif s(d(1)+1,d(2)) == -Inf
        
        N1(n1,n1)=N1(n1,n1)-1;
        N0(n1, v_to_z(d(1)+1, d(2)))=1;
        
    end
    
    if s(d(1),d(2)-1) < 10000 && s(d(1),d(2)-1) > -10000
    
        N1(n1,n1)=N1(n1,n1)-1;
       
    elseif s(d(1),d(2)-1) == -Inf
        
        N1(n1,n1)=N1(n1,n1)-1;
        N0(n1, v_to_z(d(1), d(2)-1))=1;
        
    end
    
    if s(d(1),d(2)+1) < 10000 && s(d(1),d(2)+1) > -10000
    
        N1(n1,n1)=N1(n1,n1)-1;
        
    elseif s(d(1),d(2)+1) == -Inf
        
        N1(n1,n1)=N1(n1,n1)-1;
        N0(n1, v_to_z(d(1), d(2)+1))=1;
        
    end
    
end

N=N0+N1;

%Calculation of nodes - s

%Replace 0 in corner with NaN

for i = 14:ny
    for j = 1:nx
        
        s(i,j)=0;
        
    end    
end

tol=1d-6;
erro=1;
l = waitbar(0,'s is calculating...');

while erro > tol
    
    waitbar(tol/erro)
    
    for n1=1:n
    
        d = z_to_v(:,n1);
   
        s(d(1), d(2))= -Inf;
    
    end

    f=zeros(n,1);
   
    for n1=1:n
    
        d = z_to_v(:,n1);
    
        if s(d(1)-1,d(2)) < 10000 && s(d(1)-1,d(2)) > -10000
    
            f(n1,1)=-s(d(1)-1,d(2));
        
        end
    
        if s(d(1)+1,d(2)) < 10000 && s(d(1)+1,d(2)) > -10000
    
            f(n1,1)=f(n1,1)-s(d(1)+1,d(2));
        
        end
    
        if s(d(1),d(2)-1) < 10000 && s(d(1),d(2)-1) > -10000
    
            f(n1,1)=f(n1,1)-s(d(1),d(2)-1);
        
        end
    
        if s(d(1),d(2)+1) < 10000 && s(d(1),d(2)+1) > -10000

            f(n1,1)=f(n1,1)-s(d(1),d(2)+1);
        
        end
    
    end
   
    z=N\f;

    for n1=1:n
    
        d = z_to_v(:,n1);
   
        s(d(1), d(2))=z(n1,1);
    
    end

    skp1=s;
   
    %Boundary conditions

    %Upper left boundary of the model
    
    for j = 2:24
        
        skp1(13,j) = 0.25*(s(13,j+1)+2*s(12,j)+s(13,j-1));
        
    end
    
    %Upper right boundary of the model
    
    for j = 27:nx-1
        
        skp1(13,j) = 0.25*(s(13,j+1)+2*s(12,j)+s(13,j-1));
        
    end
    
    erro = sqrt(sum(sum((skp1-s).^2)));
    
    s = skp1;
    
end

close(l)

for i = 14:ny
    for j = 1:nx
        
        s(i,j)=NaN;
        
    end  
end

s;

[X,Y] = meshgrid(0:1:Lx,0:1:Ly);
[C,h]=contour(X,Y,h,'LevelStep',0.1,'color', 'b');
clabel(C,h,'LabelSpacing',2000);
hold on
[C,s]=contour(X,Y,s,'LevelStep',0.1,'color', 'r');
clabel(C,s,'LabelSpacing',2000);

x=[24 25 25 24];
y=[Ly Ly 6 6];
fill(x,y,'k')
hold off
