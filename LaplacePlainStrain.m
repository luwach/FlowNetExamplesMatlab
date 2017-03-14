clc;
clear all;


%Parametry modelu i scianki

Lx=21;
Ly=31;

%Liczba wezlow

nx=Lx+1;
ny=Ly+1;
    
h=ones(ny,nx);

for i = 1:ny
    for j = 1:nx
        
        h(i,j)= -Inf;
        
    end    
end

%Prawe wezly

for i = 1:10
    for j = 12:nx
        
        h(i,j)=NaN;
        
    end    
end

%WARUNKI BRZEGOWE

%Cisnienie gora

h(1,1:11)=linspace(3,3,11);

%Cisnienie dol

h(11,12:nx)=linspace(0,0,11);

%Lewa strona scianki

h(1:15,11)=linspace(3,3.5,15);

%Prawa strona scianki

h(11:15,12)=linspace(0,3.5,5);

%Lewy brzeg

h(1:ny,1)=linspace(3,4.5,ny);

%Dolny brzeg

h(ny,1:nx)=linspace(4.5,4.5,nx);

%Prawy brzeg

h(11:ny,nx)=linspace(0,4.5,22);

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

%Budowanie macierzy M

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

%Iteracyjne wyznaczenie wartosci wezlow w h

%Zamiana NaN w rogu na wartosci

for i = 1:10
    for j = 12:nx
        
        h(i,j)=0;
        
    end    
end

tol=1d-6;
err=1;


while err > tol
    
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
   
    %Warunki brzegowe

    %Lewa granica modelu
    
    for i = 2:ny-1
        
        hkp1(i,1) = 0.25*(h(i+1,1)+2*h(i,2)+h(i-1,1));
        
    end
    
    %Lewa dolny rog modelu
    
    hkp1(ny,1) = 0.5*(h(ny-1,1)+h(ny,2));
    
    %Dolna granica modelu
  
    for j = 2:nx-1
        
        hkp1(ny,j) = 0.25*(h(ny,j-1)+h(ny,j+1)+2*h(ny-1,j));
        
    end
    
    %Prawy dolny rog modelu

    hkp1(ny,nx) = 0.5*(h(ny,nx-1)+h(ny-1,nx));
    
    %Prawa granica modelu

    for i = 12:ny-1
        
        hkp1(i,nx) = 0.25*(h(i+1,nx)+2*h(i,nx-1)+h(i-1,nx));
        
    end
    
    %Lewa strona scianki
    
    for i = 2:15
        
        hkp1(i,11) = 0.25*(h(i+1,11)+h(i-1,11)+2*h(i,10));
        
    end
    
    %Prawa strona scianki
    
    for i = 12:15
        
        hkp1(i,12) = 0.25*(h(i+1,12)+h(i-1,12)+2*h(i,13));
        
    end

    err = sqrt(sum(sum((hkp1-h).^2)));
    
    h = hkp1;
    
end

for i = 1:10
    for j = 12:nx
        
        h(i,j)=NaN;
        
    end    
end

h

[X,Y] = meshgrid(0:1:Lx,0:-1:-Ly);
[C,h]=contourf(X,Y,h,'LevelStep',0.2,'color', 'r');
clabel(C,h,'LabelSpacing',1000);
hold on
x=[10 11 11 10];
y=[0 0 -15 -15];
fill(x,y,'k')
