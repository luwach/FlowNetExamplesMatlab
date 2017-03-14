%Parameters of model and sheet

Lx = input('Podaj szerokosc modelu Lx [m]: ');
Ly = input('Podaj wysokosc modelu Ly [m]: ');
Ls = input('Podaj dlugosc scianki szczelnej Ls [m]: ');
d = input('Podaj oczko siatki d [m]: ');
dH = input('Podaj roznice wysokosci zwierciadla wody dH [m]: ');

%Number of nodes

nx = Lx/d+1;
ny = Ly/d+1;
ns = Ls/d+1;

%New martix

AB = linspace(dH,dH,nx);
AE = linspace(dH, dH/2+1, ny);
ED = linspace(dH/2+1, dH/2, nx);
BC = linspace(dH, dH/2, ns);
CD = linspace(dH/2, dH/2, ny-ns);
BD = [BC CD];

h = ones(ny,nx);

h(1,:)=AB;
h(:,nx)=BD;
h(:,1)=AE;
h(ny,:)=ED;

%Result of solution

tol=1d-6;
err=1;
m=0;

hkp1=h;

while err > tol
    
    m = m+1;
    
    for i = 2:ny-1
        for j = 2:nx-1
            
            hkp1(i,j) = 0.25*(h(i+1,j)+h(i,j+1)+h(i-1,j)+h(i,j-1));
            
        end
    end
    
    for i = 2:ny-1
        
        hkp1(i,1) = 0.25*(h(i+1,1)+2*h(i,2)+h(i-1,1));
        
    end
    
    hkp1(ny,1) = 0.5*(h(ny-1,1)+h(ny,2));
    
    for j = 2:nx-1
        
        hkp1(ny,j) = 0.25*(h(ny,j-1)+h(ny,j+1)+2*h(ny-1,j));
        
    end
    
    for i = 2:ns-1
        
        hkp1(i,nx) = 0.25*(h(i-1,nx)+h(i+1,nx)+2*h(i,nx-1));
        
    end
    
    err = sqrt(sum(sum((hkp1-h).^2)));
    
    h = hkp1;
    
end

h;

q1=0;

for j=2:nx-1
    
    q1 = q1 + (h(2,j)-h(4,j));
    
end

q_per_unit = (1/(2*d))*((h(2,1)-h(4,1))+d*q1+(h(2,13)-h(4,13)));

S_AE = linspace(q_per_unit,q_per_unit,ny);
S_BC = linspace(0,0,ns-1);
S_CD = linspace(0,q_per_unit,ny-(ns-1));
S_BD = [S_BC S_CD];
S_ED = linspace(q_per_unit, q_per_unit,nx);

s = zeros(ny,nx);

s(:,1)=S_AE;
s(:,nx)=S_BD;
s(ny,:)=S_ED;

erro=1;
k=0;

skp1=s;

while erro > tol
    
    k = k+1;
    
    for i = 2:ny-1
        for j = 2:nx-1
            
            skp1(i,j) = 0.25*(s(i+1,j)+s(i,j+1)+s(i-1,j)+s(i,j-1));
            
        end
    end
    
    for i = ns+1:ny-1
        
        skp1(i,nx) = 0.25*(s(i-1,nx)+s(i+1,nx)+2*s(i,nx-1));
        
    end
    
    for j = 2:nx-1
        
        skp1(1,j) = 0.25*(s(1,j-1)+s(1,j+1)+2*s(2,j));
        
    end
    
    erro = sqrt(sum(sum((skp1-s).^2)));
    
    s = skp1;
    
end

s;

[X,Y] = meshgrid(0:d:Lx,0:-d:-Ly);
ch = contour(X,Y,h,10, 'color', 'b');
hold on
cs = contour(X,Y,s,10, 'color', 'r');
line([-2 -8], [-4 -10], 'color', 'k');
hold off
