% qsho_fdrelax.m
% Solve the quantum simple harmonic oscillator for the lowest
% energy wave function using inverse power iteration relaxation
% of the linear system obtained by finite differencing

% Clear memory and show only a few digits
clear all; format short; close all;

% Define values of the independent variable
xmin = -5;
xmax = 5;
ymin = -5;
ymax = 5;
zmin = -5;
zmax = 5;
h=0.25;
x=xmin:h:xmax;
y=ymin:h:ymax;
z=zmin:h:zmax;
nx = length(x);
ny = length(y);
nz = length(z);
well_potential = 1e10;
num_states = 10;
L = nx*ny*nz;

% Indexing function
m = @(p,q,r) (p-1)*nx^2  + (q-1)*ny + r;
% Reverse Indexing
rev = @(k) [floor((k-1)/nx^2)+1 floor(mod((k-1),nx^2)/ny)+1 mod((k-1),ny)+1];


D = zeros(L,L);

e=ones(nx,1);
L=spdiags([e -2*e e], -1:1, nx, nx)/h^2;
I = speye(nx);
L3 = kron(kron(L,I),I) + kron(kron(I,L), I) + kron(kron(I,I),L);

% Derivative Matrix
for i = 1:nx
    for j = 1:ny
        for k = 1:nz
            if i > 1
                D(m(i,j,k),m(i-1,j,k)) = 1;
            end
            if i < nx
                D(m(i,j,k),m(i+1,j,k)) = 1;
            end
            if j > 1
                D(m(i,j,k),m(i,j-1,k)) = 1;
            end
            if j < ny
                D(m(i,j,k),m(i,j+1,k)) = 1;
            end
            if k > 1
                D(m(i,j,k),m(i,j,k-1)) = 1;
            end
            if k < nz
                D(m(i,j,k),m(i,j,k+1)) = 1;
            end
            D(m(i,j,k),m(i,j,k)) = -6;
        end
    end
end

D = -D/h^2;

% Add diagonal component associated with the potential
%V = zeros(length(x),length(y));
% 2D Harmonic Oscillator
%V = x'.^2 * ones(1,nx) + ones(1,ny)' * y.^2;
%psi0=cos(0.5*pi*x/5)'*cos(0.5*pi*y/5);


% 3D Infinite well
V = zeros(nx,ny,nz);
V(abs(x)>xmax/2,:,:) = well_potential;
V(:,abs(y)>ymax/2,:) = well_potential;
V(:,:,abs(z)>zmax/2) = well_potential;
psi0=zeros(nx,ny,nz);
psi0(round(nx/2),round(ny/2),round(nz/2)) = 1;

D=D+diag(reshape(V,1,nx*ny*nz));

psi=reshape(psi0,nx*ny*nz,1);

D=sparse(D);

% Analytic solution for lowest energy wave function
%xa=-5:0.1:5;
%psi_an=exp(-0.5*xa.^2);


[sol E] = eigs(D, num_states, 'sa');

for energy_level = 1:num_states
    
    figure(energy_level);
    psi = sol(:,energy_level);
    
    % Plot 95% surface
    cxy = zeros(nx,ny);
    cx = zeros(nx,1);
    cy = zeros(ny,1);
    cz = zeros(nz,1);
    summed = 0;
    [sortedValues,sortIndex] = sort(psi.^2/(sum(psi.^2)),'descend');
    for item = 1:nx*ny*nz
        if summed > 0.95
            break;
        end
        summed = summed + sortedValues(item);
    end
    
    p = patch(isosurface(x,y,z,reshape(psi,nx,ny,nz),sqrt(sortedValues(item))));
    set(p,'FaceColor','red','EdgeColor','none');
    q = patch(isosurface(x,y,z,reshape(psi,nx,ny,nz),-sqrt(sortedValues(item))));
    set(q,'FaceColor','blue','EdgeColor','none');
    camlight
    lighting gouraud
    axis equal
    axis([xmin xmax ymin ymax zmin zmax]);
    axis equal
    axis([xmin xmax ymin ymax zmin zmax]);
    
    
    % Calculating Energy of eigenstate
    energy_est = E(energy_level,energy_level);
    disp('Estimate of eigenvalue E: ');
    disp(energy_est);
    %disp(E_analy);
    
end


E_analy = zeros(5,5);
for n_x =1:5
    for n_y = 1:5
        for z = 1:5
        E_analy(n_x,n_y,nz) = (n_x*pi/(xmax-xmin) + ...
            n_y*pi/(ymax-ymin) +nz*pi/(zmax-zmin))^2;
        % E_analy(n_x,n_y) = n_x^2*pi^2/(xmax-xmin)^2 + ...
        %    2*n_x*n_y*pi^2/((xmax-xmin)*(ymax-ymin)) + ...
        %n_y^2*pi^2/(ymax-ymin)^2;
        end
    end
end

