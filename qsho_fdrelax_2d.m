% qsho_fdrelax.m
% Solve the quantum simple harmonic oscillator for the lowest
% energy wave function using inverse power iteration relaxation
% of the linear system obtained by finite differencing

% Clear memory and show only a few digits
clear all; format short; clf; close all;

% Define values of the independent variable
xmin = -1;
xmax = 1;
ymin = -1;
ymax = 1;
h=0.05;
x=xmin:h:xmax;
y=ymin:h:ymax;
nx = length(x);
ny = length(y);
well_potential = 1e10;
num_states = 1;
thresh = 1e-5;
nits = 100;

% Indexing function
m = @(p,q) (p-1)*nx + q;
% Reverse Indexing
rev = @(k) [floor((k-1)/nx)+1 mod((k-1),nx)+1];

% Size of the linear system
L=length(x)*length(y);

D = zeros(nx*ny,nx*ny);
e=zeros(nx,1);
O=spdiags([e -2*e e], -1:1, nx, nx)/h^2;
I = speye(nx);
L2 = kron(I,O) + kron(I,O);
% Derivative Matrix
for i = 1:nx
    for j = 1:nx
        if i > 1
            D(m(i,j),m(i-1,j)) = 1;
        end
        if i < nx
            D(m(i,j),m(i+1,j)) = 1;
        end
        if j > 1
            D(m(i,j),m(i,j-1)) = 1;
        end
        if j < ny
            D(m(i,j),m(i,j+1)) = 1;
        end
        D(m(i,j),m(i,j)) = -4;
    end
end

D = -D/h^2;

% Add diagonal component associated with the potential
%V = zeros(length(x),length(y));
% 2D Harmonic Oscillator
%V = x'.^2 * ones(1,nx) + ones(1,ny)' * y.^2;
%psi0=cos(0.5*pi*x/5)'*cos(0.5*pi*y/5);


% 2D Infinite well
V = zeros(nx,ny);
V(abs(x)>xmax/2,:) = well_potential;
V(:,abs(y)>ymax/2) = well_potential;
psi0=cos(pi*x/5)'*cos(pi*y/5);

D=D+diag(reshape(V,1,nx*ny));

psi0=reshape(psi0,(nx)*(ny),1);

D=sparse(D);
% Analytic solution for lowest energy wave function
%xa=-5:0.1:5;
%psi_an=exp(-0.5*xa.^2);

tic
%sol = eig(D);
% Inverse
B = inv(D);

%[L,U,P] = lu(D);

for energy_level = 1:num_states
    psi = psi0;
    psi_old = psi;
    figure(energy_level);
    for n=1:nits
        
        % Inverse iteration & normalisation
        %psi = B*psi;
        %Y = L\psi;
        %psi = U\Y;
        psi = B*psi;
        psi = psi/max(psi);
        %psi=psi_new/max(abs(psi_new(:)));
        if abs(sum(psi_old.^2-psi.^2)) < thresh
        %    break;
        end
        
        psi_old = psi;
        
    end
    
    % Plot wave function at each iteration, initial wave
    % function, and analytic solution
    subplot(1,2,1);
    surface(x,y,reshape(psi,nx,ny))
    
    %axis([xmin xmax ymin ymax -1 1]);
    xlabel('x'); ylabel('\Psi');
    title(['Iteration: ',num2str(n)]);
    drawnow;
    
    % Sort values in descending order
    cx = zeros(nx,1);
    cy = zeros(ny,1);
    [sortedValues,sortIndex] = sort(psi.^2/(sum(psi.^2)),'descend');
    summed = 0;
    for item = 1:nx*ny
        if summed < 0.95
            ind = rev(sortIndex(item));
            cx(item) = ind(1);
            cy(item) = ind(2);
        else
            break;
        end
        summed = summed + sortedValues(item);
        
    end
    subplot(1,2,2);
    scatter(cx.*h+xmin,cy.*h+ymin);
    axis([xmin xmax ymin ymax]);
    
    % Calculating Energy of eigenstate
    %Y = L\psi;
    %energy_est = psi./(L\psi);
    energy_est = psi./(B*psi);
    energy_est = reshape(energy_est,nx,ny);
    energy_est = mean(mean(energy_est(abs(x)<=xmax/2, abs(y)<=ymax/2)));
    
    disp('Estimate of eigenvalue E: ');
    disp(energy_est);
    %disp(E_analy);
    
    % Deflation
    %B = B - (psi*psi')/(energy_est*psi'*psi);
    B = B - (psi*psi')/(energy_est*psi'*psi);
    %psiA = psi/sqrt(sum(psi.^2));
    %D = D-(psiA'*D*psiA);
    %[L,U,P] = lu(D);
    %U = U-diag(energy_est*ones(nx*ny,1));
end
toc

E_analy = zeros(5,5);
for n_x =1:5
    for n_y = 1:5
        E_analy(n_x,n_y) = (n_x*pi + n_y*pi)^2/(xmax-xmin)^2;
        % E_analy(n_x,n_y) = n_x^2*pi^2/(xmax-xmin)^2 + ...
        %    2*n_x*n_y*pi^2/((xmax-xmin)*(ymax-ymin)) + ...
        %n_y^2*pi^2/(ymax-ymin)^2;
    end
end

