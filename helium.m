% helium.m
% Solves the Schrodinger equation of the helium atom using the DFT
% approximation. 
% Eigenvalues are computed using the built in matlab function eigs, finding
% the lowest eigenvalue. This is much faster than other methods.
% The Poisson equation for the electron-electron coerrelation is computed
% using a FFT implementation.

% Clear memory and figures
clear all; close all;

% Axis size, equal for x,y,z
axis_min = -5;
axis_max = 5;

% Size
L = 40;
L3 = L^3;

% Convergence tolerance
tol = 1e-6;

% Energy Correction 2*Rydberg energy
e_corr = 27.211;

% Constucting coordinates
p=linspace(axis_min,axis_max,L);
[X,Y,Z] = meshgrid(p,p,p);
X=X(:);
Y=Y(:);
Z=Z(:);
R=sqrt(X.^2+Y.^2+Z.^2);

% Stepsize
h=p(2)-p(1);

% Construct Laplacian matrix
e=ones(L,1);
D=spdiags([e -2*e e], -1:1, L, L)/h^2;
I = speye(L);
D2 = kron(kron(D,I),I) + kron(kron(I,D), I) + kron(kron(I,I),D);

% Potential Energy from nucleus
Vext = -2./R;

% Compensation Charges -> Better approximates with boundary conditions
ncomp = exp(-R.^2/2);
% Normalisation
ncomp = -2*ncomp / sum(ncomp) / h^3;
% Compensation charge solution
ncomppot = -2./R.*erf(R/sqrt(2));

% Initial Estimate of potential
Vtot = Vext;
Eold = 0;

% Self Consistent Field Iteration
for i=1:50
    % Solving Schrodinger equation with adjusted potential
    [PSI, E] = eigs(-0.5*D2+spdiags(Vtot, 0, L3, L3), 1, 'sa');
    
    %Normalising
    PSI = PSI / h^(3/2);
    
    % Density -> Electrons per unit volume
    n = 2*PSI.^2;
    
    % Exchange Potential
    Vx = -(3/pi)^(1/3)*n.^(1/3);
    
    % Hartree Potential, Solving Poission equation in atomic units
    %L3 Vh = -4*pi*n;
    rho=reshape(n+ncomp,L,L,L);
    Vh = poisson_fft3d(4*pi*rho,h,L)-ncomppot;
    
    % Total Potential
    Vtot = Vx + Vh + Vext;
    
    % Kinetic energy
    T = 2 * PSI'*(-0.5*D2)*PSI * h^3;
    
    % External Potential, Integration over all space
    Eext = sum(n.*Vext)*h^3;
    
    % Hartree Potential
    Eh = 0.5 * sum(n.*Vh)*h^3;
    
    % Exchange Potential energy -> Integrate Potential
    Ex = sum(-(3/4)*(3/pi)^(1/3)*n.^(4/3))*h^3;
    Etot = T + Eext + Eh + Ex;
    
    % Convergence criteria
    if abs(Etot - Eold) < tol
        break;
    end
    
    % Update values
    Eold = Etot;
    
end

% Display energy values
disp(['Eigenvalue ' num2str(E,5) ]);
disp(['Kinetic energy ' num2str(T,5) ]);
disp(['Exchange energy ' num2str(Ex,5) ]);
disp(['External energy ' num2str(Eext,5) ]);
disp(['Potential energy ' num2str(Eh,5) ]);
disp(['Total energy ' num2str(Etot,5) ]);

% Plotting Surface
summed=0;
% Find value on which to put isosurface
[sortedValues,sortIndex] = sort(PSI.^2/(sum(PSI.^2)),'descend');
for item = 1:D2
    if summed > 0.90
        break;
    end
    summed = summed + sortedValues(item);
end
% Plot positive surface blue
pos = patch(isosurface(p,p,p,reshape(PSI,L,L,L),sqrt(sortedValues(item))));
set(pos,'FaceColor','blue','EdgeColor','none');
% Plot negative surface red
neg = patch(isosurface(p,p,p,reshape(PSI ,L,L,L),-sqrt(sortedValues(item))));
set(neg,'FaceColor','red','EdgeColor','none');
title(['E=' num2str(E*27.21,5)],'fontsize',16);
camlight
lighting gouraud
axis equal
axis([axis_min axis_max axis_min axis_max axis_min axis_max])


