% hydrogen.m
% Solving the Schrodinger equation for the hydrogen atom with a 1/r
% potential and plotting a 95% surface of the top points.
% The eigenvalues are found using the built in Matlab function eigs which
% finds the lowest n eigenvalues and eigenvectors much faster than
% computing all of them.
% Sparse matricies are used for the Hamiltonian to save space as matricies
% get very large

% Clear memory and figures
clear all; close all;

% Axis size, equal for x,y,z
axis_min = -60;
axis_max = 60;

% Size of system
L = 80;
L3 = L^3;

% Number of states to plot
num_states = 50;

% Energy Correction 2*Rydberg energy
e_corr = 27.2;

% Constucting coordinates
p = linspace(axis_min,axis_max,L);
[X,Y,Z] = meshgrid(p,p,p);
X = X(:);
Y = Y(:);
Z = Z(:);
R = sqrt(X.^2+Y.^2+Z.^2);

% Stepsize
h = p(2)-p(1);

% Construct Laplacian matrix
e = ones(L,1);
D = spdiags([e -2*e e], -1:1, L, L)/h^2;
I = speye(L);
D2 = kron(kron(D,I),I) + kron(kron(I,D), I) + kron(kron(I,I),D);

% Coulomb potential
Vext = -1./R;

% Compute eigenvectors and values
[psis E] = eigs(-0.5*D2+spdiags(Vext, 0, L3, L3), num_states, 'sa');

% Show lowest energy level
disp(['Total energy for Hydrogen atom ' num2str(E(1)*e_corr,5) ' eV.']);

% Ploting isosurfaces
for i= 1:num_states
    psi = psis(:,i);
    figure(i)
    summed=0;
    % Find value on which to put isosurface
    [sortedValues,sortIndex] = sort((psi).^2/(sum((psi).^2)),'descend');
    for item = 1:L3
        if summed > 0.90
            break;
        end
        summed = summed + sortedValues(item);
    end
    % Plot positive surface blue
    pos = patch(isosurface(p,p,p,reshape(psi,L,L,L),sqrt(sortedValues(item))));
    set(pos,'FaceColor','blue','EdgeColor','none');
    % Plot negative surface red
    neg = patch(isosurface(p,p,p,reshape(psi ,L,L,L),-sqrt(sortedValues(item))));
    set(neg,'FaceColor','red','EdgeColor','none');
    title(['E=' num2str(E(i,i)*27.21,5)],'fontsize',18);
    xlabel('x', 'fontsize', 16)
    ylabel('y', 'fontsize', 16)
    zlabel('z', 'fontsize', 16,'Rotation',1)
    camlight
    lighting gouraud
    axis equal
    axis([axis_min axis_max axis_min axis_max axis_min axis_max])
end