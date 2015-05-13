% independant.m
% When dealing with the two electron case we make the assumption that the
% electrons are non interacting.
% Solves the problem by inverse power iteration

% Clear memory and show only a few digits
clear all; format short;

% Size of System
L = 1000;

% Define values of the independent variable
xmin = -1;
xmax = 1;
h=(xmax-xmin)/L;
x=xmin:h:xmax;
L=L+1; % Passing through zero

% Number of matrix iterations
nits=10;

% Construct matrix D2 associated with the derivative
D2=-2*eye(L);
D2=D2+diag(ones(L-1,1),+1)+diag(ones(L-1,1),-1);
D2=D2/h^2;

% Add diagonal component associated with the potential
V = zeros(L,1);
V(abs(x)>xmax/2) = 1e5;

A=-D2+diag(V);

% Compute inverse matrix
B=inv(A);

% Initial guess for the wave function
psi0=zeros(L,1);
psi0(abs(x) <= xmax/2)=cos(0.5*pi*2*x(abs(x) <= xmax/2));
Psi=zeros(2,L);

for i=1:2;
    psi=psi0;
    
    figure(1);
    for n=1:10
        
       
        
        % Inverse iteration & normalisation
        psi=B*psi;
        psi=psi/max(abs(psi));
        
    end
    Psi(i,:) = psi;
end
% Determine eigenvalue
energy_est=mean(Psi(1,:)'./(B*Psi(1,:)'))+mean(Psi(2,:)'./(B*Psi(2,:)'));

 % Plot wave function, initial wave function
        plot(x,Psi(1,:).*Psi(2,:),'ko',x,psi0.*psi0,'r*',...
            [xmax/2 xmax/2],[0 1],'k', [xmin/2 xmin/2],[0 1], 'k')
        xlabel('x'); ylabel('\Psi');
        anno=legend('Numerical solution','Initial guess');
        set (anno,'Box','off','Location','NorthWest')
        title(['Iteration: ',num2str(n)]);
        
disp('Estimate of eigenvalue E: ');
disp(energy_est);