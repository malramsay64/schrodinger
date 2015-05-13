% independant.m
% When dealing with the two electron case we make the assumption that the
% electrons are non interacting.
% Solves the problem by inverse power iteration

% Clear memory and show only a few digits
clear all; format short;

% Size of System
L = 100;

% Define values of the independent variable
xmin = -5;
xmax = 5;
h=(xmax-xmin)/L;
x=xmin:h:xmax;
L=L+1


% Number of matrix iterations
nits=10;

% Construct matrix D2 associated with the derivative
D2=-2*eye(L);
D2=D2+diag(ones(L-1,1),+1)+diag(ones(L-1,1),-1);
D2=D2/h^2;

% Add diagonal component associated with the potential
A=-D2+diag(x.^2);

% Compute inverse matrix
B=inv(A);

% Initial guess for the wave function
psi0=zeros(L,1);
psi0=cos(0.5*pi*x/5)';

Psi=zeros(2,L);

for i=1:2;
    psi=psi0;
    
    % Analytic solution for lowest energy wave function
    xa=-5:0.1:5;
    psi_an=exp(-0.5*xa.^2);
    
    figure(1);
    for n=1:10
        
        % Plot wave function at each iteration, initial wave
        % function, and analytic solution
        plot(x,psi,'ko',x,psi0,'r*',xa,psi_an,[0 0],[0 1],'k')
        xlabel('x'); ylabel('\Psi');
        anno=legend('Numerical solution','Initial guess',...
            'Analytic solution');
        set (anno,'Box','off','Location','NorthWest')
        title(['Iteration: ',num2str(n)]);
        drawnow;
        %pause(0.1)
        
        % Inverse iteration & normalisation
        psi=B*psi;
        psi=psi/max(abs(psi));
        
    end
    Psi(i,:) = psi;
end
% Determine eigenvalue
energy_est=mean(Psi(1,:)'./(B*Psi(1,:)'))+mean(Psi(2,:)'./(B*Psi(2,:)'));

disp('Estimate of eigenvalue E: ');
disp(energy_est);