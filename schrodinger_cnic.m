% schrodinger_cnic.m
% Solve the 1-D time-dependent Schrodinger equation for an initial 
% Gaussian wave packet using Crank-Nicolson

% Clear memory and show only a few digits
clear all; format short;

% Time step and spatial step
tau=0.05;
h=0.005;

% Height of potential barrier
V0=0.25e+6; 

% Parameters of initial wave function
k0=50;  % Average wavenumber
s0=0.05; % Width of Gaussian

% Total integration time and number of steps
tint=1e2*(1+h)/k0;
nsteps=floor(tint/tau)+1;

% Vector of x values
x=0:h:1;
L=length(x);

% Construct Hamiltonian matrix...
H=-2*eye(L);
H=H+diag(ones(L-1,1),+1)+diag(ones(L-1,1),-1);

% ...with periodic boundary conditions
H(1,L)=1;
H(L,1)=1;

% with Dirichlet BC
%H(1,:) = 0;
%H(end,:) = 0;

H=-0.5*H/h^2;

% Construct potential matrix corresponding to barrier
V=zeros(L,1);
ii=find((x > 0.75));
V(ii)=V0;
V(x<0.25) = V0;
V(x>0.75) = V0;

Vmat=diag(V);

% Add potential to Hamiltonian matrix
H=H+Vmat;

%H = H+diag(V);

% Construct matrix for the linear system solved at each step of 
% Crank-Nicolson. Note that i is sqrt(-1) by default
A=0.5*(eye(L)+0.5*tau*H);

% Initial wave function
C1=1./sqrt(s0*sqrt(pi));
psi=C1*exp(k0*x'); % Oscillatory part
psi=psi.*exp(-0.5*((x-0.5)/s0)'.^2); % Gaussian envelope
%psi=cos(2*pi*x-0.5)'; 

% Scale for axis
max_psi=max(abs(psi));
    
% March forwards in time
time(1)=0;
figure(1);
for n=1:nsteps

    % Update the time
    time(n+1)=time(n)+tau;   
    
    % Perform Crank-Nicolson update
    chi=A\psi;
    psi=chi-psi;
    
%     % Analytic solution for probability density
%     sig=sqrt(s0^2+time(n+1)^2/s0^2);
%     C2=1/(sqrt(pi)*sig);
%     pdens_an=C2*exp(-(x-0.5)'.^2/sig^2);
%     k=round(k0*time(n+1)/h); % k is the index the peak has reached
%     pdens_an=pdens_an([end-k+1:end 1:end-k]); % Shift to right by k
    
    % Plot psi versus position and annotate
    subplot(2,1,1),plot(x,real(psi),'b',x,imag(psi),'r');
    axis([0 1 -max_psi max_psi]);
    title(['Blue: Re(\psi),   Red: Im(\psi)      Time: ',...
        num2str(time(n+1))]);
    ylabel('\psi');
    subplot(2,1,2),plot(x,abs(psi).^2,'b')
    axis([0 1 0 max_psi^2]);
    title('Blue: Numerical solution   Green: Analytic solution');
    xlabel('Position');
    ylabel('|\psi|^2');
    drawnow;
   
end
