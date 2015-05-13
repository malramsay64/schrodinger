% independant.m
% When dealing with the two electron case we make the assumption that the
% electrons are non interacting.
% Solves the problem by inverse power iteration

% Clear memory and show only a few digits
clear all; format short;

% Size of System
L = 100;

% Define values of the independent variable
xmin = -1;
xmax = 1;
h=(xmax-xmin)/L;
x=xmin:h:xmax;
L=L+1; % Passing through zero

% Electric Constant
e_const = 1;

% Number of matrix iterations
nits=15;

% Initial guess for the wave function
psi0=zeros(L,1);
psi1=zeros(L,1);
psi0(abs(x) <= xmax/2)=cos(0.5*pi*2*x(abs(x) <= xmax/2));
psi1(abs(x) <= xmax/2)=cos(pi*2*x(abs(x) <= xmax/2));
%psi0 = psi0./(sum(psi0)*h);
Psi=zeros(2,L);
Psi(1,:) = psi1;
Psi(2,:) = psi0;
energy = zeros(2,1);

% Construct matrix D2 associated with the derivative
D2=-2*eye(L);
D2=D2+diag(ones(L-1,1),+1)+diag(ones(L-1,1),-1);
D2=D2/h^2;

V = zeros(L,1);
V0 = zeros(L,1);
V0(abs(x)>xmax/2) = 1e6;

for j = 1:nits
    for i = 1:2
        % Add diagonal component associated with the potential
        
        V = V0;
        for k = 1:L/2
            k = k+abs(xmin/h)/2;
            p = Psi(mod(i,2)+1,:).^2*h;
            o = (abs((x-xmin)-k*h)).^-1;
            o(o>2/h) = 0;
            V(k) = V(k) + e_const*sum(o.*p);
        end
        
        A=-D2+diag(V);
        psi = Psi(i,:)';
        % Compute inverse matrix
        B=inv(A);
        for n=1:nits
            % Inverse iteration & normalisation
            psi=B*psi;
            psi=psi/max(abs(psi)); 
        end
        Psi(i,:) = psi;
        energy(i) = mean(Psi(1,:)'./(A*Psi(1,:)'));
        
        % Plot Functions
        figure(1)
        plot(x,Psi(1,:),'ko',x,Psi(2,:),'r*',...
            [xmax/2 xmax/2],[0 1],'k', [xmin/2 xmin/2],[0 1], 'k')
        xlabel('x'); ylabel('\Psi');
        anno=legend('Numerical solution','Initial guess');
        set (anno,'Box','off','Location','NorthWest')
        drawnow;
        %pause(0.5)
        
    end
    
end
% Determine eigenvalue
energy_est=sum(energy);

figure(2)
% Plot wave function, initial wave function
plot(x,((Psi(1,:)).*Psi(2,:))./max(Psi(1,:).*Psi(2,:)),'ko',x,psi0.*psi0,'r*',...
    [xmax/2 xmax/2],[0 1],'k', [xmin/2 xmin/2],[0 1], 'k')
xlabel('x'); ylabel('\Psi');
anno=legend('Numerical solution','Initial guess');
set (anno,'Box','off','Location','NorthWest')

disp('Estimate of eigenvalue E: ');
disp(energy_est);