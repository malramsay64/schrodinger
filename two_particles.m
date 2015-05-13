% Solving the Schrodinger equation for two particles in an infinite well
% Clear
clear; 

% Size
L = 1;

% Number of Steps
N = 10;

% Step Size
h = L/(N+1);

% (i-1)N + j

% Energy estimate
E = 1;

% Mappping of Matrix values
n = @(p,q) (p-1)*N + q;
m = @(i,j) (i-1)*N + j;

% Delta Function
delta = @(x) 0 + (abs(x) <= h/2)*1;


% Discretised Hamiltonian
% Kinetic energy, Centered Difference approximation used for derivatives
K = zeros(N^2,N^2);
for i = 1:N
    for j = 1:N
        for p = 1:N
            for q = 1:N
                K(n(p,q),m(i,j)) = ...
                    delta(q-j)*(delta(p-i-1)-2*delta(p-i)+delta(p-i+1))+...
                    delta(p-i)*(delta(q-j-1)-2*delta(q-j)+delta(q-j+1));
            end
        end
    end
end
K = -1/(2*h^2)*K;

% Coulomb Repulsion
V = zeros(N^2, N^2);
for i = 1:N
    for j = 1:N
        % Diagonal matrix so n=m for all values
        V(n(i,j),m(i,j)) = 1/(abs(j-i)*h);
    end
end

H = K+V;

% Remove elements that correspond to singularity
index = true(1, N^2);
for i = 1:N^2
    if H(i,i) == Inf
        index(i) = false;
    end
end
H = H(index,index);


% Basis Set
C = zeros(N^2,N^2);
for i= 1:N
    for j=1:N
        for p = 1:N
            for q = 1:N
                C(m(i,j),n(p,q))= delta(p-i)*delta(q-j);
            end
        end
    end
end
C = C(index,index);

S = diag(H);

E = H*C;

