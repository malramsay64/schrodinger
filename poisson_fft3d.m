function phi = poisson_fft3d(charge_density,h,L)
% Function to solve the 3D Poisson equation using FFTs (and hence 
% implicit periodic boundary conditions) for the charge density,
% charge_density
% charge_denstiy is an array of size (L^3,1), and the resulting phi is
% returned in this format.

% Fourier transform charge density
sigma = reshape(charge_density,L,L,L);
fsigma=fftn(sigma);

% Construct arrays for FTs of phi
fphi=zeros(L,L,L);

% Loop over m and n and construct FTs of phi
for k=1:L
    for m=1:L
        for n=1:L
            % Denominator in FT of solution at grid point
            den=cos(2*pi*(k-1)/(L))+cos(2*pi*(m-1)/(L))...
                +cos(2*pi*(n-1)/(L))-3;
            if (den ~= 0)  
                % Construct FT of phi at grid point
                fphi(k,m,n)=-0.5*h^2*fsigma(k,m,n)/den;     
            end
        end
    end
end

% Invert FTs to obtain phi(x,y,z)
phi=real(ifftn(fphi));

% Reshape
phi = reshape(phi,L^3,1);

