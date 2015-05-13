% poisson_fft.m
% Solve the Poisson equation using FFTs (and hence implicit periodic
% boundary conditions)

% Clear memory and show only a few digits
clear all; format short;

% Spatial step
h=0.025;

% Vectors of x and y values
x=0:h:1; y=x; z=x;
L=length(x);

% Initial phi
phi=zeros(L);

% Charge density: a delta function at the centre of the region
sigma=zeros(L,L,L);
sigma(round(L/2),round(L/2),round(L/2))=1/h^3;
%sigma(L,L)=1/h^2;

% Fourier transform charge density
fsigma=fftn(sigma);

% Construct arrays for FTs of phi, Ex, Ey
fphi=zeros(L,L,L);
fEx=zeros(L,L,L);
fEy=zeros(L,L,L);
  
% Loop over m and n and construct FTs of phi, Ex, Ey
for m=1:L
  for n=1:L
      for k=1:L
    % Denominator in FT of solution at grid point
    den=cos(2*pi*(m-1)/L)+cos(2*pi*(n-1)/L)+cos(2*pi*(k-1)/L)-3;
    if (den ~= 0)
      
      % Construct FT of phi at grid point
      fphi(m,n,k)=-0.5*h^2*fsigma(m,n,k)/den;
      
    end
      end
  end
end

% Invert FTs to obtain phi(x,y), Ex(x,y), Ey(x,y)
phi=real(ifftn(fphi));

sigma=reshape(sigma,L^3,1);
phi2 = poisson_fft3d(sigma,h,L);
