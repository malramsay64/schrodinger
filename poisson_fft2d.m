function phi = poisson_fft2d(charge_density,h,L)


% % Solve the Poisson equation using FFTs (and hence implicit periodic
% % boundary conditions) for the charge density charge_density
% 
% % Charge density: a delta function at the centre of the region
% sigma=charge_density;
% 
% % Fourier transform charge density
fsigma=fftn(charge_density);
% 
% % Construct arrays for FTs of phi, Ex, Ey
%  fphi=zeros(L);
% % fEx=zeros(L,L);
% % fEy=zeros(L,L);
% % Initial phi
% phi=zeros(L);
% 
% % Charge density: a delta function at the centre of the region
% % sigma=zeros(L);
% % sigma(round(L/2),round(L/2))=1/h^2;
% %sigma(L,L)=1/h^2;
% 
% % Fourier transform charge density
% % fsigma=fft2(sigma);
% 
% % Construct arrays for FTs of phi, Ex, Ey
% fphi=zeros(L);
% fEx=zeros(L);
% fEy=zeros(L);
% phi =zeros(L,L);
% % Loop over m and n and construct FTs of phi, Ex, Ey
for m=1:L
  for n=1:L
      for k=1:L
      
    % Denominator in FT of solution at grid point
    den=cos(2*pi*(m-1)/L)+cos(2*pi*(n-1)/L)-2;
    if (den ~= 0)
      
      % Construct FT of phi at grid point
      fphi(m,n)=-0.5*h^2*fsigma(m,n)/den;
    end
    end
  end
end

% for k=1:L
%     for m=k:L
% 
%             den=cos(2*pi*(k-1)/L)+cos(2*pi*(m-1)/L)-2;
%     if (den ~= 0)
%       
%       % Construct FT of phi at grid point
%       fphi(k,m)=-0.5*h^2*fsigma(k,m)/den;
%             % Denominator in FT of solution at grid point
%             
%             
% %             den=cos(2*pi*(k-1)/L)+cos(2*pi*(m-1)/L)-2;
% %             if (den ~= 0)
% %                 
% %                 % Construct FT of phi at grid point
% %                 fphi(k,m)=-0.5*h^2*fsigma(k,m)/den;
% %                 
%              end
%     end
% end
%                 % Construct FTs of electric field components
%                 fEx(l,m)=-sqrt(-1)*sin(2*pi*(l-1)/L)*fphi(l,m)/h;
%                 fEy(l,m)=-sqrt(-1)*sin(2*pi*(m-1)/L)*fphi(l,m)/h;
%                 
%             end
%     end
% end
% 
% % Invert FTs to obtain phi(x,y), Ex(x,y), Ey(x,y)
% phi=real(ifft2(fphi));
% Ex=real(ifft2(fEx));
% Ey=real(ifft2(fEy));
phi = fphi;


