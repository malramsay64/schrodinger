% Clear memory and show only a few digits
clear all; format short;

% Define values of the independent variable
h=0.5; 
x=-5:h:5;
y=-5:h:5;

nx = length(x);
ny = length(y);
% Size of the linear system
L=length(x)*length(y);
x = x(2:end-1);
y = y(2:end-1);
% Number of matrix iterations
nits=1;

D=zeros(nx,ny);

% Construct 3d matrix
D=D+diag(ones(L-1,1),+1)+diag(ones(L-1,1),-1);