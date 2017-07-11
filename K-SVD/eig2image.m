function [Lambda1,Lambda2,Ix1,Ix2,Iy1,Iy2]=eig2image(Dxx,Dxy,Dyy)
% This function eig2image calculates the eigen values from the
% hessian matrix, sorted by abs value. And gives the direction
% of the ridge (eigenvector smallest eigenvalue) .
% 
% [Lambda1,Lambda2,Ix,Iy]=eig2image(Dxx,Dxy,Dyy)
%

%
% | Dxx  Dxy |
% |          |
% | Dxy  Dyy |


% Compute the eigenvectors of J, v1 and v2
tmp = sqrt((Dxx - Dyy).^2 + 4*Dxy.^2);
v2x = 2*Dxy; v2y = Dyy - Dxx + tmp;

% Normalize
mag = sqrt(v2x.^2 + v2y.^2); i = (mag ~= 0);
v2x(i) = v2x(i)./mag(i);
v2y(i) = v2y(i)./mag(i);

% The eigenvectors are orthogonal
v1x = -v2y; 
v1y = v2x;

% Compute the eigenvalues
mu1 = 0.5*(Dxx + Dyy + tmp);
mu2 = 0.5*(Dxx + Dyy - tmp);

% Sort eigen values by absolute value abs(Lambda1)<abs(Lambda2)
check=abs(mu1)>abs(mu2);

Lambda1=mu1; Lambda1(check)=mu2(check);
Lambda2=mu2; Lambda2(check)=mu1(check);


Ix=v1x; Ix(check)=v2x(check);
Iy=v1y; Iy(check)=v2y(check);

% e11 = v1x;
% e12 = v2y;
% e21 = v2x;
% e22 = v2y;

Ix1=v1x;
Ix2=v1y;

Iy1=v2x;
Iy2=v2y;





