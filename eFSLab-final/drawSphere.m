function drawSphere(x, y, z, r)

phi=linspace(0,pi,30);
theta=linspace(0,2*pi,40);
[phi,theta]=meshgrid(phi,theta);

v1=x + r*sin(phi).*cos(theta);
v2=y + r*sin(phi).*sin(theta);
v3=z + r*cos(phi); 

% [l,c]=size(v1);
% c=rand(l,c);

mesh(v1, v2, v3, c);


