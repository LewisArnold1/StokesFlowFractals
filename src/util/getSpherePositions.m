function [P rSphere] = getSpherePositions()

P=2;  % number of particles/spheres
sep=[2 0 0]; % [2 0 0]= chain of touching spheres extending in x direction
rSphere(1:P,1:3)=0;
rSphere(1,1:3)=[0 0 0];
for II=2:P
    rSphere(II,1:3)=rSphere(II-1,1:3)+sep;  %location of spheres
end