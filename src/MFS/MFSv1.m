clear all
% non-dimensionalisation with §\eq{\mu}§ and some length scale §\eq{L}§ 

rb=load('points/rb_Sphere_M.txt'); %load nodes: §\eq{r^b_i}§ (1:N,1:3)
rs=load('points/rs_Sphere_M.txt'); %load sites: §\eq{r^s_j}§ (1:M,1:3)

%remove a few sites to benefit numerics
rs(1:8:length(rs),:)=[]; 

N=length(rb); % §\eq{N}§, number of boundary nodes
M=length(rs); % §\eq{M}§, number of singularity sites

V=[1 0 0]'; % §\eq{\bm{V}}§, particle translational velocity
Om=[0 0 0]'; % §\eq{\bm{\Omega}}§, particle angular velocity (about origin)

d=constructRHS(rb,N,V,Om); % construct §\eq{d_p}§
A=matrixConstruct(rs,rb,M,N); % construct §\eq{\mathbb{A}_{pq}}§

invA=pinv(A); % calculate Moore-Penrose pseudo inverse
x=invA*d; % solve for §\eq{x_{q}}§, Equation§\reff{eq_invEqs}§
F=0;tau=0; % initialise total force (F) and torque (tau)
for s=1:M        
    for j=1:3
        q=3*(s-1)+j;    % Equation§\reff{bijection}§
        f(s,j)=x(q);    % Equation§\reff{newA}§
    end
    
    % add force from §\eq{s^{\mathrm{th}}}§ singularity to total force
    F=F+f(s,:)';
    
    % add torque (generated around origin) from §\eq{s^{\mathrm{th}}}§ 
    % force to total torque 
    tau=tau+cross(rs(s,:),f(s,:))'; % (§\eq{\bm{\tau}^s=\bm{r}^s \times \bm{f}^s}§)
end

% force/torque on particle, equal and opposite to force/torque on fluid
MFS_drag=-F  
%MFS_torque=-tau
