clear all;

% Note on non-dimensionalisation. For simplicity, all variables in this code have been
% non-dimensionalised. This means, depending on their units, they have been
% normalised with appropriate combinations of R (the radius of a single
% sphere [m]), mu (the dynamic viscosity [kg/m/s]), and V (the translational
% velocity of the particle [m/s]). That means, instead of V, R or mu, in the code,
% you will see just see 1. Also, that means the output is non-dimensional.
% For example, the non-dimensional drag on a single sphere = 6*pi. To
% 'dimensionalise' the drag, should you want to, you need to multiply it by
% mu*V*R. 



[P, rSphere]=getSpherePositions(); %subroutine that returns number and position of spheres


psi=1.35; % the under/over-relaxation factor.

N=1000; % number of boundary nodes
V=[0 1 0]'; % §\eq{\bm{V}}§, particle translational velocity

rb=pointonsphere(N);N=length(rb);
M=floor(0.8*N);
rs=0.2*pointonsphere(M); M=length(rs);  % site positions are on a scaled-down sphere.


A=matrixConstructFast(rs,rb,M,N); 
INVA=pinv(A);



AIJ=zeros(3*N,3*M,P,P);
for II=1:P
    for IJ=1:P
        if (II~=IJ)
         rRel=rSphere(IJ,:)-rSphere(II,:);
         AIJ(:,:,II,IJ)=matrixConstructFast(rs,rb+rRel,M,N);
        end  
    end
end


for II=1:P
    RHSI(:,II)=constructRHS(rb,N,V,[0 0 0]');
end


FIc=zeros(M*3,P);
FI=FIc;


d_rms=rms(rms(RHSI(:,:)));
IT=0;
zeta=1;

while(zeta>1e-4 & zeta<10)
    IT=IT+1;
    for II=1:P
        FI(:,II)=INVA*psi*RHSI(:,II);
        FIc(:,II)=FIc(:,II)+FI(:,II);
        RHSI(:,II)=(1-psi)*RHSI(:,II);
        for IJ=1:P
            if (II~=IJ)
                RHSI(:,IJ)=RHSI(:,IJ)-AIJ(:,:,II,IJ)*FI(:,II);
            end
        end
    end
    
    zeta=rms(rms(RHSI(:,:)))/d_rms;
    
end



if (zeta>10)
    disp("Did not converge, reduce under-relaxation factor")
else
    disp("Converged in " + IT + " iterations")

    dragIZ(1:P)=0;
    for IJ=1:P
        for II=1:M
            dragIZ(IJ)=dragIZ(IJ)+FIc(3+(II-1)*3,IJ);
        end
    end
    StokesDrag=6*pi;
end

dragTotal = sum(dragIZ)/(6*pi)
drag = sum(dragIZ)/(6*pi*P)