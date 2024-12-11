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

o8pi=1/(8*pi);

[P, rSphere]=getSpherePositions(); %subroutine that returns number and position of spheres


psi=1.35; % the under/over-relaxation factor.

N=50; % number of boundary nodes
V=[0 1 0]'; % §\eq{\bm{V}}§, particle translational velocity

rb=pointonsphere(N);N=length(rb);
M=floor(0.8*N);
rs=0.2*pointonsphere(M); M=length(rs);  % site positions are on a scaled-down sphere.

A=matrixConstructFast(rs,rb,M,N); 
INVA=pinv(A);

rsL=reshape(rs',3*M,1);

% AIJ=zeros(3*N,3*M,P,P);
% for II=1:P
%     for IJ=1:P
%         if (II~=IJ)
%          rRel=rSphere(IJ,:)-rSphere(II,:);
%          AIJ(:,:,II,IJ)=matrixConstructFast(rs,rb+rRel,M,N);
%         end  
%     end
% end


for II=1:P
    RHSI(:,II)=constructRHS(rb,N,V,[0 0 0]');
end


FIc=zeros(M*3,P);
FI=FIc;


d_rms=rms(rms(RHSI(:,:)));
IT=0;
zeta=1;

while(zeta>1e-4 && zeta<10)
    
    IT=IT+1;
    

    
    for II=1:P
        FI(:,II)=INVA*psi*RHSI(:,II);
        FIc(:,II)=FIc(:,II)+FI(:,II);
        Ftemp=FI(:,II);
        RHSI(:,II)=(1-psi)*RHSI(:,II);
        for IJ=1:P
            if (II~=IJ)
                 rRel=rSphere(IJ,:)-rSphere(II,:);
                 rbD=rb+rRel;
               %   At=matrixConstructFast(rs,rb+rRel,M,N); 
                 %AIJ(:,:,II,IJ)=matrixConstructFast(rs,rb+rRel,M,N);
                %RHSI(:,IJ)=RHSI(:,IJ)-AIJ(:,:,II,IJ)*FI(:,II);
                %RHSI(:,IJ)=RHSI(:,IJ)-At*FI(:,II);
 
                tmp=0;
                RHStemp=RHSI(:,IJ);
                
                for b=1:N
                    rbDx=rbD(b,1);
                    rbDy=rbD(b,2);
                    rbDz=rbD(b,3);
                    p=3*b-2;
                    for s=1:M
                        q=3*s-2;  
                  
                        Fx=Ftemp(q);Fy=Ftemp(q+1);Fz=Ftemp(q+2);
                        rbsx=rbDx-rs(s,1);
                        rbsy=rbDy-rs(s,2);
                        rbsz=rbDz-rs(s,3);
%                          rbsx=rbDx-rsL(q);
%                          rbsy=rbDy-rsL(q2);
%                          rbsz=rbDz-rsL(q3);
                        rbsn2=rbsx*rbsx+rbsy*rbsy+rbsz*rbsz;
                        rbsn=rbsn2^0.5;
                        rbsn3=rbsn*rbsn2;             
                        fdotr=Fx*rbsx+Fy*rbsy+Fz*rbsz;                       
 %                      RHSI(p,IJ)=RHSI(p,IJ)-o8pi*(Fx/rbsn+rbsx*fdotr/rbsn3);
 %                      RHSI(p+1,IJ)=RHSI(p+1,IJ)-o8pi*(Fy/rbsn+rbsy*fdotr/rbsn3);
 %                      RHSI(p+2,IJ)=RHSI(p+2,IJ)-o8pi*(Fz/rbsn+rbsz*fdotr/rbsn3);
%                         c1=-o8pi/rbsn;
%                         c2=c1*fdotr/rbsn2;
%                         RHStemp(p)=RHStemp(p)+c1*Fx+c2*rbsx;
%                         RHStemp(p+1)=RHStemp(p+1)+c1*Fy+c2*rbsy;
%                         RHStemp(p+2)=RHStemp(p+2)+c1*Fz+c2*rbsz;
        
                 
                        RHStemp(p)=RHStemp(p)-o8pi*(Fx/rbsn+rbsx*fdotr/rbsn3);
                        RHStemp(p+1)=RHStemp(p+1)-o8pi*(Fy/rbsn+rbsy*fdotr/rbsn3);
                        RHStemp(p+2)=RHStemp(p+2)-o8pi*(Fz/rbsn+rbsz*fdotr/rbsn3);
                    end
                end
                RHSI(:,IJ)=RHStemp;
            end
        end
        
    end
    
    zeta=rms(rms(RHSI(:,:)))/d_rms
    
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

sum(dragIZ)/(6*pi*P)