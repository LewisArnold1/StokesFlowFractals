function [A] = matrixConstructFast(spos, bpos, totalStokeslets, totalNodes)
t1 =pi.^(-1);
t2 =(1/8).*t1;
t3 =(1/15).*t1;

A(1:3*totalNodes,1:3*totalStokeslets)=0;
for II=1:totalNodes
    b1=bpos(II,1);b2=bpos(II,2);b3=bpos(II,3);
    IIX=1+(II-1)*3;
    for IJ=1:totalStokeslets
        r1=b1-spos(IJ,1);r2=b2-spos(IJ,2);r3=b3-spos(IJ,3);
        R2=r1*r1+r2*r2+r3*r3;
        R=sqrt(R2);
        t4 =R.^(-1).*t2;
        t6 =R2.^(-1).*t4;
        IIY=1+(IJ-1)*3;
        A(IIX,IIY)=t4+r1*r1*t6; A(IIX+1,IIY)=r1*r2*t6; A(IIX+2,IIY)=r1*r3*t6;
        A(IIX,IIY+1)=r1*r2*t6; A(IIX+1,IIY+1)=t4+r2*r2*t6; A(IIX+2,IIY+1)=r2*r3*t6;
        A(IIX,IIY+2)=r1*r3*t6; A(IIX+1,IIY+2)=r2*r3*t6; A(IIX+2,IIY+2)=t4+r3*r3*t6;
    end
end



         
