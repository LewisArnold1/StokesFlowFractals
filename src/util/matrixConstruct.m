function [A] = matrixConstruct(rs, rb, M, N)
A=zeros(3*N,3*M);rbs=zeros(3,1);  %initialize
I=eye(3);  %define identity tensor

for b=1:N
	for s=1:M
		rbs=rb(b,:)-rs(s,:);	% §\mcommentfont$r^{bs}_{i}=r^{b}_{i}-r^{s}_{i}$§
		R=norm(rbs);			% §\mcommentfont$ | r^{bs}_{i} |$§
		
		J=1/(8*pi)*(I/R+rbs'*rbs/R^3);   % §\mcommentfont $\mathbb{J}^{bs}_{ij}$\, , Eq.\,(\ref{oseen}) §
    
        for i=1:3
            for j=1:3
                p=3*(b-1)+i;   % §\mcommentfont Eq.\,(\ref{bijection})§
                q=3*(s-1)+j;   % §\mcommentfont Eq.\,(\ref{bijection})§
                A(p,q)=J(i,j); % §\mcommentfont $\mathbb{A}_{pq}$\, , Eq.\,(\ref{newA})§
            end
        end       

    end
end
return