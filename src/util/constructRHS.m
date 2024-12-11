function [d] =constructRHS(rb, N, V,Om)
% Input: §\mcommentfont $r^b_i$§, node positions (relative to centre of rotation)
% Input: §\mcommentfont $N$§, number of boundary nodes
% Input: §\mcommentfont $V$§, translational velocity of particle
% Input: §\mcommentfont $\Omega$§, angular velocity of particle
% Output: §\mcommentfont${d}_{p}$§, right-hand-side vector

d=zeros(N,1); % initialize RHS vector, §\mcommentfont $d_p$§

for b=1:N
    % calculate particle wall velocity, §\mcommentfont $\bm{v}^b_{\mathrm{wall}}=\bm{V}+\bm{\Omega}\times\bm{r}^b $§
    v_wall(b,:)=V+cross(Om,rb(b,:)'); 
    
    for i=1:3
        p=3*(b-1)+i;  		% §\mcommentfont Eq.\,(\ref{bijection})§
        d(p)=v_wall(b,i); 	% §\mcommentfont ${d}_{p}$\, , Eq.\,(\ref{newA})§
    end
end