function [Ft, Jt] = getFt(x, refTwist)
%%
global GJ nv ne voronoiRefLen

Ft = x * 0;
Jt = zeros(length(x), length(x));
for c=2:ne
    
    node0 = [x(4*(c-2)+1); x(4*(c-2)+2); x(4*(c-2)+3)]';
    node1 = [x(4*(c-1)+1); x(4*(c-1)+2); x(4*(c-1)+3)]';
    node2 = [x(4*c+1);     x(4*c+2);     x(4*c+3)]';
    
    theta_e = x(4*(c-1));
    theta_f = x(4*c);
    
    [dF, dJ] = ...
        gradEt_hessEt(node0, node1, node2, ...
        theta_e, theta_f, refTwist(c), ...
        voronoiRefLen(c), GJ);

    % Arrange in the global force vector
    ci = 4*(c-1) + 1 - 4;
    cf = 4*(c-1) + 1 + 6;

    Ft( ci: cf) = Ft( ci: cf) - dF;
    Jt( ci: cf, ci: cf) = Jt( ci: cf, ci: cf) - dJ;
end


end
