function [Fb, Jb] = getFb(x, m1, m2)
%%
global EI ne voronoiRefLen kappaBar

%% Compute grad kappa
Fb = x * 0;
Jb = zeros(length(x), length(x));

for c=2:ne
    node0 = [x(4*(c-2)+1); x(4*(c-2)+2); x(4*(c-2)+3)]';
    node1 = [x(4*(c-1)+1); x(4*(c-1)+2); x(4*(c-1)+3)]';
    node2 = [x(4*c+1);     x(4*c+2);     x(4*c+3)]';
    m1e = m1(c-1,:);
    m2e = m2(c-1,:);
    m1f = m1(c,:);
    m2f = m2(c,:);
    [dF, dJ] = ...
        gradEb_hessEb(node0, node1, node2, m1e, m2e, m1f, m2f, ...
        kappaBar(c, :), voronoiRefLen(c), EI);

    % Arrange in the global force vector
    ci = 4*(c-1) + 1 - 4;
    cf = 4*(c-1) + 1 + 6;
    Fb( ci: cf) = Fb( ci: cf) - dF;    
    Jb( ci: cf, ci: cf) = Jb( ci: cf, ci: cf) - dJ;
end


end
