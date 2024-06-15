function kappa = getkappa( x, m1, m2 )
nv = (length(x)+1)/4;
ne = nv-1;

% Computer Kappa
kappa = zeros(nv, 2);
for c=2:ne

    node0 = x(4*(c-2)+1:4*(c-2)+3);
    node1 = x(4*(c-1)+1:4*(c-1)+3);
    node2 = x(4*(c-0)+1:4*(c-0)+3);
    m1e = m1(c-1,:);
    m2e = m2(c-1,:);
    m1f = m1(c,:);
    m2f = m2(c,:);
    
    kappaL = computekappa(node0, node1, node2, m1e, m2e, m1f, m2f );

    kappa(c,1) = kappaL(1);
    kappa(c,2) = kappaL(2);
end

end
