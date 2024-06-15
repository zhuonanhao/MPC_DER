function [Fs, Js] = getFs(x)
%%
global EA nv ne refLen

%% Compute force
Fs = x * 0;
Js = zeros(length(x), length(x));
for c=1:ne
    node0 = [x(4*(c-1)+1); x(4*(c-1)+2); x(4*(c-1)+3)]';
    node1 = [x(4*c+1); x(4*c+2); x(4*c+3)]';
    [dF, dJ] = ...
        gradEs_hessEs(node0, node1, ...
        refLen(c), EA);
    
    ci = 4*(c-1)+1;
    cf = 4*(c-1)+3;
    Fs(ci:cf) = Fs(4*(c-1)+1:4*(c-1)+3) - dF(1:3);
    Fs(ci+4:cf+4) = Fs(4*c+1:4*c+3) - dF(4:6);
    
    Js(ci:cf, ci:cf) = Js(ci:cf, ci:cf) - dJ(1:3, 1:3);
    Js(ci+4:cf+4, ci+4:cf+4) = Js(ci+4:cf+4, ci+4:cf+4) - dJ(4:6, 4:6);
    Js(ci+4:cf+4, ci:cf) = Js(ci+4:cf+4, ci:cf)  - dJ(4:6, 1:3);
    Js(ci:cf, ci+4:cf+4) = Js(ci:cf, ci+4:cf+4)  - dJ(1:3, 4:6);    
end

end
