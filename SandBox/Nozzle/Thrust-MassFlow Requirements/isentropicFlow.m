function [presRatio, areaRatio] = isentropicFlow (g, M)
    [~, ~, P, ~, A] = flowisentropic(g, M);
    areaRatio = A;
    presRatio = P; 
end

    