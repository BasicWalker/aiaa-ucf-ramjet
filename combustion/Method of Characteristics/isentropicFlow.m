function areaRatio = isentropicFlow (g, M)
    [~, ~, ~, ~, A] = flowisentropic(g, M);
    areaRatio = A;
end

    