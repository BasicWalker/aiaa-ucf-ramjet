% Calculate Total Number of Nodes Created by the Number of Characteristics
function nodes = nodeCalculator(n)
    point(1) = 2;                       % Amount of nodes with 1 char line
    for i = 1:n-1
        point(i+1) = point(i) + (i+2);  %#ok<AGROW>
    end
    nodes = point(end);
end

