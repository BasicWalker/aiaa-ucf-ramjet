function [mach, machAngle] = inversePrandtlMeyer (g, nu)
    [mach, ~, machAngle] = flowprandtlmeyer(g, nu, 'nu');
end

    