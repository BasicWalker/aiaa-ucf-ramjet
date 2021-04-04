function AFT = aft_equation(mol_r, mol_p, hfMat,CpMat, T_r, T_STP)

    AFT = (sum(mol_r.*hfMat) - sum(mol_p.*hfMat) + sum(mol_r.*CpMat.*(T_r-T_STP)) +...
    sum(mol_p.*CpMat.*T_STP)) / sum(mol_p.*CpMat);

end
