function T_AFT = AFT(mol_mat, hfMat,CpMat, T_mat, T_STP) 

    T_AFT = (sum(mol_mat.*hfMat) + sum(mol_mat.*CpMat.*(T_mat-T_STP)))...
        / abs(sum(mol_mat(3:end).*CpMat(3:end)));

end
