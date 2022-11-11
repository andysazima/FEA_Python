% Pre-processing inputs to be used in FEMSolve

function [L_ele] = Preprocessing(Dimension, Total_Length, N_Elements)

switch Dimension
    case 1
    L_ele = Total_Length / N_Elements.*ones(1,N_Elements);
end
end