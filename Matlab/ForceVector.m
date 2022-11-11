function [F_global] = ForceVector(K_or_M,BCs)
% Inputs
%    K_or_M = Stiffness or Mass Matrix. dim: (n+1)x(n+1)
%    BCs    = matrix of Boundary Conditions. dim: (num BCs)x(2)
%             Left column is the DOF, right column is the BC

% Outputs
%    F_global = Global force vector

F_global = zeros(size(K_or_M,1),1);
F_global(BCs(:,1),1) = BCs(:,2);
end