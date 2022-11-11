function [K_global, K_local_out] = StiffnessMatrix(L_Elements,Nodal_Positions,D,Ngp)
n_nodes = length(L_Elements) + 1; % # of nodes
n = length(L_Elements);           % # of elements
% Gauss Integration
[GP, WT] = GaussInt(Ngp);
% Initialize global stiffness matrix and force vector
K_global = zeros(n_nodes,n_nodes);
% Create an array storing start node ID and end node ID for each element.
elements = zeros(n,2);
for iE=1:n
    st_pt=Nodal_Positions(iE,1);
    end_pt=Nodal_Positions(iE+1,1);
    elements(iE,:)=[st_pt end_pt];
end
% Assemble the global stiffness matrix and force vector
K_local = zeros(2,2);
for iK = 1:n
    st_pt = elements(iK,1);   % coordinate of starting pt of element.
    end_pt = elements(iK,2);  % coordinate of end pt.
    for i = 1:Ngp    % Loop over all Gauss points. 
        xi = GP(i);  % xi coordinate of the Gauss point number i
        X = elements(iK,:)'; % Nodal x-coordinates for element iK
        [~,B,J] = NBJvec(1,xi,X); % Evaluate B and J for element iK at Gauss point xi
        r_xi = 0.5*(st_pt + end_pt) + xi/2*(end_pt - st_pt);
        % Compute the element stiffness matrix
        K_local(:,:,i) = 4*pi* WT(i) * B'*D*B * r_xi^2 * J;
    end 
    K_local = sum(K_local,3);
    K_local_out(:,:,iK) = K_local(:,:);
    v = [iK,iK+1];
    K_global(v,v) = K_global(v,v) + K_local;
end
end