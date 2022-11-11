function [M_global, M_local_out] = MassMatrix(L_Elements,Nodal_Positions,density,Ngp,MassType)
n_nodes = length(L_Elements) + 1; % # of nodes
n = length(L_Elements);           % # of elements
rho = density;
% Gauss Integration
[GP, WT] = GaussInt(Ngp);
% Initialize global mass matrix and local mass matrix for output
M_global = zeros(n_nodes,n_nodes);
M_local_out = zeros(2,2,n);
% Create an array storing start node ID and end node ID for each element.
elements = zeros(n,2);
for iE=1:n
    st_pt=Nodal_Positions(iE,1);
    end_pt=Nodal_Positions(iE+1,1);
    elements(iE,:)=[st_pt end_pt];
end
% Assemble the global stiffness matrix and force vector
M_local = zeros(2,2);
for iK = 1:n
    st_pt = elements(iK,1);       % coordinate of starting pt of element.
    end_pt = elements(iK,2);      % coordinate of end pt.
    for i = 1:Ngp                 % Loop over all Gauss points. 
        xi = GP(i);               % xi coordinate of the Gauss point number i
        X = elements(iK,:)';      % Nodal x-coordinates for element iK
        [N,~,J] = NBJvec(1,xi,X); % Evaluate B and J for element iK at Gauss point xi
        r_xi = 0.5*(st_pt + end_pt) + xi/2*(end_pt - st_pt);
        % Compute the element stiffness matrix
        M_local(:,:,i) = 4 * pi* WT(i) * N'*rho*N * r_xi^2 * J; % Mass at the i-th Gauss Point
    end
    M_local = sum(M_local,3);
    switch MassType
        case 'consistent'
        case 'lumped'
            M_local = sum(M_local,2).*eye(size(M_local));
    end
    M_local_out(:,:,iK) = M_local(:,:); % Mass for each element (not each gauss point)
    v = [iK,iK+1];
    M_global(v,v) = M_global(v,v) + M_local;
end
end