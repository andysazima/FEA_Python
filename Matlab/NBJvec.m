% Finds the shape functions, their derivative, and the jacobian for each
% element

function [N,B,J] = NBJvec(Ord,xi,x)

% inputs:   Ord  order (n=1 for linear and n=2 for quadratic elements)
%           xi   Gauss integration point in xi
%           x    array of x-coordinates of the element

% outputs:  N    shape function vector for the element
%           B    B-matrix containing the derivatives of shape functions 
%           J    Jacobian matrix. for 1D elements, it is [dx/dxi]

    switch Ord
        case 1 % For linear elements
            r_xi = (x(2) + x(1))/2 + (x(2) - x(1))/2*xi;
            
            N = (0.5).*[ (1-xi)   (1+xi) ];
            B = [ -1/(x(2)-x(1))  1/(x(2)-x(1));
                   1/2/r_xi*(1-xi)    1/2/r_xi*(1+xi);
                   1/2/r_xi*(1-xi)    1/2/r_xi*(1+xi);];
            J = (0.5)*( x(2) - x(1) );
            
%         case 2 % For quadratic elements
%             N = [ 0.5*(xi-1)*xi    (1-xi)*(1+xi)    0.5*(xi+1)*xi ];
%             dxdxi = (0.5)*(x(3)-x(1));
%             dNdxi = [ xi-0.5    -2*xi    xi+0.5 ];
%             J = dxdxi;
%             B = dNdxi;

        otherwise
            disp('Not Available')
    end
end