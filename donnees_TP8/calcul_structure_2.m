function [u_barre_c_r] = calcul_structure_2(u_barre_c,b_c,Dx,Dy,lambda,epsilon)
% approximation de Erof pour pouvoir dériver

%number of pixels
[nb_l, nb_c] = size(u_barre_c);
N = nb_l * nb_c;

u_k = u_barre_c(:);

dxu_k = Dx * u_k;
dyu_k = Dy * u_k;

W = spdiags( 1 ./ sqrt( dxu_k .^2 + dyu_k .^2 + epsilon) , 0, N, N);

A = speye(N) - lambda .* ( - transpose(Dx) * W * Dx - transpose(Dy) * W * Dy );

u_barre_c_r = A \ b_c;
u_barre_c_r = reshape(u_barre_c_r, nb_l, nb_c);

end

