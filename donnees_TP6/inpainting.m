function [u_kp1] = inpainting(b,u_k,lambda,Dx,Dy,epsilon, D)

N = size(u_k,1);

dxu_k = Dx * u_k;
dyu_k = Dy * u_k;

W = spdiags( 1 ./ sqrt( dxu_k .^2 + dyu_k .^2 + epsilon) , 0, N, N);

W_prive = spdiags(D < 1, 0, N, N);
A = W_prive - lambda .* ( - transpose(Dx) * W * Dx - transpose(Dy) * W * Dy );

u_kp1 = A \ (W_prive * b);




end

