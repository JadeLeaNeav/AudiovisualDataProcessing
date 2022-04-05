function [U,k] = recuit_simule(U,k,AD,T,beta)

nb_lignes = size(k, 1);
nb_colones = size(k, 2);
nb_classes = size(AD,3);

new_k = randi(nb_classes - 1, nb_lignes, nb_colones);
new_k = mod(new_k + k -1, nb_classes) +1;



for i = 1:nb_lignes
    for j = 1:nb_colones

        U_ij = U(i, j);

        new_U = AD(i,j, new_k(i,j)) + beta * regularisation(i,j,k,new_k(i,j));

        p = exp( - (new_U - U_ij) / T );

        change = ( p > rand );

        if new_U < U_ij | change
            k(i, j) = new_k(i, j);
            U(i, j) = new_U;
        end

    end



end

