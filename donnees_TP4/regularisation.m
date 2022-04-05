function [r] = regularisation(i,j,k,k_ij)

nb_lignes = size(k, 1);
nb_colones = size(k, 2);

x1 = max(1, i -1);
x2 = min(nb_lignes, i +1);

y1 = max(1, j -1);
y2 = min(nb_colones, j +1);

voisins = k(x1:x2, y1:y2);

r = sum(sum((k_ij ~= voisins)));

end

