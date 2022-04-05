function [AD] = attache_donnees(I,moyennes,variances)

nb_classes = size(moyennes,2);
nb_lignes = size(I,1);
nb_colones = size(I,2);
AD = zeros(nb_lignes, nb_colones, nb_classes);

for classe = 1:nb_classes

    AD(:, :, classe) = 0.5 * ( log(variances(classe)) + (I - moyennes(classe)) .^2 / variances(classe) ); 

end

end