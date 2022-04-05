function [moyennes,variances] = estimation(I,N,couleur_classes)

fprintf('Selectionnez %d echantillons\n',N);
[nb_lignes,nb_colonnes] = size(I);
moyennes = zeros(1,N);
variances = zeros(1,N);
coordinates = [15 45 24 47 63 75 81 91; 17 36 45 68 32 67 81 91];
for k = 1:N
% 	[x1,y1] = ginput(1);
% 	while (x1<1)||(x1>nb_colonnes)||(y1<1)||(y1>nb_lignes)
% 		[x1,y1] = ginput(1);
% 	end
% 	[x2,y2] = ginput(1);
% 	while (x2<1)||(x2>nb_colonnes)||(y2<1)||(y2>nb_lignes)||(x2==x1)||(y2==y1)
% 		[x2,y2] = ginput(1);
% 	end
    x1 = coordinates(1,2 * k -1);
    x2 = coordinates(1,2 * k);
    y1 = coordinates(2,2 * k -1);
    y2 = coordinates(2,2 * k);
	line([x1,x1],[y1,y2],'Color',couleur_classes(k,:),'Linewidth',2);
	line([x1,x2],[y2,y2],'Color',couleur_classes(k,:),'Linewidth',2);
	line([x2,x2],[y2,y1],'Color',couleur_classes(k,:),'Linewidth',2);
	line([x2,x1],[y1,y1],'Color',couleur_classes(k,:),'Linewidth',2);
	drawnow;

	echantillons = [];
	for i = floor(min([y1,y2])):ceil(max([y1,y2]))
		for j = floor(min([x1,x2])):ceil(max([x1,x2]))
			echantillons = [ echantillons I(i,j) ];
		end
	end
	moyenne = mean(echantillons);
	moyennes(k) = moyenne;
	nb_donnees = length(echantillons);
	echantillons_centres = echantillons-moyenne;
	variances(k) = echantillons_centres*transpose(echantillons_centres)/nb_donnees;
end
