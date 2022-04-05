clear;
close all;
taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);

% Paramètres :
N = 0;					% Nombre de disques d'une configuration
R = 10;					% Rayon des disques
nb_points_affichage_disque = 30;
increment_angulaire = 2*pi/nb_points_affichage_disque;
theta = 0:increment_angulaire:2*pi;
rose = [253 108 158]/255;
q_max = 300;
nb_affichages = 1000;
pas_entre_affichages = floor(q_max/nb_affichages);
temps_pause = 0.0005;

beta = 1;
S = 130;
gamma = 5;
T0 = 0.1;
lambda0 = 100;
alpha = 0.99;

% Lecture et affichage de l'image :
I = imread('colonie.png');
I = rgb2gray(I);
I = double(I);
I = I(1:400,100:450);
[nb_lignes,nb_colonnes] = size(I);
figure('Name',['Detection de ' num2str(N) ' flamants roses'],'Position',[0.25*L,0,0.75*L,0.5*H]);




q = 1;
T = T0;
lambda = lambda0;
U = zeros(N,1);

%%initialisation
%tirage aléatoire du nombre de naissance
    n_new = poissrnd(lambda);
    %disp(n_new);
    
    %ajout des n_new cercles
   
    c = zeros(n_new,2);
    I_moyen =zeros(1, n_new);
    for naissance = 1:n_new
        c_n = [nb_colonnes*rand nb_lignes*rand];
        c(naissance,:) = c_n;
	    I_moyen(naissance) = calcul_I_moyen(I,c_n,R);
    end
    
    %maj du nombre de cercles
    N = N + n_new;
flag = true;

liste_q = 0;

 %tri des disques
    % lors des recouvrement : on va toujours supprimer en premier le moins
    % intéresssant
    Ui =  ones(1,N) - 2 ./ ( ones(1,N) + exp( - gamma * ( I_moyen / S - ones(1,N)))) ;
    size(Ui)
    [Ui_sorted, index] = sort(Ui,'descend');

    %calcul de U
    rayon1_x = triu(repmat(c(:,1), 1, N), 1);
    rayon1_y = triu(repmat(c(:,2), 1, N), 1);
    rayon2_x = triu(repmat(transpose(c(:,1)), N, 1), 1);
    rayon2_y = triu(repmat(transpose(c(:,2)), N, 1), 1);
    rayons = sqrt((rayon1_x - rayon2_x) .^ 2 + (rayon1_y - rayon2_y) .^ 2 );
    difference = sum(sum(triu(rayons < sqrt(2) * R , 1)));
    somme = sum(Ui) + beta * difference;

I_moyen_config = mean(somme);
liste_I_moyen_config = [I_moyen_config];

% Affichage de la configuration initiale :
subplot(1,2,1);
imagesc(I);
axis image;
axis off;
colormap gray;
hold on;
for i = 1:N
	x_affich = c(i,1)+R*cos(theta);
	y_affich = c(i,2)+R*sin(theta);
	indices = find(x_affich>0 & x_affich<nb_colonnes & y_affich>0 & y_affich<nb_lignes);
	plot(x_affich(indices),y_affich(indices),'Color',rose,'LineWidth',3);
end
pause(temps_pause);


while flag && q <= q_max


    %tirage aléatoire du nombre de naissance
    n_new = poissrnd(lambda);
    %disp(n_new);
    
    %ajout des n_new cercles
    
    c = [c ;zeros(n_new,2)];
    
    for naissance = 1:n_new
        c_n = [nb_colonnes*rand nb_lignes*rand];
        c(N + naissance,:) = c_n;
	    I_moyen(N + naissance) = calcul_I_moyen(I,c_n,R);
    end
    
    %maj du nombre de cercles
    N = N + n_new;

    %tri des disques
    % lors des recouvrement : on va toujours supprimer en premier le moins
    % intéresssant
    Ui =  ones(1,N) - 2 ./ ( ones(1,N) + exp( - gamma * ( I_moyen / S - ones(1,N)))) ;
    
    [Ui_sorted, index] = sort(Ui,'descend');
    

%     % Affichage de la configuration initiale :
%      %hold off;
%     subplot(1,2,1);
%     imagesc(I);
%     axis image;
%     axis off;
%     colormap gray;
%     hold on;
%     for i = 1:N
% 	    x_affich = c(i,1)+R*cos(theta);
% 	    y_affich = c(i,2)+R*sin(theta);
% 	    indices = find(x_affich>0 & x_affich<nb_colonnes & y_affich>0 & y_affich<nb_lignes);
% 	    plot(x_affich(indices),y_affich(indices),'Color',rose,'LineWidth',3);
%     end
%     pause(temps_pause);

    %morts
    
    %calcul de U
    rayon1_x = triu(repmat(c(:,1), 1, N), 1);
    rayon1_y = triu(repmat(c(:,2), 1, N), 1);
    rayon2_x = triu(repmat(transpose(c(:,1)), N, 1), 1);
    rayon2_y = triu(repmat(transpose(c(:,2)), N, 1), 1);
    rayons = sqrt((rayon1_x - rayon2_x) .^ 2 + (rayon1_y - rayon2_y) .^ 2 );
    difference = sum(sum(triu(rayons < sqrt(2) * R , 1)));
    somme = sum(Ui) + beta * difference;

    flag = false;
    index2= index;
    N2 = N;
    deleted_index = [];
    for i=1:N2
        U_prive_i = Ui;
        deleted_index_i = [index(i) deleted_index];
        U_prive_i(deleted_index_i) = [];
        c_prive_i = c;
        c_prive_i(deleted_index_i, :) = [];
        rayon1_x_prive_i = triu(repmat(c_prive_i(:,1), 1, N-1), 1);
        rayon1_y_prive_i = triu(repmat(c_prive_i(:,2), 1, N-1), 1);
        rayon2_x_prive_i = triu(repmat(transpose(c_prive_i(:,1)), N-1, 1), 1);
        rayon2_y_prive_i = triu(repmat(transpose(c_prive_i(:,2)), N-1, 1), 1);
        rayons_prive_i = sqrt((rayon1_x_prive_i - rayon2_x_prive_i) .^ 2 + (rayon1_y_prive_i - rayon2_y_prive_i) .^ 2 );
        difference_prive_i = sum(sum(triu(rayons_prive_i < sqrt(2) * R , 1)));
        somme_prive_i = sum(U_prive_i) + beta * difference_prive_i;
        p_i = lambda /  (lambda + exp((somme_prive_i - somme) / T));
        p = rand(1);
        if p_i > p
            
            deleted_index = [deleted_index index(i)];
            N = N-1;
            flag = true;

            cc = c;
            cc(deleted_index,:) = [];


            rayon1_x = triu(repmat(cc(:,1), 1, N), 1);
            rayon1_y = triu(repmat(cc(:,2), 1, N), 1);
            rayon2_x = triu(repmat(transpose(cc(:,1)), N, 1), 1);
            rayon2_y = triu(repmat(transpose(cc(:,2)), N, 1), 1);
            rayons = sqrt((rayon1_x - rayon2_x) .^ 2 + (rayon1_y - rayon2_y) .^ 2 );
            difference = sum(sum(triu(rayons < sqrt(2) * R , 1)));
            Ui_delete = Ui;
            Ui_delete(deleted_index) = [];
            somme = sum(Ui_delete) + beta * difference;
        end
    end

    T = alpha * T;
    lambda = alpha * lambda;
    

    c(deleted_index,:) = [];
    I_moyen(deleted_index) = [];








    
    

    %display
    hold off;
		    subplot(1,2,1);
		    imagesc(I);
		    axis image;
		    axis off;
		    colormap gray;
		    hold on;
		    for j = 1:N
			    x_affich = c(j,1)+R*cos(theta);
			    y_affich = c(j,2)+R*sin(theta);
			    indices = find(x_affich>0 & x_affich<nb_colonnes & y_affich>0 & y_affich<nb_lignes);
			    plot(x_affich(indices),y_affich(indices),'Color',rose,'LineWidth',3);
		    end
		    %pause(temps_pause);

    




    
	    % Courbe d'évolution du niveau de gris moyen :
		    liste_q = [liste_q q];
		    I_moyen_config = somme;
		    liste_I_moyen_config = [liste_I_moyen_config I_moyen_config];
		    subplot(1,2,2);
		    plot(liste_q,liste_I_moyen_config,'.-','Color',rose,'LineWidth',3);
		    axis([0 q_max -200 255]);
		    set(gca,'FontSize',20);
		    xlabel('Nombre d''iterations','FontSize',30);
		    ylabel('Niveau de gris moyen','FontSize',30);

        
        q = q+1;
        
            
    %end

    


	
end
