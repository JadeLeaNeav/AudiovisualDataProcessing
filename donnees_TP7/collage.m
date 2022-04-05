function u = collage(r,s,interieur)

    r_double = double(r);
    s_double = double(s);
    
    [nb_lignes, nb_colonnes, nb_canaux] = size(r_double);
    
    % Vectorisation de u0 :
    nb_pixels = nb_lignes * nb_colonnes;
    u0 = reshape(r_double,[nb_pixels nb_canaux]); 
    
    
    s0 = reshape(s_double,[nb_pixels nb_canaux]); 
    
    u = zeros(nb_pixels, nb_canaux);
    
    
    %%CALCUL DE A
    % Oprateur gradient :
    e = ones(nb_pixels,1);
    Dx = spdiags([-e e],[0 nb_lignes],nb_pixels,nb_pixels);
    Dx(end-nb_lignes+1:end,:) = 0;
    Dy = spdiags([-e e],[0 1],nb_pixels,nb_pixels);
    Dy(nb_lignes:nb_lignes:end,:) = 0;
    
    indices_bord_r = [2:(nb_lignes-1) 1:nb_lignes:(nb_pixels - nb_lignes +1) nb_lignes:nb_lignes:nb_pixels (nb_pixels - nb_lignes +2):(nb_pixels -1)];
    n_bord_r = size(indices_bord_r, 2);
    n_r = nb_pixels;
    
    % Matrice A du systme :
    A = -Dx'*Dx -Dy'*Dy;
    
    A(indices_bord_r,:) = sparse(1:n_bord_r,indices_bord_r,ones(n_bord_r,1),n_bord_r,n_r);
    
    
    
    for i = 1:nb_canaux
    % Second membre b du systme :
    
        %
        gx = - Dx * u0(:, i); 
        gy = - Dy * u0(:, i);
        grad_sx = -Dx * s0(:, i); 
        grad_sy = -Dy * s0(:, i); 
        gx(interieur) =  grad_sx(interieur);
        gy(interieur) = grad_sy(interieur);
        grad_g = -Dx * gx - Dy * gy;
        
    
    
    
        grad_g(indices_bord_r) = u0(indices_bord_r,i);
    
        % Rsolution du systme A*u = b :
        u(:,i) = A \ grad_g;
        
    end
    
    u = reshape(u, [nb_lignes nb_colonnes nb_canaux]);


end

