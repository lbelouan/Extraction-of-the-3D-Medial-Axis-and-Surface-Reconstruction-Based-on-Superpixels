
% Trace_Axe_Median Relie les centres des triangles qui partages un cpoté en commun
%
%   Entrées :
%   - delaunay_v : objet delaunayTri contenant la triangulation
%   - centre_final : coordonnées [y x] des centres avant filtrage
%   - centre_final_final : coordonnées [y x] des centres valides
%
%   Sorties :
%   - Gg : graphe sparse N×N où les arêtes relient des centres voisins
%   - XY : coordonnées [x y] des sommets à utiliser avec gplot ou autres

function [Gg, XY] = Trace_Axe_Median(delaunay_v, centre_final, centre_final_final)

    triangle_neighbors = neighbors(delaunay_v);  % voisins de triangles dans la triangulation

    [~, valid_idx] = ismember(round(centre_final), centre_final_final, 'rows'); % valid_idx(i) = position de centre_final(i,:) dans centre_final_final (0 si absent)
    valid_mask = valid_idx > 0;  % Crée un masque logique : true si le centre est trouvé dans centre_final_final
    map_tri_to_new = zeros(size(valid_mask)); % Initialise un vecteur de correspondance pour tous les triangles
    map_tri_to_new(valid_mask) = 1:nnz(valid_mask); % Pour chaque triangle valide, on lui attribue un nouvel index consécutif (1,2,...)

    edges = []; % Initialisation de la liste des arêtes du graphe

    for i = 1:size(triangle_neighbors,1) % Pour chaque triangle de la triangulation
        if ~valid_mask(i)
            continue; % Si le triangle i n'est pas associé à un centre valide, on le saute
        end
        for j = 1:3 % Pour les 3 voisins du triangle i
            voisin = triangle_neighbors(i,j);
            if voisin > 0 && valid_mask(voisin)
                a = map_tri_to_new(i);  % Index du centre du triangle i dans les centres valides
                b = map_tri_to_new(voisin); % Index du centre du triangle voisin
                if a < b  % éviter les doublons
                    edges = [edges; a b]; % Ajoute une arête entre a et b
                end
            end
        end
    end

    Gg = sparse(edges(:,1), edges(:,2), 1, nnz(valid_mask), nnz(valid_mask)); % Chaque arête a un poids de 1. Taille : Nvalide x Nvalide
    XY = [centre_final_final(:,2), centre_final_final(:,1)];  % Coordonnées XY des sommets du graphe, pour affichage (x = colonne, y = ligne)
end
