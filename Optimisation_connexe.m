
% OPTIMISATION_CONNEXE Fusionne les petits superpixels avec leurs voisins
%
%    Entrées :
%   - labels : matrice de segmentation 
%   - rows, cols : dimensions de l'image
%   - numSuperpixels : nombre total de superpixels
%   - seuil_connexite : pour déterminer la taille minimale acceptable
%
%   Sortie :
%   - labels : labels mis à jour après fusion des petits superpixels

function labels = Optimisation_connexe(labels, rows, cols, numSuperpixels, seuil_connexite)

    Np = rows * cols;
    taille_moyenne = Np / numSuperpixels;
    seuil_taille = seuil_connexite * taille_moyenne;

    % Taille de chaque superpixel
    superpixel_sizes = histcounts(labels(:), 0:numSuperpixels);

    % Identifier les petits superpixels
    petits_superpixels = find(superpixel_sizes < seuil_taille);

    % Fusion de chaque petit superpixel avec un voisin
    for i = 1:length(petits_superpixels)
        label_courant = petits_superpixels(i);

        % Coordonnées des pixels du superpixel
        [r, c] = find(labels == label_courant);

        % Labels voisins : haut, bas, gauche, droite
        voisins = unique([
            labels(sub2ind([rows, cols], max(r-1, 1), c));      % Haut
            labels(sub2ind([rows, cols], min(r+1, rows), c));   % Bas
            labels(sub2ind([rows, cols], r, max(c-1, 1)));      % Gauche
            labels(sub2ind([rows, cols], r, min(c+1, cols)))    % Droite
        ]);

        % Retirer le superpixel courant et nettoyer
        voisins(voisins == label_courant) = [];
        voisins = voisins(voisins > 0 & voisins <= length(superpixel_sizes));

        % Choisir le voisin le plus gros
        if ~isempty(voisins)
            [~, idx_max] = max(superpixel_sizes(voisins));
            label_fusion = voisins(idx_max);

            % Mise à jour des labels
            labels(labels == label_courant) = label_fusion;

            % Mise à jour des tailles
            superpixel_sizes(label_fusion) = superpixel_sizes(label_fusion) + superpixel_sizes(label_courant);
        end
    end
end
