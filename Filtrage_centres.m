
% FILTRAGE_CENTRES Garde les centres situés à l'intérieur du masque BW
%
%   Entrées :
%   - BW : image binaire (mask) de taille H × W
%   - centre_final : matrice N × D contenant les informations associées aux centres
%   - rayon: 
%
%   Sortie :
%   - centres_filtres : sous-ensemble de centre_final dont les coordonnées sont dans BW == 1
%   - rayons_filtres :

function [centres_filtres, rayons_filtres] = Filtrage_centres(BW, centre_final, rayon)

    [H, W] = size(BW);
    centres_filtres = [];  % initialisation
    rayons_filtres = [];

    for i = 1:length(centre_final(:,1))
        x = round(centre_final(i,1));
        y = round(centre_final(i,2));

        % Vérifie que le point est dans les dimensions de l'image
        if x >= 1 && x <= H && y >= 1 && y <= W
            if BW(x, y) == 1
                centres_filtres = [centres_filtres; centre_final(i,:)];  % ajout de la ligne complète
                rayons_filtres = [rayons_filtres; rayon(i)];
            end
        end
    end
end
