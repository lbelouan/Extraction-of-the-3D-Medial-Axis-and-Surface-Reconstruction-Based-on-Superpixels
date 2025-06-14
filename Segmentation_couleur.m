
% SEGMENTATION_COULEUR Effectue une segmentation binaire basée sur a* et b*
%
%   masque_binaire_couleur = Segmentation_couleur(a_channel, b_channel)
%
%   Entrées :
%   - a_channel : canal a* de l'image Lab
%   - b_channel : canal b* de l'image Lab
%
%   Sortie :
%   - masque_binaire_couleur : matrice logique de même taille, 1 si le pixel appartient à l'objet

function masque_binaire_couleur = Segmentation_couleur(a_channel, b_channel, L_channel)

    % Définir les seuils comme moyennes globales
    seuil_a = mean(a_channel(:));
    seuil_b = mean(b_channel(:));
    seuil_L = 1.5*mean(L_channel(:)); 
   

    % Créer le masque binaire : pixels au-dessus des deux seuils
    masque_binaire_couleur =   ((b_channel > seuil_b) & (a_channel > seuil_a)) |  ((L_channel > seuil_L) & (b_channel > seuil_b));

end
