
% SEGMENTATION_COMPACITE Réalise une segmentation par compacité des superpixels
%
%   Entrées :
%   - labels : matrice des labels des superpixels
%   - boundary : image binaire des contours des superpixels
%   - seuil_compacite : seuil numérique sur la compacité (ex : 200)
%
%   Sortie :
%   - masque_binaire_compacite : image binaire des régions jugées suffisamment compactes

function masque_binaire_compacite = Segmentation_compacite(labels, boundary)

    % Calcul de la surface (aire) de chaque superpixel
    superpixel_areas = cell2mat(struct2cell(regionprops(labels, 'Area')));

    % Calcul du périmètre à partir des pixels de contour
    perimeters = calcul_aire(labels(boundary == 1));  % À adapter si nécessaire

    % Calcul de la compacité : (périmètre²) / aire
    compacity = (perimeters'.^2) ./ superpixel_areas;
    compacity(isnan(compacity)) = 0; % Remplace tous les NaN par 0 pour eviter des erreurs
    
     % Seuillage capacité
     % seuil_compacite = median(compacity); % On prend la médiane comme seuil
     seuil_compacite = 215; % dans la pratique la médiane n'est pas optimale donc on choixit judicieusement en observant compacity

    % Génère un masque logique : 1 si la compacité est ≥ seuil
    numLabels = max(labels(:));
    superpixel_mask = ones(numLabels, 1);
    superpixel_mask(compacity < seuil_compacite) = 0;

    % Création du masque binaire global
    masque_binaire_compacite = superpixel_mask(labels);
end



function result = calcul_aire(M_labels)
    Nombre_labels = max(M_labels(:));
    [rows_label, cols_label] = size(M_labels);
    result = zeros(Nombre_labels, 1);

    for i = 1:rows_label
        for j = 1 : cols_label
             result(M_labels(i,j)) = result(M_labels(i,j)) + 1;
        end
    end
end