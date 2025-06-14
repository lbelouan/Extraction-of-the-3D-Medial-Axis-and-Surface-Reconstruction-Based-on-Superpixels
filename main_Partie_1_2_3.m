clear;
close all;

% Nombre d'images utilisées
nb_images = 36;

% Chargement des images
for i = 1:nb_images
    if i < 11
        nom = sprintf('images/viff.00%d.ppm', i-1);
    else
        nom = sprintf('images/viff.0%d.ppm', i-1);
    end
    im(:,:,:,i) = im2double(imread(nom)); % Conversion en double directement
end

% Affichage de l'image originale
figure;

subplot(1,4,1); imshow(im(:,:,:,1)); title('Image 1 originale');

% Sélection de l'image à traiter
image = im2double(im(:,:,:,1));
[rows, cols, ~] = size(image);

numSuperpixels = 200; % Nombre de superpixels souhaités

% Passage de RGB en LAB
lab_image = rgb2lab(image);

% Extraction des valeurs Lab aux positions des centres (interpolation bilinéaire)
L_channel = lab_image(:,:,1);
a_channel = lab_image(:,:,2);
b_channel = lab_image(:,:,3);

% Initialisation des centres des superpixels
S = round(sqrt((rows * cols) / numSuperpixels)); % Espacement entre centres

[x_grid, y_grid] = meshgrid(S/2:S:cols, S/2:S:rows);

centers = [x_grid(:), y_grid(:)];
centers(:,3) = interp2(L_channel, centers(:,1), centers(:,2), 'linear', 0);
centers(:,4) = interp2(a_channel, centers(:,1), centers(:,2), 'linear', 0);
centers(:,5) = interp2(b_channel, centers(:,1), centers(:,2), 'linear', 0);

numClusters = size(centers, 1); % Nombre de centres 


% Méthode du déplacement des centres vers les faibles gradients
centers = Faibles_gradients(L_channel, centers, numClusters, rows, cols);


% Algorithme SLIC 
max_iterations = 150; % Nombre d'itération de SLIC
m = 12; % coefficent de pondération distance couleur vs spatiale

[labels, centers] = Algorithme_SLIC(L_channel, a_channel, b_channel, centers, S, numClusters, rows, cols, max_iterations, m);


% Optimisation de la segmentation avec la Connexité des Superpixels
seuil_connexite = 0.5; % coefficient pour déterminer la taille minimale de Superpixel acceptable
nb_labels_avant = numel(unique(labels));
labels = Optimisation_connexe(labels, rows, cols, numSuperpixels, seuil_connexite);
nb_labels_apres = numel(unique(labels));
nb_labels_fusionne = nb_labels_avant-nb_labels_apres;

%  Recuperation bordures des superpixels
boundary = boundarymask(labels);

% Affichage des megapixels optimisés et de leurs frontières
subplot(1,4,2);
imshow(imoverlay(image, boundary,'r'));
title('Apperçu Mégapixels + frontières après optimisation connexe');


% Segmentation sur la couleur
masque_binaire_couleur = Segmentation_couleur(a_channel, b_channel, L_channel);

% Affichage de l'image segmentée sur la couleur 
subplot(1,4,3);
imshow(masque_binaire_couleur);
title('Image avec Segmentation sur la couleur');


% Segmentation sur la compacité
masque_binaire_compacite = Segmentation_compacite(labels, boundary); % Définir judicieusement seuil_compacité en observant compacity

% Affichage de l'image segmentée sur la compacité 
subplot(1,4,4);
imshow(masque_binaire_compacite);
title('Image avec Segmentation sur la compacité');

% Affichage légende
txt = sprintf(['Paramètres utilisés :\n\n', ...
               '• numSuperpixels : %d\n', ...
               '• S = %.1f\n', ...
               '• m : %d\n', ...
               '• max iterations : %d\n', ...
               '• seuil connexité : %.2f\n'], ...
               numSuperpixels, S, m, max_iterations, seuil_connexite );

annotation('textbox', [0.25, 0.12, 0.5, 0.1], ...
    'String', txt, ...
    'FitBoxToText', 'on', ...
    'EdgeColor', 'none', ...
    'FontSize', 10, ...
    'HorizontalAlignment', 'center');


% Récupération de la frontière pour les deux segmentations 
BW1 = logical(masque_binaire_couleur); % par mesure de sécurité on s'assure que l'image est logique
BW2 = logical(masque_binaire_compacite);

[startPointx1,startPointy1] = find(BW1, 1); % on trouve le premier pixel blanc de l’image (en parcourant de haut en bas et de gauche à droite).
[startPointx2,startPointy2] = find(BW2, 1);

% direction = find_initial_direction(BW, startPoint);
direction1='w'; % cette direction est la plus adaptée pour l'image 1
direction2='w';

contour1=bwtraceboundary(BW1, [startPointx1,startPointy1],direction1);
contour2=bwtraceboundary(BW2, [startPointx2,startPointy2],direction2);

% Affichage de la frontière obtenue pour la segmentation sur la couleur 
figure;
subplot(2,4,1);
imshow(BW1); % Affiche l’image binaire BW1
hold on
plot(contour1(:,2), contour1(:,1), 'r','LineWidth', 2);% Trace le contour en rouge
hold off
title('Contour détecté pour Seg Couleur ');


% Affichage de la frontière obtenue pour la segmentation sur la compacité 
subplot(2,4,5);
imshow(BW2); % Affiche l’image binaire BW1
hold on
plot(contour2(:,2), contour2(:,1), 'g','LineWidth', 2);% Trace le contour en vert
hold off
title('Contour détecté pour Seg Compacité ');


% Récupération de l'axe médian avec la méthode de Delaunay
pas_delauney = 10; % pas d'échantillonage du contour pour éviter trop de centres sur les contours

contour_echantillone1=contour1(1 :pas_delauney : end,: ); % contour échantilloné pour segmentation couleur
contour_echantillone2=contour2(1 :pas_delauney : end,: ); % contour échantilloné pour segmentation compacité

delaunay_var1 = delaunayTriangulation(contour_echantillone1(:,1),contour_echantillone1(:,2));
delaunay_var2 = delaunayTriangulation(contour_echantillone2(:,1),contour_echantillone2(:,2));

[centre_final1,rayon1]=circumcenter(delaunay_var1);
[centre_final2,rayon2]=circumcenter(delaunay_var2);

centre_final1 = round(centre_final1);
centre_final2 = round(centre_final2);

% Affichages des centres de delaunay pour Segmentation couleur
subplot(2,4,2);
imshow(BW1); % Affiche l’image binaire BW
hold on; % Permet de superposer le contour
plot(contour1(:,2), contour1(:,1), 'r','LineWidth', 2); % Trace le contour en rouge
plot(centre_final1(:,2),centre_final1(:,1), 'b.','MarkerSize',5,'LineStyle',"none") % Trace les centres en bleu
hold off
title('Contour détecté + centres non filtré pour Seg Couleur');

% Affichages des centres de delaunay pour Segmentation compacité
subplot(2,4,6);
imshow(BW2); % Affiche l’image binaire BW
hold on; % Permet de superposer le contour
plot(contour2(:,2), contour2(:,1), 'g','LineWidth', 2); % Trace le contour en vert
plot(centre_final2(:,2),centre_final2(:,1), 'b.','MarkerSize',5,'LineStyle',"none") % Trace les centres en bleu
hold off
title('Contour détecté + centres non filtré pour Seg compacité');


% Filtrage des centre en dehors du contour 
[centres_filtres1, rayon_filtre1] = Filtrage_centres(BW1, centre_final1, rayon1);
[centres_filtres2, rayon_filtre2] = Filtrage_centres(BW2, centre_final2, rayon2);
centres_xy1 = centres_filtres1(:, [2 1]);  % [colonne ligne] = [x y]
centres_xy2 = centres_filtres2(:, [2 1]);  % [colonne ligne] = [x y]

% Affichages des centres delaunay filtrés pour Segmentation Couleur 
subplot(2,4,3);
imshow(BW1); % Affiche l’image binaire BW
hold on; % Permet de superposer le contour
plot(contour1(:,2), contour1(:,1), 'r','LineWidth', 2);% Trace le contour en rouge
plot(centres_filtres1(:,2),centres_filtres1(:,1), 'b.','MarkerSize',5,'LineStyle',"none")
viscircles(centres_xy1, rayon_filtre1, 'Color', 'y');
hold off
title('Contour détecté + centres filtré pour Seg couleur');

% Affichages des centres delaunay filtrés pour Segmentation Compacité 
subplot(2,4,7);
imshow(BW2); % Affiche l’image binaire BW
hold on; % Permet de superposer le contour
plot(contour2(:,2), contour2(:,1), 'g','LineWidth', 2);% Trace le contour en rouge
plot(centres_filtres2(:,2),centres_filtres2(:,1), 'b.','MarkerSize',5,'LineStyle',"none")
viscircles(centres_xy2, rayon_filtre2, 'Color', 'y');
hold off
title('Contour détecté + centres filtré pour Seg compacité');

% Tracé des axes médians

[Gg1, XY1] = Trace_Axe_Median(delaunay_var1, centre_final1, centres_filtres1);
[Gg2, XY2] = Trace_Axe_Median(delaunay_var2, centre_final2, centres_filtres2);

%Affichage Axe médian pour la Segmentation Couleur
subplot(2,4,4);
imshow(BW1); 
hold on;
gplot(Gg1, XY1, '-r');
%plot(XY1(:,1), XY1(:,2), 'b');
title('Axe médian pour Seg Couleur');

% Axe médian pour la Segmentation Compacité
subplot(2,4,8);
imshow(BW2); 
hold on;
gplot(Gg2, XY2, '-g');
%plot(XY2(:,1), XY2(:,2), 'b');
title('Axe médian pour Seg Compacité');



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  FONCTIONS SEGMENTATION BINAIRE AVEC CRITERE DE COMPACITE%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

% Fonction pour calculer l'aire des superpixels
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


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % FONCTIONS ESTIMATION DE LA FRONTIERE %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function best_direction = find_initial_direction(BW, startPoint)
%     [rows, cols] = size(BW);
% 
%     % Définition des directions (Nord, Sud, Est, Ouest)
%     directions = {'N', 'S', 'E', 'W'};
%     offsets = [-1 0; 1 0; 0 1; 0 -1]; % Déplacements associés
% 
%     best_direction = ''; % Par défaut, aucune direction trouvée
% 
%     % Vérifier les 4 voisins
%     for d = 1:4
%         next_point = startPoint + offsets(d, :);
% 
%         % Vérifier si le point est dans l'image et fait partie du contour
%         if is_inside_image(next_point, rows, cols) && BW(next_point(1), next_point(2)) == 1
%             best_direction = directions{d}; % On choisit cette direction
%             break; % On s'arrête à la première direction trouvée
%         end
%     end
% end
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  FONCTIONS UTILES (Alternative sans Image Processing Toolbox)%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% % Fonction de conversion RGB -> Lab (remplace `rgb2lab()`)**
% function lab = rgb_to_lab(rgb)
%     % Normaliser entre 0 et 1
%     rgb = double(rgb);
% 
%     % Transformation RGB -> XYZ
%     M = [0.4124564, 0.3575761, 0.1804375;
%          0.2126729, 0.7151522, 0.0721750;
%          0.0193339, 0.1191920, 0.9503041];
% 
%     xyz = reshape(rgb, [], 3) * M';
%     xyz = reshape(xyz, size(rgb));
% 
%     % Transformation XYZ -> Lab
%     X = xyz(:,:,1) / 0.95047;
%     Y = xyz(:,:,2);
%     Z = xyz(:,:,3) / 1.08883;
% 
%     F = @(t) ((t > 0.008856) .* (t .^ (1/3)) + (t <= 0.008856) .* (7.787 * t + 16/116));
% 
%     L = 116 * F(Y) - 16;
%     a = 500 * (F(X) - F(Y));
%     b = 200 * (F(Y) - F(Z));
% 
%     lab = cat(3, L, a, b);
% end
