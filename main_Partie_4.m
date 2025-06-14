clear;
close all;

s = load('mask.mat'); % On charge les masques fourni car pas assez de temps pour enregistrer les 36 masques perso

im_mask = s.im_mask; % Extraire le champ im_mask
im_mask = ~im_mask; % Les masques sont inversés

% Affichage du 1er masque pour vérifier la segmentation
figure;
imshow(im_mask(:,:,1)); 
title('Premier masque binaire ');

[H1, W1, nb_images1] = size(im_mask);

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

load dino_Ps;
pts = load('viff.xy');

% Reconstruction des points 3D
X = []; % Contient les coordonnees des points en 3D
color = []; % Contient la couleur associee
% Pour chaque couple de points apparies
for i = 1:size(pts,1)
    % Recuperation des ensembles de points apparies
    l = find(pts(i,1:2:end)~=-1);
    % Verification qu'il existe bien des points apparies dans cette image
    if size(l,2) > 1 & max(l)-min(l) > 1 & max(l)-min(l) < 36
        A = [];
        R = 0;
        G = 0;
        B = 0;
        % Pour chaque point recupere, calcul des coordonnees en 3D
        for j = l
            A = [A;P{j}(1,:)-pts(i,(j-1)*2+1)*P{j}(3,:);
            P{j}(2,:)-pts(i,(j-1)*2+2)*P{j}(3,:)];
            R = R + double(im(int16(pts(i,(j-1)*2+1)),int16(pts(i,(j-1)*2+2)),1,j));
            G = G + double(im(int16(pts(i,(j-1)*2+1)),int16(pts(i,(j-1)*2+2)),2,j));
            B = B + double(im(int16(pts(i,(j-1)*2+1)),int16(pts(i,(j-1)*2+2)),3,j));
        end;
        [U,S,V] = svd(A);
        X = [X V(:,end)/V(end,end)];
        color = [color [R/size(l,2);G/size(l,2);B/size(l,2)]];
    end;
end;
fprintf('Calcul des points 3D termine : %d points trouves. \n',size(X,2));

% % Affichage du nuage de points 3D
% figure;
% hold on;
% for i = 1:size(X,2)
%     plot3(X(1,i),X(2,i),X(3,i),'.','col',color(:,i)/255);
% end;
% axis equal;


T=DelaunayTri(X(1,:)',X(2,:)',X(3,:)');

fprintf('Tetraedrisation terminee : %d tetraedres trouves. \n',size(T,1));

% Affichage de la tetraedrisation de Delaunay
figure;
tetramesh(T);
title ('Affichage de la tetraedrisation de Delaunay');


% Initialisation des barycentres
nb_barycentres = 5;
poids = cat(3, ...
    [1/4; 1/4; 1/4; 1/4], ...
    [0.97; 0.01; 0.01; 0.01], ...
    [0.01; 0.97; 0.01; 0.01], ...
    [0.01; 0.01; 0.97; 0.01], ...
    [0.01; 0.01; 0.01; 0.97] ...
);

tri = T.Triangulation;  % Chaque ligne contient 4 indices (sommets d'un tetra)

nb_tetra = size(tri,1);
C_g = zeros(3, nb_tetra, nb_barycentres);  % barycentres homogènes (X Y Z 1)

% Calcul des barycentres biaisés
for i = 1:nb_tetra
    verts = X(1:3, tri(i,:)); 
    for k = 1:nb_barycentres
        C_g(:,i,k) = verts * poids(:,:,k);  % barycentre biaisé 
    end
end

% A DECOMMENTER POUR VERIFICATION
% A RE-COMMENTER UNE FOIS LA VERIFICATION FAITE
% Visualisation pour vérifier le bon calcul des barycentres
% for i = 1:nb_images
%    for k = 1:nb_barycentres
%        o = P{i}*[C_g(:,i,k);1];
%        o = o./repmat(o(3,:),3,1);
%        imshow(im_mask(:,:,i));
%        hold on;
%        plot(o(2,:),o(1,:),'rx');
%        pause;
%        close;
%    end
% end

% Initialiser la liste de tétraèdres à garder
mask_keep = true(nb_tetra, 1);  


for i = 1:nb_tetra

    for k = 1:nb_barycentres
        barycentre=true;
        for j = 1:nb_images
            proj = P{j} * [C_g(:,i,k);1];
            u = round(proj(1)/proj(3));
            v = round(proj(2)/proj(3));

            if isfinite(u) && isfinite(v) && v >= 1 && v <= W1 && u >= 1 && u <= H1 
                    if im_mask(u,v,j) == 0
                         barycentre = false; 
                         break; 
                    end
            end
        end
        if barycentre==false
            mask_keep(i)=false;
            break
        end
    end   
end


% Affichage du maillage filtré
tri_bis = tri(mask_keep,:);
fprintf("Tétraèdres conservés : %d / %d\n", sum(mask_keep), nb_tetra);

figure;
trisurf(tri_bis, X(1,:), X(2,:), X(3,:), 'FaceColor','cyan','EdgeColor','k');
title('Maillage avec tétraedres filtrés ');
axis equal;

% Sauvegarde des donnees
%save donnees;

tetra_sauv=tri_bis;

% Calcul des faces du maillage à garder
faces = [
    tetra_sauv(:, [1 2 3]);
    tetra_sauv(:, [1 2 4]);
    tetra_sauv(:, [1 3 4]);
    tetra_sauv(:, [2 3 4])
];

fprintf('Nombre de faces initiales : %d faces \n',size(faces,1));

faces = sort(faces, 2);  % on trie chaque élément des ligne dans l'ordre croissant pour eviter que deux memes faces est un trie de sommet différent

faces_sorted = sortrows(faces); % on trie les lignes entières

[nombre_face,~] = size(faces_sorted);
face_mask = true(nombre_face,1);

for i=1:nombre_face-1
    if isequal(faces_sorted(i,:), faces_sorted(i+1,:)) % si deux lignes consécutives sont égales
        face_mask(i)=false; % on les retirent du masque
        face_mask(i+1)=false;
    end
end

faces_surface = faces_sorted(face_mask, :); % On enlève les faces qui apparaissent 2 fois

fprintf('Calcul du maillage final termine : %d faces. \n',size(faces_surface,1));

% Affichage du maillage
figure;
trisurf(faces_surface, X(1,:), X(2,:), X(3,:), ...
        'FaceColor', 'red', 'EdgeColor', 'k');
title('Maillage final de l’objet 3D');
axis equal;

% %Affichage du maillage final
% figure;
% hold on
% for i = 1:size(faces_surface,1)
%    plot3([X(1,faces_surface(i,1)) X(1,faces_surface(i,2))],[X(2,faces_surface(i,1)) X(2,faces_surface(i,2))],[X(3,faces_surface(i,1)) X(3,faces_surface(i,2))],'r');
%    plot3([X(1,faces_surface(i,1)) X(1,faces_surface(i,3))],[X(2,faces_surface(i,1)) X(2,faces_surface(i,3))],[X(3,faces_surface(i,1)) X(3,faces_surface(i,3))],'r');
%    plot3([X(1,faces_surface(i,3)) X(1,faces_surface(i,2))],[X(2,faces_surface(i,3)) X(2,faces_surface(i,2))],[X(3,faces_surface(i,3)) X(3,faces_surface(i,2))],'r');
% end;