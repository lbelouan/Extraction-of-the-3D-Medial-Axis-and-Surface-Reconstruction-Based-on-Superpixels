
% FAIBLES_GRADIENTS Ajuste les centres des superpixels vers les gradients les plus faibles
%
%   Entrées :
%   - L_channel : canal de luminance de l'image (matrice)
%   - centers : positions initiales des centres [x, y]
%   - numClusters : nombre de centres
%   - rows, cols : dimensions de l'image
%
%   Retourne :
%   - centers : centres mis à jour

function centers = Faibles_gradients(L_channel, centers, numClusters, rows, cols)

    % Définition des filtres Sobel
    sobel_x = [-1  0  1; -2  0  2; -1  0  1]; % bords verticaux
    sobel_y = [-1 -2 -1;  0  0  0;  1  2  1]; % bords horizontaux

    % Calcul du gradient
    Gx = conv2(L_channel, sobel_x, 'same');
    Gy = conv2(L_channel, sobel_y, 'same');
    gradMag = sqrt(Gx.^2 + Gy.^2);

    % Mise à jour des centres
    for k = 1:numClusters
        x = round(centers(k, 1));
        y = round(centers(k, 2));

        minGradient = gradMag(y, x);
        for dx = -1:1
            for dy = -1:1
                newX = min(max(x + dx, 1), cols); % pour s’assurer que l’on reste dans les limites de l’image
                newY = min(max(y + dy, 1), rows);
                if gradMag(newY, newX) < minGradient
                    minGradient = gradMag(newY, newX);
                    centers(k, 1) = newX;
                    centers(k, 2) = newY;
                end
            end
        end
    end
    centers=round(centers);
end
