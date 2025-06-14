
% Algorithme_SLIC Segmentation superpixels par l'algorithme SLIC
%
%   Entrées :
%   - L_channel, a_channel, b_channel : composantes Lab de l'image
%   - centers : tableau (numClusters x 5), avec [x y L a b] pour chaque centre
%   - S : intervalle spatial entre centres
%   - numClusters : nombre de centres = nombre de superpixels initiaux 
%   - rows, cols : dimensions de l'image
%   - max_iterations : nombre d'itérations SLIC
%   - m : pondération distance couleur vs spatiale
%
%   Sorties :
%   - labels : matrice de segmentation (rows × cols)
%    - centers : centres mis à jour

function [labels, centers] = Algorithme_SLIC(L_channel, a_channel, b_channel, centers, S, numClusters, rows, cols, max_iterations, m)

    % Initialisation
    labels = -ones(rows, cols);
    distances = inf(rows, cols);

    for iter = 1:max_iterations
        for k = 1:numClusters
            xk = round(centers(k, 1));
            yk = round(centers(k, 2));

            xmin = max(1, xk - S);
            xmax = min(cols, xk + S);
            ymin = max(1, yk - S);
            ymax = min(rows, yk + S);

            [Xg, Yg] = meshgrid(xmin:xmax, ymin:ymax);

            region_L = L_channel(ymin:ymax, xmin:xmax);
            region_a = a_channel(ymin:ymax, xmin:xmax);
            region_b = b_channel(ymin:ymax, xmin:xmax);

            Dc = sqrt((region_L - centers(k,3)).^2 + ...
                      (region_a - centers(k,4)).^2 + ...
                      (region_b - centers(k,5)).^2);

            Ds = sqrt((Xg - xk).^2 + (Yg - yk).^2);

            D = sqrt(Dc.^2 + (Ds / S * m).^2);

            update_mask = D < distances(ymin:ymax, xmin:xmax);
            distances(ymin:ymax, xmin:xmax) = min(distances(ymin:ymax, xmin:xmax), D);
            labels(ymin:ymax, xmin:xmax) = k .* update_mask + labels(ymin:ymax, xmin:xmax) .* ~update_mask;

        end

        % Mise à jour des centres
        for k = 1:numClusters
            [r, c] = find(labels == k);
            if ~isempty(r)
                centers(k, 1) = mean(c);  % x
                centers(k, 2) = mean(r);  % y
                idx = sub2ind([rows, cols], r, c);
                centers(k, 3) = mean(L_channel(idx));
                centers(k, 4) = mean(a_channel(idx));
                centers(k, 5) = mean(b_channel(idx));
            end
        end
    end

    % On repasse en valeur entière car on veut des indices entiers de lignes/clonnes
    centers(:, 1:2) = round(centers(:, 1:2));
end
