function [Gxy] = compute_network_metric(Cxy, metric)

% Read dimension of the problem 
[n_pnts, n_chans, ~, n_bands] = size(Cxy);
Gxy = zeros(n_pnts, n_chans, n_bands);

switch metric 

    % Weighted node degree
    case 'wnd'
        
        Gxy = squeeze(sum(Cxy,3));
     
    % Clustering coefficient     
    case 'cc'
        
        for t = 1 : n_pnts
            
            for b = 1 : n_bands
                    
                    graph = squeeze(Cxy(t, :, :, b));
                    Gxy(t, :, b) = ...
                        clustering_coef_wu(graph);
                    
            end
        
        end

end 

    end
