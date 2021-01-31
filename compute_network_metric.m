function [Gxy] = compute_network_metric(Cxy, metric)

% Read dimension of the problem 
[n_pnts, n_chans, ~, n_bands] = size(Cxy);
Gxy = zeros(n_pnts, n_chans, n_bands);

switch metric 
    % Weighted node degree
    case 'wnd'  
        Gxy = squeeze(sum(Cxy,3));
    % Nodal efficiency    
    case 'wne'
        for t = 1 : n_pnts
            for b = 1 : n_bands        
                graph = squeeze(Cxy(t, :, :, b));
                Gxy(t, :, b) = efficiency_wei(graph, 3);
            end
        end
    % Betweenness centrality      
    case 'bc'
        for t = 1 : n_pnts
            for b = 1 : n_bands        
                graph = squeeze(Cxy(t, :, :, b));
                Gxy(t, :, b) = betweenness_wei(graph);
            end
        end  
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
