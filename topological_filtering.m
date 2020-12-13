function [conspec_topo, decision] = topological_filtering(conspec_sig)

[n_pnts, ~, ~, n_bands] = size(conspec_sig);
conspec_topo = zeros(size(conspec_sig));

for t = 1 : n_pnts
    
    conspec_sig_par = conspec_sig; 
    
    for b = 1 : n_bands
        
        graph = squeeze(conspec_sig_par(t, ...
            :, :, b));

        % Obtain the thresholded connectivity 
        % matrix based on orthogonal MSTs
        [~, conspec_topo_par, ~, ~, ~, ~] = ...
            threshold_omst_gce_wu_very_fast(graph, 0);
        
        conspec_topo(t, :, :, b) = conspec_topo_par;
        
    end
    
end

decision = ones(size(conspec_topo));
decision(conspec_topo == 0) = 0;

