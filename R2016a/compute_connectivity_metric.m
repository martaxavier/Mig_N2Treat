function conspec = compute_connectivity_metric(cspec, metric)

%   Computes the connectivity metric specified in 'metric', from the 
%   input cross-spectrum, 'cross_spectrum'
%
%   INPUTS:
%
%           cross_spectrum  - Cross-spectrum lower triangular matrix 
%           metric          - Connectivity metric 
% 
%   OUTPUTS:
%
%           conspec         - Output connectivity spectrum 

% Compute the functional connectivity 
conspec = get_Cxy(cspec, metric);

end


    function [Cxy] = get_Cxy(cspec, metric)

    % Read dimension of the problem 
    [~,~,n_chans,~]= size(cspec);

    % Pre-allocate Cxy
    Cxy = zeros(size(cspec));

    % Go through channels 
    for c1 = 1 : n_chans

        % Lower triangular matrix 
        for c2 = 1 : c1

            switch metric 

                 % == Imaginary Coherence 
                case 'icoh'

                    % Cross-spectrum of the current pair of channels
                    Sxy = cspec(:, :, c1, c2);
                    Sxx = cspec(:, :, c1, c1);
                    Syy = cspec(:, :, c2, c2);

                    % Coherence - Cxy = Sxy/sqrt(Sxx*Syy)
                    Cxy_chan = bsxfun(@rdivide, Sxy, sqrt(Sxx));
                    Cxy_chan = bsxfun(@rdivide, Cxy_chan, sqrt(Syy));

                    % Imaginary Part of Coherency - ICxy = abs(im(Cxy))
                    Cxy_chan =  abs(imag(Cxy_chan));

                % == Phase Locking Value 
                case 'plv'


            end

            Cxy(:, :, c1, c2) = Cxy_chan;

        end 

    end

    end

    