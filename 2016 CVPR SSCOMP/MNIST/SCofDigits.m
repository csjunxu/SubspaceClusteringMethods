function [feature] = SCofDigits(data)
% Extract features by scattering transform.

% Copyright Chong You @ Johns Hopkins University, 2016
% chong.you1987@gmail.com

N = size(data, 2);

foptions.J = 3; foptions.L = 8; soptions.M = 2; soptions.oversampling = 0;
[Wop, ~] = wavelet_factory_2d([28 28], foptions, soptions);
for ii = 1:N
    if mod(ii, 100) == 0
        fprintf( '%d finished\n' , ii)
    end
    im = reshape(data(:, ii), 28, 28);
    sc_digit = scat(im, Wop);

    sc = format_scat(sc_digit);
    
    if ii == 1 % set parameters and allocate space for ``feature''
        Npath = size(sc, 1);
        feature = zeros(Npath*size(sc, 2)*size(sc, 3), N);
    end
        
    feature(:, ii) = sc(:)';
end
for iF = 1:Npath % normalize each feature image to unit norm
   feature(iF:Npath:end, :) = cnormalize( feature(iF:Npath:end, :), Inf );
end

