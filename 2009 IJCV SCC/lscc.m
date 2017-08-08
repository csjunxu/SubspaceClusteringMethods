function [sampleLabels,averageL2Error] = lscc(X,d,K,OPTIONS)
%
%   Linear Spectral Curvature Clustering (LSCC)
%
%   [sampleLabels, averageL2Error] = lscc(X,d,K) partitions the points 
%   in the N-by-D data matrix X into K clusters, 
%   each representing a d-dimensional linear subspace.
%   Rows of X correspond to points, columns correspond to variables.  
%   LSCC returns an N-by-1 vector sampleLabels containing the cluster 
%   label of each point and the averaged L2 error of the K detected 
%   clusters. Those with zero labels are detected outliers.  
%
%   [ ... ] = lscc(..., OPTIONS) 
%   allows you to specify optional parameters for the lscc algorithm. 
%   OPTIONS is a structure array with the following fields:
%
%   'n' - number of points used for computing a curvature (beside origin).
%         default = d+1, but can be larger than d+1.
%         Note that in the linear version, the origin is an additional
%         point which is automatically included for computing the polar
%         curvature.
%
%   The rest of the parameters are set in the same way as in scc.m. Please
%   type 'help scc' in the Matlab command window to see the details.
%
%   (c)2009-2011 Gilad Lerman and Guangliang Chen
%   Last updated on 08/10/2011.
%   Any questions please email glchen@math.duke.edu or lerman@umn.edu

%% set default options
ABSOLUTE_MINIMUM = 1e-15;

if nargin < 4
    OPTIONS = struct();
end

if ~isfield(OPTIONS,'n') || OPTIONS.n < d+1
    OPTIONS.n = d+1;
end

if ~isfield(OPTIONS,'c')
    OPTIONS.c = 100*K;
end
OPTIONS.c = K*ceil(OPTIONS.c/K);

if ~isfield(OPTIONS,'normalizeW')
    OPTIONS.normalizeW = 1;
end

if ~isfield(OPTIONS,'normalizeU')
    OPTIONS.normalizeU = 1;
end

if ~isfield(OPTIONS,'findOptimalSigma')
    OPTIONS.findOptimalSigma = 1;
end

if OPTIONS.findOptimalSigma && ~isfield(OPTIONS,'search_by')
    OPTIONS.search_by = 'index';
end

if ~isfield(OPTIONS,'seedType')
    OPTIONS.seedType = 'hard';
end

if ~isfield(OPTIONS,'alpha')
    OPTIONS.alpha = 0;
end

%% main body of code
[N,D] = size(X');

% take off those too close to the origin
magX = sum(X.^2,2);
pointsNearZero = find(magX<ABSOLUTE_MINIMUM);
S = X(setdiff(1:N,pointsNearZero),:);
m = size(S,1);
clear magX

% auxiliary 
averageL2Error0 = Inf;
sampleLabels0 = zeros(m,1);

% one cluster
if ~isfield(OPTIONS, 'initialLabels')
    % consider given data as one cluster
    sampleLabels1 = ones(m,1);
else
    % only for improving clusters obtained by other algorithms
    sampleLabels1 = OPTIONS.initialLabels;
end
averageL2Error1 = computing_average_L2_error(S, d, sampleLabels1);

q = OPTIONS.n; % q = d+1
while averageL2Error1 < averageL2Error0 * 0.99 || q>1

    sampleLabels0 = sampleLabels1;
    averageL2Error0 = averageL2Error1;
    
    sampledColumns = sampling_columns(sampleLabels1,K,OPTIONS);

    % use the origin as an auxiliary point to compute the curvatures
    polarCurv = computing_polar_curvatures([S;zeros(1,D)],[sampledColumns; (m+1)*ones(1,OPTIONS.c)],d);
    polarCurv = polarCurv(1:end-1,:);
    polarCurv_sorted = sort(reshape(polarCurv,1,m*OPTIONS.c));
    
    q = max(q-1,1);
    
    if OPTIONS.findOptimalSigma == 0 %% use a single sigma
        
        if isfield(OPTIONS,'sigma')
            sigma = OPTIONS.sigma;
        else
            sigma = max(polarCurv_sorted(1,ceil((m-OPTIONS.n+1)*OPTIONS.c/K)),ABSOLUTE_MINIMUM);
            %sigma = max(polarCurv_sorted(1,ceil((m-OPTIONS.n+1)*OPTIONS.c/K^(OPTIONS.n-1))),ABSOLUTE_MINIMUM);
        end

        isigma = 1/(2*sigma);
        A = exp(-polarCurv*isigma);
        sampleLabels1 = processing_affinities(A,K,OPTIONS);
        averageL2Error1 = computing_average_L2_error(S, d*ones(K,1), sampleLabels1);
        
    else %% search for optimal sigma
        
        averageL2Error1 = Inf;
        sampleLabels1 = zeros(m,1);
        
        sigma = max(polarCurv_sorted(1,ceil((m-OPTIONS.n+1)*OPTIONS.c/K)),ABSOLUTE_MINIMUM);
        
        switch OPTIONS.search_by
        
            case 'index'
    
                p = 1;
                while p <= OPTIONS.n-1 && sigma >= ABSOLUTE_MINIMUM

                    isigma = 1/(2*sigma);
                    A = exp(-polarCurv*isigma);

                    sampleLabels2 = processing_affinities(A,K,OPTIONS);
                    averageL2Error2 = computing_average_L2_error(S,d,sampleLabels2);

                    if averageL2Error1 > averageL2Error2
                        averageL2Error1 = averageL2Error2;
                        sampleLabels1 = sampleLabels2;
                    end

                    p = p+1;
                    sigma = polarCurv_sorted(1,ceil((m-OPTIONS.n+1)*OPTIONS.c/K^p));
                end % while p <= d+1 && sigma >= ABSOLUTE_MINIMUM
                
            case 'ratio'
                        
                sigmaMin = max(polarCurv_sorted(1,ceil((m-OPTIONS.n+1)*OPTIONS.c/K^q)),ABSOLUTE_MINIMUM);
                %sigma = polarCurv_sorted(1,(m-OPTIONS.n+1)*OPTIONS.c);
       
                isigma = 1/(2*sigma);
                A = exp(-polarCurv*isigma);

                while sigma >= sigmaMin

                    sampleLabels2 = processing_affinities(A,K,OPTIONS);
                    averageL2Error2 = computing_average_L2_error(S,d,sampleLabels2);

                    if averageL2Error1 > averageL2Error2
                        averageL2Error1 = averageL2Error2;
                        sampleLabels1 = sampleLabels2;
                    end

                    A = A.*A;
                    sigma = sigma/2;

                end % while sigma >= sigmaMin
      
        end % switch OPTIONS.search_type

    end % if OPTIONS.findOptimalSigma == 0
    
end % while averageL2Error1 < averageL2Error * 0.99

if averageL2Error1 < averageL2Error0
    averageL2Error0 = averageL2Error1;
    sampleLabels0 = sampleLabels1;
end

[centers,bases]= computing_centers_and_bases(S,sampleLabels0,d);

% get back those points too close to the origin
sampleLabels = zeros(N,1);
sampleLabels(setdiff(1:N,pointsNearZero),1) = sampleLabels0;
distancePoints2Flats = computing_point_to_flat_distances(X(pointsNearZero,:),centers,bases);
[ignored, sampleLabels(pointsNearZero)] = min(distancePoints2Flats,[],2);

averageL2Error = computing_average_L2_error(X, d*ones(K,1), sampleLabels0);

%do_plot_data('clusters obtained by LSCC',X,12,sampleLabels);
