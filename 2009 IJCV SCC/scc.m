function [sampleLabels,averageL2Error] = scc(X,d,K,OPTIONS)

%   Spectral Curvature Clustering (SCC)
%
%   [sampleLabels, averageL2Error] = scc(X,d,K) partitions the points 
%   in the N-by-D data matrix X into K clusters, each representing a d-flat.
%   Rows of X correspond to points, columns correspond to variables.  
%   SCC returns an N-by-1 vector sampleLabels containing the cluster 
%   label of each point and the averaged L2 error of the K detected 
%   clusters. Those with zero labels are detected as outliers.  
%
%   [ ... ] = scc(..., OPTIONS) allows you to specify optional parameters 
%   to control the algorithm.  OPTIONS is a structure array consisting of
%   the following fields:
%
%   'n' - number of points used for computing a curvature
%         default = d+2, but can be larger than d+2
%
%   'c' - number of columns sampled from the matrix A, default = K*100 
%
%   'normalizeW' - 0/1, if we normalize the matrix W. default = 1(YES)
%
%   'normalizeU' - 0/1, if we normalize the matrix U. default = 1(YES)
%
%   'findOptimalSigma' - 0, if we use a single sigma for computing A;
%                        1, we search for the best sigma (default)
%
%   'sigma' - the tuning parameter used in the affinity tensor
%             can be specified by user when findOptimalSigma = 0; 
%             otherwise the algorithm will infer its optimal value from data
%
%   'search_by' - the method of searching an interval (found by algorithm) for best sigma 
%       'index' - search by index of some discrete vector (default)
%       'ratio' - repeatedly divide the upper bound by a constant (\sqrt(2))
%
%   'alpha' - number (if >=1) or percentage (if <1) of outliers in the data; 
%             default = 0
%
%   'initialLabels' - initial labelling of data, not a required parameter for scc.
%                     This option allows scc to improve clusters obtained by 
%                     other algorithms (e.g., k-flats, GPCA)

%   'seedType' - method of selecting initial seeds for kmeans clustering
%       'hard' - the deterministic procedure presented in the IJCV paper
%       'soft' - probabilistic procedure for selecting initial seeds
%
%   (c)2009-2011 Gilad Lerman and Guangliang Chen
%   Last updated on 11/11/2011.
%   If you have any questions please email glchen@math.duke.edu or
%   lerman@umn.edu.
%
%   Most Relevant Publications:
%   1. Foundations of a Multi-way Spectral Clustering Framework for Hybrid Linear Modeling, 
%      G. Chen and G. Lerman, Foundations of Computational Mathematics, 2009.
%      DOI 10.1007/s10208-009-9043-7
%   2. Spectral Curvature Clustering (SCC), G. Chen and G. Lerman,
%      International Journal of Computer Vision, 81:317-330, 2009. 
%      DOI 10.1007/s11263-008-0178-9.

%% Set default values for OPTIONS
ABSOLUTE_MINIMUM = 1e-15;

if nargin < 4
    OPTIONS = struct();
end

if ~isfield(OPTIONS,'n') || OPTIONS.n < d+2 ... 
        || (OPTIONS.n > d+2 && d == 0)
    OPTIONS.n = d+2;
end

if ~isfield(OPTIONS,'c')
    OPTIONS.c = K*100;
else
    if mod(OPTIONS.c,K)>0
        OPTIONS.c = K*ceil(OPTIONS.c/K);
        warning('The value of the parameter c in OPTIONS has been modified to be an integer multiple of K!'); %#ok<WNTAG>
    end
end

if ~isfield(OPTIONS,'findOptimalSigma')
    OPTIONS.findOptimalSigma = 1;
end

if OPTIONS.findOptimalSigma && ~isfield(OPTIONS,'search_by')
    OPTIONS.search_by = 'index';
end

if ~isfield(OPTIONS,'normalizeW')
    OPTIONS.normalizeW = 1;
end

if ~isfield(OPTIONS,'normalizeU')
    OPTIONS.normalizeU = 1;
end

if ~isfield(OPTIONS,'seedType')
    OPTIONS.seedType = 'hard';
end

if ~isfield(OPTIONS,'alpha')
    OPTIONS.alpha = 0;
end

N = size(X,1);

%% initialization

% auxiliary 
averageL2Error = Inf;
sampleLabels = zeros(N,1);

if ~isfield(OPTIONS, 'initialLabels')
    % consider given data as one cluster
    sampleLabels1 = ones(N,1);
else
    % only for improving clusters obtained by other algorithms
    sampleLabels1 = OPTIONS.initialLabels;
end

averageL2Error1 = computing_average_L2_error(X, d, sampleLabels1);

%% main body of the code

q = OPTIONS.n; % q = d+2 by default
while averageL2Error1 < averageL2Error * 0.99 || q>2
%for iii = 1    
    sampleLabels = sampleLabels1;
    averageL2Error = averageL2Error1;
    
    %if isfield(OPTIONS,'sampledColumns')
    %    sampledColumns = OPTIONS.sampledColumns;
    %    OPTIONS = rmfield(OPTIONS,'sampledColumns');
    %else
    %   sampledColumns = sampling_columns(sampleLabels1,K,OPTIONS);
    %end

    sampledColumns = sampling_columns(sampleLabels1,K,OPTIONS);
    
    polarCurv = computing_polar_curvatures(X,sampledColumns,d);
    polarCurv_sorted = sort(reshape(polarCurv,1,N*OPTIONS.c));
    
    q = max(q-1,1);
    
    if OPTIONS.findOptimalSigma == 0 %% use a single sigma
        
        if isfield(OPTIONS,'sigma')
            sigma = OPTIONS.sigma;
        else
            sigma = max(polarCurv_sorted(1,ceil((N-OPTIONS.n+1)*OPTIONS.c/K)),ABSOLUTE_MINIMUM);
            %sigma = max(polarCurv_sorted(1,ceil((N-OPTIONS.n+1)*OPTIONS.c/K^(OPTIONS.n-1))),ABSOLUTE_MINIMUM);
        end

        isigma = 1/(2*sigma);
        A = exp(-polarCurv*isigma);
        sampleLabels1 = processing_affinities(A,K,OPTIONS);
        averageL2Error1 = computing_average_L2_error(X, d*ones(K,1), sampleLabels1);
    
    else %% search for optimal sigma
        
        averageL2Error1 = Inf;
        sampleLabels1 = zeros(N,1);
        
        sigma = max(polarCurv_sorted(1,ceil((N-OPTIONS.n+1)*OPTIONS.c/K)),ABSOLUTE_MINIMUM);
        %sigma = polarCurv_sorted(1,(N-OPTIONS.n+1)*OPTIONS.c);
        
        switch OPTIONS.search_by
        
            case 'index'
    
                p = 1;
                p_max = log((N-OPTIONS.n+1)*OPTIONS.c)/log(K);
                
                while p <= min(OPTIONS.n-1,p_max) && sigma >= ABSOLUTE_MINIMUM

                    isigma = 1/(2*sigma);
                    A = exp(-polarCurv*isigma);

                    sampleLabels2 = processing_affinities(A,K,OPTIONS);
                    averageL2Error2 = computing_average_L2_error(X,d,sampleLabels2);

                    if averageL2Error1 > averageL2Error2
                        averageL2Error1 = averageL2Error2;
                        sampleLabels1 = sampleLabels2;
                    end

                    p = p+1;
                    sigma = polarCurv_sorted(1,ceil((N-OPTIONS.n+1)*OPTIONS.c/K^p));
                end % while p <= d+1 && sigma >= ABSOLUTE_MINIMUM
                
            case 'ratio'
                        
                sigmaMin = max(polarCurv_sorted(1,ceil((N-OPTIONS.n+1)*OPTIONS.c/K^q)),ABSOLUTE_MINIMUM);

                isigma = 1/(2*sigma);
                A = exp(-polarCurv*isigma);
                    
                while sigma >= sigmaMin

                    sampleLabels2 = processing_affinities(A,K,OPTIONS);
                    averageL2Error2 = computing_average_L2_error(X,d,sampleLabels2);

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

if averageL2Error1 < averageL2Error 
    averageL2Error = averageL2Error1;
    sampleLabels = sampleLabels1;
end

% figure; do_plot_data(X,sampleLabels); title('clusters obtained by SCC');