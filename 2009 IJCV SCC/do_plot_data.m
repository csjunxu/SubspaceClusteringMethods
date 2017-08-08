function opts = do_plot_data(data,labels,opts)

[N,D] = size(data);

if N == 1 % row vector
    data = data'; % convert to column vector
    N = D; 
    D = 1; 
end

if nargin<2 || nargin==3 && isempty(labels)
    labels = ones(N,1); % one cluster
end

if nargin<3
    opts = struct();
end

if ~isfield(opts, 'markerSize')
    opts.markerSize = 12;
end

if ~isfield(opts, 'markers')
    if nargin>1 && max(labels)>1
        opts.markers ='o+v*d';
    else
        opts.markers = '.';
    end
end

if ~isfield(opts,'colors')
    if nargin>1 && max(labels)>1
        opts.colors = [1 0 0; 0 0.4 0; 0 0 1; 1 0.2 0.6; 0.8 0.4 0];
    else 
        opts.colors = [0,0,1];
    end
end

%% large D 
if D > 3
    
    if ~isfield(opts, 'view')
        opts.view = 'pca';
    end
    
    switch lower(opts.view)
        case 'first3'
            data = data(:,1:3);
        case'last3'
            data = data(:,D-2:D);
        case 'pca'
            [U,S] = svds(data - repmat(mean(data,1),N,1),6);
            %[U,S] = svds(data,6);
            data = U(:,1:3).*repmat(transpose(diag(S(1:3,1:3))), N,1);
    end
    
end

%% large K
K = max(labels); % number of clusters

if K>length(opts.markers)
    opts.markers(6:K) = randsample('o+v*dxsph', K-5, true);
    opts.colors(6:K, :) = rand(K-5,3);
else
    opts.markers = opts.markers(1:K);
    opts.colors = opts.colors(1:K,:);
end

%%
inds_0 = find(labels<1); % outliers

%figure
%axes('position',[0  0  1  1])

hold on

if D>=3
    for k = 1:K
        inds_k = find(labels==k);
        plot3(data(inds_k,1), data(inds_k,2), data(inds_k,3), opts.markers(k), 'color', opts.colors(k,:), 'MarkerSize',opts.markerSize)
    end
    plot3(data(inds_0,1),data(inds_0,2),data(inds_0,3),'k.','MarkerSize',opts.markerSize)
    %axis equal
elseif D==2
    for k = 1:K
        inds_k = find(labels==k);
        plot(data(inds_k,1), data(inds_k,2), opts.markers(k), 'color', opts.colors(k,:),'MarkerSize',opts.markerSize)
    end
    plot(data(inds_0,1),data(inds_0,2),'k.','MarkerSize',opts.markerSize)
    %axis equal
else %% D=1
    for k = 1:K
        inds_k = find(labels==k);
        plot(inds_k, data(inds_k,1), opts.markers(k), 'color', opts.colors(k,:),'MarkerSize',opts.markerSize)
    end
    plot(inds_0, data(inds_0,1),'k.','MarkerSize',opts.markerSize)
end

hold off

axis tight

% % %set(gca, 'XTick',[],'YTick', [])
% if D >= 3
%     set(gca,'XTick',[], 'XTickLabel', {},'YTick', [], 'YTickLabel', [] ,'ZTick', [], 'ZTickLabel', {} )
% else
%     set(gca,'XTick',[], 'XTickLabel', {},'YTick', [], 'YTickLabel', {})
% end
% box on
