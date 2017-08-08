function script_sc(dataset)
%SCRIPT Run spectral clustering using a sparse matrix (nearest neighbors)
%   with different parameter settings. Two data sets are used for experiments.
%
%   Input  : dataset : data set number, 1 = Corel, 2 = RCV1
%
%   Output : average and standard deviation of
%            - normalized mutual information (NMI)
%            - evd_time
%            - kmeans_time
%            - total_time
%
%   Author : Wen-Yen Chen (wychen@alumni.cs.ucsb.edu)
%			 Chih-Jen Lin (cjlin@csie.ntu.edu.tw)

%
% Parameter settings for selected data set
%
if dataset == 1 % Corel
  nn_num_array = [5 10 15 20 50 100 150 200]; % Number of nearest neighbors
  sigma = 20;
  %sigma = 0; % for self-tuning
  num_clusters = 18;
elseif dataset == 2 % RCV1
  nn_num_array = [20 50 75 100 150 200]; % Number of nearest neighbors
  sigma = 2;
  %sigma = 0; % for self-tuning
  num_clusters = 103;
end

%
% Main program
%
for j = 1:size(nn_num_array,2)
  num_neighbors = nn_num_array(1, j);
  disp(['Number of nearest neighbor: ', num2str(num_neighbors)]);

  disp('Reading data...');
  if dataset == 1 % Corel
    input_file = ['data/corel_', num2str(num_neighbors), '_NN_sym_distance.mat'];
    load(input_file, 'A'); % sparse distance matrix 'A'
    load('data/corel_label.mat', 'label');
  elseif dataset == 2 % RCV1
    input_file = ['data/rcv_', num2str(num_neighbors), '_NN_sym_distance.mat'];
    load(input_file, 'A'); % sparse distance matrix 'A'
    load('data/rcv_label.mat', 'label');
  end

  for i = 1:10 % Average over 10 runs
    disp(['Iter: ', num2str(i)]);
    [cluster_labels evd_time kmeans_time total_time] = sc(A, sigma, num_clusters);
    nmi_score = nmi(label, cluster_labels) % Calculate NMI
    accuracy_score = accuracy(label, cluster_labels) % Calculate accuracy

    all_nmi(i, j) = nmi_score;
    all_accuracy(i, j) = accuracy_score;
    all_evd_time(i, j) = evd_time;
    all_kmeans_time(i, j) = kmeans_time;
    all_total_time(i, j) = total_time;
  end

  % Output the (average/std) results to stdout
  avg_nmi = mean(all_nmi, 1)
  std_nmi = std(all_nmi, 0, 1)
  avg_accuracy = mean(all_accuracy, 1)
  std_accuracy = std(all_accuracy, 0, 1)
  avg_evd_time = mean(all_evd_time, 1)
  std_evd_time = std(all_evd_time, 0, 1)
  avg_kmeans_time = mean(all_kmeans_time, 1)
  std_kmeans_time = std(all_kmeans_time, 0, 1)
  avg_total_time = mean(all_total_time, 1)
  std_total_time = std(all_total_time, 0, 1)
end
