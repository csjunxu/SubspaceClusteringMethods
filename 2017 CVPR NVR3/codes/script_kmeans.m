function script_kmeans(dataset)
%SCRIPT Run kmeans clustering algorithm with different parameter settings.
%       Two data sets are used for experiments.
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
  num_clusters = 18;
  disp('Reading data...');
  load('data/corel_feature.mat', 'feature');
  load('data/corel_label.mat', 'label');
elseif dataset == 2 % RCV1
  num_clusters = 103;
  disp('Reading data...');
  load('data/rcv_feature.mat', 'feature');
  load('data/rcv_label.mat', 'label');
end

%
% Main program
%
for i = 1:10 % Average over 10 runs
  disp(['Iter: ', num2str(i)]);
  cluster_labels = k_means(feature, 'random', num_clusters);
  nmi_score = nmi(label, cluster_labels) % Calculate NMI
  accuracy_score = accuracy(label, cluster_labels) % Calculate accuracy

  all_nmi(i) = nmi_score;
  all_accuracy(i) = accuracy_score;
end

% Output the (average/std) results to stdout
avg_nmi = mean(all_nmi)
std_nmi = std(all_nmi)
avg_accuracy = mean(all_accuracy)
std_accuracy = std(all_accuracy)
