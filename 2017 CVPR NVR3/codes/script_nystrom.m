function script_nystrom(dataset)
%SCRIPT Run Spectral clustering using the Nystrom method with different
%   parameter settings. Two data sets are used for experiments.
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
  sample_num_array = [20 50 100 200 500 1000 1500 2000];
  sigma = 20;
  num_clusters = 18;
  disp('Reading data...');
  load('data/corel_feature.mat', 'feature');
  load('data/corel_label.mat', 'label');
elseif dataset == 2 % RCV1
  sample_num_array = [200 500 800 1000 1500 2000 2500 3000 3500];
  sigma = 2;
  num_clusters = 103;
  disp('Reading data...');
  load('data/rcv_feature.mat', 'feature');
  load('data/rcv_label.mat', 'label');
end

%
% Main program
%
for j = 1:size(sample_num_array, 2)
  num_samples = sample_num_array(1, j);
  disp(['Number of random samples: ', num2str(num_samples)]);

  for i = 1:10 % Average over 10 runs
    disp(['Iter: ', num2str(i)]);
    [cluster_labels evd_time kmeans_time total_time] = nystrom(feature, num_samples, sigma, num_clusters);
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
