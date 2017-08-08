Robust Motion Segmentation via Lossy Compression v 1.0

We have included four example motion sequences to test our motion segmentation algorithm:
1R2RC: a checkerboard sequence, 
arm: an articulated motion sequence, 
cars10: a traffic sequence, and 
oc1R2RC: a checkerboard sequence with missing entries

Before using our software, make sure that the 'helpers' directory is in the path:
> addpath 'helpers';

To use our L1-based methods for entry completion and error correction, the CVX package by Stephen Boyd (http://www.stanford.edu/~boyd/cvx/)
must be installed and in the MATLAB path.

** Trying a clean sequence (no incomplete or corrupted data) **
Here is code that will give you results for one clean sequence (The 'arm' sequence with projection down to R^5).

> epsilon = logspace(-5,3,101);
> [rawData, trueLabels] = load_sequence('arm');
> processedData = process_sequence(rawData, true);
> result = try_sequence('arm', processedData, epsilon);
> computedLabels = find_best_segmentation(result, processedData, 2, epsilon);
> err = compare_labels(trueLabels, computedLabels);

** Trying a sequence with incomplete data **
Here is code that will give you results for one incomplete sequence (The 'oc1R2RC' sequence with projection down to R^d, where d is the 'sparsity preserving' dimension).

> epsilon = logspace(-5,3,101);
> [rawData, trueLabels, mask] = load_sequence('oc1R2RC');
> processedData = process_sequence(rawData, 'sparse', 'incomplete');
> result = try_sequence('oc1R2RC', processedData, epsilon);
> computedLabels = find_best_segmentation(result, processedData, 3, epsilon);
> err = compare_labels(trueLabels, computedLabels);

** Trying a sequence with corrupted data **
Here is code that will give you results for one incomplete sequence (The 'oc1R2RC' sequence with projection down to R^d, where d is the 'sparsity preserving' dimension).

> epsilon = logspace(-5,3,101);
> [rawData, trueLabels] = load_sequence('oc1R2RC');
> processedData = process_sequence(rawData, 'sparse', 'corrupted');
> result = try_sequence('oc1R2RC', processedData, epsilon);
> computedLabels = find_best_segmentation(result, processedData, 3, epsilon);
> err = compare_labels(trueLabels, computedLabels);



