function [theta]=subspaceangle(base1, base2)
%function [theta]=subspaceangle(base1, base2)
%compute the principal angles between two subspaces
%   base1, base2    two ORTHONORMAL basis
%   theta           the subspaces angles
%
%   Example
%   subspaceangle([1 0 0 0; 0 1 0 0]',orth([0 1 1 0; 0 0 0 1]'))
%
%   ans =
%
%       1.5708 ----->(90 degrees)
%       0.7854 ----->(45 degrees)

[U,S,V]=svd(base1'*base2,'econ');
costheta=flipud(svd(base2'*base1));
theta=acos(min(1,costheta));