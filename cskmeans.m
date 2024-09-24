function [cid,nr,centers] = cskmeans(x,k,nc)
% CSKMEANS K-Means clustering - general method.
%
%   This implements the more general k-means algorithm, where 
%   HMEANS is used to find the initial partition and then each
%   observation is examined for further improvements in minimizing
%   the within-group sum of squares.
%
%   [CID,NR,CENTERS] = CSKMEANS(X,K,NC) Performs K-means
%   clustering using the data given in X. 
%   
%   INPUTS: X is the n x d matrix of data,
%   where each row indicates an observation. K indicates
%   the number of desired clusters. NC is a k x d matrix for the
%   initial cluster centers. If NC is not specified, then the
%   centers will be randomly chosen from the observations.
%
%   OUTPUTS: CID provides a set of n indexes indicating cluster
%   membership for each point. NR is the number of observations
%   in each cluster. CENTERS is a matrix, where each row
%   corresponds to a cluster center.
%
%   See also CSHMEANS


%   W. L. and A. R. Martinez, 9/15/01
%   Computational Statistics Toolbox 

warning off
[n,d] = size(x);
if nargin < 3
	% Then pick some observations to be the cluster centers.
    ind = ceil(n*rand(1,k));
	% We will add some noise to make it interesting.
	nc = x(ind,:) + randn(k,d);
end

% set up storage
% integer 1,...,k indicating cluster membership
cid = zeros(1,n);	
% Make this different to get the loop started.
oldcid = ones(1,n);
% The number in each cluster.
nr = zeros(1,k);	
% Set up maximum number of iterations.
maxiter = 100;
iter = 1;

while ~isequal(cid,oldcid)	& iter < maxiter
	% Implement the hmeans algorithm
	% For each point, find the distance to all cluster centers
	for i = 1:n
		dist = sum((repmat(x(i,:),k,1)-nc).^2,2);
		[m,ind] = min(dist);	% assign it to this cluster center
		cid(i) = ind;
	end
	% Find the new cluster centers
	for i = 1:k
		% find all points in this cluster
		ind = find(cid==i);
		% find the centroid
		nc(i,:) = mean(x(ind,:));
		% Find the number in each cluster;
		nr(i) = length(ind);
	end
 	iter = iter + 1;
end

% Now check each observation to see if the error can be minimized some more. 
% Loop through all points.
maxiter = 2;
iter = 1;
move = 1;
while iter < maxiter & move ~= 0 
	move = 0;
	% Loop through all points.
	for i = 1:n
		% find the distance to all cluster centers
		dist = sum((repmat(x(i,:),k,1)-nc).^2,2);
		r = cid(i);		% This is the cluster id for x
		dadj = nr./(nr+1).*dist';	% All adjusted distances
		[m,ind] = min(dadj);	% minimum should be the cluster it belongs to
		if ind ~= r		% if not, then move x
			cid(i) = ind;
			ic = find(cid == ind);
			nc(ind,:) = mean(x(ic,:));
			move = 1;
		end
	end
	iter = iter+1;
end
centers = nc;
if move == 0
	disp('No points were moved after the initial clustering procedure.')
else
	disp('Some points were moved after the initial clustering procedure.')
end
warning on
