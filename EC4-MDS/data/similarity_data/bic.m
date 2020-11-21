function value=bic(bictype,data,vaf,sigma,dims_nodes_clusters)

% BIC bayesian information criterion (michael.d.lee@dsto.defence.gov.au)
% value=bic(bictype,data,vaf,sigma,dims_nodes_clusters)
% 
% BICTYPE specifies the model type using the string 'mds', 'addtree', or 'adclus' (required)
% DATA is the NxN symmetric matrix of similarities or proximities (required)
% VAF is variance of the similarity values accounted for by the representation (required)
% SIGMA is the assumed level of data precision (required)
% DIMS_NODES_CLUSTERS gives the number of dimensions, internal nodes, or clusters in the representation (required)
%
% VALUE returns the bayesian information criterion value

% check the number of arguments
error(nargchk(5,5,nargin));

% check the type
validtypes=char('mds','addtree','adclus');
if isempty(strmatch(bictype,validtypes,'exact'))
   error('bictype must be specified as ''mds'', ''addtree'', or ''adclus''')
end;

% check the similarity or proximity matrix
[n check]=size(data);
if check~=n
   error('similarity or proximity matrix must be square');
end;
if ~isequal(data,data')
   error('similarity or proximity matrix must be symmetric');
end;

% check the vaf measure
if (vaf<0)|(vaf>1)
   error('vaf should be between 0 and 1');
end;

% check the sigma value
if (sigma<0)
   error('sigma should be positive');
end;

% check the dims, nodes, or clusters
if (dims_nodes_clusters<1)|(dims_nodes_clusters~=round(dims_nodes_clusters))
   error('dims, nodes, or clusters must be a positive integer');
end;

% rename variable
m=dims_nodes_clusters;

% calculate the variance of  the data matrix
bar=(sum(sum(data))-trace(data))/n/(n-1);
temp=(data-bar*ones(n)).^2;
vard=.5*(sum(sum(temp))-trace(temp));

% calculate the sum squared error
sse=(1-vaf)*vard;

switch bictype
case 'mds', value=sse/sigma^2+m*n*log(n*(n-1)/2);
case 'addtree', value=sse/sigma^2+(m+n-1)*log(n*(n-1)/2);
case 'adclus', value=sse/sigma^2+(m+1)*log(n*(n-1)/2);
end;


