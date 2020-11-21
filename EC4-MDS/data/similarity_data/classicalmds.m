function [points,vaf]=classicalmds(distance,dimensions)

% Classical multidimensional scaling, (michael.d.lee@dsto.defence.gov.au)
% [points,vaf]=classicalmds(distance,dimensions)
% 
% DISTANCE is an NxN symmetric matrix of pairwise distances or proximities (required)
% DIMENSIONS specifies the required dimensionality of the coordinate representation (required)
%
% POINTS returns an NxDIMENSIONS matrix giving the derived coordinate locations
% VAF returns the variance of the distance values accounted for by the solution

% check the number of arguments
error(nargchk(2,2,nargin));

% check the distance matrix
[n check]=size(distance);
if check~=n
   error('distance matrix must be square');
end;
if ~isequal(distance,distance')
   error('distance matrix must be symmetric');
end;

% check the number of dimensions
if (dimensions<1)|(dimensions~=round(dimensions))
   error('number of dimensions must be a positive integer');
end;

% assign shorter argument names
d=distance;
dim=dimensions;

% normalise distances to lie between 0 and 1
reshift=min(min(d));
d=d-reshift;
rescale=max(max(d));
d=d/rescale;

% calculate preliminary matrices
a=-.5*d.^2;
b=zeros(n);
for i=1:n
   for j=1:n
      b(i,j)=a(i,j)-sum(a(i,:))-sum(a(:,j))+sum(sum(a));
   end;
end;

% do ordered eigen-decomposition
[v e]=eig(b);
e=diag(e);
[val ind]=sort(e);
ind=flipud(ind);
e=diag(e(ind));
v=v(:,ind);

% return the variance accounted for
% and coordinate location of the final solution
dh=zeros(n);
for i=1:n-1;
   for j=i+1:n
      dh(i,j)=norm(v(i,1:dim)-v(j,1:dim));
   end;
end;
dh=dh+dh';
temp=corrcoef(reshape(d,n^2,1),reshape(dh,n^2,1));
vaf=temp(2,1)^2;
points=v(:,1:dim)*rescale+reshift;

