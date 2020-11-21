function [violations,variance,centrality,reciprocity]=diprox(proximity)

% DIPROX diagnostics for a proximity matrix (michael.d.lee@dsto.defence.gov.au)
% [violations,variance,centrality,reciprocity]=diprox(proximity)
% 
% PROXIMITY is an NxN symmetric matrix of pairwise proximities (required)
%
% VIOLATIONS returns the proportion of triples that violate the triangle inequality
% VARIANCE returns the variance of the proximities about their arithmetic mean
% CENTRALITY returns Tversky and Hutchinson's Centrality statistic
% RECIPROCITY returns Tversky and Hutchinson's Reciprocity statistic
%
% Histograms of the proximities, and the ratio of triples (largest edge divided
% by sum of smaller edges) are also generated

% check the number of arguments
error(nargchk(1,1,nargin));

% check the proximity matrix
[n check]=size(proximity);
if check~=n
   error('proximity matrix must be square');
end;
if ~isequal(proximity,proximity')
   error('proximity matrix must be symmetric');
end;

%rename variables
d=proximity;

% ensure zero self-distances
for i=1:n
   d(i,i)=0;
end;

% calculate the variance of  the distance matrix
dbar=(sum(sum(d))-trace(d))/n/(n-1);
temp=(d-dbar*ones(n)).^2;
variance=.5*(sum(sum(temp))-trace(temp));

flatd=tril(d,-1);
flatd=flatd(find(flatd>0));

figure;clf;
subplot(2,1,1);hist(flatd,20);
set(gca,'xlim',[0 1]);
xlabel('Proximity','fontweight','bold');
ylabel('Frequency','fontweight','bold');

x=nchoosek([1:n],3);

[num junk]=size(x);

ratio=zeros(num,1);

for i=1:num
   dists=[d(x(i,1),x(i,2)),d(x(i,1),x(i,3)),d(x(i,2),x(i,3))];
   [big ind]=max(dists);
   little=dists(setdiff([1:3],ind));
   ratio(i)=big/sum(little);
end;

violations=sum(ratio>1)/num;
subplot(2,1,2);hist(ratio,20);
xlabel('Triangle Ratio','fontweight','bold');
ylabel('Frequency','fontweight','bold');

% centrality

% temporarily set diagonal to large number
% prevents object being its own nearest neighbour
big=max(max(d))+1;
for i=1:n
   d(i,i)=big;
end;

% count nearest neighbours for each object
numnearest=zeros(n,1);
[val nearest]=min(d);
for i=1:n
   numnearest(i)=length(find(nearest==i));
end;

% calculate centrality statistic
centrality=sum(numnearest.^2)/n;


% reciprocity
reciprank=zeros(n,1);
for i=1:n
   [val ind]=sort(d(nearest(i),:));
   [val reciprank(i)]=find(ind==i);
   end;

% calculate reciprocity statistic
reciprocity=sum(reciprank)/n;

% reset diagonals
for i=1:n
   d(i,i)=0;
end;
