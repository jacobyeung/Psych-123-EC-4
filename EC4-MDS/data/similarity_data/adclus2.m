function [clusters,weights,vaf]=adclus2(similarity,numberclusters,patience)

% ADCLUS2 additive clustering, stochastic hill-climb (michael.d.lee@dsto.defence.gov.au)
% [clusters,weights,vaf]=adclus(similarity,numberclusters,patience)
% 
% SIMILARITY is an NxN symmetric matrix of pairwise similarities (required)
% NUMBERCLUSTERS specifies the number of clusters to use (required)
% PATIENCE specifies the number of local-minima restarts attempted without improvement before terminating (default=10)
%
% CLUSTERS returns an Nx(NUMBERCLUSTERS+1) matrix defining derived cluster membership plus the universal cluster
% WEIGHTS returns a vector of length (NUMBERCLUSTERS+1) containing the weights of the clusters
% VAF return the variance of the similarity values accounted for by the solution

% check the number of arguments
error(nargchk(2,3,nargin));

% check the similarity matrix
[n check]=size(similarity);
if check~=n
   error('similarity matrix must be square');
end;
if ~isequal(similarity,similarity')
   error('similarity matrix must be symmetric');
end;

% check the number of clusters
if (numberclusters<1)|(numberclusters~=round(numberclusters))
   error('number of clusters must be an integer >= 2');
end;

% set default arguments as necessary
if nargin<3, patience=10; end;

% check the maximum number of trials
if (patience<1)|(patience~=round(patience))
   error('patience must be a positive integer');
end;


%rename variables
maxpatience=patience;
s=similarity;
m=numberclusters;

% normalise similarities to lie between 0 and 1
reshift=min(min(s));
s=s-reshift;
rescale=max(max(s));
s=s/rescale;

% calculate the variance of  the similarity matrix
sbar=(sum(sum(s))-trace(s))/n/(n-1);
temp=(s-sbar*ones(n)).^2;
vard=.5*(sum(sum(temp))-trace(temp));

% express the similarity matrix as a column vector
flats=[];
for i=1:n-1
   for j=i+1:n
      flats=[flats;s(i,j)];
   end;
end;

% init other variables and constants
npairs=round(n*(n-1)/2);
besterr=vard;
veclength=n*m;
x=round(rand(veclength,1));
bestseenx=x;
bestw=zeros(m+1,1);
bestvaf=0;
patience=0;
start=1;

while patience<maxpatience
   
   patience=patience+1;
   if start==0
      msg=sprintf('%d local minima restart',patience);
      disp(msg);
      pause(.001);
   end;
   start=0;
   
   % a random 2 (or whatever) bits in the least-weighted cluster are flipped
   % to try and shake out of a local min
   x=bestseenx;
   [val ind]=sort(bestw(1:end-1));
   scrambleflip=(ind(1)-1)*n+randperm(n);
   x(scrambleflip(1:3))=1-x(scrambleflip(1:3));
   
   %  hill-climb 
   alldone=0;
   bestseen=vard;
   while alldone==0
      
      % random order in which bits will be flipped
      fliptry=randperm(veclength);
      
      % find one step
      trycount=0;
      moveon=0;
      while moveon==0
         trycount=trycount+1;
         x(fliptry(trycount))=1-x(fliptry(trycount));
         
         % augment x with universal cluster
         f=reshape(x,n,m);
         f=[f ones(n,1)];
         
         % form of cluster membership needed for non-negative least squares
         flatf=zeros(npairs,m+1);
         count=0;
         for a=1:n-1
            for b=a+1:n
               count=count+1;
               flatf(count,:)=f(a,:).*f(b,:);
            end;
         end;
         
         % find weights
         w=diag(nnls(flatf,flats));
         
         % evaluate solution
         sh=f*w*f';
         se=(s-sh).^2;
         sse=.5*(sum(sum(se))-trace(se));
         
         if sse<bestseen
            bestseen=sse;
            bestseenx=x;
            moveon=1;
         else
            x(fliptry(trycount))=1-x(fliptry(trycount));
         end;
         
         if sse<besterr
            patience=0;
            besterr=sse;
            bestvaf=1-sse/vard;
            bestx=x;
            bestw=diag(w)';
            msg=sprintf('better clustering found: accounts for %1.2f percent of the variance',bestvaf*100);
            disp(msg);
            pause(.001);
         end;
         
         if trycount==veclength
            moveon=1;
            alldone=1;
         end;
      end; % move on   
   end; % all done
end; % while patience

% return the best cluster membership matrix, its associated weights
% and percentage of variance it accounts for
clusters=reshape(bestx,n,m);
clusters=[clusters ones(n,1)];
weights=bestw*rescale+reshift;
vaf=bestvaf;