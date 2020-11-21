function [clusters,weights,vaf]=adclus(similarity,numberclusters,learnrate,maxtrials,batchsize,batchlearn,mutprob,mutshift)

% ADCLUS additive clustering (michael.d.lee@dsto.defence.gov.au)
% [clusters,weights,vaf]=adclus(similarity,numberclusters,learnrate,maxtrials,batchsize,batchlearn,mutprob,mutshift)
% 
% SIMILARITY is an NxN symmetric matrix of pairwise similarities (required)
% NUMBERCLUSTERS specifies the number of clusters to use (required)
% LEARNRATE specifies the learning rate used in optimisation (default=0.1)
% MAXTRIALS specifies the number of trees evaluated without improvement before terminating (default=3000)
% BATCHSIZE specifies the size of the batch of potential tree solutions (default=50)
% BATCHLEARN specifies the number of best solutions in the batch used in learning (default=1)
% MUTPROB specifies the probability of 'mutating' the PBIL probability vector (default=0.02)
% MUTSHIFT specifies the extent of a 'mutation' shift (default=0.05)
%
% CLUSTERS returns an Nx(NUMBERCLUSTERS+1) matrix defining derived cluster membership plus the universal cluster
% WEIGHTS returns a vector of length (NUMBERCLUSTERS+1) containing the weights of the clusters
% VAF returns the variance of the similarity values accounted for by the solution

% check the number of arguments
error(nargchk(2,8,nargin));

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
if nargin<8, mutshift=0.05; end;
if nargin<7, mutprob=0.02; end;
if nargin<6, batchlearn=1; end;
if nargin<5, batchsize=50; end;
if nargin<4, maxtrials=3000; end;
if nargin<3, learnrate=0.1; end;

% check the maximum number of trials
if (maxtrials<1)|(maxtrials~=round(maxtrials))
   error('maxtrials must be a positive integer');
end;

% check the batchsize
if (batchsize<1)|(batchsize~=round(batchsize))
   error('batchsize must be a positive integer');
end;

% check the number of batchsize used to learn
if (batchlearn<1)|(batchlearn>batchsize)|(batchlearn~=round(batchlearn))
   error('batchlearn must be a positive integer between 1 and batchsize');
end;

% learnrate, mutprob and mutshift are expected to be positive numbers
% but these constraints are not explicitly imposed

%rename variables
lr=learnrate;
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
p=.5*ones(veclength,1);
updates=50; %controls how often an update message is displayed

% main PBIL loop
tr=0;
while tr<maxtrials
   % initialise solution and evaluation storage
   sv=zeros(batchsize,veclength);
   evaluate=zeros(batchsize,1);
   % generate solutions
   for i=1:batchsize
      for j=1:veclength
         sv(i,j)=(rand<p(j));
      end;
      f=reshape(sv(i,:),n,m);
      % augment universal cluster for additive constant
      f=[f ones(n,1)];
      % form of cluster membership needed for non-negative least squares
      flatf=zeros(npairs,m+1);
      count=0;
      for x=1:n-1
         for y=x+1:n
            count=count+1;
            flatf(count,:)=f(x,:).*f(y,:);
         end;
      end;
      % find weights
      w=diag(nnls(flatf,flats));
      % evaluate solution
      sh=f*w*f';
      se=(s-sh).^2;
      evaluate(i)=.5*(sum(sum(se))-trace(se));
      
      % if solution is best one found, reset trial counter,
      % save the solution and display a message
      % otherwise increment the trial counter
      if evaluate(i)<besterr
         tr=0;
         besterr=evaluate(i);
         bestvaf=1-besterr/vard;
         best=f;
         bestw=diag(w);
         msg=sprintf('better clustering found: accounts for %1.2f percent of the variance',bestvaf*100);
         disp(msg);
         pause(.001);
      else
         tr=tr+1;
      end;
      
      % periodically display the trial counter
      if (tr>0)&(mod(tr,updates)==0)
         msg=sprintf('%d trials have elapsed without improvement',tr);
         disp(msg);
         pause(.001);
      end;
   end; %for i
   
   % sort the batch of evaluated solutions
   [evaluate order]=sort(evaluate);
   sv=sv(order,:);
   
   % use the best mupdate with competitive learning to adjust probability vectors
   for i=batchlearn:-1:1
      p=(1-lr)*p+sv(i,:)'*lr;
   end;
   
   % stochastically 'mutate' probability vector using mutshift and mutprob
   for i=1:veclength
      if rand<mutprob
         if p(i)>.5
            p(i)=p(i)*(1-mutshift);
         else
            p(i)=p(i)*(1-mutshift)+mutshift;
         end;
      end;
   end;
   
end; %while tr

% return the best cluster membership matrix, its associated weights
% and percentage of variance it accounts for
clusters=best;
weights=bestw*rescale+reshift;
vaf=bestvaf;