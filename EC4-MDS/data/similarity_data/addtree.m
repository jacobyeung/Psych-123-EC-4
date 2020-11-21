function [adjacency,lengths,vaf]=addtree(distance,nodes,learnrate,maxtrials,batchsize,batchlearn,mutprob,mutshift)

% ADDTREE additive tree (michael.d.lee@dsto.defence.gov.au)
% [adjacency,lengths,vaf]=addtree(distance,nodes,learnrate,maxtrials,batchsize,batchlearn,mutprob,mutshift)
% 
% DISTANCE is an NxN symmetric matrix of pairwise distances or proximities (required)
% NODES specifies the number of internal nodes to place in the tree (required)
% LEARNRATE specifies the learning rate used in optimisation (default=0.1)
% MAXTRIALS specifies the number of trees evaluated without improvement before terminating (default=3000)
% BATCHSIZE specifies the size of the batch of potential tree solutions (default=50)
% BATCHLEARN specifies the number of best solutions in the batch used in learning (default=1)
% MUTPROB specifies the probability of 'mutating' the PBIL probability vector (default=0.02)
% MUTSHIFT specifies the extent of a 'mutation' shift (default=0.05)
%
% ADJACENCY returns an (N+NODES)x(N+NODES) adjacency matrix defining the tree topology
% LENGTHS returns a vector of length (N+NODES) containing the arc-lengths for the tree
% VAF returns the variance of the distance values accounted for by the solution

% check the number of arguments
error(nargchk(2,8,nargin));

% check the distance matrix
[n check]=size(distance);
if check~=n
   error('distance matrix must be square');
end;
if ~isequal(distance,distance')
   error('distance matrix must be symmetric');
end;

% check the number of internal nodes
if (nodes<2)|(nodes~=round(nodes))
   error('number of internal nodes must be an integer >= 2');
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
d=distance;
m=nodes;

% normalise distances to lie between 0 and 1
reshift=min(min(d));
d=d-reshift;
rescale=max(max(d));
d=d/rescale;

% calculate the variance of  the distance matrix
dbar=(sum(sum(d))-trace(d))/n/(n-1);
temp=(d-dbar*ones(n)).^2;
vard=.5*(sum(sum(temp))-trace(temp));

% express the distance matrix as a column vector
flatd=[];
for i=1:n-1
   for j=i+1:n
      flatd=[flatd ;d(i,j)];
   end;
end;

% init other variables and constants
npairs=round(n*(n-1)/2);
besterr=vard;
updates=50; %controls how often an update message is displayed

% find length of inner-node vector and index its components
count=1;
start(1)=0;stop(1)=0;
start(2)=0;stop(2)=0;
for i=3:m
   start(i)=count;
   count=count+ceil(log2(i-1));
   stop(i)=count-1;
end;
count=count-1;

% initialise probability vectors
pprob=0.5*ones(count,1); %inner node links
numt=n*ceil(log2(m));
ptprob=0.5*ones(numt,1); %terminal node links

% main PBIL loop
tr=0;
while tr<maxtrials
   % initialise solution and evaluation storage
   svp=zeros(batchsize,count);
   svpt=zeros(batchsize,numt);
   evaluate=zeros(batchsize,1);
   % generate solutions
   for up=1:batchsize
      for j=1:count
         svp(up,j)=(rand<pprob(j));
      end;
      for j=1:numt
         svpt(up,j)=(rand<ptprob(j));
      end;
      p=svp(up,:);
      p=p';
      pt=svpt(up,:);
      pt=pt';
      
      % now evaluate solutions ...
      % first convert PBIL probability strings into adjacency matrix g
      % interpret inner node tree from p
      g=zeros(m+n);g(2,1)=1;
      for i=3:m
         j=bin2dec(int2str(p(start(i):stop(i)))')+1;
         if j>=i
            j=1;
         end;
         g(i,j)=1;
      end;
      
      % interpret terminal nodes from pt
      for i=1:n
         j=bin2dec(int2str(pt((i-1)*ceil(log2(m))+1:i*ceil(log2(m))))')+1;
         if j>m
            %j=ceil(rand*m);
            j=1;
         end;
         g(i+m,j)=1;
      end;
      
      % find edges in unique path between every pair
      % find parent list of each terminal node
      pl=zeros(n,m);
      for i=1:n
         current=find(g(m+i,:)==1);
         pl(i,current)=1;
         while current~=1
            current=find(g(current,:)==1);
            pl(i,current)=1;
         end;
      end;
      
      % object pairs x weight matrix for non-negative least squares
      wm=zeros(npairs,m+n-1);
      c=0;
      for i=1:n-1
         for j=i+1:n
            c=c+1;
            pairmeet=max(intersect(find(pl(i,:)==1),find(pl(j,:)==1)));
            current=m+i;
            wm(c,current)=1;
            while find(g(current,:)==1)~=pairmeet
               current=find(g(current,:)==1);
               wm(c,current)=1;
            end;
            current=m+j;
            wm(c,current)=1;
            while find(g(current,:)==1)~=pairmeet
               current=find(g(current,:)==1);
               wm(c,current)=1;
            end;
         end;
      end;
      
      % find arc-lengths and evaluate sum-squared error 
      w=nnls(wm,flatd);
      evaluate(up)=sum((wm*w-flatd).^2);
      
      % if solution is best one found, reset trial counter,
      % save the solution and display a message
      % otherwise increment the trial counter
      if evaluate(up)<besterr
         tr=0;
         besterr=evaluate(up);
         bestvaf=1-besterr/vard;
         best=g;
         bestw=w;
         msg=sprintf('better tree found: accounts for %1.2f percent of the variance',bestvaf*100);
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
      
   end; %for up
   
   % sort the batch of evaluated solutions
   [evaluate order]=sort(evaluate);
   svp=svp(order,:);
   svpt=svpt(order,:);
   
   % use the best mupdate with competitive learning to adjust probability vectors
   for i=batchlearn:-1:1
      pprob=(1-lr)*pprob+svp(i,:)'*lr;
      ptprob=(1-lr)*ptprob+svpt(i,:)'*lr;
   end;
   
   % stochastically 'mutate' probability vectors using mutshift and mutprob
   for i=1:count
      if rand<mutprob
         if pprob(i)>.5
            pprob(i)=pprob(i)*(1-mutshift);
         else
            pprob(i)=pprob(i)*(1-mutshift)+mutshift;
         end;
      end;
   end;
   
   for i=1:numt
      if rand<mutprob
         if ptprob(i)>.5
            ptprob(i)=ptprob(i)*(1-mutshift);
         else
            ptprob(i)=ptprob(i)*(1-mutshift)+mutshift;
         end;
      end;
   end;
   
end; %while tr

% return the best adjacency matrix, its associated arc-lengths
% and percentage of variance it accounts for
adjacency=best+best';
lengths=bestw*rescale+reshift;
vaf=bestvaf;