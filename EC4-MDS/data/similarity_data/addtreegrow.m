function [adjacency,lengths,vaf,bic]=addtreegrow(proximity,precision,labels,evidence,patience)

% ADDTREEGROW grows an additive tree model (michael.d.lee@dsto.defence.gov.au)
% [adjacency,lengths,vaf,bic]=addtreegrow(proximity,precision,labels,evidence,patience)
% 
% PROXIMITY is an NxN symmetric matrix of pairwise proximities (required)
% PRECISION specifies the mean standard error of the similarities (required)
% LABELS is a string array naming each object in the similarity matrix (numbered by default)
% EVIDENCE specifies the increase over the minimum BIC value for terminating (default=6)
% PATIENCE specifies the number of local-minima restarts attempted without improvement (default=1)
%
% ADJACENCY returns an (N+BESTNUMBERNODES)x(N+BESTNUMBERNODES) adjacency matrix
% defining the tree topology using the 'best' number of internal nodes
% LENGTHS returns a vector of length (N+BESTNUMBERNODES) containing the arc-lengths
% for the tree using the 'best' number of internal nodes
% VAF returns the variance of the proximity data accounted for by the generated models
% BIC returns the Bayesian Information Criteria for by the generated models

% check the number of arguments
error(nargchk(2,5,nargin));

% check the proximity matrix
[n check]=size(proximity);
if check~=n
   error('proximity matrix must be square');
end;
if ~isequal(proximity,proximity')
   error('proximity matrix must be symmetric');
end;

% check the precision
if precision<=0
   error('precision must be positive');
end;

% set default arguments as necessary
if nargin<3, labels='1'; for i=2:n labels=char(labels,int2str(i)); end; end;
if nargin<4, evidence=6; end;
if nargin<5, patience=1; end;

% check the object labels
[check junk]=size(labels);
if check~=n
   error('number of labels must match size of matrix');
end;

% check the evidence
if evidence<=0
   error('evidence must be positive');
end;

% check the maximum number of trials
if (patience<1)|(patience~=round(patience))
   error('patience must be a positive integer');
end;

%rename variables
maxtries=patience;
d=proximity;
sigma=precision;
labs=labels;

% normalise proximities to lie between 0 and 1, and scale precision
reshift=min(min(d));
d=d-reshift;
rescale=max(max(d));
d=d/rescale;
sig=sig/rescale;

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

% matlab options
warning off;

% init other variables and constants
npairs=round(n*(n-1)/2);
maxtries=2;
evidence=10;
adjacency=[];
lengths=[];

% fit star tree
m=1
wm=zeros(npairs,n);
c=1;
for i=1:n-1
   for j=i+1:n
      wm(c,i)=1;
      wm(c,j)=1;
      c=c+1;
   end;
end;

% find arc-lengths and evaluate sum-squared error 
w=max(wm\flatd,0);
sse=sum((wm*w-flatd).^2);

% create adjacency matrix g from x
best=zeros(m+n);
for i=2:m+n
   best(i,1)=1;
end;

bestw=[0;w];

vaf=[0 1-sse/vard];
err=[vard sum((wm*w-flatd).^2)];
bic=err/sigma^2+([0:1]+n-1)*log(npairs);

minbic=min(bic);
currentbic=bic(end);
besterr=min(err);

if bic(2)<bic(1)
   adjacency=best+best';
   lengths=bestw;
end;

msg=sprintf('best tree found: accounts for %1.2f percent of the variance',vaf(2)*100);
disp(msg);
pause(.01);

% draw results
figure(1);clf;hold on;
plot([0:m],bic,'k-','linewidth',1);
rng=get(gca,'ylim');
axis([-1 m+1 rng(1) rng(2)]);
axis([-1 m+1 minbic-2 minbic+2*evidence]);
xlabel('Number of Nodes','fontsize',14,'fontweight','bold');
ylabel('Bayesian Information Criterion','fontsize',14,'fontweight','bold');
ax1=gca;
set(ax1,'xtick',[0:m]);
ax2=axes('Position',get(ax1,'Position'),'YAxisLocation','right','Color','none');
axis([-1 m+1 0 100]);
set(ax2,'xtick',[0:m]);
p=line([0:m],100*vaf,'Color','k','LineStyle','--','Parent',ax2,'linewidth',1);
ylabel('Percentage Variance Accounted For','fontsize',14,'fontweight','bold');
pause(.001);

displaytree(d,best+best',bestw,labs,0,.5,2,.1);
pause(.01);

while currentbic-minbic<evidence
   
   if m==1
      % setup two-internal-node tree
      m=m+1
      x=zeros(m+n,1);
      x(1)=0;
      x(2)=1;
      
      %for i=m+1:m+n
      %   x(i)=ceil(rand*m);
      %end;
      
      x([1:n]+m)=1;
      
      %basis=ceil(rand*n);
      %onenode=find(d(basis,:)<mean(mean(d)));
      %twonode=find(d(basis,:)>=mean(mean(d)));
      %x(onenode+m)=1;
      %x(twonode+m)=2;
      
      candidates=[1:n];
      parents=[1 2];
   else
      % average errors per internal node
      res=zeros(n);
      c=1;
      for j=1:n-1;
         for k=j+1:n
            res(j,k)=resid(c);
            c=c+1;
         end;
      end;
      res=res+res';
      for i=1:m
         stimlist=find(bestx==i)-m;
         stimlist=stimlist(find(stimlist>0));
         avsse(i)=sum(sum(res(stimlist,:)));
      end;
      % set internal topology
      [val ind]=max(avsse);
      %candidates=find(bestx(m+1:end)==ind);
      candidates=find(bestx==ind)-m;
      candidates=candidates(find(candidates>0));
      x=[bestx(1:m);ind;bestx(m+1:end)];
      labs(candidates,:)
      m=m+1
      parents=[ind m];
   end;
   
   
   foundbetter=0;
   
   % main hillclimbing loop
   failedtries=0;
   
   while failedtries<maxtries
      
      failedthisgo=1;
      list=randperm(length(candidates));
      listind=1;
      
      while listind<length(list)
         node=candidates(list(listind))+m;
         
         % swap parent
         oldx=x;
         if x(node)==parents(1)
            x(node)=parents(2);
         else
            x(node)=parents(1);
         end;
         
         % create adjacency matrix g from x
         g=zeros(m+n);
         for i=2:m+n
            g(i,x(i))=1;
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
         %[w sse]=lsqnonneg(wm,flatd,[],options);
         w=max(wm\flatd,0);
         resid=(wm*w-flatd).^2;
         sse=sum(resid);
         
         
         if sse<besterr
            failedthisgo=0;
            foundbetter=1;
            besterr=sse;
            failedtries=0;
            bestvaf=1-besterr/vard;
            bestx=x;
            best=g;
            bestw=w;
            bestresid=resid;
            msg=sprintf('better tree found: accounts for %1.2f percent of the variance',bestvaf*100);
            disp(msg);
            pause(.01);
         else
            x=oldx;
         end;
         
         listind=listind+1;
         
      end; % while listind
      
      if failedthisgo==1
         failedtries=failedtries+1;
         shake=candidates(ceil(rand*length(candidates)))+m;
         if x(shake)==parents(1)
            x(shake)=parents(2);
         else
            x(shake)=parents(1);
         end;
      end;
   end; %while failedtries
   
   if foundbetter==0
      bestx=[bestx(1:m);1;bestx(m+1:end)];
      bestw=[bestw(1:m);0;bestw(m+1:end)];
      best=zeros(m+n);
      for i=2:m+n
         best(i,bestx(i))=1;
      end;
   end;
   
   err=[err besterr];
   currentbic=besterr/sigma^2+(m+n-1)*log(npairs);
   bic=[bic currentbic];
   vaf=[vaf 1-besterr/vard];
   if currentbic<minbic
      minbic=currentbic;
      adjacency=best+best';
      lengths=bestw;
   end;
   
   % draw results
   figure(1);clf;hold on;
   plot([0:m],bic,'k-','linewidth',1);
   rng=get(gca,'ylim');
   axis([-1 m+1 rng(1) rng(2)]);
   axis([-1 m+1 minbic-2 minbic+2*evidence]);
   xlabel('Number of Clusters','fontsize',14,'fontweight','bold');
   ylabel('Bayesian Information Criterion','fontsize',14,'fontweight','bold');
   ax1=gca;
   set(ax1,'xtick',[0:m]);
   ax2=axes('Position',get(ax1,'Position'),'YAxisLocation','right','Color','none');
   axis([-1 m+1 0 100]);
   set(ax2,'xtick',[0:m]);
   p=line([0:m],100*vaf,'Color','k','LineStyle','--','Parent',ax2,'linewidth',1);
   ylabel('Percentage Variance Accounted For','fontsize',14,'fontweight','bold');
   pause(.001);

   displaytree(d,best+best',bestw,labs,0,.5,10,.1);
   pause(.01);
   
end; % while currentbic

