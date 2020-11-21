function [clusters,weights,vaf,bic]=adclusgrow(similarity,precision,labels,evidence,patience)

% ADCLUSGROW grows an additive clustering model (michael.d.lee@dsto.defence.gov.au)
% [clusters,weights,vaf,bic]=adclusgrow(similarity,precision,labels,evidence,patience)
% 
% SIMILARITY is an NxN symmetric matrix of pairwise similarities (required)
% PRECISION specifies the mean standard error of the similarities (required)
% LABELS is a string array naming each object in the similarity matrix (numbered by default)
% EVIDENCE specifies the increase over the minimum BIC value for terminating (default=6)
% PATIENCE specifies the number of local-minima restarts attempted without improvement (default=1)
%
% CLUSTERS returns an Nx(BESTNUMBERCLUSTERS+1) matrix defining derived cluster membership plus the universal cluster
% WEIGHTS returns a vector of length (BESTNUMBERCLUSTERS+1) containing the weights of the clusters
% VAF returns the variance of the similarity data accounted for by the generated models
% BIC returns the Bayesian Information Criteria for by the generated models

% check the number of arguments
error(nargchk(2,5,nargin));

% check the similarity matrix
[n check]=size(similarity);
if check~=n
   error('similarity matrix must be square');
end;
if ~isequal(similarity,similarity')
   error('similarity matrix must be symmetric');
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

% rename variables
maxpatience=patience;
s=similarity;
sig=precision;
labs=labels;

% matlab
warning off;

% normalise similarities to lie between 0 and 1, and scale precision
reshift=min(min(s));
s=s-reshift;
rescale=max(max(s));
s=s/rescale;
sig=sig/rescale;

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
m=1;
newbic=0;
oldbic=1*log(n*(n-1)/2)+vard/sig^2;
minbic=oldbic;
fbx=[];
fbw=[];
fbm=0;
bic=[oldbic];
vaf=[0];
%matlab 5.3
%options=optimset('display','off');

thresh=.5;


% keep adding clusters while the most recent model has a BIC
% at least evidence above the lowest BIC found
while newbic-minbic<evidence
   
   if m>1
      oldbic=newbic;
   end;
   
   % choose starting vector
   % uses ADDI-S idea of adding most (averagely) similar
   % until less than *half* within-cluster similarity
   if m==1
      % first pair
      [row col]=find(s==max(max(tril(s,-1)))); %tril(s,-1) is tril without diagonal
      ind=find(row~=col);
      row=row(ind);col=col(ind);
      ind=ceil(rand*length(row));row=row(ind);col=col(ind); %break ties randomly
      new=[row col];
      ins=(sum(sum(s(new,new)))-trace(s(new,new)))/length(new); %average within new sim
      avs=zeros(n,1);
      for i=1:n
         if nnz(ismember(new,i))==0
            avs(i)=mean(s(new,i)); %average similarity to new
         end;
      end;
      while nnz(avs)>0
         [val ind]=max(avs);
         if avs(ind)>thresh*ins
            new=[new ind];
            ins=(sum(sum(s(new,new)))-trace(s(new,new)))/length(new); %average within new sim
            avs(ind)=0;
         else
            avs(ind)=0;
         end;
      end;
      x=zeros(n,1);
      x(new)=1;
   else
      tmps=s-bestseensh;
      tmps=max(tmps,0);
      [row col]=find(tmps==max(max(tril(tmps,-1)))); %tril(s,-1) is tril without diagonal
      ind=find(row~=col);
      row=row(ind);col=col(ind);
      ind=ceil(rand*length(row));row=row(ind);col=col(ind); %break ties randomly
      new=[row col];
      ins=(sum(sum(tmps(new,new)))-trace(tmps(new,new)))/length(new); %average within new sim
      avs=zeros(n,1);
      for i=1:n
         if nnz(ismember(new,i))==0
            avs(i)=mean(tmps(new,i)); %average similarity to new
         end;
      end;
      while nnz(avs)>0
         [val ind]=max(avs);
         if avs(ind)>thresh*ins
            new=[new ind];
            ins=(sum(sum(tmps(new,new)))-trace(tmps(new,new)))/length(new); %average within new sim
            avs(ind)=0;
         else
            avs(ind)=0;
         end;
      end;
      newx=zeros(n,1);
      newx(new)=1;
      x=[bestseenx;newx];
   end;
   
   % reset/recalculate for new cardinality
   patience=0;
   veclength=n*m;
   bestseenx=x;
   bestw=zeros(m+1,1);
   bestseenerr=vard;
   
   % how many local minima for fixed cardinality
   while patience<maxpatience
      
      %  hill-climb on the surface
      alldone=0;
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
            
            % find weights (nnls, lsqnonneg and slash versions)
            
            %w=diag(nnls(flatf,flats));
            
            %[w sse]=lsqnonneg(flatf,flats,[],options);
            %w=diag(w);
            
            w=diag(max(flatf\flats,0));
            
            
            % evaluate solution
            sh=f*w*f';
            se=(s-sh).^2;
            sse=.5*(sum(sum(se))-trace(se));
            
            % see if best solution
            if sse<bestseenerr
               patience=0;
               bestseenerr=sse;
               bestseenx=x;
               bestseenw=diag(w)';
               bestseensh=sh;
               moveon=1;
            else
               x(fliptry(trycount))=1-x(fliptry(trycount));
            end;
            
            if trycount==veclength
               moveon=1;
               alldone=1;
            end;
         end; % move on   
         
      end; % all done
      
      % a random selection of 3 bits are flipped to try and shake out of a local min
      % there is no principled reason for choosing 3
      x=bestseenx;
      scrambleflip=randperm(veclength);
      x(scrambleflip(1:3))=1-x(scrambleflip(1:3));
      patience=patience+1;
      
   end; % while patience
   
   % display cluster membership labels for best solution
   msg=strcat('Clusters=',num2str(m),', Variance Explained=',num2str(100*(1-bestseenerr/vard)));
   disp(msg)
   for i=1:m
      msg=strcat(num2str(bestseenw(i),3),':');
      for j=1:n
         if bestseenx((i-1)*n+j)==1
            msg=strcat(msg,' ',(labs(j,:)),',');
         end;
      end;
      disp(msg(1:end-1))
   end;
   msg=strcat(num2str(bestseenw(m+1),3),':additive constant');
   disp(msg)
   
   % calculate new BIC and VAF
   newbic=(m+1)*log(n*(n-1)/2)+bestseenerr/sig^2;
   bic=[bic newbic];
   vaf=[vaf 1-bestseenerr/vard];
   
   % see if model has best BIC
   if newbic<minbic
      minbic=newbic;
      fbx=bestseenx;
      fbw=bestseenw;
      fbm=m;
   end;
   
   % draw results
   figure(1);clf;hold on;
   plot([0:m],bic,'k-','linewidth',1);
   rng=get(gca,'ylim');
   %axis([-1 m+1 rng(1) rng(2)]);
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
   
   % now ready to add a cluster
   m=m+1;
   
end; % while newbic-oldbic

% return variables
if fbm>0
clusters=reshape(fbx,n,fbm);
weights=fbw;
else
   clusters=[];
   weights=[];
end;
