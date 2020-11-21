function [points,vaf,bic]=mdsgrow(proximity,precision,metric,labels,evidence,patience)

% MDSGROW grows a metric MDS model (michael.d.lee@dsto.defence.gov.au)
% [points,vaf,bic]=mdsgrow(proximity,precision,metric,labels,evidence,patience))
% 
% PROXIMITY is an NxN symmetric matrix of pairwise proximities (required)
% PRECISION specifies the mean standard error of the proximities (required)
% METRIC specifies the Minkowskian distance metric (default=2)
% LABELS is a string array naming each object in the similarity matrix (numbered by default)
% EVIDENCE specifies the increase over the minimum BIC value for terminating (default=6)
% PATIENCE specifies the number of attempts made at each dimensionality (default=5)
%
% POINTS returns an Nx(BESTNUMBERDIMENSIONS) matrix giving the spatial representation
% VAF returns the variance of the similarity data accounted for by the generated models
% BIC returns the Bayesian Information Criteria for by the generated models

global A flatd r npairs;

% check the number of arguments
error(nargchk(2,6,nargin));

% check the similarity matrix
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
if nargin<3, metric=2; end;
if nargin<4, labels='1'; for i=2:n labels=char(labels,int2str(i)); end; end;
if nargin<5, evidence=6; end;
if nargin<6, patience=5; end;

% check the object labels
[check maxlabc]=size(labels);
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
d=proximity;
sig=precision;
labs=labels;
r=metric;

% matlab
warning off;
% matlab 5.3
options=optimset('display','off');

% normalise proximities to lie between 0 and 1, and scale precision
reshift=min(min(d));
d=d-reshift;
rescale=max(max(d));
d=d/rescale;
sig=sig/rescale;

% calculate the variance of  the similarity matrix
dbar=(sum(sum(d))-trace(d))/n/(n-1);
temp=(d-dbar*ones(n)).^2;
vard=.5*(sum(sum(temp))-trace(temp));

% express the similarity matrix as a column vector
flatd=[];
for i=1:n-1
   for j=i+1:n
      flatd=[flatd;d(i,j)];
   end;
end;

% setup differences between coordinates
A=zeros(npairs,n);
cc=1;
for i=1:n-1
   for j=i+1:n
      A(cc,i)=1;
      A(cc,j)=-1;
      cc=cc+1;
   end;
end;

% init other variables and constants
npairs=round(n*(n-1)/2);
m=1;
newbic=0;
oldbic=1*log(n*(n-1)/2)+vard/sig^2;
minbic=oldbic;
bic=[oldbic];
vaf=[0];


% keep adding clusters while the most recent model has a BIC
% at least evidence above the lowest BIC found
while newbic-minbic<evidence
   
   patience=0;
   
   if m>1
      oldbic=newbic;
      bestprevdim=bestthisdim;
   end;
   
   bestthisdimerr=vard;
   bestthisdim=zeros(n+1,m);
   
   % how many local minima for fixed dimensionality
   while patience<maxpatience
      
      
      % minimise the residual function and find residual
      if m==1
         x0=[rand(n+1,1)*.05];
      else
         x0=[bestprevdim rand(n+1,1)*.05];
      end;
      [x sse]=lsqnonlin('mdsresidualsc',x0,[],[],options);
      x(n+1,2:end)=0;
      
      % see if best solution
      %if sse<bestseenerr
      %   bestseenerr=sse;
      %   bestseenx=x;
      %end;
      
      % see if best this dim
      if sse<bestthisdimerr
         bestthisdimerr=sse;
         bestthisdim=x;
         disp('*');
      else
         disp('.');
      end;
      
      patience=patience+1;
      
   end; % while patience
   
   % display coordinate of solution
   msg=strcat('Dimensions=',num2str(m),', Variance Explained=',num2str(100*(1-bestthisdimerr/vard)));
   disp(msg);
   disp(bestthisdim);
   %for i=1:n
   %   msg=strcat(labs(i,:),':');
   %   for j=1:m
   %         msg=strcat(msg,' ',num2str(bestthisdim(i,j),3),',');
   %      end;
   %      msg=strcat(msg(1:end-1));
   %   disp(msg);
   %end;
   %msg=strcat('additive constant:',num2str(bestthisdim(n+1,1),3));
   %disp(msg);
      
   
   % calculate new BIC and VAF
   newbic=(m*n+1)*log(n*(n-1)/2)+bestthisdimerr/sig^2;
   bic=[bic newbic];
   vaf=[vaf 1-bestthisdimerr/vard];
   
   % see if model has best BIC
   if newbic<minbic
      minbic=newbic;
      fbx=bestthisdim;
      fbm=m;
   end;
   
   % draw results
   figure(1);clf;hold on;
   plot([0:m],bic,'k-','linewidth',1);
   rng=get(gca,'ylim');
   %axis([-1 m+1 rng(1) rng(2)]);
   axis([-1 m+1 minbic-2 minbic+2*evidence]);
   xlabel('Number of Dimensions','fontsize',14,'fontweight','bold');
   ylabel('Bayesian Information Criterion','fontsize',14,'fontweight','bold');
   ax1=gca;
   set(ax1,'xtick',[0:m]);
   ax2=axes('Position',get(ax1,'Position'),'YAxisLocation','right','Color','none');
   axis([-1 m+1 0 100]);
   set(ax2,'xtick',[0:m]);
   p=line([0:m],100*vaf,'Color','k','LineStyle','--','Parent',ax2,'linewidth',1);
   ylabel('Percentage Variance Accounted For','fontsize',14,'fontweight','bold');
   pause(.01);
   
   % show configuration
   figure;clf;hold on;
   switch m
   case 1
      % draw nodes
      % may want to adjust the markersize for some labels
      labc=0;
      for i=1:n
         plot(0,bestthisdim(i,1),'ko','markersize',8+labc*4,'linewidth',1,...
            'markerfacecolor','w');
      end;
      for i=1:n
         if labc==0
            text(i,bestthisdim(i,1),deblank(labs(i,:)),...
               'horizontalalignment','left',...
               'fontname','arial','fontsize',8,'fontweight','normal');
         else
            text(rand+.1,bestthisdim(i,1),deblank(labs(i,1:labc)),...
               'horizontalalignment','center',...
               'fontname','arial','fontsize',8,'fontweight','normal');
         end;
      end;
       set(gca,'xlim',[-1 n+1]);
  case 2
      % draw nodes
      % may want to adjust the markersize for some labels
      labc=min(maxlabc,3);
      for i=1:n
         plot(bestthisdim(i,1),bestthisdim(i,2),'ko','markersize',8+labc*4,'linewidth',1,...
            'markerfacecolor','w');
      end;
      for i=1:n
         if labc==0
            text(bestthisdim(i,1),bestthisdim(i,2),deblank(labs(i,:)),...
               'horizontalalignment','center',...
               'fontname','arial','fontsize',8,'fontweight','normal');
         else
            text(bestthisdim(i,1),bestthisdim(i,2),deblank(labs(i,1:labc)),...
               'horizontalalignment','center',...
               'fontname','arial','fontsize',8,'fontweight','normal');
         end;
      end;
   otherwise
      labc=0;
         plotmatrix(bestthisdim);
      end; % switch
      
      % format axis and figure border
   set(gca,'box','on','xticklabel','none','xtick',[],...
      'yticklabel','none','ytick',[]);
   axis square;
   axis(1.1*axis);
   
   % display long label names if characters>0
   if labc>0
      text(max(get(gca,'xlim'))*1.1,max(get(gca,'ylim'))*.9,[upper(labs(:,1:labc)) labs(:,labc+1:end)],...
         'verticalalignment','top');
   end;
   
   pause(.01);
   
   % now ready to add a dimension
   m=m+1;
   
end; % while newbic-oldbic

% return variables
points=fbx*rescale+reshift;
