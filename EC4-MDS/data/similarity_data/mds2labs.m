function [points,vaf]=mds2labs(distance,dimensions,metric,labels,characters)

% MDS multidimensional scaling, version 2 with labels (michael.d.lee@dsto.defence.gov.au)
% [points,vaf]=mds2labs(distance,dimensions,metric,labels,characters)
% 
% DISTANCE is an NxN symmetric matrix of pairwise distances or proximities (required)
% DIMENSIONS specifies the required dimensionality of the coordinate representation (required)
% METRIC specifies the Minkowski distance metric operating in the space (default=2)
% LABELS is a string array of containing N labels for the terminal nodes (default uses numbers)
% CHARACTERS specifies how many characters to include in the display (default=0 displays the entire label)
%
% POINTS returns an NxDIMENSIONS matrix giving the derived coordinate locations
% VAF return the variance of the distance values accounted for by the solution

global A flatd r;

% check the number of arguments
error(nargchk(2,5,nargin));

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

% set default arguments as necessary
if nargin<3, metric=2; end;
if nargin<4, labels='1'; for i=2:n labels=char(labels,int2str(i)); end; end;
if nargin<5, characters=0; end;

% check the labels
if ~isstr(labels)
   error('labels must be a string array');
end;
[check labelwidth]=size(labels);
if check~=n
   error('labels must contain one string for each object');
end;

% check the number of characters to be displayed
if (characters<0)|(characters>labelwidth)|(characters~=round(characters))
   error('characters must be a non-negative integer, not greater than the longest label');
end;

% metric is expected to be a positive number
% but this constraints is not explicitly imposed

% assign shorter argument names
d=distance;
dim=dimensions;
r=metric;
labs=labels;
labc=characters;

% normalise distances to lie between 0 and 1
reshift=min(min(d));
d=d-reshift;
rescale=max(max(d));
d=d/rescale;

% calculate the variance of  the distance matrix
dbar=(sum(sum(d))-trace(d))/n/(n-1);
temp=(d-dbar*ones(n)).^2;
vard=.5*(sum(sum(temp))-trace(temp));

% initialise variables
options=zeros(18,1);
options(1)=-1; %display
options(2)=1e-1; %x precision
options(3)=1e-1; %residuals precision
npairs=round(n*(n-1)/2);

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

% setup target distances as vector
flatd=zeros(npairs,1);
cc=1;
for i=1:n-1
   for j=i+1:n
      flatd(cc)=d(i,j);
      cc=cc+1;
   end;
end;

% minimise the residual function and find residual
p0=rand(n,dim)*.05;
p=leastsq('mdsresiduals',p0,options);
resid=mdsresiduals(p);

% return the variance accounted for
% and coordinate location of the final solution
sse=sum(sum(resid.^2));
vaf=1-sse/vard;
points=p*rescale+reshift;

% draw the tree if it is a 2D space
if dim==2
   
   % set figure window
   figure;clf;hold on;
   
   % draw terminal nodes
   % may want to adjust the markersize for some labels
   for i=1:n
      plot(p(i,1),p(i,2),'ko','markersize',8+labc*4,'linewidth',1,...
         'markerfacecolor',[1 1 1],'markeredgecolor','w');
   end;
   
   % label terminal nodes
   for i=1:n
      if labc==0
         text(p(i,1),p(i,2),deblank(labs(i,:)),...
            'horizontalalignment','center',...
            'fontname','arial','fontsize',8,'fontweight','bold');
      else
         text(p(i,1),p(i,2),deblank(labs(i,1:labc)),...
            'horizontalalignment','center',...
            'fontname','arial','fontsize',8,'fontweight','normal');
      end;
   end;
   
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
   
end; % if dim==2

