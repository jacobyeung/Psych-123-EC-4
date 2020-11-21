function [points,corr]=displaytree(distance,adjacency,lengths,labels,characters,preserveweight,iterations,learnrate)

% DISPLAYTREE displays an additive tree (michael.d.lee@dsto.defence.gov.au)
% [points,corr]=displaytree(distance,adjacency,lengths,labels,characters,preserveweight,iterations,learnrate)
% 
% DISTANCE is an NxN symmetric matrix of pairwise distances or proximities (required)
% ADJACENCY is an (N+NODES)x(N+NODES) adjacency matrix defining the tree topology, as reutned by addtree
% LENGTHS returns a vector of length (N+NODES) containing the arc-lengths for the tree, as reutned by addtree
% LABELS is a string array of containing N labels for the terminal nodes (required)
% CHARACTERS specifies how many characters to include in the display (default=0 displays the entire label)
% PRESERVEWEIGHT specifies the relative emphasis given to preserving lengths in the tree (default=10)
% ITERATIONS specifies the number of optimisation iterations performed (default=50)
% LEARNRATE specifies the learning rate used in optimisation (default=0.05)
%
% POINTS returns the coordinate locations of each of the internal and terminal nodes (optional)
% CORR returns the correlation between the specified lengths and the displayed lengths (optional)

% check the number of arguments
error(nargchk(4,8,nargin));

% check the distance matrix
[n check]=size(distance);
if check~=n
   error('distance matrix must be square');
end;
if ~isequal(distance,distance')
   error('distance matrix must be symmetric');
end;

% check the adjacency matrix
[tot check]=size(adjacency);
if check~=tot
   error('adjacency matrix must be square');
end;
if ~isequal(adjacency,adjacency')
   error('adjacency matrix must be symmetric');
end;
if (length(find(adjacency==1))+length(find(adjacency==0)))~=prod(size(adjacency))
   error('adjacency matrix must be contain only 0 and 1');
end;

% check the lengths
check=length(lengths);
if check~=tot
   error('number of arc lengths must match adjacency matrix');
end;
if sum(lengths>=0)<tot
   error('arc lengths must be positive');
end;

% check the labels
if ~isstr(labels)
   error('labels must be a string array');
end;
[check labelwidth]=size(labels);
if check~=n
   error('labels must contain one string for each object');
end;

% set default arguments as necessary
if nargin<8, learnrate=0.05; end;
if nargin<7, iterations=50; end;
if nargin<6, preserveweight=10; end;
if nargin<5, characters=0; end;

% check the number of characters to be displayed
if (characters<0)|(characters>labelwidth)|(characters~=round(characters))
   error('characters must be a non-negative integer');
end;

% check the number of iterations
if (iterations<1)|(iterations~=round(iterations))
   error('number of iterations must be a positive integer');
end;

% learnrate and preserverate are expected to be positive numbers
% but these constraints are not explicitly imposed

% assign shorter argument names
d=distance;
lr=learnrate;
w=lengths;
labs=labels;
labc=characters;
g=adjacency;

% normalise distances to lie between 0 and 1
reshift=min(min(d));
d=d-reshift;
rescale=max(max(d));
d=d/rescale;

% and keep the lengths on the same scale
w=w/rescale;

% infer number of internal nodes
m=tot-n;

% append distance matrix to include internal nodes
% and set a minimum separation of terminals
d=max(.05,d);
d=[zeros(m),zeros(m,n);zeros(n,m),d];
[x y]=find(tril(g)==1);
for i=1:length(x)
   d(x(i),y(i))=w(x(i));
   d(y(i),x(i))=w(x(i));
end;

% initialise variables
its=0;
p=rand(m+n,2)*.01-.005;
dh=zeros(m+n);
% main optimisation loop
while (its<iterations)
   its=its+1;
   % select pinning order
   r=randperm(m+n);
   for i=1:m+n
      cc=r(i);
      for j=1:m+n
         if ~(cc==j|d(cc,j)==0)
            % update estimated distance matrix
            dh(cc,j)=norm(p(cc,:)-p(j,:));
            dh(j,cc)=dh(cc,j);
            % if not terminal-terminal pair keep the learning rate
            if ~((cc>m)&(j>m)) 
               lrt=lr;
               % otherwise lower it according to preserverate
            else
               lrt=1/preserveweight*lr;
            end;
            % update display locations
            for k=1:2
               p(j,k)=p(j,k)-lrt*(dh(cc,j)-d(cc,j))/dh(cc,j)*(p(j,k)-p(cc,k));
            end;
         end;
      end;
   end;
end;

% rescale then centre the coordinates
p=p*rescale+reshift;
p(:,1)=p(:,1)-p(1,1);
p(:,2)=p(:,2)-p(1,2);

% return points and correlation, as required
if nargout>1
   % check and checkh hold the actual and displayed arc lengths
   % used to calculate the corr measure
   check=[];
   checkh=[];
   for i=1:m+n
      for j=1:m+n
         if ~(i==j|d(i,j)==0|((i>m)&(j>m)))
            check=[check;d(i,j)];
            checkh=[checkh;dh(i,j)];
         end;
      end;
   end;
   temp=corrcoef(check,checkh);
   corr=temp(1,2);
end;
if nargout>0, points=p; end;

% draw the tree
% set figure window
figure(1);clf;hold on;

% draw arcs
for i=1:m+n
   j=find(g(i,:)==1);
   j=j(1);
   plot([p(i,1) p(j,1)],[p(i,2) p(j,2)],'k-','linewidth',1);
end;

% draw internal nodes
for i=1:m
   plot(p(i,1),p(i,2),'ko','markersize',8,'linewidth',2,...
      'markerfacecolor',[1 1 1],'markeredgecolor',[0 0 0]);
end;

% draw terminal nodes
% may want to adjust the markersize for some labels
for i=1:n
   plot(p(i+m,1),p(i+m,2),'ko','markersize',12+labc*4,'linewidth',1,...
      'markerfacecolor',[1 1 1],'markeredgecolor',[0 0 0]);
end;

% label terminal nodes
for i=1:n
   if labc==0
      text(p(i+m,1),p(i+m,2),deblank(labs(i,:)),...
         'horizontalalignment','center',...
         'fontname','arial','fontsize',8,'fontweight','bold');
   else
      text(p(i+m,1),p(i+m,2),deblank(labs(i,1:labc)),...
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

