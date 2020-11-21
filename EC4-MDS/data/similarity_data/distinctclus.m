function [clusters,weights,vaf,gcc,complexity,rho]=distinctclus(similarity,precision,labels,evidence,patience)

% [clusters,weights,vaf,gcc,gc,rho]=distinctclus(similarity,precision,labels,evidence,patience)
 
% rename variables
maxpatience=patience;
s=similarity;
sig=precision;
labs=labels;

n=size(s,1);

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
newgcc=0;
oldgcc=.5*(vard/sig^2 + 1*log(n*(n-1)/(4*pi*sig^2)) + log(npairs));
mingcc=oldgcc;
fbx=[];
fbw=[];
fbm=0;
bic=[oldbic];
gcc=[oldgcc];
fform=[log(npairs)];
complexity=[.5*(1*log(n*(n-1)/(4*pi*sig^2)) + log(npairs))];
vaf=[0];
%matlab 5.3
%options=optimset('display','off');

thresh=.5;


% keep adding clusters while the most recent model has a GCC
% at least evidence above the lowest GCC found
%while newbic-minbic<evidence
while newgcc-mingcc<evidence
   
   if m>1
      oldgcc=newgcc;
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
      x=[x;round(rand)]; % add rho
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
      x=[bestseenx(1:n*(m-1));newx;bestseenx(n*(m-1)+1:end);round(rand)]; % add new cluster & rhos
   end;
      
   % randomly choose starting vector
   %if m==1
   %   x=round(rand(n+1,1));
   %else
   %   newx=round(rand(n,1));
   %   x=[bestseenx(1:n*(m-1));newx;bestseenx(n*(m-1)+1:end);round(rand)]; % add new cluster & rhos
   %end
   
   % reset/recalculate for new cardinality
   patience=0;
   veclength=(n+1)*m;
   bestseenx=x;
   bestw=zeros(m+1,1);
   bestseenw=zeros(m,1);
   bestseenf=zeros(n,m);
   bestseenerr=inf;
   bestseengcc=inf;
   bestseenc=inf;
   
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
            f=reshape(x(1:n*m),n,m);
            f=[f ones(n,1)];
            
            % obtain rho
            rho=x(n*m+1:end)';
            
            % HACK FOR ADCLUS OR DFCLUS 
            %rho=ones(1,m); % adclus
            rho=zeros(1,m); % DFclus
            
            % form of cluster membership needed for non-negative least squares
            flatf=zeros(npairs,m+1);
            count=0;
            for a=1:n-1
              	for b=a+1:n
              		count=count+1;
                	flatf(count,:)=f(a,:).*f(b,:).*[rho 1] - f(a,:).*(1-f(b,:)).*(1-[rho 1]).*.5 - (1-f(a,:)).*f(b,:).*(1-[rho 1]).*.5;
              	end;
            end;
            
            % FIND WEIGHTS (nnls, lsqnonneg and slash versions)
            
            %w=diag(nnls(flatf,flats));
            
            %[w sse]=lsqnonneg(flatf,flats,[],options);
            %w=diag(w);
            
            [state,freq]=warning;
         	warning off
            w=diag(max(flatf\flats,0));
            wCF=w*diag([rho 1]);
            wDF=w*diag(1-[rho 1]);
          	eval(['warning ' state]);eval(['warning ' freq])
            
            % FIND SSE
            sh=f*wCF*f'-.5*f*wDF*(1-f)'-.5*(1-f)*wDF*f';
            se=(s-sh).^2;
            sse=.5*(sum(sum(se))-trace(se));
            
            % FIND G
            
            r=[rho 1]'; % we need rho as a column vector, and need to add additive constant
            
            % If both feaures x & y are CF (given by the .*(r*r') term at the end)...
				a1=f'*f; % a1 counts the number of stimuli that have both x and y 
				g1=a1.*(a1-1)/2; % how many pairs of stimuli have both?
				g1=g1.*(r*r');

				% If x is CF and y is DF (given by the .*(r*(1-r)) term at the end)...
				a2=f'*(1-f); % a2 counts the number of stimuli that have x and not y
				g2=a1.*a2; % how many pairs are there such that both have x and only one has y?
				g2=-.5*g2.*(r*(1-r'));
	
				% If x is DF and y is CF (given by the .*((1-r)*r') term at the end)...
				a3=(1-f')*f; % a3 counts the number of stimuli that have y and not x
				g3=a1.*a3; % how many pairs are there such that both have y and only one has x? 
				g3=-.5*g3.*((1-r)*r');
	
				% If both x and y are DF (given by the .*((1-r)*(1-r')) term at the end)...
				a4=(1-f')*(1-f); % a4 counts the number of stimuli that have neither x nor y
				g4=a1.*a4 + a2.*a3; % how many pairs s.t. either (a) one has both and the other neither, or (b) each has one?
				g4=.25*g4.*((1-r)*(1-r'));

				% find complexity matrix
				G=g1+g2+g3+g4;
            
         	% avoid degenerate F
         	temp=det(G);
         	if temp<=0 % matrices that definitely correspond to degenerate solutions
         	   fform=inf; 
         	else
         	   fform=log(temp);
         	   if temp<1e-1; % matrices that probably correspond to degenerate solutions
         	      fform=inf;
         	      warning('Potential rank deficiency. Skip.')
         	   end
            end
                     
         	% find GCC
         	gc=.5*(sse/sig^2 + (m+1)*log(n*(n-1)/(4*pi*sig^2)) + fform);
            c=.5*((m+1)*log(n*(n-1)/(4*pi*sig^2)) + fform);
            
            % see if best solution
            %if sse<bestseenerr %& gc<10e10
            if gc<bestseengcc
               patience=0;
               bestseenerr=sse;
               bestseengcc=gc;
               bestseenc=c;
               bestseenvaf=1-sse/vard;
               bestseenx=x;
               bestseenf=f;
               bestseenfform=fform;
               bestseenrho=rho;
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
      nflip=min(length(x),5);
      %nflip=round(length(x)/3);
      x(scrambleflip(1:nflip))=1-x(scrambleflip(1:nflip));
      patience=patience+1;
      
   end; % while patience
      
   % calculate new GCC and VAF
   newbic=(m+1)*log(n*(n-1)/2)+bestseenerr/sig^2;
   bic=[bic newbic];
	newgcc=bestseengcc;
   gcc=[gcc bestseengcc];
   complexity=[complexity bestseenc];
   vaf=[vaf 1-bestseenerr/vard];
   fform=[fform bestseenfform];
   
   % see if model has best BIC
   if newbic<minbic
      minbic=newbic;
   %   fbx=bestseenx;
   %   fbw=bestseenw;
   %   fbm=m;
   end;
   
   % see if model has best GCC
   if newgcc<mingcc
      mingcc=newgcc;
      fbx=bestseenx;
      fbw=bestseenw;
      fbm=m;
   end;
   
   % draw results
   figure(1);clf;hold on;
   plot([0:m],gcc,'k-','linewidth',1);
   rng=get(gca,'ylim');
   %axis([-1 m+1 rng(1) rng(2)]);
   axis([-1 m+1 mingcc-2 mingcc+2*evidence]);
   xlabel('Number of Clusters','fontsize',14,'fontweight','bold');
   ylabel('Geometric Complexity Criterion','fontsize',14,'fontweight','bold');
   ax1=gca;
   set(ax1,'xtick',[0:m]);
   
   ax2=axes('Position',get(ax1,'Position'),'YAxisLocation','right','Color','none');
   axis([-1 m+1 0 100]);
   set(ax2,'xtick',[0:m]);
   p=line([0:m],100*vaf,'Color','k','LineStyle','--','Parent',ax2,'linewidth',1);
   ylabel('Percentage Variance Accounted For','fontsize',14,'fontweight','bold');
   
   ax3=axes('Position',get(ax1,'Position'),'YAxisLocation','left','Color','none');
   p=line([0:m],bic,'Color','k','LineStyle',':','Parent',ax3,'linewidth',1);
   axis([-1 m+1 minbic-2 minbic+2*evidence]);
   set(gca,'ytick',[]);
   %ylabel('Bayesian Information Criterion','fontsize',14,'fontweight','bold');
   
   % results in command window
   displayrepn_d(bestseenf,bestseenw,bestseenvaf,bestseenrho,labs,m,bestseengcc,bestseenc)
   drawnow

   pause(.001);
 
   % now ready to add a cluster
   m=m+1;
   
end; % while new-old

% return variables
if fbm>0
clusters=[reshape(fbx(1:fbm*n),n,fbm) ones(n,1)];
weights=fbw;
rho=fbx(fbm*n+1:end);
else
   clusters=[];
   weights=[];
   rho=[];
end;
