afunction [order,bayesfactors,bic,parameters]=bayesgroups(data1,data2)

% [order,bayesfactors,bic,parameters]=bayesgroups(data1,data2)
% version 1.0, 31 March 2002 (michael.lee@psychology.adelaide.edu.au)
%
% BAYESGROUPS uses Bayesian Information Criterion to estimate
% the Bayes Factors for four Gaussian models on two groups of data
% Model 1: Same Mean, Same Variance (1m1v)
% Model 2: Different Means, Same Variance (2m1v)
% Model 3: Same Mean, Different Variances (1m2v)
% Model 4: Different Means, Different Variances (2m2v)
% 
% DATA1 is a vector with the data values from the first group
% DATA2 is a vector with the data values from the second group
%
% The Figure window shows the best-fitting distributions under each model
% together with the raw data, and the models ordered in terms of their
% Bayes Factors
% ORDER list the models from most likely to least likely
% BAYESFACTORS is a 4x1 vector of the Bayes Factors for each model (in the
% order given by ORDER), in relation to the most likely model
% BIC is a 4x1 vector with the BIC values for each model, in order of likelihood
% PARAMETERS is structure giving the best-fitting parameters for each model

% check the number of arguments
error(nargchk(2,2,nargin));

% two data sets here
d1=data1;
d2=data2;

% get into columns, if needed
[r c]=size(d1);
if c==1
    d1=d1';
end;
[r c]=size(d2);
if c==1
    d2=d2';
end;

% set axis width, height
w=.4;
h=.3;

% squeeze data into middle x in plots, and plot curves with res points
x=3/4;
res=50;

% don't try to show bayesfactors above this
maxbf=100;

% calculate other bounds for symmetry
b=(1-2*w)/3;

% create figures, axis handles
figure(1);clf
ax(1)=axes('Units','normalized','Position',[b 1-h-b w h]);
ax(2)=axes('Units','normalized','Position',[2*b+w 1-h-b w h]);
ax(3)=axes('Units','normalized','Position',[b 1-2*h-2*b w h]);
ax(4)=axes('Units','normalized','Position',[2*b+w 1-2*h-2*b w h]);
axbic=axes('Units','normalized','Position',[b b 1-2*b 1-2*h-4*b]);

% best-fitting parameters, fit index, and bic for 4 models

% calculate lengths
n1=length(d1);
n2=length(d2);

% MODEL 1 (1M1V)
m1=mean([d1 d2]);
s1=std([d1 d2],1);
fit1=(n1+n2)*log(2*pi*s1^2)+1/s1^2*sum(([d1 d2]-m1).^2);
bic1=fit1+2*log(n1+n2);

% MODEL 2 (2M1V)
m21=mean(d1);
m22=mean(d2);
s2=sqrt((n1*std(d1,1)^2+n2*std(d2,1)^2)/(n1+n2));
fit2=(n1+n2)*log(2*pi*s2^2)+1/s2^2*sum((d1-m21).^2)+1/s2^2*sum((d2-m22).^2);
bic2=fit2+3*log(n1+n2);

%MODEL 3 (1M2V)
sx=sum(d1);
sy=sum(d2);
sx2=sum(d1.^2);
sy2=sum(d2.^2);
c=[(-n1^2*n2-n2^2*n1) (n1*sx*n2+2*n1^2*sy+n2*sy*n1+2*n2^2*sx) ...
        (-2*n1*sx*sy-n1^2*sy2-2*n2*sx*sy-n2^2*sx2) n1*sx*sy2+n2*sy*sx2];
r=roots(c);
for i=1:3
    m=r(i);
    if isreal(m)
    ss1=sum((d1-m).^2);
    ss2=sum((d2-m).^2);
    s31=sqrt(ss1/n1);
    s32=sqrt(ss2/n2);
    tryfit(i)=n1*log(2*pi*s31^2)+n2*log(2*pi*s32^2)+1/s31^2*sum((d1-m).^2)+1/s32^2*sum((d2-m).^2);
else
    tryfit(i)=inf;
end; % if
end; % for i
[val ind]=min(tryfit);
m3=r(ind);
ss1=sum((d1-m3).^2);
ss2=sum((d2-m3).^2);
s31=sqrt(ss1/n1);
s32=sqrt(ss2/n2);
fit3=n1*log(2*pi*s31^2)+n2*log(2*pi*s32^2)+1/s31^2*sum((d1-m3).^2)+1/s32^2*sum((d2-m3).^2);
bic3=fit3+3*log(n1+n2);

% MODEL 4 (2M2V)
m41=mean(d1);
m42=mean(d2);
s41=std(d1,1);
s42=std(d2,1);
fit4=n1*log(2*pi*s41^2)+n2*log(2*pi*s42^2)+1/s41^2*sum((d1-m41).^2)+1/s42^2*sum((d2-m42).^2);
bic4=fit4+4*log(n1+n2);

% draw the figure

% draw the data
% range and domain
r1=max([d1 d2])-min([d1 d2]);
lo=min([d1 d2])-(r1-r1*x)/2/x;
hi=max([d1 d2])+(r1-r1*x)/2/x;
dom=[lo:(hi-lo)/(res-1):hi];

% now the points
for i=1:4
    ph=plot(d1,0,'ko');
    set(ph,'markerfacecolor','w','markersize',6,'linewidth',1,'parent',ax(i));
    ph=plot(d2,-.1,'ko');
    set(ph,'markerfacecolor','k','markersize',6,'linewidth',1,'parent',ax(i));
set(ax(i),'xlim',[lo hi]);
set(ax(i),'ylim',[-.2 1.2],'ytick',[],'yticklabel','none');
set(ax(i),'box','on');
set(ax(i),'fontsize',14);
end;

% plot the fits
% model 1
c=1/sqrt(2*pi*s1^2)*exp(-((dom-m1).^2)/2/s1^2);
c=c/max(c);
ph=plot(dom,c,'k-');
set(ph,'linewidth',2,'parent',ax(1));
%set(ax(1),'xtick',[min([d1 d2]) m1 max([d1 d2])],'xticklabel',num2str([min([d1 d2]);m1;max([d1 d2])],'%1.2f'));
%set(ax(1),'xtick',[min([d1 d2])max([d1 d2])],'xticklabel',num2str([min([d1 d2]);max([d1 d2])],'%1.2f'));
set(ax(1),'xtick',[],'xticklabel','none');
set(get(ax(1),'title'),'string','1m1v','fontsize',14,'fontweight','bold','verticalalignment','middle');
% model 2
c1=1/sqrt(2*pi*s2^2)*exp(-((dom-m21).^2)/2/s2^2);
c2=1/sqrt(2*pi*s2^2)*exp(-((dom-m22).^2)/2/s2^2);
den=max([c1 c2]);
c1=c1/den;
c2=c2/den;
ph=plot(dom,c1,'k:');
set(ph,'linewidth',2,'parent',ax(2));
ph=plot(dom,c2,'k--');
set(ph,'linewidth',2,'parent',ax(2));
%set(ax(2),'xtick',[min([d1 d2]) m21 m22 max([d1 d2])],'xticklabel',num2str([min([d1 d2]);m21;m22;max([d1 d2])],'%1.2f'));
%set(ax(2),'xtick',[min([d1 d2])max([d1 d2])],'xticklabel',num2str([min([d1 d2]);max([d1 d2])],'%1.2f'));
set(ax(2),'xtick',[],'xticklabel','none');
set(get(ax(2),'title'),'string','2m1v','fontsize',14,'fontweight','bold','verticalalignment','middle');
% model 3
c1=1/sqrt(2*pi*s31^2)*exp(-((dom-m3).^2)/2/s31^2);
c2=1/sqrt(2*pi*s32^2)*exp(-((dom-m3).^2)/2/s32^2);
den=max([c1 c2]);
c1=c1/den;
c2=c2/den;
ph=plot(dom,c1,'k:');
set(ph,'linewidth',2,'parent',ax(3));
ph=plot(dom,c2,'k--');
set(ph,'linewidth',2,'parent',ax(3));
%set(ax(3),'xtick',[min([d1 d2]) m3 max([d1 d2])],'xticklabel',num2str([min([d1 d2]);m3;max([d1 d2])],'%1.2f'));
%set(ax(3),'xtick',[min([d1 d2])max([d1 d2])],'xticklabel',num2str([min([d1 d2]);max([d1 d2])],'%1.2f'));
set(ax(3),'xtick',[],'xticklabel','none');
set(get(ax(3),'title'),'string','1m2v','fontsize',14,'fontweight','bold','verticalalignment','middle');
% model 4
c1=1/sqrt(2*pi*s41^2)*exp(-((dom-m41).^2)/2/s41^2);
c2=1/sqrt(2*pi*s42^2)*exp(-((dom-m42).^2)/2/s42^2);
den=max([c1 c2]);
c1=c1/den;
c2=c2/den;
ph=plot(dom,c1,'k:');
set(ph,'linewidth',2,'parent',ax(4));
ph=plot(dom,c2,'k--');
set(ph,'linewidth',2,'parent',ax(4));
%set(ax(4),'xtick',[min([d1 d2]) m41 m42 max([d1 d2])],'xticklabel',num2str([min([d1 d2]);m41;m42;max([d1 d2])],'%1.2f'));
%set(ax(4),'xtick',[min([d1 d2])max([d1 d2])],'xticklabel',num2str([min([d1 d2]);max([d1 d2])],'%1.2f'));
set(ax(4),'xtick',[],'xticklabel','none');
set(get(ax(4),'title'),'string','2m2v','fontsize',14,'fontweight','bold','verticalalignment','middle');

% plot the odds ratios according to the bic
% sort the bic and draw the squares
bic=[bic1 bic2 bic3 bic4];
[val ind]=sort(bic);
odds=exp(diff(val)/2);
locs=[1;cumprod(odds)'];
bayesfactors=locs;
locs=locs(find(locs<=maxbf))
ph=plot(locs,.3,'ks');
set(ph,'markerfacecolor','k','linewidth',2,'markersize',8,'parent',axbic);
set(axbic,'xlim',[0 locs(end)+1],'ylim',[0 1]);
set(axbic,'xtick',locs,'xticklabel',num2str(locs,'%1.1f'),'ytick',[],'yticklabel','none');
% now the text labels
labs=char('1m1v','2m1v','1m2v','2m2v');
th=text(locs,.6*ones(length(locs),1),labs(ind(1:length(locs)),:));
set(th,'fontsize',14,'horizontalalignment','center','verticalalignment','middle');
set(axbic,'fontsize',14);
set(get(axbic,'title'),'string','Bayes Factors','fontsize',14,'fontweight','bold','verticalalignment','middle');


% stuff to return
order=labs(ind,:);
bic=val;
parameters.mean_1m1v=m1;
parameters.var_1m1v=s1^2;
parameters.mean1_2m1v=m21;
parameters.mean2_2m1v=m22;
parameters.var_2m1v=s2^2;
parameters.mean_1m2v=m3;
parameters.var1_1m2v=s31^2;
parameters.var2_1m2v=s32^2;
parameters.mean1_2m2v=m41;
parameters.mean2_2m2v=m42;
parameters.var1_2m2v=s41^2;
parameters.var2_2m2v=s42^2;

