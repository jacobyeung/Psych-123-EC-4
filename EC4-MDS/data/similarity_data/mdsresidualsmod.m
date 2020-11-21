function residuals = mds(x);
global A flatd r n dim;
x=[zeros(3,1);x];
x=reshape(x,dim,n);
x=x';
residuals=(sum(abs((A*x)).^r,2)).^(1/r)-flatd;
