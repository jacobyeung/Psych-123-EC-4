function residuals= mdsresiduals3(x);
global A flatd r n dim;
x=reshape(x,n,dim);
residuals=(sum(abs((A*x)).^r,2)).^(1/r)-flatd;
      
