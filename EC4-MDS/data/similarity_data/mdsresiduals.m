function residuals = mds(x);
global A flatd r;
residuals=(sum(abs((A*x)).^r,2)).^(1/r)-flatd;
