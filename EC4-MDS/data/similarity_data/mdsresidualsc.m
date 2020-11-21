function residuals = mds(x);
global A flatd r npairs;
residuals=(sum(abs(A*x(1:end-1,:)).^r,2)).^(1/r)+ones(npairs,1)*abs(x(end,1))-flatd;
