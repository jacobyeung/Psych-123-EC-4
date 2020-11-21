function displayrepn_b(clusters,weights,vaf,rho,labs,m,gcc,gc)

%displayrepn(clusters,weights,vaf,rho,labs,m,gcc,gc)

n=size(labs,1);
%m=size(clusters,2);

X=fopen('contrast_out.txt','a+');


% display cluster membership labels for best solution
disp(sprintf('\n'));
fprintf(X,'\n');
msg=strcat('Clusters=',num2str(m),', Variance Explained=',num2str(100*vaf),...
   ', GCC=',num2str(gcc),', GC=',num2str(gc));
disp(msg)
fwrite(X,msg);
fprintf(X,'\n');


[w,ind]=sort(weights(1:m));
w=flipud(w);
ind=flipud(ind);

for x=0:1
	for i=1:m
   	if rho(ind(i))==1;msg=['CF, ']; else; msg=['DF, ']; end
   	msg=[msg num2str(w(i),'%0.3f') ': '];
   	for j=1:n
   	   if clusters(j,ind(i))==1
   	      msg=[msg ' ' deblank(labs(j,:)) ','];
   	   end;
      end;
      if rho(ind(i))==x
         disp(msg(1:end-1))
         fwrite(X,msg);
			fprintf(X,'\n');
      end
	end
end

msg=['CF, ' num2str(weights(m+1),'%0.3f') ':  additive constant'];
disp(msg)
%disp(sprintf('\n'));
fwrite(X,msg);
fprintf(X,'\n');

fclose(X);


