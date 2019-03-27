% full(gallery('tridiag',5,-1,1,1));
% full(gallery('tridiag',5,1,2,1));

%  x=-2*pi:0.1:2*pi;
% f1=sin(x);
% f2=sin(x)./x;
% plotyy(x,f1,x,f2)

% x1=[1:49]
% f1=1+5.*x1.^(-1)-3.*x1.^(-4)
% x2=[50:100]
% f2=x2.^(1/2)
% plotyy(x1,f1,x2,f2);

% function nr=divide(x,y)
% nr=0;
% for i=1:length(x)
%     if nnz(~mod(y,x(i)))~=0
%         nr=nr+1;
%     end
% end

% x=[3 4 6 9 8]
% y=[8 6 4]
% dim=0
% numel( intersect( x, y ) )
% for(i=1:length(x))
% if(numel(find(x(i)==y))~=0)
% dim=dim+1
% end
% end

%  elem=0;
% for i=1:10000
%  e1=randi(6,1,24)
%  e2=randi(6,1,24)
% % if(sum((e1==6).*(e2==6))~=0)
% %     elem=elem+1;
% % end
% % elem
% % end
% % frecv=elem/10000
% 
% 
% if(nnz(~(mod(e2,e1)))>0)
%     elem=elem+1;
% end
% end
% 
% frecv=elem/10000
% 

% x=input('');
% y=input('');
% if (sum(y)==1 && sum(y>=0)==length(x))
%     disp('media este: ')
%     sum(x.*y)
%     disp('dispersia este: ')
%     sum((x.^2.*y))-sum(x.*y)
% else
%     disp('eroare, nu este tabel de dispersie');
% end

% a=binornd(1,19,0.2,[1 19])
% sum(a==3)/length(a)
% sum(a<=3)/length(a)
% sum(a>=3)/length(a)

% p=poissrnd(17,1000,1)
% sum(p==5)/1000
% sum(p==10)/1000
% sum(p==15)/1000
% poisspdf([5 10 15],17)

% close all
% a=-2;
% b=6;
% m=2;
% sigma=1;
% x=(m-4*sigma):0.01:(4*sigma+m);
% lmb=0.5;
% 
% 
% plot(x,normpdf(x,m,sigma),'r',x,exppdf(x,lmb),'b',x,unifpdf(x,a,b),'g',x,normcdf(x,m,sigma),'--r',x,expcdf(x,lmb),'--b',x,unifcdf(x,a,b),'--g')
% 
% legend('normpdf','exppdf','unifpdf','normcdf','expcdf','unifcdf')


% m1=1;
% m2=2;
% s1=1;
% s2=1;
% r=0.3;
% 
% x1=(m1-4*s1):0.01:(m1+4*s1);
% x2=(m2-4*s2):0.01:(m2+4*s2);
% X1=repmat(x1,length(x1),1);
% X2=repmat(x2,length(x2),1);
% 
% f=(1/(1*pi*s1*s2)).*exp(((-1/(2*(1-r^2))).*(((X1-m1).^2)/s1^2)-(2*r.*((X1-m1).*(X2-m2))/(s1*s2))+((X2-m2).^2/s2^2)))
% mesh(f)


