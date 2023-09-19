clc;
clear;
load('Sample_Process.mat');

%%%plot%%%%%%%%%%%%%%
M=input("enter number of sample functions wanted");
for i=1:M
    figure;
    plot(t,X(i,:));
    title(['Sample Function ',num2str(i)]);
end
%%plot(t,X(1:M,:));

%%ensemble mean%%%
enmean=mean(X);
figure;
plot(t,enmean);
title('Ensemble Mean');

%%%%time mean%%%%
N=input("enter the number of sample function you want to calculate time average to");
while N>M
fprintf('error, number must be less than %d \n',M);
N=input("enter the number of sample function you want to calculate time average to");
end
tmean=mean(transpose(X(N,:)))/0.1;
disp(tmean);

%%ensemble ACF%%%%%%%%%%
i=input('inset i');
j=input('inset j');
autoc=X(:,i).*X(:,j);
autocmean=mean(autoc);
figure;
plot(t,autocmean,".");

%%%time ACF%%%%
NN=input("enter the number of sample function you want to calculate time average to");

multi=0;
for tao=1:101
    for q=1:101
        if (q+tao)<102
        multi=multi+X(NN,q)*X(NN,q+tao);
       end
    end
    TimeACF(tao)=multi;
    multi=0;
end
disp(TimeACF);
figure;
plot(-t,TimeACF,t,TimeACF);

%%autocorr(X(NN),10);

%%%total average power%%%
% mult=0;
% for w=1:100
%     for q=1:101
%         mult=mult+X(w,q)*X(w,q);
%     end
%     M(w)=mult;
%     mult=0;
% end
% timeAVG=mean(transpose(M));
% disp("this is A");
% disp(timeAVG);

for v=1:100
    tvp(v)=sum(X(v,:).*X(v,:));
end
%%disp("this is B");
%%disp(mean(tvp))
pAVG=mean(tvp);


%%PSD%%%
psd=abs(fft(TimeACF,100)).^2/101;
figure,bar(psd);

%%%%%%%%%%
for i=1:100
    for j=1:100
        ACF(i,j)=mean(X(:,i).*X(:,j));
    end
end
disp (ACF)
figure;
x1=1:100;
x2=1:100;
surf(ACF)









