clc;
clear;
theta = unifrnd(0,2 *pi,1000);
omega_c = 3/4*pi;
A = 4;
t=[-4.99:0.01:5];
for i=1:length(t)
    Y=A*sin(omega_c *t(1,i)+theta)+0.5*A*cos(2*omega_c *t(1,i)+theta/3);
end

%%%plot%%%%%%%%%%%%%%
M=input("enter number of sample functions wanted");
for i=1:M
    figure;
    plot(t,Y(i,:));
    title(['Sample Function ',num2str(i)]);
end
%%plot(t,Y(1:M,:));

%%ensemble mean%%%
enmean=mean(Y);
figure;
plot(t,enmean);
title('Ensemble Mean');

%%%%time mean%%%%
N=input("enter the number of sample function you want to calculate time average to");
while N>M
fprintf('error, number must be less than %d \n',M);
N=input("enter the number of sample function you want to calculate time average to");
end
tmean=mean(transpose(Y(N,:)))/0.01;
disp(tmean);

%%ensemble ACF%%%%%%%%%%
i=input('inset i');
j=input('inset j');
autoc=Y(:,i).*Y(:,j);
autocmean=mean(autoc);
figure;
plot(t,autocmean,".");

%%%time ACF%%%%
NN=input("enter the number of sample function you want to calculate time average to");

multi=0;
for tao=1:length(t)
    for q=1:length(t)
        if (q+tao)<1001
        multi=multi+Y(NN,q)*Y(NN,q+tao);
       end
    end
    TimeACF(tao)=multi;
    multi=0;
end
%disp(TimeACF);
figure;
plot(t,TimeACF);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%
%%%%%%%
%%%%
%


%%autocorr(Y(NN),10);

%%%total average power%%%
% mult=0;
% for w=1:100
%     for q=1:101
%         mult=mult+Y(w,q)*Y(w,q);
%     end
%     M(w)=mult;
%     mult=0;
% end
% timeAVG=mean(transpose(M));
% disp("this is A");
% disp(timeAVG);

for v=1:1000
    tvp(v)=sum(Y(v,:).*Y(v,:));
end
%%disp("this is B");
%%disp(mean(tvp))
timeAVG=mean(tvp);


%%PSD%%%
psd=abs(fft(TimeACF,100)).^2/length(t);
figure,bar(psd);

%%%%%%%%%%
for i=1:1000
    for j=1:1000
        ACF(i,j)=mean(Y(:,i).*Y(:,j));
    end
end
%figure;
surf(ACF);
colorbar;









