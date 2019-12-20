clear all;clc;
load work.txt;
l_s=length(work);
n=ceil(l_s/2);
for i=1:l_s
    t1(i)=work(i,1);
    f1(i)=work(i,2);
end
Fs=1/(t1(2)-t1(1));%采样频率
%subplot(3,1,1);
%plot(t1,f1);
%xlabel('Theta'); 
%ylabel('Voltage');
grid
ps=fft(f1,l_s);
A=2*abs(ps)/l_s;
A(1)=A(1)/2;
%画出频域图
tt=Fs*(0:n-1)/l_s*1000;
tt1=tt/tt(2);
%subplot(3,1,2);
%%subplot(3,1,3);
%bar(tt1,A(1:n),0.5);
%pu=A/max(A);
bar(A(2:28));
%set(gca,'xtick',[0:1:30]);
title('空载反电势');
xlabel('次数'); 
ylabel('幅值 T');
grid
%可以将tt和A的数据保存，他们就是频谱分析的结果，分别是频率和幅值
sum = 0;
for n=1:25
    sum = sum+A(n)*A(n);
end
q =sqrt(sum/2);

sum1 = 0;
for n=1:25
    sum1 = sum1+A(n)*A(n);
end
    sum1 = sum1-A(4)*A(4);
q1 =sqrt(sum1/2);