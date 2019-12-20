function [Spectrum]=fun_fft1d(X,time)
%% 数据处理
Spectrum.xdata=time';
Spectrum.ydata=X';
N=length(Spectrum.xdata);% 采样点数
Fs=1/(Spectrum.xdata(2)-Spectrum.xdata(1));% 采样频率
%% FFT 分解
Spectrum.ComplexData=fft(Spectrum.ydata,[],1);
Spectrum.Amplitude=2*abs(Spectrum.ComplexData)/N;% 幅值转化
Spectrum.Amplitude(1,:)=Spectrum.Amplitude(1,:)/2;% 直流分量幅值转化
Spectrum.pu(:,1)=Spectrum.Amplitude(:,1)/Spectrum.Amplitude(1,1);
Spectrum.Frequency=(0:1:N/2)*Fs/N;% 频率转化,双边谱转化为单边谱
Spectrum.Phase=angle(Spectrum.ComplexData)*180/pi;% 相位求解
Spectrum.Order=Spectrum.Frequency/Spectrum.Frequency(2);% 阶次求解
end
