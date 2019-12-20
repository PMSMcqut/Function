function [Fourier]=fun_fft2d(X,f0,T,TimeDomain)
%% 采样数据
N.Time=length(TimeDomain.Time);% 时间采样点数
Fs.Time=1/(T/N.Time);% 时间采样频率
N.Space=length(TimeDomain.Space);
Fs.Space=1/(TimeDomain.Space(2)-TimeDomain.Space(1));% 空间采样频率
%% 电磁力信号FFT分析
% 径向电磁力分解
Fourier.ComplexData=fft2(X,N.Space,N.Time);
Fourier.Amplitude=4*abs(Fourier.ComplexData)/(N.Space*N.Time);% 幅值转化
Fourier.Amplitude(1,:)=Fourier.Amplitude(1,:)/2;% 直流分量幅值转化
Fourier.Amplitude(:,1)=Fourier.Amplitude(:,1)/2;% 直流分量幅值转化
Fourier.Amplitude=fftshift(Fourier.Amplitude);
Fourier.ComplexData=fftshift(Fourier.ComplexData);
Fourier.Phase=angle(Fourier.ComplexData)*180/pi;% 相位求解

Fourier.Force.Frequency=(0:1:N.Time/2)*Fs.Time/N.Time;% 时间频率转化
Fourier.TimeOrder=Fourier.Force.Frequency/Fourier.Force.Frequency(2);% 时间阶次求解
Fourier.Force.SpaceOrder=(0:1:N.Space/2)*Fs.Space/N.Space;
Fourier.SpaceOrder=Fourier.Force.SpaceOrder/Fourier.Force.SpaceOrder(2);% 空间阶次求解

%% 图形处理
% FourierTemp.ComplexData=Fourier.ComplexData(1:2:size(Fourier.ComplexData,1),1:2:size(Fourier.ComplexData,2));
% 取出偶数阶次与偶数倍频
FourierTemp.Amplitude=Fourier.Amplitude(1:2:size(Fourier.Amplitude,1),1:2:size(Fourier.Amplitude,2));
FourierTemp.Phase=Fourier.Phase(1:2:size(Fourier.Phase,1),1:2:size(Fourier.Phase,2));
Fourier.Real.Frequency=Fourier.Force.Frequency(1:2:size(FourierTemp.Amplitude,2));
% Fourier.Real.ComplexData=FourierTemp.ComplexData(:,size(FourierTemp.ComplexData,2)/2+1:size(FourierTemp.ComplexData,2));
% 双边谱变成单边谱
Fourier.Real.Amplitude=FourierTemp.Amplitude(:,size(FourierTemp.Amplitude,2)/2+1:size(FourierTemp.Amplitude,2));
Fourier.Real.Phase=FourierTemp.Phase(:,size(FourierTemp.Phase,2)/2+1:size(FourierTemp.Phase,2));
Fourier.Real.TimeOrder=Fourier.Real.Frequency/f0;
Fourier.Real.SpaceOrder=(2*([1:size(Fourier.Real.Amplitude,1)]-(size(FourierTemp.Amplitude,1)/2+1)));
end