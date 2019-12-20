function [Fourier]=fun_fft2d(X,N,Fs,f)
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
FourierTemp.Amplitude=Fourier.Amplitude(1:2:size(Fourier.Amplitude,1),1:2:size(Fourier.Amplitude,2));
FourierTemp.Phase=Fourier.Phase(1:2:size(Fourier.Phase,1),1:2:size(Fourier.Phase,2));
Fourier.Real.Amplitude=FourierTemp.Amplitude(:,size(FourierTemp.Amplitude,2)/2+1:size(FourierTemp.Amplitude,2));
Fourier.Real.Phase=FourierTemp.Phase(:,size(FourierTemp.Phase,2)/2+1:size(FourierTemp.Phase,2));
Fourier.Real.Frequency=Fourier.Force.Frequency(1:2:size(FourierTemp.Amplitude,2));
Fourier.Real.TimeOrder=Fourier.Real.Frequency/f;
Fourier.Real.SpaceOrder=(2*([1:size(Fourier.Real.Amplitude,1)]-(size(FourierTemp.Amplitude,1)/2+1)));
end