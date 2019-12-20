function [Fourier]=fun_fft2d(X,f0,T,TimeDomain)
%% ��������
N.Time=length(TimeDomain.Time);% ʱ���������
Fs.Time=1/(T/N.Time);% ʱ�����Ƶ��
N.Space=length(TimeDomain.Space);
Fs.Space=1/(TimeDomain.Space(2)-TimeDomain.Space(1));% �ռ����Ƶ��
%% ������ź�FFT����
% ���������ֽ�
Fourier.ComplexData=fft2(X,N.Space,N.Time);
Fourier.Amplitude=4*abs(Fourier.ComplexData)/(N.Space*N.Time);% ��ֵת��
Fourier.Amplitude(1,:)=Fourier.Amplitude(1,:)/2;% ֱ��������ֵת��
Fourier.Amplitude(:,1)=Fourier.Amplitude(:,1)/2;% ֱ��������ֵת��
Fourier.Amplitude=fftshift(Fourier.Amplitude);
Fourier.ComplexData=fftshift(Fourier.ComplexData);
Fourier.Phase=angle(Fourier.ComplexData)*180/pi;% ��λ���

Fourier.Force.Frequency=(0:1:N.Time/2)*Fs.Time/N.Time;% ʱ��Ƶ��ת��
Fourier.TimeOrder=Fourier.Force.Frequency/Fourier.Force.Frequency(2);% ʱ��״����
Fourier.Force.SpaceOrder=(0:1:N.Space/2)*Fs.Space/N.Space;
Fourier.SpaceOrder=Fourier.Force.SpaceOrder/Fourier.Force.SpaceOrder(2);% �ռ�״����

%% ͼ�δ���
% FourierTemp.ComplexData=Fourier.ComplexData(1:2:size(Fourier.ComplexData,1),1:2:size(Fourier.ComplexData,2));
% ȡ��ż���״���ż����Ƶ
FourierTemp.Amplitude=Fourier.Amplitude(1:2:size(Fourier.Amplitude,1),1:2:size(Fourier.Amplitude,2));
FourierTemp.Phase=Fourier.Phase(1:2:size(Fourier.Phase,1),1:2:size(Fourier.Phase,2));
Fourier.Real.Frequency=Fourier.Force.Frequency(1:2:size(FourierTemp.Amplitude,2));
% Fourier.Real.ComplexData=FourierTemp.ComplexData(:,size(FourierTemp.ComplexData,2)/2+1:size(FourierTemp.ComplexData,2));
% ˫���ױ�ɵ�����
Fourier.Real.Amplitude=FourierTemp.Amplitude(:,size(FourierTemp.Amplitude,2)/2+1:size(FourierTemp.Amplitude,2));
Fourier.Real.Phase=FourierTemp.Phase(:,size(FourierTemp.Phase,2)/2+1:size(FourierTemp.Phase,2));
Fourier.Real.TimeOrder=Fourier.Real.Frequency/f0;
Fourier.Real.SpaceOrder=(2*([1:size(Fourier.Real.Amplitude,1)]-(size(FourierTemp.Amplitude,1)/2+1)));
end