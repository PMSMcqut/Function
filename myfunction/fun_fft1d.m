function [Spectrum]=fun_fft1d(X,time)
%% ���ݴ���
Spectrum.xdata=time';
Spectrum.ydata=X';
N=length(Spectrum.xdata);% ��������
Fs=1/(Spectrum.xdata(2)-Spectrum.xdata(1));% ����Ƶ��
%% FFT �ֽ�
Spectrum.ComplexData=fft(Spectrum.ydata,[],1);
Spectrum.Amplitude=2*abs(Spectrum.ComplexData)/N;% ��ֵת��
Spectrum.Amplitude(1,:)=Spectrum.Amplitude(1,:)/2;% ֱ��������ֵת��
Spectrum.pu(:,1)=Spectrum.Amplitude(:,1)/Spectrum.Amplitude(1,1);
Spectrum.Frequency=(0:1:N/2)*Fs/N;% Ƶ��ת��,˫����ת��Ϊ������
Spectrum.Phase=angle(Spectrum.ComplexData)*180/pi;% ��λ���
Spectrum.Order=Spectrum.Frequency/Spectrum.Frequency(2);% �״����
end
