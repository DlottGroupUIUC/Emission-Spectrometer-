function [ velocity, unc , out_velocities,out_FFT_intensity ] = LTFT_error_fn_v4( time,intensity,time_bound1,time_bound2,initial_speed )
%LTFT Provides the velocity and uncertainty from an interferogram over a
%longer time window

warning off

pi=find(time>time_bound1,1,'first');
pf=find(time>time_bound2,1,'first');

L=pf-pi;
NFFT=2^nextpow2(L*10);
Y=fft(intensity(pi:pf,1),NFFT)/L;
Yshort=2*abs(Y(1:NFFT/2+1));
YL=fft(intensity(pi:(pf-3*L/4)),NFFT)/(L/4);
YR=fft(intensity((pi+3*L/4):pf),NFFT)/(L/4);
f = (((.04)*10^-9)^-1)/2*linspace(0,1,NFFT/2+1)*0.775*10^-9;

yfilt=Yshort(f>(initial_speed/2));
ffilt=f(f>(initial_speed/2));
yfilt=yfilt(ffilt<(initial_speed*1.5));
ffilt=ffilt(ffilt<(initial_speed*1.5));
gaus_fit=fit(ffilt',yfilt,'Gauss1');
parameters=coeffvalues(gaus_fit);

velocity=parameters(2);
unc=parameters(3)*(2^0.5);

out_velocities=ffilt;
out_FFT_intensity=yfilt;

warning on

end

