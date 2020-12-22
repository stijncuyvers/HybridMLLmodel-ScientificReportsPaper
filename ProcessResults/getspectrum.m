% ===============================================================
%  FUNCTION: Get spectrum of pulse train
%  Release year: 2020
%  Author: Stijn Cuyvers
%  Affiliation: Photonics Research Group, Ghent University - imec, Belgium
% ===============================================================

function [  validdata,lambdas, spectrum, lambdas_envelope, envelope, lambdas_phase, phases ] = getspectrum( outputmonitor, lambda0, dt )
% based on time-domain trace (outputmonitor), calculate the spectrum 

outputmonitor_secondhalf=outputmonitor(round(0.5*length(outputmonitor)):length(outputmonitor));
minpeakheight=max(abs(outputmonitor_secondhalf))*0.3;
[peaks,peakpos]=findpeaks(abs(outputmonitor),'MinPeakHeight',minpeakheight);
n_o_peaks=length(peaks);


if n_o_peaks<40
    validdata=0;
    lambdas=-1;
    phases=-1;
    spectrum=-1;
    envelope=-1;
    lambdas_envelope=-1;
    lambdas_phase=-1;
else
    validdata=1;

pos_2=round(0.5*(peakpos(n_o_peaks)+peakpos(n_o_peaks-1)));
pos_1=round(0.5*(peakpos(round(0.5*n_o_peaks))+peakpos(round(0.5*n_o_peaks)-1)));
OUTPUTmonitor_selection=outputmonitor(pos_1:pos_2);        

% add noise, if desired
%OUTPUTmonitor_selection=OUTPUTmonitor_selection+0.01.*( randn(length(OUTPUTmonitor_selection),1) + 1i*randn(length(OUTPUTmonitor_selection),1)   );

c=3*1e8;
V=linspace(-0.5*2*pi/dt,0.5*2*pi/dt,length(OUTPUTmonitor_selection));
W=V+2*pi*c/lambda0;
f=(W./(2*pi)).*1e-12;
spectrum=(dt/length(V))*abs(fftshift(fft(OUTPUTmonitor_selection))).^2;

envelope_n_o_points=5000;
spectrum_envelope=zeros(1,envelope_n_o_points);
f_envelope=zeros(1,envelope_n_o_points);
n_o_samp=floor(length(spectrum)/envelope_n_o_points);
for i=1:envelope_n_o_points
    if i==length(spectrum)
        [temp_max,pos_max]=max(10*log10(spectrum(1+(i-1)*n_o_samp:length(spectrum))/1e-3));
    else
        [temp_max,pos_max]=max(10*log10(spectrum(1+(i-1)*n_o_samp:i*n_o_samp)/1e-3));
    end
    spectrum_envelope(i)=temp_max;
    f_envelope(i)=f(pos_max+(i-1)*n_o_samp);
end

lambdas=(c./f).*1e-3;
spectrum=10*log10(spectrum/1e-3);
lambdas_envelope=(c./f_envelope).*1e-3;
envelope=spectrum_envelope;

[peaks,peakpos]=findpeaks(spectrum,'MinPeakHeight',max(spectrum)-20);
phases=angle(fftshift(fft(OUTPUTmonitor_selection)));
phases=phases(peakpos);
lambdas_phase=lambdas(peakpos);
%end
end

