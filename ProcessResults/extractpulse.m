function [ selection, fit, t ] = extractpulse( signal )
    % function to extract a pulse from the datatrace for visualization

    OUTPUTmonitor=signal;
    time_span=2000*1e-12; % should be modified
    n_o_timesteps=length(OUTPUTmonitor);
    dt=20e-15; % should be modified
    OUTPUTmonitor_selection=OUTPUTmonitor(n_o_timesteps-round(time_span/dt):n_o_timesteps);
    [peaks,peakpositions]=findpeaks(abs(OUTPUTmonitor_selection).^2, 'MinPeakHeight',max(abs(OUTPUTmonitor_selection).^2)/10);

    pos_1=0.5*(peakpositions(length(peaks)-1)+peakpositions(length(peaks)-2));   
    pos_2=0.5*(peakpositions(length(peaks)-1)+peakpositions(length(peaks)));   

    steadystatepulsenergy=sum(abs(OUTPUTmonitor_selection(pos_1:pos_2)).^2)*dt*1e12;
    pulsereprate=1e-9/((pos_2-pos_1)*dt);

    % determine pulsewidth of output pulses
    pulsewidth=-1;
    startpeakpos=peakpositions(length(peaks)-1);
    peakpower=abs(OUTPUTmonitor_selection(startpeakpos)).^2;
    for i=startpeakpos:length(OUTPUTmonitor_selection)
        if abs(OUTPUTmonitor_selection(i))^2<=peakpower/2
            pulsewidth=(2*(i-startpeakpos)*dt)*1e12;
            break;
        end
    end


    span=20e-12;
    nr_samples_selection=ceil(span/dt); %ceil((10*pulsewidth*1e-12)/dt);
    lower_bound=max(1,round(startpeakpos-nr_samples_selection/2));
    upper_bound=min(length(OUTPUTmonitor_selection),round(startpeakpos+nr_samples_selection/2));
    toffset=(length(OUTPUTmonitor)-length(OUTPUTmonitor_selection))*dt*1e12;
    taxis=linspace(toffset+lower_bound*dt*1e12,toffset+upper_bound*dt*1e12,nr_samples_selection+1);

    ysamples=transpose(abs(OUTPUTmonitor_selection(lower_bound:upper_bound)).^2);
    xdata=taxis;
    length(ysamples)
    length(xdata)
    t_center=taxis(find(ysamples==max(ysamples)));
    P0=max(ysamples);
    ti=t_center;
    fun=@(x,xdata) (P0*(sech((xdata-x(2))./x(1)).^2)); %sech^2 fit
    x0=[0.5,ti];
    x=lsqcurvefit(fun,x0,xdata,ysamples);
    pulsewidth=1.763*x(1);
    
    selection=abs(OUTPUTmonitor_selection(lower_bound:upper_bound)).^2;
    fit=fun(x,xdata);
    t=taxis;

end

