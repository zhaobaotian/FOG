% for visulization and mannual check
for i = 1:D_Theta.nchannels
    timetamps_temp = ThetaBurstTimestamps{i};
    Signal_temp = D_Theta(i,:,1);
    TimeInterval = zeros(D_Theta.nsamples,1);
    ThetaPowerEnvelope_temp = ThetaPowerEnvelope(i,:)';
    for j = 1:size(timetamps_temp,1)
        TimeInterval(int16((timetamps_temp(j,1)/1000)*500):int16((timetamps_temp(j,3)/1000)*500)) = max(Signal_temp);
    end
        plotECG(D_Theta.time,[Signal_temp' repmat(Threshold75th(i),[D_Theta.nsamples 1]) TimeInterval ThetaPowerEnvelope_temp])
end