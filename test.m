BetaPeakCounter = 0;
for i = 1:D_tf.nchannels    
    if BetaPeakLocMatrix(i) > 0
        BetaPeakCounter = BetaPeakCounter + 1;
        BetaPeakMatrix(ThetaPeakCounter,:) = PowerSpectrumNormalized(i,BetaPeakLocMatrix(i) - 5:BetaPeakLocMatrix(i) + 5);
    end
end