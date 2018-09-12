function [ave_coh_otime,coh_otime_map,cuttime] = time_coherence(powspctrm_map,crsspctrm_map,cyclenum,freq,timevec,time_interval)



winsize = round((cyclenum./freq)./(2*time_interval));

maxhfwin = max(winsize);

padded_crsspctrm_map = cat(3,NaN(size(crsspctrm_map,1),size(crsspctrm_map,2),maxhfwin),crsspctrm_map, ...
    NaN(size(crsspctrm_map,1),size(crsspctrm_map,2),maxhfwin));

padded_powspctrm_map = cat(4,NaN(size(powspctrm_map,1),size(powspctrm_map,2),size(powspctrm_map,3),maxhfwin), ...
    powspctrm_map,NaN(size(powspctrm_map,1),size(powspctrm_map,2),size(powspctrm_map,3),maxhfwin));



coh_otime_map = zeros(size(crsspctrm_map));


for freq_ind = 1:length(freq)
    
    if freq_ind == 1
        fprintf('%03i%% ', floor((freq_ind/length(freq))*100));
    else
        fprintf('\b\b\b\b\b%03i%% ', floor((freq_ind/length(freq))*100));
    end
    if freq_ind == length(freq)
        fprintf('\n');
    end
 
%   disp(num2str(freq_ind));  
    
  winsize_tmp = winsize(freq_ind);
 

  for time_ind = 1: size(coh_otime_map,3)
      
  centerpos = time_ind+maxhfwin;    
     
   numerator = nanmean(padded_crsspctrm_map(:,freq_ind,centerpos-winsize_tmp:centerpos+winsize_tmp),3); 
    
%    numerator = numerator.^2;
   
   denominator_a = nanmean(padded_powspctrm_map(:,1,freq_ind,centerpos-winsize_tmp:centerpos+winsize_tmp),4);
   
   denominator_a = squeeze(denominator_a);
   
   denominator_b = nanmean(padded_powspctrm_map(:,2,freq_ind,centerpos-winsize_tmp:centerpos+winsize_tmp),4);
   
   denominator_b = squeeze(denominator_b);
   
   coh_otime_map(:,freq_ind,time_ind) = numerator./sqrt(denominator_a.*denominator_b);
   
  end

end


coh_otime_map = abs(imag(coh_otime_map));

coh_otime_map = coh_otime_map.^2;

% coh_otime_map = abs(coh_otime_map);

ave_coh_otime = squeeze(nanmean(coh_otime_map,1));

ave_coh_otime = ave_coh_otime(:,maxhfwin+1:end-maxhfwin);

coh_otime_map = coh_otime_map(:,:,maxhfwin+1:end-maxhfwin);

cuttime       = timevec(maxhfwin+1:end-maxhfwin);

end



