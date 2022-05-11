function  [fr,t,fano,mean_FR] = firing_rate(spike_ts,sr,bin_sz,rec_dur)

if nargin<4
    rec_dur = spike_ts(end)/sr + bin_sz;
else
    
end

fr=zeros(1,round((rec_dur)/(bin_sz)));
t = linspace(0,rec_dur,round((rec_dur)/(bin_sz)));
% Binarize the data
Spike_binarized=zeros(1,round(rec_dur*sr)+(bin_sz*sr)+1);
for i=1:numel(spike_ts)
    idx=round(spike_ts(i));
    if idx<0
    else
 Spike_binarized(1,idx)=1;
    end
end  

% firing rate
for i=1:size(fr,2)
    win_in=(i-1)*(bin_sz*sr)+1;
    win_out=i*(bin_sz*sr)+1;
    spike_temp = Spike_binarized(1,win_in:win_out);
    fr(1,i) = sum(spike_temp)/(bin_sz);
    isi = [];
    if numel(diff(find(spike_temp==1)))==1
    isi = diff(find(spike_temp==1))./sr;
    elseif isempty(diff(find(spike_temp==1)))
    isi = NaN;    
    else
    isi(1,:) = diff(find(spike_temp==1))./sr; 
    end
    fano(i) = nanstd(isi)^(2)/nanmean(isi); % fano factor within bin
end

% Mean fr
mean_FR = mean(fr,2);

end