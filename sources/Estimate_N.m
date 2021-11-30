function N=Estimate_N(data,fs)
% it estimates the period length of the data using auto-correlation function estimates  
    N=0;
    d_found=[];
    data_length=size(data,2);          
    data_channels=size(data,1); 
    data=data-repmat(mean(data,2),1,data_length);      
   
    % we go over all data channels and check the individual periodicities
    for index_c=1:data_channels          
        x=(abs(xcorr(data(index_c,1:data_length),'biased'))); 
        x=x/max(x);
        try
            [~, pindecies]=findpeaks(x(floor(length(x)/2)-1:end),'MinPeakHeight',0.3,'MinPeakDistance',fs);              
        catch
            [~, pindecies]=findpeaks(x(floor(length(x)/2)-1:end),'MinPeakHeight',0.3,'MinPeakDistance',0.1*fs);              
        end
        d=diff(pindecies);
        d_found=[d_found d(find(d>fs/3))];
    end

    if isempty(d_found) return; end
    
    % we take the mode of the detected periodicities
    if mode(d_found)==mean(d_found) d_found=mean(d_found); end
    
    if length(d_found)>2
        [n,bin] = hist(d_found,unique(d_found));
        [~,idx] = sort(-n);   
    else
       n=1; bin=d_found; idx=1; 
    end
    rms_differences=zeros(1,length(idx));
    N_test=zeros(1,length(idx));

    % we test if this periodicity (full band excitation) or the double of
    % this length has to be used (sparse excitation)
    index=1;
    for testN=bin(idx)

        Nest=testN;
        if data_length>4*Nest
            if(mean(rms(data(:,1:2*Nest)'-data(:,2*Nest+1:4*Nest)'))<0.5*mean(rms(data(:,1:Nest)'-data(:,Nest+1:2*Nest)')))          
                Nest=Nest*2;
            end
        end
        if data_length>2*Nest        
            rms_differences(index)=mean(rms(data(:,1:Nest)'-data(:,Nest+1:2*Nest)'));
        else
            rms_differences(index)=inf;
        end
        N_test(index)=Nest;
        index=index+1;
    end
    [~, idx]=min(rms_differences);
    
    N=N_test(idx);
    % if everything fails, we return with 0 period length
    if isnan(N)
        N=0;
    end
end