function Ptr = Estimate_Ptr(data,P,N)
    % estimates the number of transient distrubed periods
    data_length=size(data,2);          
    data_channels=size(data,1);       
    Ptr=0;
    
    if P==1 return; end    
    % if there are multiple periods, for each channel and period the rms
    % values are calculated as it is detailed in the MSSP paper
    for i = 1:data_channels   
        ylast=data(i,(P-1)*N+1:P*N);
        for px=0:P-1
            index=px*N+1:(px+1)*N;       
            rms_y(i,px+1)=rms(data(i,index)-ylast,2);          
        end

        % the treshold number
        rmslimit=3*std(rms_y(i,end-1:end));            

        for px=1:P    
            if(rms_y(i,px)>rmslimit)
                Ptr=max(Ptr,px);        
            else
                break;
            end
        end    
    end

