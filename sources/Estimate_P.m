function P=Estimate_P(data,N)
    
    P=0;
    data_length=size(data,2);          
    data_channels=size(data,1);             
    data=data-repmat(mean(data,2),1,data_length);    
    
    if data_length<2*N 
        P=1;     
    else
        C_channels=zeros(data_channels,1);

        for index=1:data_channels
            shift=1;
            C=(corrcoef(data(index,1:N),data(index,N+1:(shift+1)*N)));
            while C(1,2)>0.9 && shift<floor((data_length-N)/N)
              shift=shift+1;
              C=(corrcoef(data(index,1:N),data(index,shift*N+1:(shift+1)*N)));
            end
            if shift~=floor((data_length-N)/N) 
                C_channels(index)=shift;          
            end      
            if shift==floor((data_length-N)/N) && C(1,2)>0.9
                C_channels(index)=shift+1;          
            end
        end
        P=mode(C_channels);
    end
end