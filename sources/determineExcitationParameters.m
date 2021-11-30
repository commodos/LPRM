function options = determineExcitationParameters(options)
% it determines the excitation frequency lines
    
    % treshold settings are loaded
    Test_settings;

    % we consider the channel with highest energy level 
    [~,highest_u]=max(rms(options.u,2));       
    x=options.u(highest_u,:); 
    
    % determine which lines are excited
    X=db(fft(x((options.P-1)*options.N+1:options.P*options.N))/sqrt(options.N));          
    if(~isfield(options,'RMSthreshold'))
        options.RMSthreshold=floor(max(X(2:ceil(options.N/2)))-10);
    end
    
    % the excited indicies 
    if(~isfield(options,'ind'))
       options.ind.exc=find(X(1:floor(length(X)/2))>options.RMSthreshold);            
    end
    if(~isfield(options.ind,'exc'))
       options.ind.exc=find(X(1:floor(length(X)/2))>options.RMSthreshold);            
    end

    % lowest and highest frequency lines
    if(~isfield(options,'fmin')) options.fmin=options.fs/options.N*(options.ind.exc(1)-1); end
    if(~isfield(options,'fmax')) options.fmax=options.fs/options.N*(options.ind.exc(end)-1); end
    
end