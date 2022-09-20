function [ options ] = checkDefaultValuesAIO( options )
%CheckDefaultValuesAIO is to check the data segmentation and to create
% default values for the missing parameters
% options variable is an All-In-One formated structure

    % it loads the treshold variables    
    Test_settings;
    
    % rotating the input-output data into the same dimensions
    [Ni,li] = size(options.u); if li < Ni, options.u = options.u'; [Ni,li] = size(options.u); end
    [No,lo] = size(options.y); if lo < No, options.y = options.y'; [No,lo]  = size(options.y); end
    options.NInputchannels=Ni;
    options.NOutputchannels=No;
    
    % length is the shortest data of input-output
    l=min(li,lo);

    
    % initializing number of realizations variable
    if(~isfield(options,'R')) options.R=0; end 
    if(isempty(options.R)) options.R=0; end 
    if(isempty(options.fs)) options.fs=1; end 
    
    % estimate the period length and compare it with the manually given one
    if options.N>0 % comparing the automated-manual ones
        data=options.u;
        testN=Estimate_N(data,options.fs);
        if testN>0 && abs(options.N-Estimate_N(data,options.fs))>options.N/10 
            warning('Manual period length differs from the automated period length.');
        end
        clear data;
    else % estimating the period length 
        if options.N==0 options.N=Estimate_N(options.u,options.fs); end
        if options.N==0 options.N=l; warning('Nonperiodic data detected.'); end        
    end
    
   
    % estimating the number of periods
    if options.N>0 && options.P==0 options.P=Estimate_P(options.u,options.N); end
    % estimating the period length if the pervious estimate failed
    if options.N==0 options.N=min(li,lo); end
    % estimating the number of periods if the pervious estimate failed
    if options.P==0 options.P=max(1,floor(min(li,lo)./options.N)); end
    
    % estimating the number of realizations if the pervious estimate failed
    if options.N>0 && options.P>0 && options.R==0 options.R=floor(l/options.P/options.N); end
    % estimating the number of samples in a period if the pervious estimate failed
    if options.P>0 && options.R>0 && options.N==0 options.N=floor(l/options.P/options.R); end
    % estimating the number of periods if the pervious estimate failed
    if options.P==0 && options.R>0 && options.N>0 options.P=floor(l/options.N/options.R); end
    
    % if nothing works, then error
    if options.N==0 || options.P==0 || options.R==0 error('Automated data seqmentation failed'); end

    if(isfield(options,'estimateTransient'))
       if options.estimateTransient
           if ~isfield(options,'Ptr')
               options.Ptr=0;
           end
       end
    end
    
    if options.N>0 && options.P>2 && ~isfield(options,'Ptr') options.Ptr=Estimate_Ptr(options.y,options.P,options.N); end
    if options.N>0 && options.P==1 && ~isfield(options,'Ptr') options.Ptr=0; end
    if ~isfield(options,'Ptr') options.Ptr=0; end

    if options.N*options.P*options.R>l error('Seqmentation parameters are incompatible with data'); end

    % frequency scale
    f0 = options.fs/options.N;
    options.f=0:f0:floor((options.N-1)/2)*f0;

    % determine which lines are excited
    options = determineExcitationParameters(options);
    
    % determine which solver has to be used
    if(~isfield(options,'solver')) options = determineSolver(options); end
    if isempty(options.solver) options = determineSolver(options); end

    % initialize estimateTransient options
    if(~isfield(options,'estimateTransient')) 
         if options.P<3 
             options.estimateTransient=1
         elseif(options.P-options.Ptr)<2 
             options.estimateTransient=1;
         else 
             options.estimateTransient=0; 
         end
    end
    

    % determine the degree and bandwidth parameters
    if(~isfield(options,'degree')) options.degree=0; end %degree of poly       
    if(options.degree==0) options.degree=2; end %degree
    if(~isfield(options,'bw')) options.bw=0; end %degree of poly    
    if(strcmp(upper(options.solver),'LPM') && options.bw==0)         
       options.bw=(options.degree+1)*(options.NInputchannels+options.estimateTransient)+1;         
    end
    if(strcmp(upper(options.solver),'LRM') && options.bw==0)         
        options.bw=(options.degree+1)*(options.NInputchannels+options.estimateTransient)+options.degree+1; 
    end       
    

    % setting up all frequency indicies
    options.ind.all = 1:floor(options.N/2); % all frequencies
    
    
    end
