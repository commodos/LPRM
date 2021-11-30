    %% Segment data into realizations, periods
    y_t=reshape(options.y(:,1:options.N*options.P*options.R)',options.N*options.P,options.R,options.NOutputchannels);
    u_t=reshape(options.u(:,1:options.N*options.P*options.R)',options.N*options.P,options.R,options.NInputchannels);
    
    
    %% Remove transients term 
    u_t_full=u_t;
    y_t_full=y_t;    
    u_t=u_t(options.N*options.Ptr+1:end,:,:);      
    y_t=y_t(options.N*options.Ptr+1:end,:,:);
    options.P_used=options.P-options.Ptr;


    %% initialize variables
    
    options.G           = zeros(floor(options.N/2),options.NOutputchannels,options.NInputchannels);
    options.std_G_n     = zeros(floor(options.N/2),options.NOutputchannels,options.NInputchannels);
    options.std_G       = zeros(floor(options.N/2),options.NOutputchannels,options.NInputchannels);
    options.gammaSquare = zeros(floor(options.N/2),options.NOutputchannels,options.NInputchannels);     
   
    options.U_est       = zeros(floor(options.N/2),options.NInputchannels);
    options.Y_est       = zeros(floor(options.N/2),options.NOutputchannels);    
    options.T_est       = zeros(floor(options.N/2),options.NOutputchannels);    

    options.std_Y_est_n = zeros(floor(options.N/2),options.NOutputchannels);  
    options.std_U_est_n = zeros(floor(options.N/2),options.NInputchannels);
    options.SNR_U       = zeros(floor(options.N/2),options.NInputchannels);
    options.SNR_Y       = zeros(floor(options.N/2),options.NOutputchannels);
    options.SNR_FRF     = zeros(floor(options.N/2),options.NOutputchannels,options.NInputchannels);
    
    
    var_G_n             = zeros(floor(options.N/2),options.NOutputchannels,options.NInputchannels);
    var_G               = zeros(floor(options.N/2),options.NOutputchannels,options.NInputchannels);
    var_G_s             = zeros(floor(options.N/2),options.NOutputchannels,options.NInputchannels); 
    
    G_m                 = zeros(options.R,floor(options.N/2),options.NOutputchannels,options.NInputchannels); 
    var_G_m             = zeros(options.R,floor(options.N/2),options.NOutputchannels,options.NInputchannels); 
    var_U_m             = zeros(options.R,floor(options.N/2),options.NInputchannels);
    var_Y_m             = zeros(options.R,floor(options.N/2),options.NOutputchannels);
    U_m                 = zeros(options.R,floor(options.N/2),options.P_used,options.NInputchannels);
    Y_m                 = zeros(options.R,floor(options.N/2),options.P_used,options.NOutputchannels);
    T_m                 = zeros(options.R,floor(options.N/2),options.P_used,options.NOutputchannels);
    U_mean              = zeros(options.R,floor(options.N/2),options.NInputchannels);
    Y_mean              = zeros(options.R,floor(options.N/2),options.NOutputchannels);
    var_U_n             = zeros(floor(options.N/2),options.NInputchannels);
    var_Y_n             = zeros(floor(options.N/2),options.NOutputchannels);      
    
    var_U_est_n          = zeros(floor(options.N/2),options.NInputchannels);
    var_Y_est_n          = zeros(floor(options.N/2),options.NOutputchannels);
    
    %% estimate input-output spectra
    
    for m = 1:options.R        
        U=fft(reshape(u_t(:,m,:),options.N,options.P_used,options.NInputchannels))/sqrt(options.N);
        U=U(options.ind.all,:,:);
        U_mean(m,:,:)  = mean(U,2);   
        var_U_m(m,:,:)=var(U,0,2);
        U_m(m,:,:,:)=U; 
 

        Y=fft(reshape(y_t(:,m,:),options.N,options.P_used,options.NOutputchannels))/sqrt(options.N);
        Y=Y(options.ind.all,:,:);     
        Y_mean(m,:,:)  = mean(Y,2);     
        var_Y_m(m,:,:)=var(Y,0,2);
        Y_m(m,:,:,:)=Y;
        
        U=fft(reshape(u_t_full(:,m,:),options.N,options.P,options.NInputchannels))/sqrt(options.N);
        U=U(options.ind.all,:,:);
        U_m_full(m,:,:,:)=U;  
        
        Y=fft(reshape(y_t_full(:,m,:),options.N,options.P,options.NOutputchannels))/sqrt(options.N);
        Y=Y(options.ind.all,:,:);
        Y_m_full(m,:,:,:)=Y;          
        
    end
    
    