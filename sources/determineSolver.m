function options = determineSolver(options)
% it determines that LPM or LRM has to be used

% loading treshol settings
Test_settings;
l=size(options.u,1);
[~,lowest_y]=min(rms(options.y,2));
[~,lowest_u]=min(rms(options.u,2));


SNR_U=0; SNR_Y=0;

% we estimate the SNRs
if options.P>1
    U1=abs(fft(options.u(lowest_u,(options.P-1)*options.N+1:(options.P-0)*options.N)/sqrt(options.N))); U1=U1(options.ind.exc);          
    U2=abs(fft(options.u(lowest_u,(options.P-2)*options.N+1:(options.P-1)*options.N)/sqrt(options.N))); U2=U2(options.ind.exc);
    SNR_U=20 * log10(rms(U1) / rms(U1-U2));

    Y1=abs(fft(options.y(lowest_y,(options.P-1)*options.N+1:(options.P-0)*options.N)/sqrt(options.N))); Y1=Y1(options.ind.exc);          
    Y2=abs(fft(options.y(lowest_y,(options.P-2)*options.N+1:(options.P-1)*options.N)/sqrt(options.N))); Y2=Y2(options.ind.exc);
    SNR_Y=20 * log10(rms(Y1) / rms(Y1-Y2));
end

    % if there are at least two steady states periods then we don't have to
    % estimate transient
    if (options.P-options.Ptr)<2
        if(~isfield(options,'estimateTransient')) options.estimateTransient=1; end
    else
        if(~isfield(options,'estimateTransient')) options.estimateTransient=0; end
    end
    % if the SNR is band then we go foor LPM
    if SNR_Y>Test.TresholdforLPM || options.P==1
        options.solver='LRM';
    else
        options.solver='LPM'; 
    end
    
    