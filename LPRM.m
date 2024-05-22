function [options] = LPRM(varargin)
% LPRM: iteration free Local Polynomial/Rational Method
% LPRM Create and uses an All-In-One (options structure) format 
% supported by MUMI, SAMI toolboxes. 
%
% usage:
% [options]=LPRM(options), options is in AIO format 
% OR
% [options]=LPRM(u,y,fs,N,P,R) 
%   u: input signal, time-domain data
%   y: output signal, time-domain data
%   fs: sampling frequency, optional, default is 1
%   N: number of samples in a (processing) block, optional
%   P: number of periods per realization, optional
%   R: number of realization, optional
% OTHER INPUT PARAMETERS
%   P_tr: number of periods to be discarded
%   RMSthreshold: used to determine the lines parameters (on r or u)
%   fmin: min frequency of interest, default: automatically determined
%   fmin: max frequency of interest, default: automatically determined
%   Ptr: number of transient blocks, default: automatically determined
%   solver: default: automatically determined 
%     LPM: Local Polynomial Method is used
%     LRM: Local Rational Method is used
%   bw and degree parameters: bandwith and degree (default automatically set)
%   estimateTransient: default 1 (1=estimate the transient)
%   
%   Output paramters in the AIO structure:
%   NInputchannels: number of input channels
%   NOutputchannels: number of output channels
%   G: The estimated averaged FRF, size floor(N/2)*NOutputchannels*NInputchannels
%   std_G_n: The noise FRF estimate, size floor(N/2)*NOutputchannels*NInputchannels
%   std_G_BLA: The std estimate of the FRF, size floor(N/2)*NOutputchannels*NInputchannels
%   SNR_FRF: The SNR estimate of the FRF, size floor(N/2)*NOutputchannels*NInputchannels
%   U_est: The estimated averaged input spectrum, size floor(N/2)*NInputchannels
%   std_U_est_n: The estimated noise of the input, size floor(N/2)*NInputchannels
%   SNR_U: The estimated SNR of the input, size floor(N/2)*NInputchannels
%   Y_est: The estimated averaged output spectrum, size floor(N/2)*NOutputchannels
%   std_Y_est_n: The estimated noise of the output, size floor(N/2)*NOutputchannels
%   SNR_Y: The estimated SNR of the output, size floor(N/2)*NOutputchannels
%   T_est: The estimated highest level of transient spectrum, size floor(N/2)*NOutputchannels
%   gammaSquare: multiple coherence of the processed data, size floor(N/2)*NOutputchannels*NInputchannels
%
%
%   version 1.3
% 	Dr. Péter Zoltán CSURCSIA, 2010
%   Last modified September 2022 

% add the foler of the auxiliary source codes
addpath sources;

%% process the input paramters
Nvar=length(varargin);
if Nvar==0; error('No function parameters'); end
if Nvar==1 && ~isnumeric(varargin{1})
    options=varargin{1};
else
    for i = 1:Nvar
        if isnumeric(varargin{i})
            if i==1
                options.u=varargin{i};
            elseif isnumeric(varargin{i-1})                
               switch i
                   case 2, options.y=varargin{i};
                   case 3, options.fs=varargin{i};
                   case 4, options.N=varargin{i};
                   case 5, options.P=varargin{i};
                   case 6, options.R=varargin{i};                   
               end
            end
        elseif i>1 && i<=Nvar-1
            if ~isnumeric(varargin{i})
                options=setfield(options,varargin{i},varargin{i+1});
                i=i+1;
            end
        end        
    end    
end
 
%% pre-process the input parameters
if ~isfield(options,'u') || ~isfield(options,'y') error('Input/output data required'); end
if ~isfield(options,'fs') options.fs=1; end
if ~isfield(options,'N') options.N=0; end
if ~isfield(options,'P') options.P=0; end
if ~isfield(options,'R') options.R=0; end

% rotate everything in the same order
[Ni,li] = size(options.u); if li < Ni, options.u = options.u'; [Ni,li] = size(options.u); end
[No,lo] = size(options.y); if lo < No, options.y = options.y'; [No,lo]  = size(options.y); end

% length of the measurement
l=min([li lo]); 

if l<li warning('Input signal is longer than the output signal, it will be truncated.'); options.u=options.u(:,1:l); end
if l<lo warning('Output signal is longer than the input signal, it will be truncated.'); options.y=options.y(:,1:l); end

options.NInputchannels=Ni;
options.NOutputchannels=No;

% perform input correlation check-up
if li>0 
    checkCorrelation(options.u); 
end

% determine default values
options=checkDefaultValuesAIO(options);
       
% segment data
segmentData;

%% estimating FRM
hwaitbar=waitbar(0,'Estimating FRM...'); 

if (strcmp(upper(options.solver),'LPM') || strcmp(upper(options.solver),'LRM') || strcmp(upper(options.solver),'H1'))
    % if LPM or LRM is used
    if (strcmp(upper(options.solver),'LPM') || strcmp(upper(options.solver),'LRM'))
        if options.bw>length(options.ind.exc)
            warning('The bandwith is larger than the number of excitation points, consider using more data samples or decrease the level of RMSthreshold')
        end
        warning off;
        % go over different realizations
        for m = 1:options.R   
            waitbar(0.05+m/options.R*0.6,hwaitbar,'Estimating FRM...'); pause(0.0001);
            Gk=zeros(options.P_used,options.NOutputchannels,options.NInputchannels,floor(options.N/2));
            % go over different periods
            for index_p=1:options.P_used
                for ry=1:options.NOutputchannels
                    [Gk(index_p,ry,:,:) T_m(m,:,index_p,ry)]=lrm_fd(U_m(m,:,index_p,:),Y_m(m,:,index_p,ry),options.ind.exc,strcmp(upper(options.solver),'LRM'),options.degree,options.degree_tr,options.bw,options.estimateTransient,options.fs);
                end
            end
            
           % process the FRM and variance estimates frequency-wise
           for k = options.ind.exc
                G_m(m,k,:,:)=squeeze(mean(Gk(:,:,:,k),1));
                var_G_m(m,k,:,:)=var(Gk(:,:,:,k),0,1)/options.P_used;  
           end  
        end
        warning on;
    else
        % If H1 is used
        Gk=zeros(length(options.ind.all),options.P_used,options.NOutputchannels,options.NInputchannels);
        % go over different realizations
        for m = 1:options.R   
            waitbar(0.05+m/options.R*0.6,hwaitbar,'Estimating FRM...'); pause(0.0001);
            % go over the excited frequencies
            for k = options.ind.exc
                % go over the periods
                for index_p=1:options.P_used
                    Uk=squeeze(U_m(m:ceil(options.R/options.NInputchannels):options.R,k,index_p,:)).';
                    for ry=1:options.NOutputchannels
                        Yk=squeeze(Y_m(m:ceil(options.R/options.NInputchannels):options.R,k,index_p,ry)).';
                        Gk(k,index_p,ry,:)=(Yk*Uk')/(Uk*Uk');
                    end
                end
            end
            % process the FRM and variance estimates
            G_m(m,options.ind.exc,:,:)=squeeze(mean(Gk(options.ind.exc,:,:,:),2));
            var_G_m(m,options.ind.exc,:,:)=var(Gk(options.ind.exc,:,:,:),0,2)/options.P_used;
        end    
    end
else
    error('unsupported solver')
end

    waitbar(0.70,hwaitbar,'Processing FRM...'); pause(0.01);

    % initialize the transient estimate (if it was not estimated)
    if(options.estimateTransient==0) T_m(:,:,:,:)=0;   end

    %% input-output spectra estimation
    options.U_est(options.ind.all,:)=squeeze(mean(U_mean(:,options.ind.all,:),1));
    options.Y_est(options.ind.all,:)=squeeze(mean(Y_mean(:,options.ind.all,:),1));
    if(options.estimateTransient)
        for ry=1:options.NOutputchannels
            options.T_est(:,ry)=squeeze(max(max(T_m(:,:,:,ry),[],3),[],1));
        end
    end

    
    % store the corresponding variables per realization
    options.G_m=G_m;
    options.U_m=U_m;
    options.Y_m=Y_m;
    options.T_m=T_m;


    % calculate the FRM and total/noise std's
    for ru = 1:options.NInputchannels
        for ry = 1:options.NOutputchannels
            options.G(:,ry,ru)       = mean(G_m(:,:,ry,ru),1);
            var_G_BLA(:,ry,ru)       = var(G_m(:,:,ry,ru),0,1)/options.R;              
            var_G_BLA_n(:,ry,ru)     = mean(var_G_m(:,:,ry,ru),1)/options.R;
            options.std_G_n(:,ry,ru) = sqrt(var_G_BLA_n(:,ry,ru));
            options.std_G(:,ry,ru)   = sqrt(var_G_BLA(:,ry,ru));            
        end
    end
        
    waitbar(0.90,hwaitbar,'Processing FRM...'); pause(0.01);

    % interpolate the FRM over nonexcited frequencies
    options.G(:,:,:) = interp1(options.ind.exc,options.G(options.ind.exc,:,:),options.ind.all,'linear',0); %linear interpolation over the whole frequency domain
    % Sample noise variances
    options.std_U_est_n(:,:) = sqrt(squeeze(mean((var_U_m)./options.P_used,1))/options.R);
    options.std_Y_est_n(:,:) = sqrt(squeeze(mean((var_Y_m)./options.P_used,1))/options.R);   
     
    %% estimating multiple coherence function
    for ry=1:options.NOutputchannels
            for k = options.ind.exc

                   % to estimate the coherence function
                   U=reshape(squeeze(U_m(:,k,:,:)),options.R*options.P_used,options.NInputchannels).';
                   Y=reshape(squeeze(Y_m(:,k,:,ry)),options.R*options.P_used,1).';
                   SUU=U*U'/(options.P_used*options.R); 
                   SYY=Y*Y'/(options.P_used*options.R); 
                   SUY=U*Y'/(options.P_used*options.R);
                   % coherence from data     
                   diagonal=diag(real(SUY'*pinv(SUU)*SUY./SYY));
                   options.gammaSquare(k,ry,:)=min(diagonal,ones(size(diagonal)));%diag(real(SYU*pinv(SUU)*SYU'./SYY));                                                
            end
    end
    
    options.gammaSquare(:,:,:) = interp1(options.ind.exc,options.gammaSquare(options.ind.exc,:,:),options.ind.all,'linear',0); %linear interpolation over the whole frequency domain
 
    %% calculating the Singal-To-Noise ratios
    for rx=1:options.NInputchannels                      
        options.SNR_U(options.ind.exc,rx)=(abs(options.U_est(options.ind.exc,rx))./(options.std_U_est_n(options.ind.exc,rx)));
    end  
    for rx=1:options.NOutputchannels
        options.SNR_Y(options.ind.exc,rx)=(abs(options.Y_est(options.ind.exc,rx))./(options.std_Y_est_n(options.ind.exc,rx)));
    end
    for ry = 1:options.NOutputchannels
        for ru = 1:options.NInputchannels
            options.SNR_FRF(options.ind.exc,ry,ru) =(abs(options.G(options.ind.exc,ry,ru)./(options.std_G_n(options.ind.exc,ry,ru))));
        end
    end
    


close(hwaitbar)
% return the structure
AIO=options;
