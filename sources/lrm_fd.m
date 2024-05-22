function [G,T] = lrm_fd(U,Y,lines,LRM,d,d_tr,bw,Tr,fs)
% lrm_fd Multivariable Transfer Function Estimate using Local Polynomial or Local Rational Method.
% [G,T]=lrm_fd(U,Y,lines,LRM,d,bw,Tr) estimates the transfer
%   function G at frequencies f of the system with time-domain input u and output y
%   U: input signal, frequency-domain data, one block
%   Y: output signal, frequency-domain data, one block, one channel
%   LRM: optional, if 1 LRM is used, if 0 then LPM is used, default is 1
%   D: optional, degree of the polynomials in transfer function, default is D
%   D_TR: optional, degree of the polynomials in transient terms , default is 2
%   BW: optional, bandwidth of the sliding window, default (D+1)*(#Input+LRM+TR)+LRM*d
%   TR: optional, estimate transient if 1
%
%   version 1.3
% 	Dr. Péter Zoltán CSURCSIA, September 2022

U=squeeze(U);
Y=squeeze(Y);
[li,Ni] = size(U); if li < Ni, U = U.'; [li,Ni] = size(U); end
[lo,No] = size(Y); if lo < No, Y = Y.'; [lo,No] = size(Y); end
N=li;

G = zeros(Ni,N);
T = zeros(N,1);


center=0;

for k = lines    
    if 1>k-floor(bw/2)
        [r POLY  POLY_T POLY_Y] = lpm_create_polynomials(k,lines,LRM,d,d_tr,bw,Tr,Ni);
    elseif 1<=k-floor(bw/2) && ~center
        center=1;
        [r POLY  POLY_T POLY_Y] = lpm_create_polynomials(k,lines,LRM,d,d_tr,bw,Tr,Ni);
    elseif((k>floor(N/2)-(ceil(bw/2)-1)))        
        [r POLY  POLY_T POLY_Y] = lpm_create_polynomials(k,lines,LRM,d,d_tr,bw,Tr,Ni);
    end
        
    f_selected=k+r;

    
    K=[kron(squeeze(U(f_selected,:)),ones(1,d+1)).*POLY];
    L=Y(f_selected);                                        
    
    if(LRM)
        K=[K -repmat(L,1,d).*POLY_Y];
    end

    K=[K POLY_T];
        
    theta=K\L;
    G(:,k)=theta(1:(d+1):(d+1)*(Ni));   
    if(Tr) T(k)=theta(end-d); end     

end


function [r POLY  POLY_T POLY_Y] = lpm_create_polynomials(k,lines,LRM,d,d_tr,bw,Tr,Ni)
% lpm_create_polynomials 
%   D: polynomials in transfer function
%   D_TR: polynomials in transient term
%   BW: bandwidth of the sliding window
%   TR: estimate transient if 1
%   LRM: LRM if 1, LPM is 0
%
%   version 1.2
% 	Dr. Péter Zoltán CSURCSIA, September 2022


POLY_T=[];
POLY_Y=[];

r=-floor(bw/2):ceil(bw/2)-1; % central frequencies
k_last=lines(end);

if r(1)+k<1
    r=r-(r(1)+k)+1;
elseif k_last-k<r(end)
    r=r-(r(end)-(k_last-k));
end

polynomial_t=[]; 
polynomial=[]; 
polynomial_Y=[]; 

for order=0:d polynomial(:,order+1)=r.^order; end
for order=0:d_tr polynomial_t(:,order+1)=r.^order; end

POLY=repmat(polynomial,1,Ni); 

if Tr
    POLY_T=polynomial_t;
end

if LRM
    for order=1:d polynomial_Y(:,order)=r.^order; end
    POLY_Y=polynomial_Y;     
end
