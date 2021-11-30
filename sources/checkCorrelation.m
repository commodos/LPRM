function checkCorrelation(u)
% this function performs correlation check-up for the input signal

% treshold variables are in this file:
Test_settings

if(size(u,1)>size(u,2)) u=u.'; end
Ni=size(u,1);
N=size(u,2);

% gives warning if the signals are correlated
if N>0
    r=abs(corrcoef(u'));
    for rx=1:Ni 
        for rx2=rx:Ni
            if(r(rx,rx2)>Test.TresholdforCorrelation && rx2~=rx)
                warning(sprintf('input %i is CORRELATED with input %i',rx,rx2));
            end
        end
    end
end