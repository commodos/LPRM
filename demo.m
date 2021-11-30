clear all;
close all;

% loading the test data (car frame testing)
load test_data.mat;

%% suggestions to try out:
% this estimate uses relative frequency index
options=LPRM(u,y); 

% using explicitly LPM method
options=LPRM(u,y,fs,'solver','LPM'); 

% gives correlation warning
options=LPRM([u(1,:);u(1,:)],y) 

% gives warning because the automated period length is
% not the same as the manually given one
options=LPRM(u,y,[],2048) 

% using LRM with bandwith of 10 
options=LPRM(u,y,'bw',10) 
% using LRM with bandwith of 40
options=LPRM(u,y,'bw',40) 

%%
% These results are shown in the article 

% classical H1 estimate using 1 block
options_H1_1block=LPRM(u,y,fs,'solver','H1','R',1,'P',1);

% LRM estimate using 1 block
options_LRM_1block=LPRM(u,y,fs,'R',1,'P',1);

% classical H1 estimate using all data
options_H1_all_data=LPRM(u,y,fs,'solver','H1');

% LRM estimate using all data
options_LRM_all_data=LPRM(u,y,fs,'estimateTransient',1);


%% compare cases with all data
figure
index_fig=0;
for index_1=1:2
    for index_2=1:2
        index_fig=index_fig+1;
        subplot(2,2,index_fig); hold on; grid on;
        plot(options_H1_all_data.f,20*log10(abs(squeeze(options_H1_all_data.G(:,index_1,index_2)))));
        plot(options_LRM_all_data.f,20*log10(abs(squeeze(options_LRM_all_data.G(:,index_1,index_2)))));
        legend('H_1 using all data','LRM using all data')
        xlabel('frequency [Hz]'); ylabel('Magnitude [dB]');
    end
end

%% compare cases with 1 block
figure;
index_fig=0;
for index_1=1:2
    for index_2=1:2
        index_fig=index_fig+1;
        subplot(2,2,index_fig); hold on; grid on;
        plot(options_H1_1block.f,20*log10(abs(squeeze(options_H1_1block.G(:,index_1,index_2)))));
        plot(options_LRM_1block.f,20*log10(abs(squeeze(options_LRM_1block.G(:,index_1,index_2)))));
        legend('H_1 using 1 block','LRM using 1 block')
        xlabel('frequency [Hz]'); ylabel('Magnitude [dB]');
    end
end



%% compare cases for the first FRF - this is in the article
figure;
index_1=1;
index_2=2;

subplot(1,3,1); hold on; grid on;
title('H_1 estimates')
plot(options_H1_all_data.f,20*log10(abs(squeeze(options_H1_all_data.G(:,index_1,index_2)))));
plot(options_H1_1block.f,20*log10(abs(squeeze(options_H1_1block.G(:,index_1,index_2)))));
xlim([options_H1_1block.fmin options_H1_1block.fmax]);
legend('using all data','using 1 block');
xlabel('frequency [Hz]'); ylabel('Magnitude [dB]');

subplot(1,3,2); hold on; grid on;
title('LRM estimates')
plot(options_LRM_all_data.f,20*log10(abs(squeeze(options_LRM_all_data.G(:,index_1,index_2)))));
plot(options_LRM_1block.f,20*log10(abs(squeeze(options_LRM_1block.G(:,index_1,index_2)))));
xlim([options_LRM_1block.fmin options_LRM_1block.fmax]);
legend('using all data','using 1 block');
xlabel('frequency [Hz]'); ylabel('Magnitude [dB]');

subplot(1,3,3); hold on; grid on;
title('H_1 and LRM estimates for 1 block')
plot(options_H1_1block.f,20*log10(abs(squeeze(options_H1_1block.G(:,index_1,index_2)))));
plot(options_LRM_1block.f,20*log10(abs(squeeze(options_LRM_1block.G(:,index_1,index_2)))));
legend('H_1','LRM');
xlabel('frequency [Hz]'); ylabel('Magnitude [dB]');


