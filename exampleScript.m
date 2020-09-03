clc;clear;close all
addpath(genpath(fileparts(which('exampleScript.m'))))% please do not change the file's name

%% %%%%%%%%%%%%%%%%%%%%%% SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%

%%%% Generate phase  %%%%%
[exp_phase,phase]=WrapNoisyGaussPhase(480,640,[120,160],5,6,2); % simulated noisy phase


%%%%%% Run SPAF  %%%%%%%%% (our method)
[SPAFphasor,unwrappedPhase]=SPAFUnwrap(exp_phase,10,10);        % SPAF + unwrap


%%%%%% Run MNF  %%%%%%%%% (for comparison) 
PhiUnwrap = cunwrap(phase);                                     % unwrapped by smartPAF


%%%%%%  Display  %%%%%%%%%

figure(100),mesh(angle(SPAFphasor));colormap(jet)               % Filtered (SPAF) phase is noisy but it does not have residue
title('wrapped SPAF Phase');axis off;colorbar;

figure(101),
subplot(1,2,1);imagesc(PhiUnwrap);caxis([-5 40]);
title('MNF');axis off
h=colorbar;ylabel(h, 'Phase (rad)', 'FontSize', 16);
subplot(1,2,2);imagesc(unwrappedPhase);caxis([-5 40]);      
title('SPAF');axis off
h=colorbar;ylabel(h, 'Phase (rad)', 'FontSize', 16);
colormap(jet)

%% %%%%%%%%%%%%%%%%%%%%%% NEURON %%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% read data %%%%%%%%%
exampleNeuron=[7514,7293]; %Neurons are extracted from the multplexed hologram using croping method
input=2; % select 1 or 2
phasor = exp(-1i.*double(imread([num2str(exampleNeuron(input)),'.png'])).*2*pi/255);


%%%%%% Run SPAF  %%%%%%%%%
[SPAFphasor,unwrappedPhase]=SPAFUnwrap(phasor,6,7);             % unwrapped by SPAF



%%%%%%  Display  %%%%%%%%%
displayPhase=imresize(abs(zerolevel(unwrappedPhase,5)),3);      %background zeroleveling and interpolation for display
figure(200),
mesh(displayPhase);colorbar;colormap(jet);zlim([-.5 20])%;caxis([0 11])
title('SPAF for Neuron');axis off 
h=colorbar;ylabel(h, 'Phase (rad)', 'FontSize', 16)




%%



%  Function #1 generate a wrapped noisy guassian phase

function [ex_phase,phase,gaussPhase]=WrapNoisyGaussPhase(size_x,size_y,width_gaus,wrap_order,snr_gaus,snr_speckle)

    if nargin<4;size_x=480;size_y=640;end               % image size
    if nargin<3;width_gaus=[size_x/4,size_y/4];end      % gaussian width
    if nargin<4;wrap_order=5;end                        % wrap value
    if nargin<5;snr_gaus=6;end                          % white gaussian SNR 
    if nargin<6;snr_speckle=2;end                       % speckle noise value

    [x,y]=meshgrid(1:size_x,1:size_y); % Creat space

    gaussPhase=exp(-(((x-size_x/2)./width_gaus(1)).^2+...
        ((y-size_y/2)./width_gaus(2)).^2)).*2.*pi.*wrap_order;  %creat Gaussian phase

    for jj=1:4          %phase shifting
        
        ima=cos(gaussPhase+(jj-1).*pi./2)+1;                    %interferogram
        ima_speck= imnoise(ima,'speckle', snr_speckle);         %add speckle noise
        N_phasor(:,:,jj)=awgn(ima_speck,snr_gaus,'measured');   %add Gaussian noise
        
    end
    
    ex_phase(:,:)=N_phasor(:,:,1)-N_phasor(:,:,3)+...
        1i.*(N_phasor(:,:,4)-N_phasor(:,:,2));                  %extract phasor 
    
    phase=angle(ex_phase);                                      %extract phase
    
end


%  Function #2 for displa

function zeroleve=zerolevel(imagex,m)
    n=m/100;% m=percentage of background;
    u=sort(imagex(:));
    v=u(1:round(length(u)*n));
    meanmin=mean(v);
    zeroleve=imagex-meanmin;
end

