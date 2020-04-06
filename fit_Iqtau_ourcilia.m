function [ Frequency, Damping, Amplitude, Offset, GOF] = fit_Iqtau_ourcilia( Iqtau, max_q, max_tau )
%fit_Iqtau_cilia Fits every row of the Iqtau matrix with a dampened
%oscillatory function


%{
% Version 1.0
% © Luigi Feriani, Maya Juenet, Pietro Cicuta, 2017 (lf352@cam.ac.uk) 
% 
% fit_Iqtau_ourcilia.m is licensed under a Creative Commons 
% Attribution-NonCommercial-NoDerivatives 4.0 International License.
% 
% Original work
% Feriani, L., et al., Biohpysical Journal 2017
% "Assessing the Collective Dynamics of Motile Cilia in Cultures of Human
% Airway Cells by Multiscale DDM"
% http://dx.doi.org/10.1016/j.bpj.2017.05.028
%
%}

if nargin < 3 || isempty(max_tau)
    max_tau = size(Iqtau,2);
end
if nargin < 2 || isempty(max_q)
    max_q = size(Iqtau,1);
end


Amplitude = nan(max_q, 1, 'double');
Frequency = nan(max_q, 1, 'double');
Damping   = nan(max_q, 1, 'double');
Offset    = nan(max_q, 1, 'double');
GOF       = nan(max_q, 1, 'double');

%% preliminary search for qstart
%|    _
% \  / \
%  \/   '.__
% Since the shape of the <Iqtau>tau is as above, looking for the negative
% peak and then mediating over q only from then on

qstart = [];
if max_q > 2 %otherwise error in findpeaks
    [~,qstart]=findpeaks(-1*mean(Iqtau(1:max_q,1:max_tau),2),'NPeaks',1);
end
if isempty(qstart) || qstart > max_q/2
    qstart = 1;
end


%% preliminary fit of the <Iqtau>q to get good parameters for the fit of each I(q',tau)
xx = (1:max_tau)';
yy = mean(Iqtau(qstart:max_q,1:max_tau))';

init_Offset = mean(yy(round(end/2):end));
init_Amplitude = abs(-yy(1)+init_Offset);
init_DecayTime = length(yy)/8;

autocorr = xcorr(yy-init_Offset);
dummy_freq = fftshift(fft( autocorr(length(yy):end) ));
dummy_freq = dummy_freq(ceil(length(dummy_freq)/2)+1:end);
[pk, init_f]=findpeaks(abs(dummy_freq));
init_f = init_f(pk == max(pk));
init_Freq = 2*pi*init_f/length(autocorr(length(yy):end));

% [fit_out_test, gof]=fminsearch('funz_fittest', [init_Amplitude init_Freq init_DecayTime init_Offset],...
%     optimset('MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',1e-11*init_Amplitude), yy,  xx, xx);
% % init_Amplitude   = fit_out_test(1);
% init_Frequency     = fit_out_test(2);
% init_DecayTime     = fit_out_test(3);
% % init_Offset      = fit_out_test(4);
%
%
% if init_Frequency <= 0
%     init_Frequency = init_Freq;
% end

ft = fittype('Ampl*(1-exp(cos(Freq*xx)))*exp(-Damp*xx)+Offset','Independent','xx');
fo = fitoptions('Method','NonLinearLeastSquare',...
    'StartPoint',[init_Amplitude, 1/init_DecayTime, init_Freq, init_Offset],...
    'Lower',[0 0 abs(init_Freq)/2 0],...
    'Upper',[max(yy)+init_Offset, Inf, abs(init_Freq)*2, max(yy)],...
    'MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',1e-9*init_Amplitude,'Weights',1./sqrt(xx));

alt_fo = fitoptions('Method','NonLinearLeastSquare',... %alternative fitoptions in case the first ones yield error
    'StartPoint',[init_Amplitude, 1/init_DecayTime, init_Freq, init_Offset],...
    'MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',1e-9*init_Amplitude,'Weights',1./sqrt(xx));

try
    [fit_output, gof] = fit(xx,yy,ft,fo);
catch
    
    if isempty(init_Freq)
        init_Freq = 5;      % absolutely random
    end
    alt_fo.StartPoint = [init_Amplitude, 1/init_DecayTime, init_Freq, init_Offset];
    [fit_output, gof] = fit(xx,yy,ft,alt_fo);
end

init_Frequency = fit_output.Freq;
if init_Frequency <= 0, init_Frequency = init_Freq; end
init_Damping = fit_output.Damp;

%     figure
%     plot(xx,yy,'.');hold on
%     plot(xx,fit_output(xx),'r');
%     pause

%% fit of the I(q',tau) curves
hpool = gcp; %if no parpool opened, opens a new one
parfor qq = 1:max_q
% for qq = 1:max_q
    
    xx = (1:max_tau)';
    yy = Iqtau(qq,1:max_tau)';
    
    %% initialisation of the fit parameters
    
    init_Offset = mean(yy(round(end/2):end));
    init_Amplitude = abs(-yy(1)+init_Offset);
    %     init_Damping = length(yy)/4;
    
    
    if init_Amplitude < 0.03*init_Offset, continue; end;       %if there's any signal
    
    
    fo = fitoptions('Method','NonLinearLeastSquare',...
        'StartPoint',[init_Amplitude, init_Damping, init_Frequency, init_Offset],...
        'Lower',[0 0 abs(init_Frequency)/2 0],...
        'Upper',[max(yy)+init_Offset, Inf, abs(init_Frequency)*2, max(yy)],...
        'MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',1e-9*init_Amplitude,'Weights',1./sqrt(xx));
    
    alt_fo = fitoptions('Method','NonLinearLeastSquare',... %alternative fitoptions in case the first ones yield error
        'StartPoint',[init_Amplitude, init_Damping, init_Frequency, init_Offset],...
        'MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',1e-9*init_Amplitude,'Weights',1./sqrt(xx));
    
    try
        [fit_output, gof] = fit(xx,yy,ft,fo);
    catch
        [fit_output, gof] = fit(xx,yy,ft,alt_fo);
    end
    
    Amplitude(qq) = fit_output.Ampl;
    Frequency(qq) = abs(fit_output.Freq);
    Damping(qq)   = fit_output.Damp;
    Offset(qq)    = fit_output.Offset;
        GOF(qq) = gof.rsquare;
     
    
    
end
Frequency = Frequency/(2*pi); 

end

