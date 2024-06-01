function [input_t, input] = generate_OUinput_rectifiedsine(t, intervals, ninstances, inputDC, inputpk2base, inputfreq, inputphaseshift, thetamult, sigmamult)

nstates = numel(inputDC);

if ~exist('thetas','var')
    outhetas = {2*pi*16, 2*pi*1}; % Kk-SR 2023 OU noise inputs
end
if ~exist('sigmamult','var')
    sigmamult = {1,1};
end
if ~exist('thetamult','var')
    thetamult = {1,1};
end


%% generating input to OU_lowpass 

intervalindices = max((1:length(diff(intervals)))*(abs(diff(t(:)'>intervals(:),[],1))==1), 1); %@(t) find(abs(diff(t(:)'>intervals(:)))==1);
if isinf(intervals(1))
    intervals(1) = -1e100+1e10*rand(1);
end


meanip = [];
thetas = [];
cov_uncorrelatednoise = [];
for i=1:nstates
    thetas = [thetas; outhetas{i}(:).*thetamult{i}(:).*ones(ninstances,1)];
    cov_uncorrelatednoise = [cov_uncorrelatednoise; sigmamult{i}(:).*ones(ninstances,1)];
    ipDC = inputDC{i};
    ipf1 = inputpk2base{i};
    ipfreq = inputfreq{i};
    ipphase0 = inputphaseshift{i};
    meanip = [ ...
            meanip; 
            ones(ninstances,1) .*...
            ( ipDC( :, intervalindices)  ...
              +  ipf1( :, intervalindices) ...
                 .*  abs( ...
                         cos(ipphase0(:,intervalindices) + 2*pi*ipfreq(:,intervalindices).*(t-intervals(intervalindices))) ...
                        ) ...
            ) 
        ];
end
wnoise = cov_uncorrelatednoise.*meanip.*randn(size(meanip));
infind = isinf(thetas(:));
thetas_woinf = thetas(:); thetas_woinf(infind) = 1;

% setup and run OU_noise input generation model
OUpop = OU_noise_lopass_form(ninstances*nstates, {wnoise, meanip, t}, thetas_woinf(:));
OUpop.input([], t);
solver = @(updateFn, tlist, y0) eulerMethod(updateFn, tlist, y0);

OUpop.ode(solver);



%%
input_t = OUpop.EIpairs.t;
input   = OUpop.EIpairs.R;
input(infind,:) = meanip(infind,:) + wnoise(infind,:);
end