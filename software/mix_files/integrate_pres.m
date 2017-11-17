function [disp,vel]=integrate_pres(cal,head,hpf_cutoff,number)
% function [disp,vel]=integrate_pres(cal,head,hpf_cutoff,number)
% integrate pressure data to get
% vertical cable velocities & displacements
% hpf_cutoff (optional) - hpf filter cutoff in Hz. Should be adjusted for every cruise!
% For scs07 0,02Hz worked well, for tao data - 0.04 seems to do better
% For shorter data (less than 2 min) the default value (0.04Hz) should be increased!
% number (optional) - number of points over which acceleration data is extrapolated
% on both sides to get rid of boundary problems with filtering
% Default value is 2 min worth of datapoints

if nargin<3
    hpf_cutoff=0.04;
end
if nargin<4
    number=head.samplerate(head.sensor_index.P)*120;
end

% dt=1./head.primary_sample_rate;sr=1/dt;
sr=head.samplerate(head.sensor_index.P); dt=1/sr;

[b,a]=butter(2,2*hpf_cutoff/sr,'high');

% extrapolate acc data
% extrapolate a VARIABLE over the
% extra NUMBER of points
% on both sides
len=length(cal.P);
flipsize=min(len,number);

in = cal.P * 6894.757 / 1027/ 9.81 ; % in m
if size(in,1)>1
    in=in';
end
fin=fliplr(in(1:flipsize)); fin=fin(1:end-1);
bin=fliplr(in(end-flipsize+1:end)); bin=bin(2:end);
in=[fin in  bin];

pf = filtfilt(b, a, in);
nsmooth = round(0.5/dt); % smooth over 1/2 sec (bit noise)
vel = filtfilt(b,a,smooth(diff(pf,1), nsmooth)/dt);
disp = filtfilt(b,a,cumtrapz(vel) * dt);

vel = vel(flipsize-1:flipsize+len-2);
disp = disp(flipsize-1:flipsize+len-2);

if size(cal.P) ~= size(vel)
    vel = vel'; disp = disp';
end
