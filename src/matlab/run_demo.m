function records = run_demo()
% Synthetic public-safe cylindrical sector model; not operational radar data.
c.frequency_hz=2.4e9; c.c0=299792458; c.panel_count=20; c.row_count=12;
c.az_spacing_wavelength=0.52; c.vertical_spacing_wavelength=0.48;
c.seed=2026; c.trials=100; c.panel_correlation=0.30;
c.within_panel_correlation=0.96; c.displacement_std_wavelength=0.012;
c.mean_displacement_wavelength=0; c.radial_coupling=0.0002;
c.vertical_coupling=0.015; c.steer_az_deg=0; c.target_az_deg=1.5;
c.beam_el_deg=3; c.analog_el_deg=0; c.calibration_half_span_deg=5;
c.peak_half_span_deg=6;
g=geometry(c); rng(c.seed,'twister');
wa=weights(c,g,'steered_analog'); wb=weights(c,g,'fixed_analog');
cal=calibration(c,g.xyz,wb);
[reference_peak,reference_power]=peak(c,g.xyz,wb);
assert(isfinite(reference_peak));
trial=[]; scenario={}; beam_error_deg=[]; monopulse_error_deg=[]; relative_peak_db=[];
names={'nominal_steered','nominal_fixed','deformed_steered','deformed_fixed'};
for k=1:c.trials
    delta=displacement(c,g); actual=g.xyz+delta;
    for j=1:4
        if j==1, xyz=g.xyz; w=wa;
        elseif j==2, xyz=g.xyz; w=wb;
        elseif j==3, xyz=actual; w=wa;
        else, xyz=actual; w=wb; end
        [bp,power]=peak(c,xyz,w); est=estimate(c,xyz,w,cal);
        trial(end+1,1)=k-1; %#ok<AGROW>
        scenario{end+1,1}=names{j}; %#ok<AGROW>
        beam_error_deg(end+1,1)=bp-reference_peak; %#ok<AGROW>
        monopulse_error_deg(end+1,1)=est-c.target_az_deg; %#ok<AGROW>
        relative_peak_db(end+1,1)=10*log10(power/reference_power); %#ok<AGROW>
    end
end
records=table(trial,scenario,beam_error_deg,monopulse_error_deg,relative_peak_db);
if ~exist('results','dir'), mkdir('results'); end
writetable(records,fullfile('results','synthetic_demo_matlab.csv'));
fprintf('Synthetic demo: %d rows; radius %.4f m\n',height(records),g.radius);
end

function g=geometry(c)
assert(mod(c.panel_count,2)==0 && c.panel_count>=4 && c.row_count>=3);
g.lambda=c.c0/c.frequency_hz;
g.radius=c.panel_count*c.az_spacing_wavelength*g.lambda/(2*pi);
g.dz=c.vertical_spacing_wavelength*g.lambda;
ph=((0:c.panel_count-1)'-(c.panel_count-1)/2)*(2*pi/c.panel_count);
g.panel=repmat((1:c.panel_count)',c.row_count,1);
g.row=kron((1:c.row_count)',ones(c.panel_count,1));
g.phi=ph(g.panel); z=(g.row-(c.row_count+1)/2)*g.dz;
g.xyz=[g.radius*cos(g.phi),g.radius*sin(g.phi),z];
end

function delta=displacement(c,g)
n=c.panel_count; idx=(1:n)';
K=c.panel_correlation.^abs(idx-idx');
[V,D]=eig((K+K')/2);
shared=V*(sqrt(max(diag(D),0)).*randn(n,1));
q=sqrt(c.within_panel_correlation)*shared(g.panel)+sqrt(1-c.within_panel_correlation)*randn(numel(g.panel),1);
s=g.lambda*(c.mean_displacement_wavelength+c.displacement_std_wavelength*q);
r=c.radial_coupling*s; z=c.vertical_coupling*s;
delta=[-sin(g.phi).*s+cos(g.phi).*r,cos(g.phi).*s+sin(g.phi).*r,z];
end

function a=steering(c,xyz,az,el)
az=az(:)'; if isscalar(el), el=repmat(el,size(az)); end
el=el(:)'; u=[cosd(el).*cosd(az);cosd(el).*sind(az);sind(el)];
a=exp(1i*2*pi*c.frequency_hz/c.c0*(xyz*u));
end

function w=weights(c,g,architecture)
n=c.panel_count; m=c.row_count;
azwin=0.54-0.46*cos(2*pi*(0:n-1)'/(n-1));
% MATLAB taylorwin needs Signal Processing Toolbox; use normalized Hamming
% in BOTH implementations only if Taylor parity is explicitly configured.
elwin=taylor_local(m,min(4,m),30); elwin=elwin/sum(elwin);
if strcmp(architecture,'steered_analog'), analog_el=c.beam_el_deg;
elseif strcmp(architecture,'fixed_analog'), analog_el=c.analog_el_deg;
else, error('Unknown architecture'); end
k=2*pi/g.lambda;
horizontal=g.xyz(:,1)*cosd(c.steer_az_deg)+g.xyz(:,2)*sind(c.steer_az_deg);
analog=azwin(g.panel).*exp(-1i*k*cosd(analog_el)*horizontal);
zrow=((0:m-1)'-(m-1)/2)*g.dz;
digital=elwin.*exp(-1i*k*zrow*sind(c.beam_el_deg))/sum(azwin);
w.sum=analog.*digital(g.row);
signs=ones(size(g.panel)); signs(g.panel<=n/2)=-1;
w.diff=w.sum.*signs;
end

function v=taylor_local(N,nbar,sll)
A=acosh(10^(sll/20))/pi;
s2=nbar^2/(A^2+(nbar-0.5)^2); mm=1:nbar-1; F=zeros(size(mm));
for j=1:numel(mm)
    other=mm(mm~=mm(j));
    F(j)=(-1)^(j+1)*prod(1-mm(j)^2./(s2*(A^2+(mm-0.5).^2)))/(2*prod(1-mm(j)^2./other.^2));
end
nn=(0:N-1)'; v=1+2*cos(2*pi*(nn-N/2+0.5)*mm/N)*F';
end

function [s,d]=response(c,xyz,w,az,el)
a=steering(c,xyz,az,el); s=w.sum.'*a; d=w.diff.'*a;
end

function cal=calibration(c,xyz,w)
az=c.steer_az_deg+linspace(-c.calibration_half_span_deg,c.calibration_half_span_deg,1001);
[s,d]=response(c,xyz,w,az,c.beam_el_deg);
ratio=imag(d./s); valid=abs(s)>.1*max(abs(s)) & isfinite(ratio);
[~,center]=min(abs(az-c.steer_az_deg));
slope=sign(ratio(center+1)-ratio(center)); lo=center; hi=center;
while lo>1 && valid(lo-1) && slope*(ratio(lo)-ratio(lo-1))>1e-12, lo=lo-1; end
while hi<numel(az) && valid(hi+1) && slope*(ratio(hi+1)-ratio(hi))>1e-12, hi=hi+1; end
assert(hi-lo>=10,'No monotonic calibration region');
[cal.ratio,order]=sort(ratio(lo:hi)); subset=az(lo:hi); cal.az=subset(order);
cal.min_sum=.1*max(abs(s));
end

function est=estimate(c,xyz,w,cal)
[s,d]=response(c,xyz,w,c.target_az_deg,c.beam_el_deg);
v=imag(d/s);
if abs(s)<=cal.min_sum || ~isfinite(v) || v<cal.ratio(1) || v>cal.ratio(end)
    est=NaN;
else
    est=interp1(cal.ratio,cal.az,v,'linear');
end
end

function [bp,power]=peak(c,xyz,w)
grid=c.steer_az_deg+(-c.peak_half_span_deg:.25:c.peak_half_span_deg);
[s,~]=response(c,xyz,w,grid,c.beam_el_deg); [~,idx]=max(abs(s).^2);
if idx==1 || idx==numel(grid), bp=NaN; power=NaN; return; end
fun=@(az) power_neg(c,xyz,w,az);
bp=fminbnd(fun,grid(idx-1),grid(idx+1),optimset('TolX',1e-8));
power=-fun(bp);
end

function p=power_neg(c,xyz,w,az)
[s,~]=response(c,xyz,w,az,c.beam_el_deg); p=-abs(s).^2;
end
