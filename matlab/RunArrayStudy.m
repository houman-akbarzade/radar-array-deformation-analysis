function results = RunArrayStudy(mode,c)
% شبیه‌سازی مستقل آزیموت؛ اجرا با حالت پایه یا تحلیل حساسیت
if nargin<1, mode='baseline'; end
if nargin<2, c=ArrayStudyConfig(); end
assert(any(strcmpi(mode,{'baseline','sensitivity','study','all'})),'Invalid mode.');
assert(c.selected_beam>=1 && c.selected_beam<=numel(c.stack_elevations));
c.beam_el=c.stack_elevations(c.selected_beam);
assert(c.nmc>=2 && c.nstats>=1000);
if ~exist(c.output_dir,'dir'), mkdir(c.output_dir); end
rng(c.seed,'twister');
g=create_array_geometry(c);
if strcmpi(mode,'study')
    results.studies=run_parametric_studies(c,g);
    results.strong=run_strong_case(c,g);
    return;
end
validate_generator(c,g);
results=run_monte_carlo(c,g);
plot_results(c,g,results);
if c.run_frequency, results.frequency=run_frequency_analysis(c,g); end
if strcmpi(mode,'all')
    results.studies=run_parametric_studies(c,g);
    results.strong=run_strong_case(c,g);
end
save(fullfile(c.output_dir,'baseline_results.mat'),'c','g','results');
writetable(results.summary,fullfile(c.output_dir,'scenario_summary.csv'));
disp(results.summary);
if strcmpi(mode,'sensitivity')
    sweep=run_sensitivity(c,g);
    save(fullfile(c.output_dir,'sensitivity_results.mat'),'sweep');
    writetable(sweep,fullfile(c.output_dir,'sensitivity_summary.csv'));
elseif ~any(strcmpi(mode,{'baseline','all'}))
    error('Invalid mode.');
end
fprintf('Results saved in %s\n',c.output_dir);
end

function g=create_array_geometry(c)
% Each row lists the 20 active panels; the 60-column cylinder remains implicit.
assert(numel(c.panel_az_deg)==c.panels);
[ph,z]=meshgrid(c.panel_az_deg*pi/180,((0:c.rows-1)-(c.rows-1)/2)*c.dz);
g.phi=reshape(ph.',[],1); g.z=reshape(z.',[],1);
g.xyz=[c.radius*cos(g.phi),c.radius*sin(g.phi),g.z];
g.panel=repmat((1:c.panels)',c.rows,1);
g.row=kron((1:c.rows)',ones(c.panels,1));
g.n=size(g.xyz,1);
assert(g.n==c.rows*c.panels && max(abs(hypot(g.xyz(:,1),g.xyz(:,2))-c.radius))<1e-12);
assert(isequal(g.panel(1:c.panels),(1:c.panels)'));
end

function m=create_deformation_model(c,g)
% همبستگی نمایی روی مسیر باز قطاع؛ بدون همسایگی دوری
assert(c.f>=c.f_min && c.f<=c.f_max);
assert(c.rho_panel>=0 && c.rho_panel<=1);
assert(c.rho_within>=0 && c.rho_within<=1);
assert(c.sigma_tip>=0 && c.phase_center_power>0);
arc=c.radius*c.panel_az_deg(:)*pi/180;
dist=abs(arc-arc');
if c.rho_panel==0
    K=eye(c.panels);
elseif c.rho_panel==1
    K=ones(c.panels);
else
    ell=median(diff(arc))/(-log(c.rho_panel)); K=exp(-dist/ell);
end
% تجزیه طیفی برای پشتیبانی از کوواریانس مثبت‌نیمه‌معین
[V,D]=eig((K+K')/2); d=diag(D);
assert(min(d)>-1e-10,'Invalid panel covariance.');
m.L=V*diag(sqrt(max(d,0))); m.K=K;
xi=(c.f-c.f_min)/(c.f_max-c.f_min);
m.xi=xi^c.phase_center_power;
if strcmp(c.bending_shape,'uniform')
    m.gain=m.xi^2*(6-4*m.xi+m.xi^2)/3;
elseif strcmp(c.bending_shape,'tip_force')
    m.gain=m.xi^2*(3-m.xi)/2;
else
    error('Unknown bending shape.');
end
m.scale=ones(g.n,1);
if strcmp(c.mechanical_model,'wind_parametric')
    % تصویر میدان مشترک روی جهت مماسی؛ مدل تجربی، نه پیش‌بینی فشار واقعی
    h=1+c.height_gradient*g.z/max(abs(g.z));
    m.scale=sin(c.wind_az_deg*pi/180-g.phi).*h;
elseif ~strcmp(c.mechanical_model,'base')
    error('Unknown mechanical model.');
end
m.c=c; m.g=g;
end

function d=generate_displacement(m)
% سه مؤلفه از یک میدان بار مشترک؛ نسبت‌ها فرض پارامتری هستند
c=m.c; g=m.g;
p=m.L*randn(c.panels,1);
q=sqrt(c.rho_within)*p(g.panel)+sqrt(1-c.rho_within)*randn(g.n,1);
s=m.gain*m.scale.*(c.mean_tip+c.sigma_tip*q);
r=c.radial_ratio*s; z=c.vertical_ratio*s;
d.s=s; d.r_local=r; d.z=z;
d.xyz=[-sin(g.phi).*s+cos(g.phi).*r,cos(g.phi).*s+sin(g.phi).*r,z];
new=g.xyz+d.xyz;
d.delta_radius=hypot(new(:,1),new(:,2))-c.radius;
d.delta_phi=angle(exp(1i*(atan2(new(:,2),new(:,1))-g.phi)));
end

function a=steering_response(c,g,xyz,az,el)
% پاسخ میدان برای جهت‌های متناظر؛ زاویه‌ها برحسب درجه‌اند
if nargin<5, el=c.target_el; end
az=az(:)'; el=el(:)';
if isscalar(el), el=repmat(el,size(az)); end
assert(numel(el)==numel(az));
u=[cosd(el).*cosd(az);cosd(el).*sind(az);sind(el)];
a=exp(1i*2*pi*c.f/c.c0*(xyz*u));
if strcmp(c.element_pattern,'cosine_front')
    ef=max(0,cos(g.phi-az*pi/180)).*cosd(el); a=a.*ef;
elseif ~strcmp(c.element_pattern,'isotropic')
    error('Unknown element pattern.');
end
end

function w=make_weights(c,g,architecture,beam_el)
% دو خروجی چپ و راست هر سطر؛ ضرایب مستقیم کامباینر با علامت منفی
if nargin<4, beam_el=c.beam_el; end
k=2*pi*c.f/c.c0;
azwin=0.54-0.46*cos(2*pi*(0:c.panels-1)'/(c.panels-1));
if strcmp(c.az_taper,'uniform'), azwin=ones(c.panels,1);
elseif ~strcmp(c.az_taper,'hamming'), error('Unknown azimuth taper.'); end
if strcmp(c.el_taper,'taylor'), elwin=taylor_local(c.rows,c.taylor_nbar,c.taylor_sll_db);
elseif strcmp(c.el_taper,'uniform'), elwin=ones(c.rows,1);
else, error('Unknown elevation taper.'); end
if strcmp(architecture,'A'), analog_el=beam_el;
elseif strcmp(architecture,'B'), analog_el=c.analog_el;
else, error('Unknown architecture.'); end
horizontal=g.xyz(:,1)*cosd(c.steer_az)+g.xyz(:,2)*sind(c.steer_az);
b=azwin(g.panel).*exp(-1i*k*cosd(analog_el)*horizontal);
% Two subarray channels per row: left, right; 32 channels total.
channel=2*(g.row-1)+1+(g.panel>c.panels/2);
w.combiner=sparse(channel,(1:g.n)',b,2*c.rows,g.n);
zrow=((0:c.rows-1)'-(c.rows-1)/2)*c.dz;
d=elwin.*exp(-1i*k*zrow*sind(beam_el));
w.digital_sum=kron(d,ones(2,1))/(sum(azwin)*sum(elwin));
w.digital_diff=w.digital_sum.*repmat([-1;1],c.rows,1);
% وزن معادل برای تطبیق خروجی معماری صریح با مدل المانی
w.sum=conj((w.digital_sum.'*w.combiner).');
w.diff=conj((w.digital_diff.'*w.combiner).');
w.architecture=architecture; w.beam_el=beam_el;
end

function v=taylor_local(N,nbar,sll)
% پنجره تیلور متقارن؛ sll عدد مثبتِ تضعیف لوب فرعی است
assert(nbar>=2 && nbar<=N && nbar==round(nbar) && sll>0);
A=acosh(10^(sll/20))/pi;
s2=nbar^2/(A^2+(nbar-0.5)^2); m=1:nbar-1; F=zeros(size(m));
for j=1:numel(m)
    other=m(m~=m(j));
    F(j)=(-1)^(j+1)*prod(1-m(j)^2./(s2*(A^2+(m-0.5).^2)))/(2*prod(1-m(j)^2./other.^2));
end
n=(0:N-1)'; v=1+2*cos(2*pi*(n-N/2+0.5)*m/N)*F';
v=v/(1+2*sum(F));
end

function [s,t]=calculate_array_response(c,g,xyz,w,az,el)
if nargin<7, el=c.target_el; end
a=steering_response(c,g,xyz,az,el);
channels=w.combiner*a;
s=w.digital_sum.'*channels; t=w.digital_diff.'*channels;
end

function cal=make_calibration(c,g,w)
az=c.steer_az+c.cal_offsets;
[s,t]=calculate_array_response(c,g,g.xyz,w,az,c.calibration_el);
v=imag(t./s);
valid=abs(s)>0.1*max(abs(s)) & isfinite(v);
dv=diff(v); k=find(abs(az-c.steer_az)==min(abs(az-c.steer_az)),1);
sg=sign(dv(min(k,numel(dv))));
assert(sg~=0 && valid(k),'Invalid monopulse calibration near boresight.');
lo=k; hi=k;
while lo>1 && valid(lo-1) && sg*dv(lo-1)>1e-12, lo=lo-1; end
while hi<numel(az) && valid(hi+1) && sg*dv(hi)>1e-12, hi=hi+1; end
assert(hi-lo>=10,'Insufficient monotone calibration interval.');
cal.az=az(lo:hi); cal.ratio=v(lo:hi);
cal.min_sum=0.1*max(abs(s)); cal.target_el=c.calibration_el;
[cal.ratio,ix]=sort(cal.ratio); cal.az=cal.az(ix);
assert(c.target_az>=min(cal.az) && c.target_az<=max(cal.az),'Target outside calibration interval.');
end

function [estimate,valid]=calculate_monopulse_estimation(c,g,xyz,w,cal)
[s,t]=calculate_array_response(c,g,xyz,w,c.target_az);
v=imag(t/s);
valid=abs(s)>cal.min_sum && isfinite(v) && v>=cal.ratio(1) && v<=cal.ratio(end);
if valid, estimate=interp1(cal.ratio,cal.az,v,'linear'); else, estimate=NaN; end
end

function [peak,power,boundary]=calculate_beam_peak(c,g,xyz,w)
grid=c.steer_az+(-c.peak_half_span:c.peak_grid_step:c.peak_half_span);
[s,~]=calculate_array_response(c,g,xyz,w,grid,c.beam_el); [~,i]=max(abs(s).^2);
boundary=i==1 || i==numel(grid);
if boundary, peak=NaN; power=NaN; return; end
objective=@(az) -sum_power(c,g,xyz,w,az);
peak=fminbnd(objective,grid(i-1),grid(i+1),optimset('TolX',1e-8,'Display','off'));
power=-objective(peak);
end

function p=sum_power(c,g,xyz,w,az)
[s,~]=calculate_array_response(c,g,xyz,w,az,c.beam_el); p=abs(s).^2;
end

function out=run_monte_carlo(c,g)
m=create_deformation_model(c,g);
wA=make_weights(c,g,'A'); wB=make_weights(c,g,'B');
calB=make_calibration(c,g,wB);
[bB,pB]=calculate_beam_peak(c,g,g.xyz,wB);
assert(isfinite(bB),'Reference B peak outside search interval.');
% Check channel combiner against element-domain weights.
av=steering_response(c,g,g.xyz,[c.steer_az c.target_az],c.target_el);
[scheck,tcheck]=calculate_array_response(c,g,g.xyz,wB,[c.steer_az c.target_az]);
assert(max(abs(scheck-wB.sum'*av))<1e-12 && max(abs(tcheck-wB.diff'*av))<1e-12,'Channel equivalence failed.');
% وزن‌های A و B در ارتفاع صفر باید یکسان باشند
cz=c; cz.analog_el=0; wa0=make_weights(cz,g,'A',0); wb0=make_weights(cz,g,'B',0);
assert(max(abs(wa0.sum-wb0.sum))<1e-12,'Zero-elevation equivalence failed.');
names={'A_IdealSingle','B_NominalStack','C_DeformedSingle','D_DeformedStack'};
out.beam=nan(c.nmc,4); out.mono=nan(c.nmc,4); out.gain_db=nan(c.nmc,4);
out.invalid=false(c.nmc,4); out.boundary=false(c.nmc,4); out.s=zeros(g.n,c.nmc);
out.element_az_shift_deg=zeros(g.n,c.nmc);
% حالت‌های بدون تغییر شکل ثابت‌اند و فقط یک بار محاسبه می‌شوند
for j=1:2
    if j==1, w=wA; else, w=wB; end
    [b,p,edge]=calculate_beam_peak(c,g,g.xyz,w);
    [est,ok]=calculate_monopulse_estimation(c,g,g.xyz,w,calB);
    out.beam(:,j)=b-bB; out.mono(:,j)=est-c.target_az;
    out.gain_db(:,j)=10*log10(p/pB); out.invalid(:,j)=~ok; out.boundary(:,j)=edge;
end
for imc=1:c.nmc
    d=generate_displacement(m); out.s(:,imc)=d.s; xyz=g.xyz+d.xyz;
    out.element_az_shift_deg(:,imc)=d.delta_phi*180/pi;
    for j=3:4
        if j==3, w=wA; else, w=wB; end
        [b,p,edge]=calculate_beam_peak(c,g,xyz,w);
        [est,ok]=calculate_monopulse_estimation(c,g,xyz,w,calB);
        out.beam(imc,j)=b-bB; out.mono(imc,j)=est-c.target_az;
        out.gain_db(imc,j)=10*log10(p/pB);
        out.invalid(imc,j)=~ok; out.boundary(imc,j)=edge;
    end
    if imc==1, out.example=d; end
end
out.names=names; out.reference_peak_B=bB; out.reference_power_B=pB; out.cal=calB;
out.summary=summary_table(out);
out.correlation_CD=corr_pair(out.beam(:,4),out.mono(:,4));
fprintf('Target: azimuth %.3f deg, elevation %.3f deg\n',c.target_az,c.target_el);
fprintf('Selected digital beam: %.3f deg; analog elevation: %.3f deg\n',c.beam_el,c.analog_el);
fprintf('Calibration B elevation: %.3f deg; interval [%.3f, %.3f] deg\n',c.calibration_el,min(calB.az),max(calB.az));
fprintf('Reference B beam azimuth: %.9f deg\n',bB);
fprintf('Scenario D signed error correlation: %.6f\n',out.correlation_CD);
fprintf('D invalid estimates: %d; boundary peaks: %d\n',sum(out.invalid(:,4)),sum(out.boundary(:,4)));
end

function t=summary_table(o)
Scenario={}; Metric={}; Bias=[]; StdPopulation=[]; RMSE=[]; P025=[]; P975=[];
MeanCI_Low=[]; MeanCI_High=[]; InvalidCount=[];
for j=1:4
    for k=1:2
        if k==1, x=o.beam(:,j); label='Beam_deg'; else, x=o.mono(:,j); label='Monopulse_deg'; end
        nbad=sum(~isfinite(x)); x=x(isfinite(x)); n=numel(x);
        if n==0, mu=NaN; sd=NaN; rm=NaN; pc=[NaN NaN]; ci=[NaN NaN];
        else
            mu=mean(x); sd=std(x,1); rm=sqrt(mean(x.^2)); pc=percentiles(x,[.025 .975]);
            if n>1, se=std(x,0)/sqrt(n); else, se=NaN; end
            ci=mu+[-1 1]*1.96*se;
            assert(abs(rm^2-mu^2-sd^2)<1e-9);
        end
        Scenario{end+1,1}=o.names{j}; Metric{end+1,1}=label;
        Bias(end+1,1)=mu; StdPopulation(end+1,1)=sd; RMSE(end+1,1)=rm;
        P025(end+1,1)=pc(1); P975(end+1,1)=pc(2);
        MeanCI_Low(end+1,1)=ci(1); MeanCI_High(end+1,1)=ci(2); InvalidCount(end+1,1)=nbad;
    end
end
t=table(Scenario,Metric,Bias,StdPopulation,RMSE,P025,P975,MeanCI_Low,MeanCI_High,InvalidCount);
end

function p=percentiles(x,q)
x=sort(x(:));
if numel(x)==1, p=repmat(x,size(q));
else, p=interp1(1:numel(x),x,1+(numel(x)-1)*q,'linear'); end
end

function r=corr_pair(x,y)
ok=isfinite(x)&isfinite(y); x=x(ok); y=y(ok);
if numel(x)<3 || std(x)==0 || std(y)==0, r=NaN; else, R=corrcoef(x,y); r=R(1,2); end
end

function validate_generator(c,g)
% این آزمون در فرکانس بالا انجام می‌شود تا حالت جابه‌جایی صفر مانع آزمون نباشد
v=c; v.f=c.f_max; v.mechanical_model='base'; v.mean_tip=0;
if v.sigma_tip==0, v.sigma_tip=0.02*v.c0/v.f_max; end
m=create_deformation_model(v,g); rng(c.seed+1,'twister');
x=zeros(g.n,c.nstats);
for i=1:c.nstats, d=generate_displacement(m); x(:,i)=d.s; end
sd=std(x,1,2); ratio=sd/v.sigma_tip;
R=corrcoef(x'); expected=c.rho_within*m.K(g.panel,g.panel)+(1-c.rho_within)*eye(g.n);
fprintf('Generator mean SD / expected SD: %.6f\n',mean(ratio));
fprintf('Generator max correlation error: %.6f\n',max(abs(R(:)-expected(:))));
assert(abs(mean(ratio)-1)<0.04,'Generator variance check failed.');
assert(max(abs(R(:)-expected(:)))<0.08,'Generator correlation check failed.');
% بررسی استقلال نمونه‌های متوالی به‌صورت تشخیصی
fprintf('Generator lag-one correlation: %.6f\n',corr_pair(x(1,1:end-1)',x(1,2:end)'));
f=figure('Name','Generator validation');
subplot(1,2,1); imagesc(expected); axis image; colorbar; title('Expected element correlation');
subplot(1,2,2); imagesc(R); axis image; colorbar; title('Empirical element correlation');
save_figure(c,f,'01_generator_validation');
save(fullfile(c.output_dir,'generator_validation.mat'),'ratio','R','expected');
rng(c.seed,'twister');
end

function plot_results(c,g,o)
d=o.example; xyz=g.xyz+d.xyz; wA=make_weights(c,g,'A'); wB=make_weights(c,g,'B');
f=figure('Name','Geometry');
scatter3(g.xyz(:,1),g.xyz(:,2),g.xyz(:,3),12,'b','filled'); hold on;
scatter3(xyz(:,1),xyz(:,2),xyz(:,3),12,'r');
quiver3(g.xyz(:,1),g.xyz(:,2),g.xyz(:,3),d.xyz(:,1),d.xyz(:,2),d.xyz(:,3),0,'k');
axis equal; grid on; xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
legend('Nominal','Deformed','Displacement, true scale'); title('Array geometry'); save_figure(c,f,'02_geometry');
f=figure('Name','Displacements');
subplot(2,3,1); plot(c.panel_az_deg,mean(reshape(d.s,c.panels,c.rows),2)*1e3,'o-'); xlabel('Panel azimuth (deg)'); ylabel('Mean dS (mm)'); grid on;
labels={'Tangential (mm)','Radial local (mm)','Vertical (mm)','Radius change (mm)','Azimuth change (deg)'};
data={d.s*1e3,d.r_local*1e3,d.z*1e3,d.delta_radius*1e3,d.delta_phi*180/pi};
for j=1:5, subplot(2,3,j+1); histogram(data{j},20); xlabel(labels{j}); ylabel('Element count'); end
save_figure(c,f,'03_displacement_one_realization');
az=c.steer_az+c.pattern_offsets;
ws={wA,wB,wA,wB}; positions={g.xyz,g.xyz,xyz,xyz};
S=zeros(4,numel(az)); T=S;
for j=1:4, [S(j,:),T(j,:)]=calculate_array_response(c,g,positions{j},ws{j},az); end
f=figure('Name','Pattern');
plot(az,20*log10(max(abs(S.'),1e-8)/sqrt(o.reference_power_B))); grid on;
xlabel('Azimuth (deg)'); ylabel('Array response (dB relative to nominal B peak)'); ylim([-60 5]);
legend('A','B','C','D'); title(sprintf('Azimuth cut at target elevation %.1f deg',c.target_el)); save_figure(c,f,'04_pattern');
f=figure('Name','Monopulse'); hold on;
for j=1:4
    valid=abs(S(j,:))>o.cal.min_sum & az>=min(o.cal.az) & az<=max(o.cal.az);
    plot(az(valid),imag(T(j,valid)./S(j,valid)));
end
grid on; xlabel('Azimuth (deg)'); ylabel('Imaginary difference / sum'); legend('A','B reference','C','D'); save_figure(c,f,'05_monopulse_curve');
f=figure('Name','Errors');
subplot(1,2,1); plot_error_histograms(o.beam(:,3:4)); xlabel('Beam pointing error (deg)'); ylabel('Count'); legend('C','D');
subplot(1,2,2); plot_error_histograms(o.mono(:,3:4)); xlabel('Monopulse error (deg)'); ylabel('Count'); legend('C','D'); save_figure(c,f,'06_errors');
f=figure('Name','Error relationship');
scatter(o.beam(:,4),o.mono(:,4),10,'filled'); grid on;
xlabel('Beam pointing error (deg)'); ylabel('Monopulse estimation error (deg)');
title(sprintf('Scenario D relative to B, correlation %.3f',o.correlation_CD)); save_figure(c,f,'07_error_relationship');
f=figure('Name','Convergence');
for j=3:4
    x=o.mono(:,j); ok=isfinite(x); x(~ok)=0;
    den=cumsum(ok); y=sqrt(cumsum(x.^2)./max(den,1)); plot(1:c.nmc,y); hold on;
end
xlabel('Realization count'); ylabel('Conditional monopulse RMSE (deg)'); grid on; legend('C','D'); save_figure(c,f,'08_convergence');
f=figure('Name','Scenarios');
bar(reshape(o.summary.RMSE,2,4)'); set(gca,'XTick',1:4,'XTickLabel',{'A','B','C','D'});
ylabel('RMSE relative to reference B (deg)'); legend('Beam','Monopulse'); grid on; save_figure(c,f,'09_scenarios');
plot_beams_3d(c,g,o);
end

function sweep=run_sensitivity(c,g)
% تثبیت بذر میان نقاط برای مقایسه زوجی؛ در هر نقطه تحقق‌ها مستقل‌اند
axes_names={'sigma_tip','f','rho_within','rho_panel','target_az','target_el'};
values={c.c0/c.f_max*[0 .005 .01 .02 .03 .05 .07 .10],linspace(c.f_min,c.f_max,5),[0 .5 .9 .995 1],[0 .35 .7 .95 1],[-1.5 -1 0 0.5 1 1.5],[0 2 4 6 10]};
Parameter={}; Value=[]; Scenario={}; BeamRMSE=[]; MonoRMSE=[]; InvalidCount=[];
for ia=1:numel(axes_names)
    vb=values{ia}; rb=zeros(numel(vb),4); rm=rb;
    for iv=1:numel(vb)
        v=c; v.(axes_names{ia})=vb(iv);
        if strcmp(axes_names{ia},'target_el'), v.beam_el=vb(iv); v.calibration_el=vb(iv); end
        rng(c.seed,'twister'); o=run_monte_carlo(v,g);
        rb(iv,:)=o.summary.RMSE(1:2:end)'; rm(iv,:)=o.summary.RMSE(2:2:end)';
        for j=1:4
            Parameter{end+1,1}=axes_names{ia}; Value(end+1,1)=vb(iv); Scenario{end+1,1}=o.names{j};
            BeamRMSE(end+1,1)=rb(iv,j); MonoRMSE(end+1,1)=rm(iv,j); InvalidCount(end+1,1)=sum(o.invalid(:,j));
        end
    end
    f=figure('Name',axes_names{ia});
    subplot(1,2,1); plot(vb,rb,'o-'); xlabel(axes_names{ia},'Interpreter','none'); ylabel('Beam RMSE (deg)'); grid on; legend('A','B','C','D');
    subplot(1,2,2); plot(vb,rm,'o-'); xlabel(axes_names{ia},'Interpreter','none'); ylabel('Monopulse RMSE (deg)'); grid on; legend('A','B','C','D');
    save_figure(c,f,['sensitivity_' axes_names{ia}]);
end
sweep=table(Parameter,Value,Scenario,BeamRMSE,MonoRMSE,InvalidCount);
end

function save_figure(c,f,name)
saveas(f,fullfile(c.output_dir,[name '.png']));
savefig(f,fullfile(c.output_dir,[name '.fig']));
end

function plot_error_histograms(x)
% مرز بازه‌ها مشترک و رسم خطی است تا رنگ سوم حاصل از هم‌پوشانی ایجاد نشود
v=x(isfinite(x));
if isempty(v), text(0.1,0.5,'No valid estimates'); return; end
lo=min(v); hi=max(v);
if hi==lo, lo=lo-1e-6; hi=hi+1e-6; end
edges=linspace(lo,hi,26);
histogram(x(:,1),edges,'DisplayStyle','stairs','LineWidth',1.5); hold on;
histogram(x(:,2),edges,'DisplayStyle','stairs','LineWidth',1.5);
end

function plot_beams_3d(c,g,o)
% نمودار سه‌بعدی پاسخ روی آزیموت و ارتفاع؛ شعاع فیزیکی یا گین مطلق نیست
[AZ,EL]=meshgrid(c.plot_az,c.plot_el); az=AZ(:)'; el=EL(:)';
xyz1=g.xyz+o.example.xyz; xyzs={g.xyz,g.xyz,xyz1,xyz1};
patterns=cell(1,4); envelopes=cell(1,4); centres=cell(1,4);
ref=sqrt(o.reference_power_B);
f=figure('Name','3D receive patterns A B C D','Position',[50 50 1400 900]);
for j=1:4
    if j==1 || j==3, arch='A'; beams=c.beam_el;
    else, arch='B'; beams=c.stack_elevations; end
    % میدان المان‌ها یک‌بار برای همه پرتوهای این سناریو ساخته می‌شود
    a=steering_response(c,g,xyzs{j},az,el);
    db=zeros(numel(beams),numel(az));
    for ib=1:numel(beams)
        w=make_weights(c,g,arch,beams(ib));
        response=w.digital_sum.'*(w.combiner*a);
        db(ib,:)=20*log10(max(abs(response)/ref,1e-10));
    end
    patterns{j}=reshape(db.',size(AZ,1),size(AZ,2),numel(beams));
    envelope=reshape(max(db,[],1),size(AZ)); envelopes{j}=envelope; centres{j}=beams;
    subplot(2,2,j);
    surf(AZ,EL,max(envelope,c.plot_floor_db),'EdgeColor','none'); hold on;
    % مرز سه‌دسی‌بل هر پرتو جداگانه؛ پوش، جمع همدوس یا توان کل نیست
    for ib=1:numel(beams)
        zbeam=patterns{j}(:,:,ib); level=max(zbeam(:))-3;
        C=contourc(c.plot_az,c.plot_el,zbeam,[level level]); pos=1;
        while pos<size(C,2)
            n=C(2,pos); xy=C(:,pos+1:pos+n);
            plot3(xy(1,:),xy(2,:),repmat(level+0.05,1,n),'k-','LineWidth',1);
            pos=pos+n+1;
        end
        [val,ip]=max(zbeam(:));
        text(AZ(ip),EL(ip),val+1,sprintf('%g deg',beams(ib)),'FontSize',8);
    end
    plot3(c.target_az,c.target_el,2,'rp','MarkerFaceColor','r','MarkerSize',9);
    xlabel('Azimuth (deg)'); ylabel('Elevation (deg)'); zlabel('Response relative to B peak (dB)');
    if j==1, title('A: ideal single receive beam');
    elseif j==2, title('B: nominal receive stack, strongest-beam envelope');
    elseif j==3, title('C: deformed single receive beam');
    else, title('D: deformed receive stack, strongest-beam envelope'); end
    caxis([c.plot_floor_db 2]); zlim([c.plot_floor_db 3]); colorbar; view(40,35); grid on;
end
save_figure(c,f,'10_receive_patterns_3d');
reference='All dB referenced to nominal B selected-beam peak; stack envelope is max across separate beams, not their sum.';
save(fullfile(c.output_dir,'receive_patterns_3d.mat'),'AZ','EL','patterns','envelopes','centres','reference');
end

function fresult=run_frequency_analysis(c,g)
% شدت مکانیکی نوک ثابت است؛ بذر مشترک، همان تحقق باد را در فرکانس‌های مختلف مقایسه می‌کند
freq=c.frequency_list(:); assert(all(freq>=c.f_min & freq<=c.f_max));
assert(c.nmc_frequency>=500,'Frequency analysis requires at least 500 realizations.');
previous_rng=rng; cleanup=onCleanup(@() rng(previous_rng)); %#ok<NASGU>
fresult.frequency_hz=freq; fresult.runs=cell(numel(freq),1);
meanS=zeros(size(freq)); stdS=meanS; meanPhi=meanS; stdPhi=meanS;
gain=meanS; theoryMean=meanS; theoryStd=meanS;
beam=zeros(numel(freq),4); mono=beam; bias=beam; gainMean=beam;
F_Hz=[]; summary=[];
for i=1:numel(freq)
    v=c; v.f=freq(i); v.nmc=c.nmc_frequency;
    rng(c.seed,'twister');
    fprintf('\nFrequency %d / %d: %.4f GHz\n',i,numel(freq),freq(i)/1e9);
    o=run_monte_carlo(v,g); m=create_deformation_model(v,g);
    % ابتدا آمار هر المان در مجموعه تحقق‌ها، سپس میانگین بین المان‌ها
    meanS(i)=mean(mean(o.s,2)); stdS(i)=mean(std(o.s,1,2));
    meanPhi(i)=mean(mean(o.element_az_shift_deg,2));
    stdPhi(i)=mean(std(o.element_az_shift_deg,1,2));
    gain(i)=m.gain;
    theoryMean(i)=mean(m.gain*m.scale*c.mean_tip);
    theoryStd(i)=mean(abs(m.gain*m.scale)*c.sigma_tip);
    beam(i,:)=o.summary.RMSE(1:2:end)'; mono(i,:)=o.summary.RMSE(2:2:end)';
    bias(i,:)=o.summary.Bias(2:2:end)';
    for j=1:4
        x=o.gain_db(:,j); x=x(isfinite(x));
        if isempty(x), gainMean(i,j)=NaN; else, gainMean(i,j)=mean(x); end
    end
    fresult.runs{i}=o;
    if i==1, summary=o.summary; else, summary=[summary;o.summary]; end %#ok<AGROW>
    F_Hz=[F_Hz;repmat(freq(i),height(o.summary),1)]; %#ok<AGROW>
end
summary=[table(F_Hz),summary];
EffectiveGain=gain; MeanTangential_m=meanS; MeanElementSDTangential_m=stdS;
TheoryMeanTangential_m=theoryMean; TheoryMeanElementSDTangential_m=theoryStd;
MeanElementAzShift_deg=meanPhi; MeanElementSDAzShift_deg=stdPhi;
Frequency_Hz=freq;
fresult.displacement_summary=table(Frequency_Hz,EffectiveGain,MeanTangential_m,MeanElementSDTangential_m,TheoryMeanTangential_m,TheoryMeanElementSDTangential_m,MeanElementAzShift_deg,MeanElementSDAzShift_deg);
fresult.error_summary=summary;
fresult.beam_rmse=beam; fresult.mono_rmse=mono; fresult.mono_bias=bias;
fresult.mean_peak_response_db_relative_B=gainMean;
writetable(summary,fullfile(c.output_dir,'frequency_error_summary.csv'));
writetable(fresult.displacement_summary,fullfile(c.output_dir,'frequency_displacement_summary.csv'));
save(fullfile(c.output_dir,'frequency_results.mat'),'c','fresult');
f=figure('Name','Frequency-dependent effective displacement','Position',[80 80 1100 750]);
subplot(2,2,1); plot(freq/1e9,meanS*1e3,'o-',freq/1e9,theoryMean*1e3,'--');
xlabel('Frequency (GHz)'); ylabel('Mean tangential displacement (mm)'); legend('Empirical','Model'); grid on;
subplot(2,2,2); plot(freq/1e9,stdS*1e3,'o-',freq/1e9,theoryStd*1e3,'--');
xlabel('Frequency (GHz)'); ylabel('Mean per-element SD of dS (mm)'); legend('Empirical','Model'); grid on;
subplot(2,2,3); plot(freq/1e9,meanPhi,'o-'); xlabel('Frequency (GHz)'); ylabel('Mean element azimuth shift (deg)'); grid on;
subplot(2,2,4); plot(freq/1e9,stdPhi,'o-'); xlabel('Frequency (GHz)'); ylabel('Mean per-element SD of azimuth shift (deg)'); grid on;
save_figure(c,f,'11_frequency_displacement_statistics');
f=figure('Name','Frequency-dependent radar errors','Position',[80 80 1100 750]);
subplot(2,2,1); plot(freq/1e9,beam,'o-'); xlabel('Frequency (GHz)'); ylabel('Azimuth beam RMSE relative to B (deg)'); legend('A','B','C','D'); grid on;
subplot(2,2,2); plot(freq/1e9,mono,'o-'); xlabel('Frequency (GHz)'); ylabel('Monopulse RMSE (deg)'); legend('A','B','C','D'); grid on;
subplot(2,2,3); plot(freq/1e9,bias,'o-'); xlabel('Frequency (GHz)'); ylabel('Monopulse bias (deg)'); legend('A','B','C','D'); grid on;
subplot(2,2,4); plot(freq/1e9,gainMean,'o-'); xlabel('Frequency (GHz)'); ylabel('Mean peak response relative to B (dB)'); legend('A','B','C','D'); grid on;
save_figure(c,f,'12_frequency_error_comparison');
disp(fresult.displacement_summary);
end

function study=run_parametric_studies(c,g)
% جاروب تک‌پارامتری در چند فرکانس؛ پارامترهای دیگر در مقدار پایه ثابت‌اند
assert(c.nmc_study>=500,'Study requires at least 500 realizations.');
assert(all(c.study_frequencies>=c.f_min & c.study_frequencies<=c.f_max));
assert(all(c.sigma_tip_list>=0));
state=rng; cleanup=onCleanup(@() rng(state)); %#ok<NASGU>
root=fullfile(c.output_dir,'parametric_studies');
if ~exist(root,'dir'), mkdir(root); end
params={'sigma_tip','mean_tip','target_el'};
lists={c.sigma_tip_list,c.mean_tip_list,c.target_el_list};
study=struct();
for ip=1:3
    list=lists{ip}; nf=numel(c.study_frequencies); nv=numel(list);
    beam=zeros(nv,nf,4); mono=beam; bias=beam; invalid=beam; edges=beam;
    gain=beam; beamBias=beam; summary=[]; frequency=[]; value=[];
    folder=fullfile(root,params{ip}); if ~exist(folder,'dir'), mkdir(folder); end
    for jf=1:nf
        for iv=1:nv
            v=c; v.f=c.study_frequencies(jf); v.nmc=c.nmc_study;
            v.(params{ip})=list(iv);
            if strcmp(params{ip},'target_el')
                % آزمایش هدف در مرکز پرتو: وزن A و وزن دیجیتال B با ارتفاع عوض می‌شوند
                v.beam_el=v.target_el; v.calibration_el=v.target_el;
            end
            rng(c.seed,'twister');
            fprintf('\nStudy %s: point %d/%d, frequency %.3f GHz, value %.6g\n',params{ip},iv,nv,v.f/1e9,list(iv));
            o=run_monte_carlo(v,g);
            beam(iv,jf,:)=reshape(o.summary.RMSE(1:2:end),1,1,4);
            mono(iv,jf,:)=reshape(o.summary.RMSE(2:2:end),1,1,4);
            bias(iv,jf,:)=reshape(o.summary.Bias(2:2:end),1,1,4);
            beamBias(iv,jf,:)=reshape(o.summary.Bias(1:2:end),1,1,4);
            for j=1:4
                invalid(iv,jf,j)=mean(o.invalid(:,j))*100;
                edges(iv,jf,j)=mean(o.boundary(:,j))*100;
                x=o.gain_db(:,j); x=x(isfinite(x));
                if isempty(x), gain(iv,jf,j)=NaN; else, gain(iv,jf,j)=mean(x); end
            end
            % داده هر نقطه فوراً ذخیره می‌شود تا اجرای طولانی قابل بررسی باشد
            save(fullfile(folder,sprintf('point_f%02d_v%02d.mat',jf,iv)),'v','o');
            if isempty(summary), summary=o.summary; else, summary=[summary;o.summary]; end %#ok<AGROW>
            frequency=[frequency;repmat(v.f,height(o.summary),1)]; %#ok<AGROW>
            value=[value;repmat(list(iv),height(o.summary),1)]; %#ok<AGROW>
        end
    end
    ParameterValue=value; Frequency_Hz=frequency;
    summary=[table(ParameterValue,Frequency_Hz),summary];
    writetable(summary,fullfile(folder,'statistics.csv'));
    data.parameter=params{ip}; data.values=list; data.frequency_hz=c.study_frequencies;
    data.beam_rmse=beam; data.mono_rmse=mono; data.mono_bias=bias;
    data.beam_bias=beamBias; data.invalid_percent=invalid; data.boundary_peak_percent=edges;
    data.mean_peak_response_db=gain; data.summary=summary;
    study.(params{ip})=data;
    if ip<3, x=list*1e3; unit=[params{ip} ' (mm)']; else, x=list; unit='Target and beam elevation (deg)'; end
    f=figure('Name',['Study ' params{ip}],'Position',[70 70 1200 850]);
    % نمودار اصلی B و D؛ برای هر فرکانس رنگ مشترک و برای معماری نوع خط متفاوت
    metrics={beam,mono,bias,invalid}; labels={'Beam RMSE relative to B (deg)','Monopulse RMSE (deg)','Monopulse bias (deg)','Invalid monopulse estimates (%)'};
    cols=lines(nf);
    for k=1:4
        subplot(2,2,k); hold on; legends={};
        for jf=1:nf
            plot(x,metrics{k}(:,jf,2),'--','Color',cols(jf,:),'LineWidth',1.2);
            plot(x,metrics{k}(:,jf,4),'-o','Color',cols(jf,:),'LineWidth',1.4);
            legends{end+1}=sprintf('B %.3f GHz',c.study_frequencies(jf)/1e9); %#ok<AGROW>
            legends{end+1}=sprintf('D %.3f GHz',c.study_frequencies(jf)/1e9); %#ok<AGROW>
        end
        xlabel(unit,'Interpreter','none'); ylabel(labels{k}); grid on; legend(legends,'Location','best');
    end
    save_figure(c,f,sprintf('%02d_study_%s_multifrequency',12+ip,params{ip}));
    % مقایسه A تا D در بالای باند؛ اختلاف معماری و پاسخ قله جدا نمایش داده می‌شوند
    [~,jf]=max(c.study_frequencies);
    f=figure('Name',['Architectures ' params{ip}],'Position',[80 80 1100 750]);
    subplot(2,2,1); plot(x,reshape(beam(:,jf,:),nv,4),'o-'); xlabel(unit,'Interpreter','none'); ylabel('Beam RMSE relative to B (deg)'); legend('A','B','C','D'); grid on;
    subplot(2,2,2); plot(x,reshape(mono(:,jf,:),nv,4),'o-'); xlabel(unit,'Interpreter','none'); ylabel('Monopulse RMSE (deg)'); legend('A','B','C','D'); grid on;
    subplot(2,2,3); plot(x,reshape(gain(:,jf,:),nv,4),'o-'); xlabel(unit,'Interpreter','none'); ylabel('Mean peak response relative to B (dB)'); legend('A','B','C','D'); grid on;
    subplot(2,2,4); plot(x,reshape(edges(:,jf,:),nv,4),'o-'); xlabel(unit,'Interpreter','none'); ylabel('Boundary peaks (%)'); legend('A','B','C','D'); grid on;
    save_figure(c,f,sprintf('%02d_study_%s_architectures',15+ip,params{ip}));
end
save(fullfile(root,'study_results.mat'),'c','study');
fprintf('\nParametric studies saved in %s\n',root);
end

function strong=run_strong_case(c,g)
% تکرار شکل‌های پایه در پوشه مستقل؛ شدت فرض آزمایش است، نه باد با سرعت معلوم
v=c; v.f=c.f_max; v.target_el=c.strong_target_el;
v.beam_el=v.target_el; v.calibration_el=v.target_el;
v.stack_elevations=v.target_el+[-8 -4 0 4 8]; v.selected_beam=3;
v.sigma_tip=c.strong_sigma_tip; v.mean_tip=c.strong_mean_tip;
v.plot_el=min(v.stack_elevations)-8:0.5:max(v.stack_elevations)+8;
v.output_dir=fullfile(c.output_dir,'strong_elevation_case');
if ~exist(v.output_dir,'dir'), mkdir(v.output_dir); end
state=rng; cleanup=onCleanup(@() rng(state)); %#ok<NASGU>
fprintf('\nStrong case: target elevation %.2f deg, tip mean %.2f mm, tip SD %.2f mm\n',v.target_el,v.mean_tip*1e3,v.sigma_tip*1e3);
validate_generator(v,g); rng(c.seed,'twister');
strong=run_monte_carlo(v,g); plot_results(v,g,strong);
save(fullfile(v.output_dir,'strong_results.mat'),'v','g','strong');
writetable(strong.summary,fullfile(v.output_dir,'scenario_summary.csv'));
disp(strong.summary);
end
