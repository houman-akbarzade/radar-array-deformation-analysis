function c=ArrayStudyConfig()
% Public, illustrative parameters. Distances in metres, angles in degrees.
c.c0=299792458;
c.f_min=5.2e9; c.f_max=5.8e9; c.f=5.8e9;
c.total_panels=60; c.panels=20; c.rows=16;
c.panel_pitch_deg=360/c.total_panels;
c.active_panel_ids=(21:40)';
c.panel_az_deg=((1:c.panels)-(c.panels+1)/2)*c.panel_pitch_deg;
c.dz=0.5*c.c0/c.f_max;
c.radius=c.dz/(2*sin(pi/c.total_panels));
c.seed=3719; c.nmc=500; c.nstats=10000;
c.run_frequency=true; c.frequency_list=linspace(c.f_min,c.f_max,5);
c.nmc_frequency=500;
c.study_frequencies=[5.2 5.5 5.8]*1e9;
lambda_ref=c.c0/c.f_max;
c.sigma_tip_list=lambda_ref*[0 .005 .01 .02 .04 .06 .08 .10];
c.mean_tip_list=lambda_ref*[-.06 -.03 -.01 0 .01 .03 .06];
c.target_el_list=[0 2 4 6 8]; c.nmc_study=500;
c.strong_sigma_tip=.08*lambda_ref; c.strong_mean_tip=.03*lambda_ref;
c.strong_target_el=6;
c.steer_az=0; c.analog_el=0;
c.target_az=1; c.target_el=4;
c.stack_elevations=[-4 0 4 8 12]; c.selected_beam=3;
c.calibration_el=4;
c.az_taper='hamming'; c.el_taper='taylor';
c.taylor_nbar=5; c.taylor_sll_db=30;
c.cal_offsets=-3:0.005:3; c.pattern_offsets=-25:0.05:25;
c.peak_half_span=4; c.peak_grid_step=0.1;
% Same correlated panel load and shared tangential/radial/vertical field.
c.sigma_tip=.02*lambda_ref;
c.rho_panel=.35; c.rho_within=.995;
c.radial_ratio=1e-4; c.vertical_ratio=.02;
c.mechanical_model='base'; c.bending_shape='uniform';
c.phase_center_power=1; c.wind_az_deg=0; c.mean_tip=0;
c.height_gradient=0; c.element_pattern='isotropic';
c.plot_az=-12:0.5:12; c.plot_el=-12:0.5:18; c.plot_floor_db=-35;
c.output_dir=fullfile(pwd,'output_matlab');
end
