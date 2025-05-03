clc;close;clear;
lat_sat=0;
long_sat=83;
eirp_sat=51.5;
G_T_sat=3;
sat_ant_size=2.2;
sat_backoff=3;
twta=140;
T_sys=500;
twta_dB=10*log10(twta);
sat_ant_gain=eirp_sat-twta_dB;
downlink=10.95;
uplink=14;
wvl_up=physconst('LightSpeed')/(uplink*1e9);
wvl_down=physconst('LightSpeed')/(downlink*1e9);
band_noise=43.2;
earth_input=input('Location of earth station = ','s');
part_earth = strsplit(earth_input,';');
lat_str=strtrim(part_earth{1});
lon_str=strtrim(part_earth{2});
lat_e=str2double(lat_str(1:end-1));
lon_e=convert_longitude(lon_str);
ht=36000e3;
ro=6378e3;
rs=42164e3;
central_angle=acosd(cosd(lat_e)*cosd(long_sat-lon_e));
ro_rs=ro/rs;
rs_ro=rs/ro;
el_angle=atand((rs/ro-cosd(central_angle))/sind(central_angle))-central_angle;
fprintf('The elevation angle of earth station to satellite = %0.2f\n',el_angle)
d=rs*sqrt(1.02288235-0.30253825*cosd(central_angle));
Range=d/1000;
fprintf('The range between the station and the satellite = %0.2f km\n',Range);
ant_size=5;
ap_eff=0.68;
loss_pol=0.5;
loss_point=0.5;
loss_atm=0.3;
trans_hpa_ob=1;
pow_otpt=10*log10(600);
R=3.7e7;
FSL_up=free_space(R,wvl_up);
Gt_up=10*log10(ap_eff*(pi*ant_size/wvl_up)^2);
Pow_rcv=pow_otpt+Gt_up+sat_ant_gain-FSL_up-loss_atm-loss_point-loss_pol;
noise_up=10*log10(physconst('boltzmann')*T_sys*band_noise*1e6);
cn_up=Pow_rcv-noise_up;
fprintf('Isotropic power at satellite receive end = %0.2f dB\n',Pow_rcv)
fprintf('C/N at uplink end in clear air = %0.2f dB\n',cn_up)
ant_noise=30;
lna_noise=110;
loss_loc=3;
dn_ant_size=0.67;
G_dn=10*log10(ap_eff*(pi*dn_ant_size/wvl_down)^2);
FSL_down=FSL_up*(uplink/downlink);
noise_down=10*log10(physconst('Boltzmann')*(ant_noise+lna_noise)*band_noise*1e6);
P_dn_rcv=eirp_sat+G_dn-FSL_down-loss_loc-loss_pol-loss_point-loss_atm;
fprintf('User terminal isotropic power requirement = %0.2f dB\n',P_dn_rcv);
cn_down=P_dn_rcv-noise_down;
fprintf('C/N for user terminal at downlink end in clear air = %0.2f dB\n',cn_down)
c_no=-10*log10(10^(-cn_up/10)+10^(-cn_down/10));
fprintf('Overall C/No in clear air = %0.2f dB\n',c_no);
R_001=42;
hs=0.9;
theta=el_angle;
phi=lat_e;
a=0.0367;b=1.154;
loss_rain_up1=rain_attn(R_001,hs,theta,phi,uplink,0.0367,1.154);
fprintf('Rain attenuation in uplink = %0.2f db\n',loss_rain_up1);
trans_otpt_pow=twta_dB-sat_backoff-loss_rain_up1;
Pow_rcv_rain=Pow_rcv-loss_rain_up1;
fprintf('Isotropic power at satellite receive end in rain = %0.2f dB\n',Pow_rcv_rain);
cn_up_rain=cn_up-loss_rain_up1;
fprintf('C/n at satellite receive end in rain = %0.2f dB\n',cn_up_rain);
loss_rain_down=rain_attn(R_001,hs,theta,phi,downlink,0.0101,1.276);
fprintf('Rain attenuation in downlink = %0.2f dB\n',loss_rain_down);
P_dn_rcv_rain=P_dn_rcv-loss_rain_down;
fprintf('Isotropic power at user terminal end in rain = %0.2f dB\n',P_dn_rcv_rain);
t_sky_rain=t_sky(loss_rain_down);
Ts_rain=lna_noise+t_sky_rain;
delta_noise=10*log10(Ts_rain/(ant_noise+lna_noise));
cn_down_rain=cn_down-loss_rain_down-delta_noise;
fprintf('C/N at user terminal end in rain = %0.2f dN\n',cn_down_rain);
c_no_rain = -10*log10(10^(-cn_up_rain/10)+10^(-cn_down_rain/10));
fprintf('Overall C/No i nrain = %0.2f dB\n',c_no_rain);
f=[1,2,4,6,7,8,10,12,15,20,25,30,35,40,45,50,60,70,80,90,100,120,150,200,300,400];
a=[0.0000387,0.000154,0.000650,0.000175,0.00301,0.00454,0.0101,0.0188,0.0367,0.0751,0.124,0.187,0.263,0.350,0.443,0.536,...
    0.707,0.851,0.975,1.06,1.12,1.18,1.31,1.45,1.36,1.32];
b=[0.912,0.963,1.121,1.308,1.332,1.327,1.276,1.217,1.154,1.099,1.061,1.021,0.979,0.939,0.903,0.873,0.826,0.793,0.769,...
    0.753,0.743,0.731,0.710,0.689,0.688,0.683];
pl_rain=rain_attn(R_001,hs,theta,phi,f,a,b);
loglog(f,pl_rain);
grid on;hold on;
title('Rain attenuation as function of frequency');
xlabel('Frequency(GHz)');ylabel('Rain attenuation (dB)');
scatter(uplink,loss_rain_up1,'g*');scatter(downlink,loss_rain_down,'r*');
legend('Rain attenuation','Uplink attenuation','Downlink attenuation','Location','southeast');

function fsl = free_space(R,wvl)
loss=(4*pi*R)/wvl;
fsl=20*log10(loss);
end

function east_long = convert_longitude(longitude)
lon_deg=str2double(longitude(1:end-1));
if endsWith(longitude,'W')
    east_long = 360-lon_deg;
else
    east_long=lon_deg;
end
end

function rain_db = rain_attn(R_001,hs,theta,phi,freq,a,b)
    hr=3+0.028*phi;
    Ls=(hr-hs)./sind(theta);
    Lg=Ls.*cosd(theta);
    Lo=35.*exp(-0.015.*R_001);
    r_001=1./(1+(Lg./Lo));
    gamma_r=a.*((R_001).^b);
    rh_001=1./(1+0.78.*sqrt(Lg.*gamma_r./freq)-0.38.*(1-exp(-2.*Lg)));
    tau=atand((hr-hs)./(Lg.*rh_001));
    if tau>theta
        Lr=(Lg.*rh_001)./cosd(theta);
    else
        Lr=(hr-hs)./sind(theta);
    end
    if abs(phi)<36
        sci=36-abs(phi);
    else
        sci=0;
    end
    v_001= 1./ (1 + sqrt(sind(theta)) .* (31 .* (1- exp(-1.*(theta./(1+sci)))) .* sqrt(Lr.*gamma_r)./(freq.^2) - 0.45));
    Le=Lr.*v_001;
    rain_db=gamma_r.*Le;
end

function tsk=t_sky(attn)
G=10^(-attn/10);
tsk=270*(1-G);
end