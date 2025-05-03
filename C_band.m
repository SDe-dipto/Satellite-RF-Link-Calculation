clear;clc;close
lat_sat = 0;
sat_long=input('Longitude of GEO satellite in question = ','s');
t_sys_sat = 75;
sat_ant_size = 2.2;
sat_ant_gain = 28;
eirp_sat = 40;
IF_bw=27;
earth_input=input('Location of earth station = ','s');
part_earth = strsplit(earth_input,';');
lat_str=strtrim(part_earth{1});
lon_str=strtrim(part_earth{2});
lat_e=str2double(lat_str(1:end-1));
lon_e=convert_longitude(lon_str);
uplink=6.75e9;
downlink=4.2e9;
lamb_up=3e8/uplink;
l_atm=0.2;
hub_size=7;
ap_eff=0.55;
hpa_tx=400;
hpbw_tx=70*lamb_up/hub_size;
hpa_obo=0.5;
misc_losses=0.3;
lon_s=str2double(sat_long(1:end-1));
ht=36000e3;
ro=6378e3;
rs=42164e3;
central_angle=acosd(cosd(lat_e)*cosd(lon_s-lon_e));
ro_rs=ro/rs;
rs_ro=rs/ro;
el_angle=atand((rs/ro-cosd(central_angle))/sind(central_angle))-central_angle;
fprintf('The elevation angle of earth station to satellite = %0.2f\n',el_angle)
d=rs*sqrt(1.02288235-0.30253825*cosd(central_angle));
R=d/1000;
fprintf('The range between the station and the satellite = %0.2f km\n',R);
Gt_dB = antenna_gain(hub_size,uplink,ap_eff);
hpa_dBW=10*log10(hpa_tx);
hpa_dBm=hpa_dBW + 30;
eirp_tx=Gt_dB+hpa_dBW-hpa_obo-misc_losses;
fprintf('Effective Isotropic Radiative power of Tx antenna = %0.1f dBW\n',eirp_tx);
x=linspace(0,6e4);
y=fsl(uplink,x);
plot(x,y,'-g',x,fsl(downlink,x),'-black');
grid on;hold on;
L_fsl=fsl(uplink,R);
scatter(R,L_fsl,'r*');
xlabel('Distance (in km)');ylabel('Free Space Path Loss (dB)');title('Free Space Loss');

Pow_rcv = eirp_tx+sat_ant_gain-L_fsl-l_atm;
fprintf('The power required at receiver antenna at satellite end = 50.2f dB \n',Pow_rcv)
k=1.38e-23;
up_ns_pow=10*log10(k*1e6*t_sys_sat);
cn_up=Pow_rcv-up_ns_pow;
fprintf('The required C/N in clear air = %0.3f dB-Hz\n',cn_up);
vsat_ant_size=1.2;
vsat_fsl=fsl(downlink,R);
vsat_ap_eff=0.63;
Grx_vsat=antenna_gain(vsat_ant_size,downlink,vsat_ap_eff);
Lrx_vsat=1;
T_sys_vsat=500;
scatter(R,vsat_fsl,'b*');
legend('Uplink loss','Downlink loss','Uplink operating point','Downlink operating point','Location','southeast');

lamb_down=physconst('LightSpeed')/downlink;
transponder_tx_size=(lamb_down/pi)*sqrt((10^(sat_ant_gain/10))/vsat_ap_eff);
fprintf('Downlink - Size of Satellite transponder = %0.2f m\n',transponder_tx_size);
Pr_down = eirp_sat + Grx_vsat+sat_ant_gain-vsat_fsl-l_atm;
fprintf('Received Power at VSAT end = %0.2f dBW\n',Pr_down)
down_ns_pow=10*log10(k*T_sys_vsat);
cn_down=Pr_down-down_ns_pow;
fprintf('Downlink C/No in clear air = %0.2f dB-Hz\n',cn_down);
C_No_total=-10*log10(10^(-cn_up/10)+10^(-cn_down/10));
fprintf('Overall C/No for C-band GSO satellite = %0.3f dB-Hz\n',C_No_total)

function east_long = convert_longitude(longitude)
lon_deg=str2double(longitude(1:end-1));
if endsWith(longitude,'W')
    east_long = 360-lon_deg;
else
    east_long=lon_deg;
end
end

function free_space=fsl(freq,dist)
wvl=3e8/freq;
free_space=abs(20*log10(3e8/(4*pi))-20*log10(1e9)-20*log10(freq/1e6)-20*log10(dist));
end

function Gain=antenna_gain(ant_size,freq,eff)
wvl=3e8/freq;
gn=eff*power(pi*ant_size/wvl,2);
Gain=10*log10(gn);
end