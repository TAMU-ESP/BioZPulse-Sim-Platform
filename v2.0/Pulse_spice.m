clc
clear
close all
fclose('all');
set(0,'DefaultFigureWindowStyle','docked')

% General settings
net.param.write_netlist='pulse_test';  % Name of results folder
net.param.write_figures=0;                  % Save Figures into results folder
net.param.write_mat=1;                      % Save Matlab Worksapce into .mat file results folder
net.param.plot_eit=0;                       
net.param.debug=1;                          % Run/Skip SPICE Simulation 0=>Run Simulation, 1=> Skip simulation and plot only conductivity images for debugging purposes

% Electrode Configuration

ac_freq=[1 10e3 10e3];  % Frequency Analysis [number of freq. pts, fstart, fend]
z_en=0;                 % Enable RC model (0=> R model, 1=> RC model)
se_en=1;                % Enable Skin-electrode impedance RC model (0=> S.C., 1=> RC model)
% Tissue paramters
%Tissue  Conduct(S/m)Permittivity
%Skin    2.04E-4     1.13E+3
%Bone    2.0E-2      5.0E+2
%Fat     4.3E-2      9.12E+2
%Muscles 3.4E-1      2.59E+4
%Blood   7.0E-1      5.25E+3


% skin,fat,muscle,blood,blood
tissue_arr={'Skin','Fat','Muscle','Blood','Blood'};

center_z_abs_arr=[0.5,0.6,0.7,0.7,0.7];
center_x_abs_arr=[0.5,0.5125,0.5125,0.275,0.775];
center_y_abs_arr=[0.5,0.5,0.5,0.5,0.5];
len_z_abs_arr=[1,0.8,0.6,0.05,0.05];
len_x_abs_arr=[1,0.875,0.775,0.025,0.025];
len_y_abs_arr=[1,1,1,1,1];


center_z_arr=[0.5,0.6,0.7,0.7,0.7];
center_x_arr=[0.5,0.5125,0.5125,0.275,0.775];
center_y_arr=[0.5,0.5,0.5,0.5,0.5];
len_z_arr=[1,0.8,0.6,0.05,0.05];
len_x_arr=[1,0.875,0.775,0.025,0.025];
len_y_arr=[1,1,1,1,1];
art_en=[0 0 0 0 0];
r_ratio=1; 
r_art_ratio=1e-2;

net.param.date=datestr(now,'yyyy-mm-dd'); 
net.param.path=strcat('./results/',net.param.date,'_',net.param.write_netlist,'/');
net.param.ac_freq=ac_freq;
net.param.z_en=z_en;
net.param.se_en=se_en;

for i=1:length(center_x_arr)
    net.param.obj(i).tissue=tissue_arr(i);
    net.param.obj(i).r_ratio=r_ratio;
    net.param.obj(i).center.x=center_x_arr(i);
    net.param.obj(i).center.y=center_y_arr(i);
    net.param.obj(i).center.z=center_z_arr(i);
    net.param.obj(i).len.x=len_x_arr(i);
    net.param.obj(i).len.y=len_y_arr(i);
    net.param.obj(i).len.z=len_z_arr(i);
    net.param.obj(i).art_en=art_en(i);  
    net.param.obj(i).r_art_ratio=r_art_ratio;
end

% Reference
%Create electrodes
[net]=pulse_elec_array_fn(net);
%Create Grid
net=img_gen_fn(net);
net=grid_gen_fn(net);
% Simulate SPICE Netlist
net_ref{1}=run_sim_fn(net);

% Plot Impedance/Conductivity Image
plot_single_img_fn(net_ref{1}(1),2); % 1=> electrodes, 2=> cond+elec

% Plot SPICE Simulation Results V/I
if net.param.debug==0
    plot_spice_arr_fn(net_ref,[1],[1],round(net_ref{1}(1).Ni/2)+1,2);   
end


if net.param.write_mat==1 && net.param.debug==0
    if ~exist(net.param.path,'dir')
        mkdir(net.param.path);
    end
    save(strcat(net.param.path,'netlist.mat'));
end

