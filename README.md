# BioZPulse-Sim-Platform

[BioZPulse-Sim-Platform](./figures/grid.jpg)

BioZPulse Simulation Platform for Arterial Pulse Wave Modeling by Bassem Ibrahim, Drew A. Hall and Roozbeh Jafari.
This repository includes the source code of the bio-impedance simulation platform used in the following paper:

Bassem Ibrahim, Drew A. Hall, Roozbeh Jafari, Bio-impedance Simulation Platform using 3D Time-Varying Impedance Grid for Arterial Pulse Wave Modeling, IEEE Biomedical Circuits and Systems Conference (BioCAS), October 17-19, 2019, Nara, Japan. [(Paper)][1].

BioZPulse simulation platform source code files:

MATLAB Core Functions:
- spice_netlist_3d_gen_fn.m  		Generate SPICE netlist file (netlsit.cir)
- spice_netlist_3d_run_fn.m 		Run LTSPICE Simulator
- spice_netlist_3d_pp_fn.m  		Post-processing for voltage and PTT calculations
- spice_netlist_3d_plot_fn.m 		Plot Output Results
- spice_netlist_3d_plot_fn.m 		Plot Output Results in Frequency domain
- xyz2node_fn.m                		Convert from X,Y,Z coordinates to node index
- find_delay_fn.m					Calculate time delay of sensed pulse signal (Sinewave)	
- sinefit_fn.m                		Sinewave paramters extraction by fitting to a regression model 
- read_spice_out_fast_gnd_fn.m  	Read and sparse SPICE output file
- importfile_ltspice_raw_fast_fn.m	Import LTSPICE output file to MATLAB
- importfile_ltspice_header.m       Import the header of the LTSPICE output file to MATLAB
- imp_image_gen_fn.m                Adjust bio-impoedance values to MATLAB


MATLAB Scripts:
- ltspice_3d_model_dc_spacing.m     Simulate DC Bio-Z(V_DC) for different electrode spacing 
- ltspice_3d_model_dc_Y.m           Simulate DC Bio-Z(V_DC) for different vertical distance of sensing electrodes 
- ltspice_3d_model_deltaR1.m		Simulate arterial pulse amplitude(deltaV_PP) when sensing location is aligned with the artery
- ltspice_3d_model_deltaR2.m 		Simulate arterial pulse amplitude(deltaV_PP) when sensing location is 1.5cm away from the artery
- ltspice_3d_model_freq.m			Simulate DC Bio-Z(V_DC) and arterial pulse amplitude(deltaV_PP) versus frequency
- ltspice_3d_model_pene.m			Simulate arterial pulse amplitude(deltaV_PP) versus different artery's depth (Z_A)
- ltspice_3d_model_ptt.m			Simulate PTT and pulse amplitude(deltaV_PP) for different sensing locations relative to the artery

Getting Started:
- Install MATLAB. This is commercial software available from The MathWorks. For system requirements and installation instructions, please refer to their documentation.
- Install LTSPICE Simulator. It is a free SPICE simulator from Analog Devices. Download link: https://www.analog.com/en/design-center/design-tools-and-calculators/ltspice-simulator.html
- Update the variable "net.param.ltspice_path" in MATLAB function "spice_netlist_3d_run_fn.m" with your LTSPICE.exe path
- Run the required MATLAB Script


Notes:
- These MATLAB files were tested using MATLAB 2018b version

## Contact us
Bassem Ibrahim: bassem@tamu.edu

Drew A. Hall: drewhall@ucsd.edu

Roozbeh Jafari: rjafari@tamu.edu

[1]: https://www.dropbox.com/s/v1uuoazmgx41jrz/BioCAS2019_final.pdf?dl=0
