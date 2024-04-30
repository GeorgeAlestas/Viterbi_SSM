#library imports
import numpy as np
import pycbc.psd
import pycbc.noise
from pycbc import frame
from pycbc.detector import Detector

#function to compute the optimum snr given a timeseries and a psd
def optimum_snr(ts, psd, flow_opt_snr, fhigh_opt_snr):

	#convert timeseries to frequencyseries
	fs = ts.to_frequencyseries()

	#find frequencies to be integrated
	idxs = (fs.sample_frequencies>flow_opt_snr) & (fs.sample_frequencies<fhigh_opt_snr) 
	freqs = fs.sample_frequencies[idxs] 
	htilde = np.array(fs)[idxs]

	#interpolate psd to the same frequencies as fs
	interp_psd = np.interp(freqs, psd.sample_frequencies, np.array(psd))

	#return optimum snr
	return np.sqrt(4*np.trapz(np.abs(htilde)**2/interp_psd, dx=freqs[1]-freqs[0]))

#also import myTaylorT3
import sys
sys.path.insert(1, './myTaylorT3')
from myTaylorT3 import myTaylorT3

import time
start_runtime = time.time()

distance_values = list(float(i)/1000 for i in range(5,151,2))

for distance in distance_values:
	#inputs about the frames
	chunck_duration = 4096       #Duration of chunks
	n_chunck = 5                #Number of chunks
	sampling_rate = 1024         #sampling rate
	flow_noise = 10              #lower noise to include 
	t_start = 1000000000         #start GPS time
	frame_folder = f'./Frames{distance}/'   #folder to save the frames
	frame_name = '_SSSM_m1_%.2g_m2_%.2g_dL_%.3g_tc_%.f_%.f-%.f.gwf' #name of the frames
	channel_name = 'STRAIN_SSSM'  #channel name

	#input on the injection
	m1 = 1.2e-2              #primary mass
	m2 = 1.2e-2              #secondary mass
	#distance = 0.05         #distance in Mpc
	ra = 1.7               #right ascension
	dec = 1.7              #declination
	pol = 0.2              #polarization
	inc = 0                #inclination
	ifos = ['H1', 'L1']    #interferometers in which to inject
	#coal_time = t_start + chunck_duration*n_chunck + 100 #coalescence time
	coal_time = t_start + 23915.024 #coalescence time
	ifo_psd_funcs = {'H1': pycbc.psd.aLIGOAdVO3LowT1800545, 'L1': pycbc.psd.aLIGOAdVO3LowT1800545}

	#information to compute the optimum SNR
	flow_snr = 20
	fhigh_snr = sampling_rate/2

	#initialize PSDs
	ifo_psds = {}
	delta_f = 1/chunck_duration
	flen = int(sampling_rate*chunck_duration)+1
	for ifo in ifos:
		ifo_psds[ifo] = ifo_psd_funcs[ifo](flen, delta_f, flow_noise)

	#initialize myTaylorT3 class
	wf_generator = myTaylorT3(m1=m1, m2=m2, distance=distance, inclination=0, sampling_rate=sampling_rate, coal_time=coal_time)

	#if frame folder does not exist, create it
	import os
	if not os.path.exists(frame_folder): os.makedirs(frame_folder)

	#loop over chunks
	snr_opt2_ifo = {}
	for i_chunck in range(n_chunck):
		#initial and final times
		t0 = t_start + i_chunck*chunck_duration
		tf = t0 + chunck_duration
		#compute the polarizations
		hp, hc  = wf_generator.tdstrain(t0, tf, PyCBC_TimeSeries=True)
		#loop over interferometers
		for ifo in ifos:
			#project GWs into the detector
			h_ifo = Detector(ifo).project_wave(hp, hc, ra, dec, pol, method="lal")
			
			#generate the time domain noise
			noise_ifo = pycbc.noise.noise_from_psd(sampling_rate*chunck_duration, 1/sampling_rate, ifo_psds[ifo])
			noise_ifo.start_time = noise_ifo.start_time + t0
		
			#inject the GW template in the noise
			strain_ifo = noise_ifo.inject(h_ifo)
			
			#save the frame files
			output_frame_file = frame_folder+ifo+frame_name%(m1, m2, distance, coal_time, t0, chunck_duration)
			channel_ifo = ifo+':'+channel_name
			frame.write_frame(output_frame_file, channel_ifo, strain_ifo)
			
			#compute the SNR in this chunck
			snr_opt2_chunck = optimum_snr(h_ifo, ifo_psds[ifo], flow_snr, fhigh_snr)**2
			#add it to the total
			if i_chunck == 0:
				snr_opt2_ifo[ifo] = snr_opt2_chunck
			else:
				snr_opt2_ifo[ifo] = snr_opt2_ifo[ifo] + snr_opt2_chunck
			
			#print information
			#print('Injected waveform in %s between t0=%.fs and tf=%.fs, snr in this chunck=%.2f, total acumulated snr in this ifo = %.2f'%(ifo, t0, tf, np.sqrt(snr_opt2_chunck), np.sqrt(snr_opt2_ifo[ifo])))

	#log the total snr
	snr_opt2_tot = 0
	output_string = 'Optimum SNR'
	for ifo in ifos:
		#add the snr
		snr_opt2_tot = snr_opt2_tot + snr_opt2_ifo[ifo]
		#add ifo to string
		output_string = output_string + ' in %s: %.2f,'%(ifo, np.sqrt(snr_opt2_ifo[ifo]))

	#print('\n',output_string,' Total = %.2f'%(np.sqrt(snr_opt2_tot)))	

		
	#Runtime
	#print("\nRuntime: %s seconds" % (time.time() - start_runtime))

