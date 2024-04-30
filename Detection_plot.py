import soapcw as soap
import numpy as np
import matplotlib.pyplot as plt
import pycbc
import lalpulsar
import lal
import time

class SFT:
    
    
    def __init__(self,nsft=None,tsft=None,delta_f=None,fmin=None,fmax=None,det_name=None,tstart=None):
        """
        Class to contain SFTs
        """
        self.sft = None
        self.norm_sft_power = None
        self.tstart = tstart
        self.nsft = nsft
        self.tsft = tsft
        self.epochs = None
        self.delta_f = delta_f
        self.det_name = det_name
        self.fmin = fmin
        self.fmax = fmax
        
        if self.delta_f is None and self.tsft is not None:
            self.delta_f = 1./self.tsft

class LoadSFT:
  def __init__(self,sftpath,fmin=None,fmax=None,norm=False,summed=False,filled=False,remove_sft=True,save_rngmed=False,tmin=None,tmax=None,vetolist = None):
        """
        Load an SFT from multiple detectors
        args
        -----------------
        sftpath: str or np.ndarray
            path to the sft files, for multiple files separate by semicolon, 'filename1;filename2;.....' can input from multiple detectors
            or input an np.ndarray of shape (Ndetectors, Ntimesteps, Nfreqbins)
        fmin: float
            minimum frequency to load from sft
        fmax: float
            maximum frequency to load from sft
        norm: bool or int
            normalise sft to running median, if integer running median of that many bins
        summed: bool or int
            sum normalised spectrograms over number of bins (default 48 (1day of 1800s)) or set to value
        filled: bool
            fill the gaps in the sfts with the expected value (2 for normalised sfts, nan for sfts)
        remove_sft: bool
            remove original sft after normalising to running median
        save_rndmed: bool
            save the running median as an estimate of noise floor
        tmin: float
            start time for sfts
        tmax: float
            end time for sfts
        vetolist: list
            list of frequency bins to set to expected value (2 or nan) (and +/- 3 bins of selected bin)
        det_names: list
            list of detector names (only necessary when using np.ndarray as input to sftpath)
        """


        if type(sftpath) is str:
            # load the sft
            self.get_sft(sftpath,fmin=fmin,fmax=fmax,tmin=tmin,tmax=tmax)
        elif type(sftpath) is dict:
            self.get_sft_from_array(sftpath, fmin=fmin, fmax=fmax, tmin=tmin, tmax=tmax)
        else:
            raise Exception(f"Please use a type str or np.ndarray for the input of sftpath, you used type {type(sftpath)}")

  
  def get_sft(self,sftpath,fmin=None,fmax=None,tmin=None,tmax=None,vetolist = None):
        '''
        load an sft to a numpy array, 
        args
        -------
        sftpath : string
            path to the sft file
        detector: string
            which detector i.e 'H1'
        fmin    : float
            minimum frequency
        fmax    : float
            maximum frequency
        tmin    : float
            min time
        tmax    : float
            max time
        '''
        
        # set up contraints for data
        constraints = lalpulsar.SFTConstraints()

        # set fmin and fmax to all sfts as default
        if fmin is None:
            fmin = -1
        if fmax is None:
            fmax = -1

        # only select sfts with start time (epoch) in range tmin-tmax
        if tmin is not None:
            self.tmin_gps = lal.LIGOTimeGPS(int(tmin),0)
            constraints.minStartTime=self.tmin_gps
        if tmax is not None:
            self.tmax_gps = lal.LIGOTimeGPS(int(tmax),0)
            constraints.maxStartTime=self.tmax_gps

        # create catalogue of sfts under specified constraints
        catalogue = lalpulsar.SFTdataFind(sftpath,constraints)

        # load the sfts from the catalogue above (currently uses multi SFTs as there is a memory leak in LoadSFTs)
        # sfts = lalpulsar.LoadSFTs(catalogue,fmin,fmax)
        sfts = lalpulsar.LoadMultiSFTs(catalogue,fmin,fmax)
        # keep sft set
        #self.sfts_old = sfts

        # define length of sfts in sft index and number of frequency bins
        # N = sfts.length
        self.det_names = []
        tsft = []
        nbins = []
        for det in sfts.data:
            tsft.append(1.0/det.data[0].deltaF)
            nbins.append(det.data[0].data.length)
        
        if len(set(tsft)) > 1:
            print("Warning: tsft not the same between detectors.")
            
        if len(set(nbins)) > 1:
            print("Warning: different detectors do not have the same number of frequency bins")
        
        for det in sfts.data:
            # get detectors name
            detname = det.data[0].name
            self.det_names.append(detname)

            # initialise SFT for detector
            data = SFT()

            # set parameters of sft
            data.nsft = det.length
            data.nbins = det.data[0].data.length
            data.delta_f = det.data[0].deltaF
            data.f0 = det.data[0].f0
            data.tsft = 1.0/det.data[0].deltaF

            #define range of frequency bin centers
            data.frequencies = np.arange(data.nbins)*data.delta_f + data.f0
            
            # create empty arrays for likleyhoods and epochs
            data.sft = np.zeros((data.nsft,data.nbins)).astype(np.complex_)
            data.epochs = np.zeros(data.nsft)
            # set fmin and fmax
            data.fmin = det.data[0].f0
            data.fmax = data.fmin + data.nbins/data.tsft

            #for i,sft in enumerate(sfts.data):
            for i,sft in enumerate(det.data):
                # fill sft with data
                data.sft[i,:] = sft.data.data
                # record epoch of each sft
                data.epochs[i] = sft.epoch
                
            # save sft for detector
            setattr(self,detname,data)

def find_rms_n(pul_track,vit_track,ref_CSh):
    diff = []
    for elem in range(len(pul_track)):
        weight = ref_CSh[elem]/np.sum(ref_CSh[elem])
        pathsqs = weight*((np.array(pul_track[elem])-np.array(vit_track[elem]))**2)
        #diff.append(np.sum(np.median(np.array(pathsqs))))
        diff.append(1./len(pathsqs)*np.sum(np.array(pathsqs)))

    return np.sqrt(1./(len(diff))*np.sum(diff))


sftdir0H1 = "./Frames0Noise/sfts_8.5sec/H-1_H1_8.5SFT_MSFT-1000"
gpstime_start = 1000000000
fband = 512
nsft = 2560
tsft = 8.5
data0H1 = np.zeros((int(tsft) * 512, nsft), dtype=np.complex128)
for k in range(nsft):
    sftfile0H1 = f"{sftdir0H1}/H-1_H1_8.5SFT_MSFT-{gpstime_start + k * int(tsft):d}-8.5.sft"
    data0H1[:,k] = LoadSFT(sftfile0H1).H1.sft[0]

rms_dl = []
distance_values = list(float(i)/1000 for i in range(5,151,2))

for distance in distance_values:
	sftdirH1 = f"./Frames{distance}/sfts_8_5sec/H-1_H1_8.5SFT_MSFT-1000"
	gpstime_start = 1000000000
	fband = 512
	nsft = 2560
	tsft = 8.5
	dataH1 = np.zeros((int(tsft) * 512, nsft), dtype=np.complex128)
	for k in range(nsft):
	    sftfileH1 = f"{sftdirH1}/H-1_H1_8.5SFT_MSFT-{gpstime_start + k * int(tsft):d}-8.5.sft"
	    dataH1[:,k] = LoadSFT(sftfileH1).H1.sft[0]
	    

	tr_1 = soap.tools.tr_p(1.0)

	PSDH1 = np.median(np.abs(dataH1)**2, axis=1)/(2*np.log(2))

	CShusterH1 = np.transpose((np.abs(dataH1)**2)/PSDH1[:, np.newaxis])

	fmin = 61
	delta_f = 1/int(tsft)
	N0 = int(fmin/delta_f)
	fmax = 130
	Nmax = int(fmax/delta_f)
	CShusterH1=CShusterH1[:, N0:Nmax]

	one_tracks_ng = soap.single_detector(tr_1, CShusterH1)

	PSD0H1 = np.median(np.abs(dataH1)**2, axis=1)/(2*np.log(2)) #np.median(np.abs(data0H1)**2, axis=1)/(2*np.log(2))

	CShuster0H1 = np.transpose((np.abs(data0H1)**2)/PSD0H1[:, np.newaxis])

	fmin = 61
	delta_f = 1/int(tsft)
	N0 = int(fmin/delta_f)
	fmax = 130
	Nmax = int(fmax/delta_f)
	CShuster0H1=CShuster0H1[:, N0:Nmax]

	one_tracks_ng0 = soap.single_detector(tr_1, CShuster0H1)

	CShuster0H1_track=np.transpose(np.take_along_axis(CShuster0H1, one_tracks_ng0.vit_track[:,np.newaxis], axis=1))

	rms = find_rms_n([one_tracks_ng0.vit_track],[one_tracks_ng.vit_track],CShuster0H1_track)
	rms_dl.append((distance, rms))

np.savetxt('rms_dl_new.txt', rms_dl)

distances = [result[0] for result in rms_dl]
rms_values = [result[1] for result in rms_dl]

# Plot the RMS values vs. distances
plt.plot(distances, rms_values, marker='o', linestyle='None', markersize=3)
plt.xlabel("d$_{L}$")
plt.ylabel("RMS")
plt.savefig("rms_vs_dL_new.pdf")
