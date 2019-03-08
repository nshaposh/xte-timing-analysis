import os
import numpy as np
import astropy.io.fits as pyfits
#from xtefilt import xtefilt
from utils import *
from copy import copy
                
def pds(
    info,
    modes=['auto'],
    outfile='pds.fps',
    rebin=1.01,
    fnyq=256.0,
    fft_int=256.0,
    expr="ELV > 5.0 and OFFSET < 0.1 and NUM_PCU_ON > 0",
    dtcor=True,
    chans=[0,255],
    dtmethod='vikhl',
    norm='rms',
    verbosity=1
):

    """
    xtepds - procedure to calculate the power density spectrum (PDS) for XTE 
    observation which is specified by the info dictionary parameter and is either. 

    Parameter[default]:

    modes['auto'] - XTE modes to use in PDS calculation. 'auto' means modes
                    are calculated using pyxte.xteobs.suggest_fft_modes() function


    outfile['pds.fits'] - name of output file containing power spectrum.
                          Empty string '' or 'none' means do not produce any output file.

    rebin[1.02]  -   log rebinning factor (similar to POWSPEC rebin parameter)

    fnyq [256.0] -   Nyquist frequency of the resulting PDS in Hertz.

    fft_int[256.0] - Interval for individual FFTs (sec).

    expr["ELV > 5.0 and OFFSET < 0.1 and NUM_PCU_ON > 0"]
                   - data screening expression. If empty, no screening will be done.

    chans[0,255] -   lowest and highest PCA channels (0-255) to include in analysis.
                     If data mode does not allow the exact channel range to be
                     extracted, only mode bins which are completely covered by
                     the specified PCA range will be used.
    
    gtifile[''] -    Good Time Intervals (GTI) file. 
                     Empty string '' or 'none' means do not use GTI file.

    norm['rms'] -    PDS normalization (!!!to be explained!!!)

    dtcor [True] -   
                     if True, then perform dead time correction
    dtmethod['zhang'] -
                     method for deadtime correction. Can also be "vikhl".
                     '' or 'none' means no dead time cerrection
    
    """

    if verbosity>0:
        print("=============================================================")
        print("XTEPDS (XTE Power spectrum extractor)")
        print("Author: Nikolai Shaposhnikov, shaposhnikov@gmail.com")
        print("=============================================================")

    bin_size = 0.5/fnyq
    fft_arr_size = int(fft_int/bin_size)
    xav = np.zeros(fft_arr_size//2-1)
    x2av = np.zeros(fft_arr_size//2-1)
    fft_arr = np.zeros(fft_arr_size,dtype=int)
    if modes.__class__.__name__ == 'str': modes = [modes]

    if verbosity>0:
        print("       Bin Size",bin_size," seconds")
        print("       Nyquist frequency:",fnyq," Hz")
        print("       FFTs will be calculated for %f sec. intervals and powers"%fft_int)
        print("       averaged. Resuting PDS will be rebinned with %f rebinning factor."%rebin)
        print("------------------------------------------------------------------------")

    freqt = np.fft.fftfreq(int(fft_int/bin_size),bin_size)
    freq= freqt[1:fft_arr_size//2]
#    freq = fnyq*array(range(1,fft_arr_size/2))/(fft_arr_size/2.0)
    df = freq[2]-freq[1]
    nph = 0
    ngamma = 0
    nrange = 0

# Making GTI from filter file

    filter_file = info['filter']
    gti_file = info['gti']

    #print filter_file

    if not os.path.exists(filter_file):
        #filter_file = xtefilt(xobs)
        print(f"Filter file {filter_file} is not found!")
        print("No deadtime correction or good time intervals filtering will be done.")
        dtcor = False
    else:
        filtgti = makegti(filter_file,expr)
        xfl = filter_file
        
    if expr.strip() =="":
        expr = "True"

#    filtgti = makegti(filter_file,expr)
#    print "GTI:",filtgti

#    if os.path.exists(filt): 
#        xfl = filt
#    elif os.path.exists(filter_file): 
#        xfl = filter_file
#    else:
#        if verbosity>0:

    
    if modes[0] == 'auto':
        
        if verbosity>0: print("Selecting PCA modes for PDS calculation.")
        fft_modes = suggest_fft_modes(info,verbosity=0)
        if fft_modes == -1:
            print("Failed to find proper modes.")
            return -1
        if verbosity > 0:
            print("Using following PCA modes:")
            for m in fft_modes:
                print("      %s"%m)
        print("------------------------------------------------------------------------")

    else:
#        print "modes",modes
        fft_modes = modes

    
    maxres = info['tres'][fft_modes[0]]
    for m in fft_modes:
        maxres = max(maxres,info['tres'][m])

    if maxres > bin_size:
        bin_size = maxres
        if verbosity>0:

            print("WARNING: Time resolution of the data modes %f sec. is lower than needed"%maxres)
            print("for the Nyquist frequency of %f Hz. Will calculate PDS with Nyquist "%fnyq)
            print("frequency of %f Hz."%(0.5/bin_size))
            print("------------------------------------------------------------------------")

            
        fnyq = 0.5/bin_size

    if verbosity>0:
        print("Maximum time resolution: %f"%maxres)

    if dtcor:
#    if True:
        filt = pyfits.open(xfl)
        fdata = filt[1].data.field('NUM_PCU_ON')
        ftime = filt[1].data.field('Time')
        pcu2hk = pyfits.open(info['fh5c'][0])
        pcu2hk_ind = int(0)
        dsvle = pcu2hk[1].data.field('dsvle')
        std1f = pyfits.open(info['pcadata']['Standard1'][0])
        std1_file = int(0)
        std1t = std1f[1].data.field('Time')
        std1f_len = len(std1t)
        xepcu0 = std1f[1].data.field('XeCntPcu0')
        xepcu1 = std1f[1].data.field('XeCntPcu1')
        xepcu2 = std1f[1].data.field('XeCntPcu2')
        xepcu3 = std1f[1].data.field('XeCntPcu3')
        xepcu4 = std1f[1].data.field('XeCntPcu4')
        vpcnt = std1f[1].data.field('VpCnt')
        remcnt = std1f[1].data.field('RemainingCnt')
        vlecnt = std1f[1].data.field('VLECnt')
        std1_arr_ind = int(0)
        filt_ind = int(0)
        std1f_ind = int(0)
        dsvle_value = 0.0
        rxenon = 0.0
        rvle = 0.0
        pcuon = 0.0
        npcuon = int(0)
        td = 1.0e-5
        n_dsvle = int(0)
        cor_time = 0.0
        if dtmethod == 'vikhl':
            def dtvik(td,tb,r,f):
                dt = 0.0
                r0 = r/(1-r*td)
                for i in range(21):
                    k = i - 10
                    ff = 2.0*np.pi*(f + k/tb)
                    x1 = r0*(1.0 - np.cos(ff*td))
                    x2 = r0*np.sin(ff*td)
                    dt += (x1+ff*x2)/(x1*x1+(x2+ff)**2)/(ff*tb)**2
                dt = -16.0*np.sin(np.pi*f*tb)**2*dt/r
                return dt

#    if os.path.exists(pth+"/"+params['gtifile']):



    try:

        gtif = pyfits.open(gtifile)
        gti = True
        gtistart = gtif[1].data.field('START')
        gtistop = gtif[1].data.field('STOP')

    except:

        gti = False


    fft_file_number = np.zeros(len(fft_modes),dtype=int)
    files = []
    fnames = []
    dsize = []

    for m in fft_modes:
        if verbosity>0: 
            print("Opening file:",os.path.basename(info['pcadata'][m][0])," mode", m)
        files.append(pyfits.open(info['pcadata'][m][0],memmap=True))
        fnames.append(info['pcadata'][m][0])

    gti_ind = int(0)

#    tstart = filtgti[0][gti_ind]
#    tstop = filtgti[1][gti_ind]

#    print tstart,tstop
    
    tstart = []
    tstop = []
    sci_mode = []
    numchan = []
    chns = []

    for f in files:
        datamode = f[1].header['DATAMODE']
#        tstart.append(f[0].header['TSTART'])
    #    tstop.append(f[0].header['TSTOP'])
        tstart.append(f[1].data.field(0)[0])
        tstop.append(f[1].data.field(0)[-1])
        sci_mode.append(f[1].header['EXTNAME'])
        nnn = info['nchan'][datamode]
        numchan.append(nnn)
        iii = []
        for ich in range(nnn):
            ch1 = info['pcachan_low'][datamode][ich]
            ch2 = info['pcachan_high'][datamode][ich]
#            print ch1,ch2,iii,chans[0],chans[1]

            if ch1>=chans[0] and ch2<=chans[1]: 
                
#                print ch1,ch2,iii,chans[0],chans[1]
                iii.append(ich)
        if iii == []:
            chns.append((-1,-1))            
        else:
            chns.append((min(iii),max(iii)))
#        print datamode,chns[0][0]
            

#    print chns
    time = filtgti[0][0]
    if dtcor: time = max(time,std1t[std1f_ind])
    dataind = np.zeros(len(files),dtype=int)

    stopfft = 0

#    for i in range(len(filtgti[0])):
#        print filtgti[0][i],filtgti[1][i]


    while not stopfft:

        skip = int(0)

        if time+fft_int > filtgti[1][gti_ind]:
            gti_ind += 1
            if gti_ind >= len(filtgti[1]):
                stopfft = True
            else:
                time = filtgti[0][gti_ind]
            skip = 4
            skip_message = "GTI changeover"

    #Making sure all datamodes are on during the interval

        for i in range(len(files)):
        
            if ( time <= tstart[i] ) or (time+fft_int >= tstop[i]+files[i][1].header['TIMEDEL']) :
                if verbosity>1: 
                    print("Mode",fft_modes[i],",file",os.path.basename(fnames[i])," does not cover interval.")
                skip = 1
                skip_message = "Data not available."

# makeing sure that STD data is available through the interval
# if not that skip the interval
        if dtcor:

  # Setting time pointer for Standard1 data
            if (time + fft_int >= std1t[-1]+128.0):
                skip = 2
                skip_message = "Std 1 data is not available"

            while std1t[std1f_ind]+0.125*float(std1_arr_ind) < time - 0.125 and skip == 0:
                std1_arr_ind += 1
                if std1_arr_ind == 1023:
                    std1f_ind += 1
                    std1_arr_ind = 0
                        #        else:
                        #            print "Std1 stop",std1_file, std1f_ind
                        #            std1f.close()
                        #           if std1_file < len(pcadata['Standard1'])-1:
                        #                print "Switch Std1 file."
                        #                std1_file += 1
                        #                        #                std1f = pyfits.open(pcadata['Standard1'][std1_file]+".gz")
                        #                std1f_ind = int(0)
                        #               std1_arr_ind = int(0)
                        #                std1t = std1f[1].data.field('Time')
                        #                std1f_len = len(std1t)
                        #                xepcu0 = std1f[1].data.field('XeCntPcu0')
                        #                xepcu1 = std1f[1].data.field('XeCntPcu1')
                        #                xepcu2 = std1f[1].data.field('XeCntPcu2')
                        #                xepcu3 = std1f[1].data.field('XeCntPcu3')
                        #                xepcu4 = std1f[1].data.field('XeCntPcu4')
                        #                vpcnt = std1f[1].data.field('VpCnt')
                        #                remcnt = std1f[1].data.field('RemainingCnt')
                        #                vlecnt = std1f[1].data.field('VLECnt')
    #                time = max(time,std1t[std1f_ind])
    #    else:
    #        print "POWSP Warning: No Std1 mode. Interval will be skipped."
#        skip = 1

#   settting HK pointer
    
            while pcu2hk[1].data.field('Time')[pcu2hk_ind] < time-8.0: pcu2hk_ind += 1


# Setting time pointer to filter file 

	
        while ftime[filt_ind]+0.125*float(std1_arr_ind) < time: filt_ind += 1

        ngamma = 0
        fft_arr = np.zeros(fft_arr_size,dtype=float)

#    GTI block, checking if the data are entirely within  GTIs

        if gti:
            ingti = int(0)
            for j in range(len(gtistart)):
                if time >= gtistart[j] and time+fft_int <= gtistop[j]:
                    ingti += 1
            if ingti == 0:
                skip = 3
                skip_message = "Rejected due to GTIs."
                continue
        else:
            ingti = 1
        
    
        if not skip:

            for i in range(len(files)):
            
                times = files[i][1].data.field(0)

                while times[dataind[i]] < time: dataind[i] += 1
                while times[dataind[i]] > time: dataind[i] -= 1

                tsize = len(times)
            
# Accumulating data Science Array formats (Binned and Single Bit)

                if sci_mode[i] == 'XTE_SA':
                    events = files[i][1].data.field(1)
                    ngr = int(files[i][1].header['TIMEDEL']/files[i][1].header['1CDLT2'])
                    b_ind = int(0)
                    st = max(1,int(bin_size/files[i][1].header['1CDLT2']))
                    barr = np.zeros(ngr,dtype=float)

                    if numchan[i]>1 and chns[i] != (-1,-1):
                        barr = np.sum(events[dataind[i]][chns[i][0]:chns[i][1]+1,:],axis=0)
                    elif chns[i] != (-1,-1):
                        barr = events[dataind[i]]

                    for k in range(fft_arr_size):
                        fft_arr[k] += np.sum(barr[b_ind:b_ind+st])
                        b_ind += st
                        if b_ind >= ngr:
                            b_ind = 0
                            dataind[i] += 1
                            barr = np.zeros(ngr,dtype=float)
                            if numchan[i]>1 and chns[i] != (-1,-1):
                                barr = sum(events[dataind[i]][chns[i][0]:chns[i][1]+1,:],axis=0)
                            elif chns[i] != (-1,-1):
                                barr = events[dataind[i]]

# Accumulating data in Science Event format (Event)

                if sci_mode[i] == 'XTE_SE':
                    events = files[i][1].data.field('PHA')
                    for k in range(fft_arr_size):
                        t1 = time + float(k+1)*bin_size
                        while times[dataind[i]] < t1:
                            if events[dataind[i]]>=chns[i][0] and events[dataind[i]]<=chns[i][1]: fft_arr[k] += 1
                            dataind[i] += 1

            ngamma = np.sum(fft_arr)
            nph += ngamma

#  Accumulating data for deatime correction  from HK and filter file

            if dtcor:
    #            cor_time = 0.0

                npcu = fdata[filt_ind]
                while std1t[std1f_ind]+0.125*float(std1_arr_ind) < time+fft_int:
                

                    while pcu2hk[1].data.field('Time')[pcu2hk_ind] < std1t[std1f_ind]+0.125*float(std1_arr_ind)-pcu2hk[1].header['TIMEDEL']:
                        dsvle_value += dsvle[pcu2hk_ind]
                        pcu2hk_ind += 1
                        n_dsvle += 1

                    while ftime[filt_ind] < std1t[std1f_ind]+0.125*float(std1_arr_ind)-filt[0].header['TIMEDEL']:
                        filt_ind += 1
                        npcu = float(fdata[filt_ind])
                        pcuon += npcu
                        npcuon += 1
#                	print pcuon, npcuon

                    rtmp = 0.0
                    if npcu > 0:
                        rtmp += xepcu0[std1f_ind][std1_arr_ind]
                        rtmp += xepcu1[std1f_ind][std1_arr_ind]
                        rtmp += xepcu2[std1f_ind][std1_arr_ind]
                        rtmp += xepcu3[std1f_ind][std1_arr_ind]
                        rtmp += xepcu4[std1f_ind][std1_arr_ind]
                        rtmp += vpcnt[std1f_ind][std1_arr_ind]
                        rtmp += remcnt[std1f_ind][std1_arr_ind]
                        rxenon += rtmp
                        rvle += vlecnt[std1f_ind][std1_arr_ind]

                    std1_arr_ind += 1
                    cor_time += 0.125
                    if std1_arr_ind == 1023:
#                	print std1f_ind,pcu2hk[1].data.field('Time')[pcu2hk_ind],pcu2hk_ind,n_dsvle
                        std1f_ind += 1
                        std1_arr_ind = 0

#        print "Stop HK."

            nrange += 1
            ftrans = np.fft.fft(fft_arr -float(ngamma) * np.ones(fft_arr_size,dtype=float)/float(fft_arr_size))
            power = 2.0*abs(ftrans[1:(fft_arr_size//2)])**2/float(ngamma)-2.0
            xav = xav + power
            x2av = x2av + power*power

#        End interval loop

        if verbosity > 0:
            if skip == 0:
                print(" -- START = %f, Rate = %f"%(time-fft_int,ngamma/fft_int))
            else:
                print(" -- START = %f, Skipped (code %i). %s"%(time-fft_int,skip,skip_message))


        time += fft_int

        for i in range(len(files)):
        
            if time >= tstop[i]:

                if verbosity>1: print("Closing file:",os.path.basename(fnames[i]),". Mode:", fft_modes[i],".")
                files[i].close()
                modefiles = info['pcadata'][fft_modes[i]]
                ind = modefiles.index(fnames[i])
#                print ind,len(modefiles)
                if ind == len(modefiles) - 1:
                    if verbosity>1: print("Mode",fft_modes[i]," exhausted. Stoping!")
                    stopfft = 1
                else:
                    fnames[i] = modefiles[ind+1]
#                    print modefiles
                    if verbosity>1: print("Opening file:",os.path.basename(fnames[i]),". Mode:", fft_modes[i],".")
                    files[i] = pyfits.open(fnames[i],memmap=True)
                    tstart[i] = files[i][1].data.field(0)[0]
                    tstop[i] = files[i][1].data.field(0)[-1]
                #                dataind[i] = 0
                    time = max(max(tstart),time)
                    dataind = np.zeros(len(files),dtype=int)

        if (time + fft_int >= std1t[-1]) & (std1_file < len(info['pcadata']['Standard1'])-1) & dtcor:
            std1f.close()
            std1_file += 1
            if verbosity>1: print("Switch Std1 file.")
            std1f = pyfits.open(info['pcadata']['Standard1'][std1_file])
            std1f_ind = int(0)
            std1_arr_ind = int(0)
            std1t = std1f[1].data.field('Time')
            std1f_len = len(std1t)
            xepcu0 = std1f[1].data.field('XeCntPcu0')
            xepcu1 = std1f[1].data.field('XeCntPcu1')
            xepcu2 = std1f[1].data.field('XeCntPcu2')
            xepcu3 = std1f[1].data.field('XeCntPcu3')
            xepcu4 = std1f[1].data.field('XeCntPcu4')
            vpcnt = std1f[1].data.field('VpCnt')
            remcnt = std1f[1].data.field('RemainingCnt')
            vlecnt = std1f[1].data.field('VLECnt')
            time = max(time,std1t[std1f_ind])


                    
    if verbosity>0:
        print('-----XTEPDS finished data pass. Finalizing PDS ---------')

#print "Total intervals:", nrange

    xav = xav/float(nrange)
    x2av = x2av/float(nrange)
    sigma = np.zeros(fft_arr_size//2-1)
    sigma = np.sqrt(abs(x2av - xav*xav)/nrange)
#print sigma
    xax_e = np.zeros(fft_arr_size//2-1)
    xax_e[1:] = (freq[1:]-freq[0:len(freq)-1])*0.5
    xax_e[0] = (freq[1]-freq[0])*0.5


#   REBINNING

    freq_reb = []
    xax_e_reb = []
    power_reb = []
    error_reb = []

    i = int(0)
    co = 1.0

    while i < len(freq):
        
        f1 = freq[i]-xax_e[i]
        f2 = freq[i]+xax_e[i]
        pow = 0.0
        err = 0.0
        l = 0
        nn = int(co)

        while l <= nn and i < len(freq):

            f2 = max(f2,freq[i]+xax_e[i])
            pow += xav[i]
            err += sigma[i]
            i += 1
            l += 1
        
        freq_reb.append(0.5*(f1+f2))
        xax_e_reb.append(0.5*(f2-f1))
        power_reb.append(pow/float(l))
        error_reb.append(err/float(l)/np.sqrt(float(l)))

        co = co*rebin
        
    if verbosity>0: 
        print("%i frequency bins are rebinned into %i."%(len(freq),len(freq_reb)))
        print("-------------------------------------------------------------------")

# Dead time correction

    power_reb_cor = copy(power_reb)

    if dtcor:

        if verbosity>0: 
            if dtmethod == 'zhang':
                print("PCA Dead time correction per Zhang et al. (1995), The Astrophysical Journal, 449:930-935")
            if dtmethod == 'vikhl':
                print("PCA Dead time correction per Vikhlinin et al. (1994), A&A, 287:73-79")
        print("------------------------------------------------------------------")

        td = 10.0e-6
        pcuon = pcuon/float(npcuon)
        rxenon = rxenon/cor_time/pcuon
        rvle = rvle/cor_time/pcuon
#        rvle = rvle/fft_int/float(nrange)/pcuon
        dsvle_value = int(dsvle_value/float(n_dsvle))
        
        if int(dsvle_value) == 0: vledt =12.0e-6
        if int(dsvle_value) == 1: vledt =76.0e-6
        if int(dsvle_value) == 2: vledt =150.0e-6
        if int(dsvle_value) == 3: vledt =500.0e-6

        print("VLE settings: DSVLE = %d, VLEDT = %f, Rvle = %f, PCUon = %f, Rxenon = %f"%(dsvle_value,vledt,rvle,pcuon,rxenon))

        mu = 1.0/(1.0-rvle*vledt)
#    print rxenon,rvle,dsvle_value,vledt,pcuon,npcuon


#power_reb_cor = power_reb

#            print "Start dead time correction"

        if dtmethod == 'zhang':
            for i in range(len(power_reb_cor)):
                power_reb_cor[i] = power_reb[i] + 4.0*rxenon*td*(1.0-td/bin_size)/pcuon
                power_reb_cor[i] += \
                    2.0*rxenon*td*(td/bin_size)*(len(xav)-1.0)/len(xav)*np.cos(np.pi*freq[i]/fnyq)
                power_reb_cor[i] +=  \
                    - 2.0*rvle*rxenon*(np.sin(np.pi*vledt*freq[i])/(3.14159*freq[i]))**2

        if dtmethod == 'vikhl':
            coeff = float(nph)/pcuon/fft_int/float(nrange)
            for i in range(len(power_reb_cor)):
                power_reb_cor[i] = power_reb[i] - coeff*dtvik(td,bin_size,rxenon*mu,freq[i])*mu \
                    - 2.0*rvle*coeff*(np.sin(np.pi*vledt*freq[i])/(np.pi*freq[i]))**2

#            print "End dead time corr"

#    power_reb_cor = zeros(len(power_reb),dtype=float)
    if norm == 'rms':
        coeff = fft_int*float(nrange)/float(nph)
        for i in range(len(power_reb_cor)):
            power_reb_cor[i] = coeff*power_reb_cor[i]
            power_reb[i] = coeff*power_reb[i]
            error_reb[i] = coeff*error_reb[i]



    if (outfile != '' and outfile != 'none'):

        fr_col = pyfits.Column(name='FREQUENCY', format='E', array = freq )
        xax_e_col = pyfits.Column(name='XAX_E', format='E', array = xax_e )
        pow_col = pyfits.Column(name='POWER', format='E', array = xav )
        error_col = pyfits.Column(name='ERROR', format='E', array = sigma )

        fr_col_reb = pyfits.Column(name='FREQUENCY', format='E', array = freq_reb )
        xax_e_col_reb = pyfits.Column(name='XAX_E', format='E', array = xax_e_reb )
        pow_col_reb = pyfits.Column(name='POWER_NODTCOR', format='E', array = power_reb )
        pow_col_reb_cor = pyfits.Column(name='POWER', format='E', array = power_reb_cor )
        error_col_reb = pyfits.Column(name='ERROR', format='E', array = error_reb )


        cols = pyfits.ColDefs([fr_col, xax_e_col, pow_col, error_col])
        cols_reb  = pyfits.ColDefs([fr_col_reb, xax_e_col_reb, pow_col_reb, pow_col_reb_cor, error_col_reb])
        
        hdr = pyfits.Header()
        hdr['EXTNAME'] = 'XTE_PDS'
        hdr['OBSID'] = info['id']
        hdr['TARGET'] = info['target']
        primary_hdu = pyfits.PrimaryHDU(header=hdr)
        
        tbhdu = pyfits.BinTableHDU.from_columns(cols)

        tbhdu1 = pyfits.BinTableHDU.from_columns(cols_reb)
#        tbhdu1.header[''] = 


        fps_hdu_list = pyfits.HDUList([primary_hdu,tbhdu1])

        if os.path.exists(outfile):

            print("File ",outfile," exists. Writing into",str(os.getpid())+outfile)
            outfile = str(os.getpid())+outfile

        fps_hdu_list.writeto(outfile)    

    if dtcor:
        std1f.close()
        filt.close()
        pcu2hk.close()

    print("---XTEPDS Finished -------------------------------------------------------")

    return power_reb_cor,error_reb,freq_reb,xax_e_reb,outfile


#  END of the xtepds function
