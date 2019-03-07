import astropy.io.fits as pyfits
from astropy.table import Table
import numpy as np
import os
import sys
from bs4 import BeautifulSoup
from urllib.request import urlretrieve
import requests

def listFD(url, ext=''):
    page = requests.get(url).text
    #print(page)
    soup = BeautifulSoup(page, 'html.parser')
    return [url + '/' + node.get('href') for node in soup.find_all('a') if node.get('href').endswith(ext)]


def get_xte_data_http(obsid,folders=['obs','pca','acs','stdprod'],download_path = '.'):
    
    cur_path = os.getcwd()
    if download_path == ".":
        download_path == os.getcwd()


    legacy_url = "http://legacy.gsfc.nasa.gov/FTP/xte/data/archive"
    
    idsplit = obsid.split("-")
    if len(idsplit) != 4:
       print("Wrong obsid format!")
       return -1
    
    pid = idsplit[0]
    ao1 = pid[1]
    ao0 = pid[0]

    if (ao1 == "0"): ao = "AO" + ao0
    if (ao0 == "9" and ao1 != "0"): ao = "AO1" + str(int(ao1)-1).strip()
    dirl = "/".join([ao,"P"+pid, obsid])
    base_url = legacy_url+'/'+dirl

#   Making Obsid directory
 
    if (not os.path.exists(download_path)):  
        os.mkdir(download_path)
        
    obs_path = os.path.join(download_path,obsid)
    if (not os.path.exists(obs_path)): 
        os.mkdir(obs_path)


    for folder in folders:

        folder_url = base_url +'/'+folder
        folder_path = os.path.join(obs_path,folder)
        if folder == 'obs': 
            folder_url = base_url
            folder_path = obs_path
        if (not os.path.exists(folder_path)):  
            os.mkdir(folder_path)
        print("Getting "+folder)

        for f in listFD(folder_url):        
            file_name = f.split('/')[-1]
            try:
                if file_name[0] == '?': continue
            except:
                continue
            print("    Getting "+file_name)
            try:
                urlretrieve(f,os.path.join(folder_path,file_name))
            except:
                print(f"Failed to get {f}")

    print("Done")
    
def get_xte_data(obsid,folders=['obs','pca','acs','stdprod'],download_path = '.'):

    import ftplib
    
    cur_path = os.getcwd()
    if download_path == ".":
        download_path == os.getcwd()

    try:
        legacy = ftplib.FTP("legacy.gsfc.nasa.gov")
    except:
        print("Error connecting to legacy")
        
    legacy.login()
    legacy.cwd("xte/data/archive")

    idsplit = obsid.split("-")
    if len(idsplit) != 4:
       print("Wrong obsid format!")
       return -1
    
    pid = idsplit[0]
    ao1 = pid[1]
    ao0 = pid[0]
    
    if (ao1 == "0"): ao = "AO" + ao0
    if (ao0 == "9" and ao1 != "0"): ao = "AO1" + str(int(ao1)-1).strip()
    dirl = "/".join([ao,"P"+pid, obsid])
    
    try:
        legacy.dir(dirl)
        
    except Exception as e:
        print(e)
        return -1
    
    try:
       legacy.cwd(dirl)
      
    except:
        return -1
        
    
#   Making Obsid directory
 
    if (not os.path.exists(download_path)):  
        
        os.system("mkdir "+download_path)
    obs_path = download_path+'/'+obsid
    if (not os.path.exists(obs_path)): os.system("mkdir "+obs_path)
#        os.chdir(download_path+'/'+obsid)
        
    if "obs" in folders:

        list = legacy.nlst()
            
        for f in list:

            if (f[0] == "F"): 
                print("Getting "+f)
                #if os.path.exists(f) and clobber == False: 
                fil = open(obs_path+'/'+f,'wb')
                legacy.retrbinary("RETR "+f,fil.write) 
                fil.close()
        
    for dir in folders:

        if dir == 'obs': continue

        if (not os.path.exists(dir)):  
            os.system("mkdir "+obs_path+'/'+dir)
        print("Getting "+dir)
        #os.chdir(dir)
        #updir = legacy.pwd()
#        print updir
        legacy.cwd(dir)
        list = legacy.nlst()

        for f in list:
            print("    Getting "+f)
            try:
                with open(obs_path+'/'+dir+'/'+f,'wb') as fil: legacy.retrbinary("RETR "+f,fil.write) 
            except:
                print(f"Failed to get {f}")

        legacy.cwd("..")
        #os.chdir("..")

    
#    print list

    legacy.close()
    os.chdir(cur_path)
    #self.init_pcadata()
    print( "Download finished")


def get_init_data(obs_path):
    
    fmi_file = os.path.join(obs_path,"FMI")
    fmi = pyfits.open(fmi_file)
    info = {}
            
    info['id'] = fmi[1].data.field('ObsID')[0]
    info['target'] =  fmi[1].data.field('Source')[0]

    info['startdate'] = fmi[1].data.field('StartDate')[0]
    info['stopdate'] =  fmi[1].data.field('StopDate')[0]
    info['starttime'] = fmi[1].data.field('StartTime')[0]
    info['stoptime'] =  fmi[1].data.field('StopTime')[0]
    stdprod_path  = os.path.join(obs_path,"stdprod")
    
    if os.path.exists(stdprod_path):
        
        idnum = ''.join(info['id'].split('-'))
        filter_file = os.path.join(stdprod_path,f'x{idnum}.xfl.gz')
        gti_file = os.path.join(stdprod_path,f'x{idnum}.gti.gz')
        if os.path.exists(filter_file): 
            info['filter'] = filter_file
        if os.path.exists(gti_file): 
            info['gti'] = gti_file
            
    pca_idx = fmi[1].data.field('PCA_Index_File')[0].strip()
    tmp_path = obs_path+'/'+pca_idx
    if os.path.exists(tmp_path):
        fipc_file = tmp_path
    elif  os.path.exists(tmp_path+'.gz'):
        fipc_file = tmp_path+'.gz'
#            print fipc_file
    fipc = pyfits.open(fipc_file)
#            print fipc
    fi = fipc[1].data
#            print fi.field('PCU2Hk')
    fisize = len(fi.field('ObsID'))

    fh5c_set  = set()
    for f in fi.field('PCU2Hk'):
        fh = f.strip()
        pth = obs_path +'/'+fh
        if (fh != '') and os.path.exists(pth):
            fh5c_set.add(pth)
        elif (f.strip() != '') and os.path.exists(pth+'.gz'):
            fh5c_set.add(pth+'.gz')
    fh5c = list(fh5c_set)
#            for j in range(fisize):
#               if fi.field('PCU2Hk')[j].strip() != '' and self.obsid +'/'+fi.field('PCU2Hk')[j].strip() not in self.fh5c: self.fh5c.append(self.obspath+'/'+fi.field('PCU2Hk')[j].strip())
    info['have_obs'] = True
    info['fh5c'] = fh5c
    
    
    modes = []
    sbmodes = []
    emodes = []
    chpos = []
    bmodes = []
    cbmodes = []
    tlamodes = []
    gxmodes = []
    pfmodes = []
    phase = []
    exposure = dict(Standard1=0.0,Standard2=0.0)
    tres = dict(Standard1=[],Standard2=[])
    channels = dict(Standard1=[],Standard2=[])
    nchan = dict(Standard1=[],Standard2=[])
    pcachan_low = dict(Standard1=[],Standard2=[])
    pcachan_high = dict(Standard1=[],Standard2=[])
    pcadata = dict(Standard1=[],Standard2=[])
        
    if info['have_obs']:
        try:
            for k in [4,7,10,13,16,19]:
                for m in range(fisize):
                    mn = fi.field(k)[m].strip()
                    mf = fi.field(k+1)[m].strip()
                    if mn == '' or mf == '': continue
#                        print mn,mf
                    if mn[:len(mn)-1] == 'Standard1': mn = 'Standard1'
                    if mn[:len(mn)-1] == 'Standard2': mn = 'Standard2'
                    if mn == '': continue
                    if mn[:2] == 'SB' and mn not in sbmodes: 
                        sbmodes.append(mn)
                    if mn[:2] == 'E_' and mn not in emodes: emodes.append(mn)
                    if mn[:2] == 'B_' and mn not in bmodes: bmodes.append(mn)
                    if mn[:3] == 'CB_' and mn not in cbmodes: cbmodes.append(mn)
                    if mn[:4] == 'TLA_' and mn not in tlamodes: tlamodes.append(mn)
                    if mn[:3] == 'PF_' and mn not in pfmodes: pfmodes.append(mn)
#        if mn[:9] == 'GoodXenon' and mn not in gxmodes: gxmodes.append(mn)
                    if mn[:9] == 'GoodXenon': continue
#                    print mn

                    if mn not in modes:
                        modes.append(mn)
                        pcadata[mn] = []
                        pcachan_low[mn] = []
                        pcachan_high[mn] = []
                        exposure[mn] = 0.0
                    
                    fname = os.path.join(obs_path,mf)
                    if not os.path.exists(fname):
                        fname = os.path.join(obs_path,mf+'.gz')
                    fid = pyfits.open(fname)
                    tddes = fid[1].header['TDDES2']
                    #print(tddes)
            #            xte_s = fid[1].header['EXTNAME']
            #            print tddes
                    chan_desc = tddes.split('C[')[1].split(']')[0]
                    
                    ch1,*tmp,ch2 = chan_desc.split('~')
                    #print(ch1,ch2)
                    ch1 = ch1.strip('(')
                    ch2 = ch2.strip(')')
                    #if zzz[0] == '(': zzz = zzz[1:]
                    #www = xxx.split('~')[1]
                    #if www[-1] == ')': www = www[:len(www)-1]
    
                    # getting PCA channel range for a mode from TDDES keyword 
                    try:
                        chans = [int(ch1),int(ch2)]
                    except Exception as e:
                        print("Warning: Unable to get channels for mode:"+mn)
                        print("Using the mode name.")
                        nchan[mn] = mn.split('_')[2][0:-1]
                    else:
                        channels[mn] = chans
                            
                    #   if we have event mode, use TEVTB2 to get detailed channel information
                    if mn[:2] == 'E_':

                        tevtb2 = fid[1].header['TEVTB2']
                        #print(tevtb2)
                        chan_desc = tevtb2.split('C[')[1].split(']')[0]
        
                    chan_split = chan_desc.split(',')
#                        print xxx_split,mn
                    nch = int(0)
                    
                    # bin channels
                    low_list = []
                    high_list = []
                    for s in chan_split:
                        ss = s.strip(')')
                        sss = ss.strip('(')

                        if len(sss.split(';')) == 2:
                            ssss = sss.split(';')
                            n1 = int(ssss[0].split('~')[0])
                            n2 = int(ssss[0].split('~')[1])
                            n3 = int(ssss[1])
                            nch += int((n2-n1+1)/n3)

                            for l in range(int((n2-n1+1)/n3)):
                                low_list.append(int(n1+l*n3))
                                high_list.append(int(n1+(l+1)*n3-1))
#                                    print "1",self.pcachan_low[mn][-1],self.pcachan_high[mn][-1]

                        elif len(sss.split(':')) == 2:
                            ssss = sss.split(':')
                            n1 = int(ssss[0])
                            n2 = int(ssss[1])
                            nch += n2-n1+1

                            for l in range(n2-n1+1):
                                low_list.append(int(n1+l))
                                high_list.append(int(n1+l))
#                                    print "2",self.pcachan_low[mn][-1],self.pcachan_high[mn][-1]  
                        else:
                            nch +=1
                            if sss.find('~') == -1:
                                low_list.append(int(sss))                                    
                                high_list.append(int(sss)) 
#                                    print "3",self.pcachan_low[mn][-1],self.pcachan_high[mn][-1]   
                            else:
                                ssss = sss.split('~')
                                low_list.append(int(ssss[0]))
                                high_list.append(int(ssss[1]))                            
#                                    print "4",self.pcachan_low[mn][-1],self.pcachan_high[mn][-1]   

                    pcachan_low[mn] = low_list
                    pcachan_high[mn] = high_list
                    nchan[mn] = nch
#                        print self.pcachan_low[mn][-1],self.pcachan_high[mn][-1]

                    if mn[:3] == 'PF_':
                            
                        xxx,yyy,*tmp = tddes.split('P[')[1].split(']')

                        for mm in xxx.split(','):
                            for mmm in mm.split(';'):
                                phase.append(mmm)

                    if mn.strip() == "Standard2": 
                        tres[mn] = 16.0
                    elif mn.strip() == "Standard1": 
                        tres[mn] = 0.125
                    elif mn[:2] == 'E_' or mn[:2] == 'D_' or mn[:4] == 'TLA_': 
                        tres[mn] = float(fid[1].header['TIMEDEL'])
                    elif mn[:2] == 'B_' or mn[:2] == 'SB' or mn[:2] == 'PF' or mn[:3] == 'CB_': 
                        tres[mn] = float(tddes.split('T[')[1].split(']')[0].split(';')[1])
    
                    fid.close()
                    exposure[mn] += fi.field(2)[m]-fi.field(1)[m]
                    if mn[:2] == 'E_':
                        coln = fi.columns[k].name[:3]
                        #print(coln[:3])
                        tcode = mf[9:24]
                        #print(tcode)
                        mf = 'pca/SE_'+tcode+'.evt'
                        fname = os.path.join(obs_path,mf)
                        if not os.path.exists(fname):
                            fname = os.path.join(obs_path,mf+'.gz')


                        #mf = fi.field('PCA_'+coln+'_Event')[m]
                        #print(mf,mn)
                    if fname not in pcadata[mn]: 
                        pcadata[mn].append(fname)

            havedata = True
            have_pca = True
            info['modes']   = modes 
            info['sbmodes'] = sbmodes
            info['emodes']  = emodes
            info['chpos']   = chpos
            info['bmodes']  = bmodes
            info['cbmodes']  = cbmodes
            info['tlamodes'] = tlamodes
            info['gxmodes'] = gxmodes
            info['pfmodes'] = pfmodes
            info['phase']   = phase
            info['exposure'] = exposure
            info['tres']     = tres
            info['channels'] = channels
            info['nchan']    = nchan
            info['pcachan_low'] = pcachan_low
            info['pcachan_high'] = pcachan_high
            info['pcadata']  = pcadata
            info['have_pca'] = True
            info['havedata'] = True


        except Exception as e:
            print("Error initializing the observation.")
            #print("Error was: "+str(e))
            print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e).__name__, e)
            info['have_pca'] = False
            info['havedata'] = False

    else:
        info['have_pca'] = False
        info['havedata'] = False

    return info


def suggest_fft_modes(info,verbosity=1):

        if not info['havedata']: return -1
        first_chan = int(8)
        last_chan = first_chan
        fft_modes = []

        if len(info['emodes']) == 0 and len(info['sbmodes']) == 0:
#        print "POWSP: No Event or Single Bit mode found."
            if len(info['bmodes']) > 0:
                if verbosity>0: print("Choosing Binned modes")
                for mm in info['bmodes']: fft_modes.append(mm)

        elif len(info['emodes']) == 1 and info['channels'][info['emodes'][0]][0] == 0 and info['channels'][info['emodes'][0]][1] >= 249:
            #        print "Event covering the whole PCA range is found."
            for mm in info['emodes']: fft_modes.append(mm)

        elif len(info['sbmodes']) > 0:
            sb_max_chan = info['channels'][info['sbmodes'][0]][1]
            for sbm in info['sbmodes']: sb_max_chan = max(info['channels'][sbm][1],sb_max_chan)
            sb_min_chan = info['channels'][info['sbmodes'][0]][0]
            for sbm in info['sbmodes']: sb_min_chan = min(info['channels'][sbm][0],sb_min_chan)
#        if sb_min_chan == 0 and sb_max_chan == 249:
            for mm in info['sbmodes']: fft_modes.append(mm)
            for em in info['emodes']:
            #            print em,channels[em]
                if  info['channels'][em][0] == sb_max_chan +1:
                #                print em,channels[em][0],sb_max_chan
                    fft_modes.append(em)

        elif len(info['bmodes']) > 0:
            b_max_chan = info['channels'][info['bmodes'][0]][1]
            for bm in info['bmodes']: b_max_chan = max(info['channels'][bm][1],b_max_chan)
            b_min_chan = info['channels'][info['bmodes'][0]][0]
            for bm in info['bmodes']: b_min_chan = min(info['channels'][bm][0],b_min_chan)
#        if b_min_chan == 0 and b_max_chan == 249:
            for mm in info['bmodes']: fft_modes.append(mm)
            for em in info['emodes']:
                if  info['channels'][em][0] == b_max_chan +1: 
                    fft_modes.append(info['emodes'][0])
        elif len(info['emodes']) > 1:
            for m in info['emodes']:
                if info['channels'][m][0] == 0 and info['channels'][m][1] >= 249: fft_modes = [m]
    #        print "Event covering the whole PCA range is found."
#          fft_modes = emodes

        if verbosity>0:
            print("Selected modes")
            for m in fft_modes: print(m)

        return fft_modes


def makegti(ffile,expr = "ELV > 5.0 and OFFSET < 0.1 and NUM_PCU_ON > 0",verbosity=1):

    """ Creates a Good Time Intervals (GTI) using condition specified by the cond parameter.
    Returns two lists of GTI start and stop times correspondingly.
    
    Parameters:
    ffile - XTE observation filter file. If filter file is not found, function returns -1
    cond["ELV > 5.0 and OFFSET < 0.1 and NUM_PCU_ON > 0"] 
          - expression to be used to calculate GTIs. Python syntaxis should be used NOT FTOOLS expression.
          Namely "NUM_PCU_ON == 5 and ELV > 5.0" is valid condition expression, while 
          "NUM_PCU_ON.EQ.5.AND.ELV.GT.5.0" is valid a will result in error ar garbage results.
    
    gtifile - name of the resulting GTI file. Default is '' or no GTI file.

    """
 

    if not os.path.exists(ffile):
        print("Can't find  filter file %s"%ffile)
        return -1

    else:
        try:
            filtered_df = Table.read(ffile).to_pandas().query(expr)
            ff = pyfits.open(ffile)
            time = ff[1].data.field('Time') 
            #elv =  ff[1].data.field('ELV')
            #offset = ff[1].data.field('Offset')
        except:
            if (verbosity > 0): print("File %s is not valid XTE filter file. Exit.")
            return -1
        

    tstart = time[0]
    tstop = time[-1]
    step = ff[1].header['TIMEDEL'] + 0.1

    gti_start = filtered_df[filtered_df['Time'] - filtered_df['Time'].shift(1) > step].Time.values
    gti_end = filtered_df[filtered_df['Time'] - filtered_df['Time'].shift(-1) < -step].Time.values
    
    i = 0
    gti_start = np.concatenate([[tstart],gti_start])
    gti_stop = np.concatenate([[tstop],gti_end])

    return gti_start,gti_stop

    #cols =  ff[1].columns.names
    #cond_split = expr.split()
    #i = 0
    #ingti = False

    #while i<(len(time)):


#Substitue values from a filter file into a condition string
     #   cond_split = expr.split()        
     #   for k in range(len(cond_split)):

     #       if cond_split[k] in cols:
      #          col = ff[1].data.field(cond_split[k])
     #           cond_split[k] = str(ff[1].data.field(cond_split[k])[i])

      #  c = " ".join(cond_split)
      #  print(c)
      #  exec("result = "+c)

       # if not ingti and result:
       #     gti_start.append(time[i])
       #     ingti = True
       # 
       # if ingti and not result:
       #     gti_stop.append(time[i])
       #     ingti = False

            
        #i += 1

    #if ingti:
    #    gti_stop.append(time[i-1])      

#    print gti_start,gti_stop

    #return gti_start,gti_stop

