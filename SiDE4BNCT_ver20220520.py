#!/usr/bin/env python
# ##### import libraries ######
from pickle import TRUE
import tkinter as tk
from tkinter import filedialog
import tkinter.ttk as ttk
import os
import math
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import scipy as sci
from scipy import interpolate
import struct
from struct import unpack
import time
import sys
import gc
from PIL import Image, ImageTk

##### version #####
VERSION = "ver20220520 (GUI)"

################################
########## functions ###########
################################

# acquire file path
def FilePathFunc(iDir, inBox):
    fTyp = [("all files", "*.dat; *.out; *.txt")]
    fPath = filedialog.askopenfilename(initialdir = iDir, filetypes = fTyp)
    inBox.delete(0, tk.END)
    inBox.insert(tk.END, os.path.split(fPath)[-1])
    if inBox.get() != '':
        statusBar["text"] = 'Specified file!'
    else:
        statusBar["text"] = 'ERROR: Failed to specify file!!'
# acquire directory path
def DirectoryPathFunc(inBox):
    iDir = r"./"
    dirPath = filedialog.askdirectory(initialdir = iDir)
    inBox.delete(0, tk.END)
    inBox.insert(tk.END, dirPath)
    if inBox.get() != '':
        statusBar["text"] = 'Specified directory!'
    else:
        statusBar["text"] = 'ERROR: Failed to specify directory!!'
# acquire final directory ID
def FinalDirectoryNameFunc(iDir, inBox):
    dirPath = filedialog.askdirectory(initialdir = iDir)
    inBox.delete(0, tk.END)
    inBox.insert(tk.END, os.path.split(dirPath)[-1])
    if inBox.get() != '':
        statusBar["text"] = 'Specified directory!'
    else:
        statusBar["text"] = 'ERROR: Failed to specify directory!!'
# spe2dpf for neutron DPFs
def spe2dpf4neutronDPFsFunc(inBox_workDir,
                            inBox_tabN, inBox_tabGpri, inBox_tabGsec,
                            inBox_engListN, inBox_engListG,
                            inBox_gkerma, inBox_nkerma, inBox_bkerma,
                            inBox_CbN, inBox_TNR,
                            inBox_CBEt, inBox_CBEm, inBox_CBEs, inBox_CBEn,
                            inBox_RBEgamma, inBox_RBEthermal, inBox_RBEfast,
                            inBox_BINS_MONO_N, inBox_BINS_MONO_G,
                            inBox_BINS_SPE_N, inBox_BINS_SPE_G_PRI, inBox_BINS_SPE_G_SEC,
                            inBox_EMIN_SPE_N, inBox_EMIN_SPE_G_PRI, inBox_EMIN_SPE_G_SEC,
                            inBox_EMAX_SPE_N, inBox_EMAX_SPE_G_PRI, inBox_EMAX_SPE_G_SEC,
                            inBox_TALLIES,
                            inBox_nDPFName):
    # check PATH
    if inBox_workDir.get() != '' and os.path.isdir(inBox_workDir.get()) == True:
        pass
    else:
        print(" ERROR: Working directory does not exist!")
        statusBar["text"] = 'ERROR: Working directory does not exist!'
        return
    if inBox_tabN.get() != '' and os.path.isdir(inBox_tabN.get()) == True:
        pass
    else:
        print(" ERROR: Directory of neutron table data does not exist!")
        statusBar["text"] = 'ERROR: Directory of neutron table data does not exist!'
        return
    if inBox_nkerma.get() != '' and os.path.isfile(os.path.join(r"kerma", inBox_nkerma.get())) == True:
        pass
    else:
        print(" ERROR: Neutron kerma file does not exist!")
        statusBar["text"] = 'ERROR: Neutron kerma file does not exist!'
        return
    if (inBox_BINS_MONO_N.get() != ''
        and inBox_BINS_MONO_G.get() != ''
        and inBox_BINS_SPE_N.get() != ''
        and inBox_BINS_SPE_G_PRI.get() != ''
        and inBox_BINS_SPE_G_SEC.get() != ''
        and inBox_EMIN_SPE_N.get() != ''
        and inBox_EMIN_SPE_G_PRI.get() != ''
        and inBox_EMIN_SPE_G_SEC.get() != ''
        and inBox_EMAX_SPE_N.get() != ''
        and inBox_EMAX_SPE_G_PRI.get() != ''
        and inBox_EMAX_SPE_G_SEC.get() != ''
        and inBox_TALLIES.get() != ''):
        pass
    else:
        print(" ERROR: Some or all of parameters in table data are missing!")
        statusBar["text"] = 'ERROR: Some or all of parameters in table data are missing!'
        return
    if (inBox_CbN.get() != ''
        and inBox_TNR.get() != ''
        and inBox_CBEt.get() != ''
        and inBox_CBEm.get() != ''
        and inBox_CBEs.get() != ''
        and inBox_CBEn.get() != ''
        and inBox_RBEgamma.get() != ''
        and inBox_RBEthermal.get() != ''
        and inBox_RBEfast.get() != ''):
        pass
    else:
        print(" ERROR: Some or all of radiobiological coefficients are missing!")
        statusBar["text"] = 'ERROR: Some or all of radiobiological coefficients are missing!'
        return
    if inBox_engListN.get() != '' and os.path.isfile(os.path.join(r"phits", inBox_engListN.get())):
        pass
    else:
        print(" ERROR: Monoenergy list of neutron does not exist!")
        statusBar["text"] = 'ERROR: Monoenergy list of neuton does not exist!'
        return
    if inBox_nDPFName.get() != '':
        pass
    else:
        print(" ERROR: Directory name to generate DPF files is missing!")
        statusBar["text"] = 'ERROR: Directory name to generate DPF files is missing!'
        return
    # main
    start_time = time.time()
    print(' ############################################')
    print(' ########## spe2dpf (neutron DPFs) ##########')
    print(' ############################################')
    print(' Computing & outputting...')
    # generate DPFs storaging directory
    dpf_path = os.path.join(inBox_workDir.get(), inBox_nDPFName.get())
    os.makedirs(dpf_path, exist_ok=True)
    # constants
    PI = PI = 3.141592653589793
    # read parameters
    CbN = float(inBox_CbN.get())
    TNR = float(inBox_TNR.get())
    CBEt = float(inBox_CBEt.get())
    CBEm = float(inBox_CBEm.get())
    CBEs = float(inBox_CBEs.get())
    CBEn = float(inBox_CBEn.get())
    RBEgamma = float(inBox_RBEgamma.get())
    RBEthermal = float(inBox_RBEthermal.get())
    RBEfast = float(inBox_RBEfast.get())
    BINS_MONO_N = int(inBox_BINS_MONO_N.get())
    BINS_MONO_G = int(inBox_BINS_MONO_G.get())
    BINS_SPE_N = int(inBox_BINS_SPE_N.get())
    BINS_SPE_G_PRI = int(inBox_BINS_SPE_G_PRI.get())
    BINS_SPE_G_SEC = int(inBox_BINS_SPE_G_SEC.get())
    EMIN_SPE_N = float(inBox_EMIN_SPE_N.get())
    EMIN_SPE_G_PRI = float(inBox_EMIN_SPE_G_PRI.get())
    EMIN_SPE_G_SEC = float(inBox_EMIN_SPE_G_SEC.get())
    EMAX_SPE_N = float(inBox_EMAX_SPE_N.get())
    EMAX_SPE_G_PRI = float(inBox_EMAX_SPE_G_PRI.get())
    EMAX_SPE_G_SEC = float(inBox_EMAX_SPE_G_SEC.get())
    TALLIES = int(inBox_TALLIES.get())
    # loading values of monoenergy (neutron)
    LIST_E_mono_n = []
    LIST_E_mono_n_flt = []
    data_path = os.path.join(r"phits", inBox_engListN.get())
    with open(data_path, mode = 'r') as f:
        for i in range(BINS_MONO_N):
            dum = f.readline()
            dum = dum.replace("\n", "")
            LIST_E_mono_n.append(dum)
            LIST_E_mono_n_flt.append(float(dum))
    # loading and inter/extrapolating neutron kerma coefficients
    LIST_eng_nkerma = []
    LIST_val_nkerma_hyd = []
    LIST_val_nkerma_nit = []
    LIST_val_nkerma_oth = []
    LIST_val_nkerma_n = []
    LIST_nkerma_hyd = []
    LIST_nkerma_nit = []
    LIST_nkerma_oth = []
    LIST_nkerma_n = []
    with open(os.path.join(r"kerma", inBox_nkerma.get()), mode = 'r') as f:
        dum = f.readline() # skipping header part
        for i in range(134):
            dum = f.readline()
            dum = dum.replace("\n", "")
            dum = dum.split()
            LIST_eng_nkerma.append(math.log(float(dum[0])))
            LIST_val_nkerma_hyd.append(math.log(float(dum[1])))
            LIST_val_nkerma_nit.append(math.log(float(dum[2])))
            LIST_val_nkerma_oth.append(math.log(float(dum[3])))
            LIST_val_nkerma_n.append(math.log(float(dum[4])))
        nkerma_hyd = interpolate.interp1d(LIST_eng_nkerma, LIST_val_nkerma_hyd, kind='linear', fill_value='extrapolate') # usage: math.exp(kerma(math.log(ENG)))
        nkerma_nit = interpolate.interp1d(LIST_eng_nkerma, LIST_val_nkerma_nit, kind='linear', fill_value='extrapolate') # usage: math.exp(kerma(math.log(ENG)))
        nkerma_oth = interpolate.interp1d(LIST_eng_nkerma, LIST_val_nkerma_oth, kind='linear', fill_value='extrapolate') # usage: math.exp(kerma(math.log(ENG)))
        nkerma_n = interpolate.interp1d(LIST_eng_nkerma, LIST_val_nkerma_n, kind='linear', fill_value='extrapolate') # usage: math.exp(kerma(math.log(ENG)))
    for i in range(BINS_MONO_N):
        LIST_nkerma_hyd.append(math.exp(nkerma_hyd(math.log(LIST_E_mono_n_flt[i]))))
        LIST_nkerma_nit.append(math.exp(nkerma_nit(math.log(LIST_E_mono_n_flt[i]))))
        LIST_nkerma_oth.append(math.exp(nkerma_oth(math.log(LIST_E_mono_n_flt[i]))))
        LIST_nkerma_n.append(math.exp(nkerma_n(math.log(LIST_E_mono_n_flt[i]))))
    plt.plot(LIST_E_mono_n_flt, LIST_nkerma_hyd, linestyle="dotted", color="cyan", label="hydrogen")
    plt.plot(LIST_E_mono_n_flt, LIST_nkerma_nit, linestyle="dotted", color="yellowgreen", label="nitrogen")
    plt.plot(LIST_E_mono_n_flt, LIST_nkerma_oth, linestyle="dotted", color="magenta", label="other nuclei")
    plt.plot(LIST_E_mono_n_flt, LIST_nkerma_n, linestyle="solid", color="red", label="ICRU soft tissue")
    plt.legend()
    plt.title('neutron kerma coeff. ')
    plt.xlabel('neutron energy [MeV]')
    plt.ylabel('kerma coeff. [pGy $\mathregular{cm^2}$]')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(color='r', linestyle='dotted', linewidth=1)
    plt.savefig(inBox_workDir.get() + r'\nkerma.png')
    plt.close()
    # generating neutron DPFs
    data_path = inBox_tabN.get()
    out_path = os.path.join(dpf_path, r"hyd.dpf")
    with open(out_path, mode = 'w')  as f:
        f.write('#constants: borofalan concentartion in normal tissue %.1f, tumor/normal tissue ratio of borofalan concentration %.1f, RBE for gamma-ray %.1f, RBE for nitrogen dose %.1f, RBE for hydrogen dose and other neutron dose %.1f, CBE for tumor %.1f, CBE for mucosa %.1f, CBE for skin %.1f, CBE for other normal tissue %.2f\n' %(CbN, TNR, RBEgamma, RBEthermal, RBEfast, CBEt, CBEm, CBEs, CBEn))
    out_path = os.path.join(dpf_path, r"nit.dpf")
    with open(out_path, mode = 'w')  as f:
        f.write('#constants: borofalan concentartion in normal tissue %.1f, tumor/normal tissue ratio of borofalan concentration %.1f, RBE for gamma-ray %.1f, RBE for nitrogen dose %.1f, RBE for hydrogen dose and other neutron dose %.1f, CBE for tumor %.1f, CBE for mucosa %.1f, CBE for skin %.1f, CBE for other normal tissue %.2f\n' %(CbN, TNR, RBEgamma, RBEthermal, RBEfast, CBEt, CBEm, CBEs, CBEn))
    out_path = os.path.join(dpf_path, r"oth.dpf")
    with open(out_path, mode = 'w')  as f:
        f.write('#constants: borofalan concentartion in normal tissue %.1f, tumor/normal tissue ratio of borofalan concentration %.1f, RBE for gamma-ray %.1f, RBE for nitrogen dose %.1f, RBE for hydrogen dose and other neutron dose %.1f, CBE for tumor %.1f, CBE for mucosa %.1f, CBE for skin %.1f, CBE for other normal tissue %.2f\n' %(CbN, TNR, RBEgamma, RBEthermal, RBEfast, CBEt, CBEm, CBEs, CBEn))
    out_path = os.path.join(dpf_path, r"hyd_err.dpf")
    with open(out_path, mode = 'w')  as f:
        f.write('#constants: borofalan concentartion in normal tissue %.1f, tumor/normal tissue ratio of borofalan concentration %.1f, RBE for gamma-ray %.1f, RBE for nitrogen dose %.1f, RBE for hydrogen dose and other neutron dose %.1f, CBE for tumor %.1f, CBE for mucosa %.1f, CBE for skin %.1f, CBE for other normal tissue %.2f\n' %(CbN, TNR, RBEgamma, RBEthermal, RBEfast, CBEt, CBEm, CBEs, CBEn))
    out_path = os.path.join(dpf_path, r"nit_err.dpf")
    with open(out_path, mode = 'w')  as f:
        f.write('#constants: borofalan concentartion in normal tissue %.1f, tumor/normal tissue ratio of borofalan concentration %.1f, RBE for gamma-ray %.1f, RBE for nitrogen dose %.1f, RBE for hydrogen dose and other neutron dose %.1f, CBE for tumor %.1f, CBE for mucosa %.1f, CBE for skin %.1f, CBE for other normal tissue %.2f\n' %(CbN, TNR, RBEgamma, RBEthermal, RBEfast, CBEt, CBEm, CBEs, CBEn))
    out_path = os.path.join(dpf_path, r"oth_err.dpf")
    with open(out_path, mode = 'w')  as f:
        f.write('#constants: borofalan concentartion in normal tissue %.1f, tumor/normal tissue ratio of borofalan concentration %.1f, RBE for gamma-ray %.1f, RBE for nitrogen dose %.1f, RBE for hydrogen dose and other neutron dose %.1f, CBE for tumor %.1f, CBE for mucosa %.1f, CBE for skin %.1f, CBE for other normal tissue %.2f\n' %(CbN, TNR, RBEgamma, RBEthermal, RBEfast, CBEt, CBEm, CBEs, CBEn))
    for i in range(BINS_MONO_N):
        data_name = LIST_E_mono_n[i]+'MeV-nspe.dat'
        depth = 0
        with open(os.path.join(data_path, data_name), mode = 'r') as f:
            DPF_hyd = []
            DPF_nit = []
            DPF_oth = []
            DPF_hyd_err = []
            DPF_nit_err = []
            DPF_oth_err = []
            while True:
                if depth >= TALLIES:
                    break
                dum = f.readline()
                dum = dum.replace("\n", "")
                dum = dum.split()
                if len(dum) < 4:
                    continue
                if dum[1] == 'e-lower':
                    dum_hyd = 0
                    dum_nit = 0
                    dum_oth = 0
                    dum_hyd_err = 0
                    dum_nit_err = 0
                    dum_oth_err = 0
                    for j in range(BINS_SPE_N):
                        dum = f.readline()
                        dum = dum.replace("\n", "")
                        dum = dum.split()
                        dum_hyd += float(dum[2])*math.exp(nkerma_hyd(math.log((float(dum[0])*float(dum[1]))**0.5)))*RBEfast*PI*5*5
                        dum_nit += float(dum[2])*math.exp(nkerma_nit(math.log((float(dum[0])*float(dum[1]))**0.5)))*RBEthermal*PI*5*5
                        dum_oth += float(dum[2])*math.exp(nkerma_oth(math.log((float(dum[0])*float(dum[1]))**0.5)))*RBEfast*PI*5*5
                        dum_hyd_err += (float(dum[2])*float(dum[3])*math.exp(nkerma_hyd(math.log((float(dum[0])*float(dum[1]))**0.5)))*RBEfast*PI*5*5)**2
                        dum_nit_err += (float(dum[2])*float(dum[3])*math.exp(nkerma_nit(math.log((float(dum[0])*float(dum[1]))**0.5)))*RBEthermal*PI*5*5)**2
                        dum_oth_err += (float(dum[2])*float(dum[3])*math.exp(nkerma_oth(math.log((float(dum[0])*float(dum[1]))**0.5)))*RBEfast*PI*5*5)**2
                    DPF_hyd.append(dum_hyd)
                    DPF_nit.append(dum_nit)
                    DPF_oth.append(dum_oth)
                    DPF_hyd_err.append(dum_hyd_err**0.5)
                    DPF_nit_err.append(dum_nit_err**0.5)
                    DPF_oth_err.append(dum_oth_err**0.5)
                    depth += 1
        out_path = os.path.join(dpf_path, r"hyd.dpf")
        with open(out_path, mode = 'a')  as f:
            for depth in range(TALLIES):
                f.write('%10.4e\t' %(DPF_hyd[depth]))
            f.write('\n')
        out_path = os.path.join(dpf_path, r"nit.dpf")
        with open(out_path, mode = 'a')  as f:
            for depth in range(TALLIES):
                f.write('%10.4e\t' %(DPF_nit[depth]))
            f.write('\n')
        out_path = os.path.join(dpf_path, r"oth.dpf")
        with open(out_path, mode = 'a')  as f:
            for depth in range(TALLIES):
                f.write('%10.4e\t' %(DPF_oth[depth]))
            f.write('\n')
        out_path = os.path.join(dpf_path, r"hyd_err.dpf")
        with open(out_path, mode = 'a')  as f:
            for depth in range(TALLIES):
                f.write('%10.4e\t' %(DPF_hyd_err[depth]))
            f.write('\n')
        out_path = os.path.join(dpf_path, r"nit_err.dpf")
        with open(out_path, mode = 'a')  as f:
            for depth in range(TALLIES):
                f.write('%10.4e\t' %(DPF_nit_err[depth]))
            f.write('\n')
        out_path = os.path.join(dpf_path, r"oth_err.dpf")
        with open(out_path, mode = 'a')  as f:
            for depth in range(TALLIES):
                f.write('%10.4e\t' %(DPF_oth_err[depth]))
            f.write('\n')
    print(' spe2dpf (neutron DPFs) successfully finished!')
    statusBar["text"] = 'Finished computing neutron DPFs! (Processing time: %.4f min)' %((time.time()-start_time)/60)
# spe2dpf for boron DPFs
def spe2dpf4boronDPFsFunc(inBox_workDir,
                          inBox_tabN, inBox_tabGpri, inBox_tabGsec,
                          inBox_engListN, inBox_engListG,
                          inBox_gkerma, inBox_nkerma, inBox_bkerma,
                          inBox_CbN, inBox_TNR,
                          inBox_CBEt, inBox_CBEm, inBox_CBEs, inBox_CBEn,
                          inBox_RBEgamma, inBox_RBEthermal, inBox_RBEfast,
                          inBox_BINS_MONO_N, inBox_BINS_MONO_G,
                          inBox_BINS_SPE_N, inBox_BINS_SPE_G_PRI, inBox_BINS_SPE_G_SEC,
                          inBox_EMIN_SPE_N, inBox_EMIN_SPE_G_PRI, inBox_EMIN_SPE_G_SEC,
                          inBox_EMAX_SPE_N, inBox_EMAX_SPE_G_PRI, inBox_EMAX_SPE_G_SEC,
                          inBox_TALLIES,
                          inBox_bDPFName):
    # check PATH
    if inBox_workDir.get() != '' and os.path.isdir(inBox_workDir.get()) == True:
        pass
    else:
        print(" ERROR: Working directory does not exist!")
        statusBar["text"] = 'ERROR: Working directory does not exist!'
        return
    if inBox_tabN.get() != '' and os.path.isdir(inBox_tabN.get()) == True:
        pass
    else:
        print(" ERROR: Directory of neutron table data does not exist!")
        statusBar["text"] = 'ERROR: Directory of neutron table data does not exist!'
        return
    if inBox_bkerma.get() != '' and os.path.isfile(os.path.join(r"kerma", inBox_bkerma.get())) == True:
        pass
    else:
        print(" ERROR: Boron kerma file does not exist!")
        statusBar["text"] = 'ERROR: Boron kerma file does not exist!'
        return
    if (inBox_BINS_MONO_N.get() != ''
        and inBox_BINS_MONO_G.get() != ''
        and inBox_BINS_SPE_N.get() != ''
        and inBox_BINS_SPE_G_PRI.get() != ''
        and inBox_BINS_SPE_G_SEC.get() != ''
        and inBox_EMIN_SPE_N.get() != ''
        and inBox_EMIN_SPE_G_PRI.get() != ''
        and inBox_EMIN_SPE_G_SEC.get() != ''
        and inBox_EMAX_SPE_N.get() != ''
        and inBox_EMAX_SPE_G_PRI.get() != ''
        and inBox_EMAX_SPE_G_SEC.get() != ''
        and inBox_TALLIES.get() != ''):
        pass
    else:
        print(" ERROR: Some or all of parameters in table data are missing!")
        statusBar["text"] = 'ERROR: Some or all of parameters in table data are missing!'
        return
    if (inBox_CbN.get() != ''
        and inBox_TNR.get() != ''
        and inBox_CBEt.get() != ''
        and inBox_CBEm.get() != ''
        and inBox_CBEs.get() != ''
        and inBox_CBEn.get() != ''
        and inBox_RBEgamma.get() != ''
        and inBox_RBEthermal.get() != ''
        and inBox_RBEfast.get() != ''):
        pass
    else:
        print(" ERROR: Some or all of radiobiological coefficients are missing!")
        statusBar["text"] = 'ERROR: Some or all of radiobiological coefficients are missing!'
        return
    if inBox_engListN.get() != '' and os.path.isfile(os.path.join(r"phits", inBox_engListN.get())):
        pass
    else:
        print(" ERROR: Monoenergy list of neutron does not exist!")
        statusBar["text"] = 'ERROR: Monoenergy list of neuton does not exist!'
        return
    if inBox_bDPFName.get() != '':
        pass
    else:
        print(" ERROR: Directory name to generate DPF files is missing!")
        statusBar["text"] = 'ERROR: Directory name to generate DPF files is missing!'
        return
    # main
    start_time = time.time()
    print('##########################################')
    print('########## spe2dpf (boron DPFs) ##########')
    print('##########################################')
    print(' Computing & outputting...')
    # generate DPFs storaging directory
    dpf_path = os.path.join(inBox_workDir.get(), inBox_bDPFName.get())
    os.makedirs(dpf_path, exist_ok=True)
    # constants
    PI = PI = 3.141592653589793
    # read parameters
    CbN = float(inBox_CbN.get())
    TNR = float(inBox_TNR.get())
    CBEt = float(inBox_CBEt.get())
    CBEm = float(inBox_CBEm.get())
    CBEs = float(inBox_CBEs.get())
    CBEn = float(inBox_CBEn.get())
    RBEgamma = float(inBox_RBEgamma.get())
    RBEthermal = float(inBox_RBEthermal.get())
    RBEfast = float(inBox_RBEfast.get())
    BINS_MONO_N = int(inBox_BINS_MONO_N.get())
    BINS_MONO_G = int(inBox_BINS_MONO_G.get())
    BINS_SPE_N = int(inBox_BINS_SPE_N.get())
    BINS_SPE_G_PRI = int(inBox_BINS_SPE_G_PRI.get())
    BINS_SPE_G_SEC = int(inBox_BINS_SPE_G_SEC.get())
    EMIN_SPE_N = float(inBox_EMIN_SPE_N.get())
    EMIN_SPE_G_PRI = float(inBox_EMIN_SPE_G_PRI.get())
    EMIN_SPE_G_SEC = float(inBox_EMIN_SPE_G_SEC.get())
    EMAX_SPE_N = float(inBox_EMAX_SPE_N.get())
    EMAX_SPE_G_PRI = float(inBox_EMAX_SPE_G_PRI.get())
    EMAX_SPE_G_SEC = float(inBox_EMAX_SPE_G_SEC.get())
    TALLIES = int(inBox_TALLIES.get())
    # loading values of monoenergy (neutron)
    LIST_E_mono_n = []
    LIST_E_mono_n_flt = []
    data_path = os.path.join(r"phits", inBox_engListN.get())
    with open(data_path, mode = 'r') as f:
        for i in range(BINS_MONO_N):
            dum = f.readline()
            dum = dum.replace("\n", "")
            LIST_E_mono_n.append(dum)
            LIST_E_mono_n_flt.append(float(dum))
    # loading and inter/extrapolating boron kerma coefficients
    LIST_eng_bkerma = []
    LIST_val_bkerma = []
    LIST_bkerma = []
    with open(os.path.join(r"kerma", inBox_bkerma.get()), mode = 'r') as f:
        dum = f.readline() # skipping header part
        for i in range(11):
            dum = f.readline()
            dum = dum.replace("\n", "")
            dum = dum.split()
            LIST_eng_bkerma.append(math.log(float(dum[0])*1e-6)) # convert eV -> MeV
            LIST_val_bkerma.append(math.log(float(dum[1])))
        bkerma = interpolate.interp1d(LIST_eng_bkerma, LIST_val_bkerma, kind='linear', fill_value='extrapolate') # usage: math.exp(kerma(math.log(ENG)))
    for i in range(BINS_MONO_N):
        LIST_bkerma.append(math.exp(bkerma(math.log(LIST_E_mono_n_flt[i]))))
    plt.plot(LIST_E_mono_n_flt, LIST_bkerma, linestyle="solid", color="blue", label="ICRU soft tissue")
    plt.title('boron kerma coeff. ')
    plt.xlabel('neutron energy [MeV]')
    plt.ylabel('kerma coeff. [pGy $\mathregular{cm^2}$/\u03bcg$\mathregular{(^{10}}$B)]')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(color='r', linestyle='dotted', linewidth=1)
    plt.savefig(inBox_workDir.get() + r"\bkerma.png")
    plt.close()
    # generating boron DPFs
    data_path = inBox_tabN.get()
    out_path = os.path.join(dpf_path, r"boron_normal.dpf")
    with open(out_path, mode = 'w')  as f:
        f.write('#constants: borofalan concentartion in normal tissue %.1f, tumor/normal tissue ratio of borofalan concentration %.1f, RBE for gamma-ray %.1f, RBE for nitrogen dose %.1f, RBE for hydrogen dose and other neutron dose %.1f, CBE for tumor %.1f, CBE for mucosa %.1f, CBE for skin %.1f, CBE for other normal tissue %.2f\n' %(CbN, TNR, RBEgamma, RBEthermal, RBEfast, CBEt, CBEm, CBEs, CBEn))
    out_path = os.path.join(dpf_path, r"boron_skin.dpf")
    with open(out_path, mode = 'w')  as f:
        f.write('#constants: borofalan concentartion in normal tissue %.1f, tumor/normal tissue ratio of borofalan concentration %.1f, RBE for gamma-ray %.1f, RBE for nitrogen dose %.1f, RBE for hydrogen dose and other neutron dose %.1f, CBE for tumor %.1f, CBE for mucosa %.1f, CBE for skin %.1f, CBE for other normal tissue %.2f\n' %(CbN, TNR, RBEgamma, RBEthermal, RBEfast, CBEt, CBEm, CBEs, CBEn))
    out_path = os.path.join(dpf_path, r"boron_mucosa.dpf")
    with open(out_path, mode = 'w')  as f:
        f.write('#constants: borofalan concentartion in normal tissue %.1f, tumor/normal tissue ratio of borofalan concentration %.1f, RBE for gamma-ray %.1f, RBE for nitrogen dose %.1f, RBE for hydrogen dose and other neutron dose %.1f, CBE for tumor %.1f, CBE for mucosa %.1f, CBE for skin %.1f, CBE for other normal tissue %.2f\n' %(CbN, TNR, RBEgamma, RBEthermal, RBEfast, CBEt, CBEm, CBEs, CBEn))
    out_path = os.path.join(dpf_path, r"boron_tumor.dpf")
    with open(out_path, mode = 'w')  as f:
        f.write('#constants: borofalan concentartion in normal tissue %.1f, tumor/normal tissue ratio of borofalan concentration %.1f, RBE for gamma-ray %.1f, RBE for nitrogen dose %.1f, RBE for hydrogen dose and other neutron dose %.1f, CBE for tumor %.1f, CBE for mucosa %.1f, CBE for skin %.1f, CBE for other normal tissue %.2f\n' %(CbN, TNR, RBEgamma, RBEthermal, RBEfast, CBEt, CBEm, CBEs, CBEn))
    out_path = os.path.join(dpf_path, r"boron_normal_err.dpf")
    with open(out_path, mode = 'w')  as f:
        f.write('#constants: borofalan concentartion in normal tissue %.1f, tumor/normal tissue ratio of borofalan concentration %.1f, RBE for gamma-ray %.1f, RBE for nitrogen dose %.1f, RBE for hydrogen dose and other neutron dose %.1f, CBE for tumor %.1f, CBE for mucosa %.1f, CBE for skin %.1f, CBE for other normal tissue %.2f\n' %(CbN, TNR, RBEgamma, RBEthermal, RBEfast, CBEt, CBEm, CBEs, CBEn))
    out_path = os.path.join(dpf_path, r"boron_skin_err.dpf")
    with open(out_path, mode = 'w')  as f:
        f.write('#constants: borofalan concentartion in normal tissue %.1f, tumor/normal tissue ratio of borofalan concentration %.1f, RBE for gamma-ray %.1f, RBE for nitrogen dose %.1f, RBE for hydrogen dose and other neutron dose %.1f, CBE for tumor %.1f, CBE for mucosa %.1f, CBE for skin %.1f, CBE for other normal tissue %.2f\n' %(CbN, TNR, RBEgamma, RBEthermal, RBEfast, CBEt, CBEm, CBEs, CBEn))
    out_path = os.path.join(dpf_path, r"boron_mucosa_err.dpf")
    with open(out_path, mode = 'w')  as f:
        f.write('#constants: borofalan concentartion in normal tissue %.1f, tumor/normal tissue ratio of borofalan concentration %.1f, RBE for gamma-ray %.1f, RBE for nitrogen dose %.1f, RBE for hydrogen dose and other neutron dose %.1f, CBE for tumor %.1f, CBE for mucosa %.1f, CBE for skin %.1f, CBE for other normal tissue %.2f\n' %(CbN, TNR, RBEgamma, RBEthermal, RBEfast, CBEt, CBEm, CBEs, CBEn))
    out_path = os.path.join(dpf_path, r"boron_tumor_err.dpf")
    with open(out_path, mode = 'w')  as f:
        f.write('#constants: borofalan concentartion in normal tissue %.1f, tumor/normal tissue ratio of borofalan concentration %.1f, RBE for gamma-ray %.1f, RBE for nitrogen dose %.1f, RBE for hydrogen dose and other neutron dose %.1f, CBE for tumor %.1f, CBE for mucosa %.1f, CBE for skin %.1f, CBE for other normal tissue %.2f\n' %(CbN, TNR, RBEgamma, RBEthermal, RBEfast, CBEt, CBEm, CBEs, CBEn))
    for i in range(BINS_MONO_N):
        data_name = LIST_E_mono_n[i]+'MeV-nspe.dat'
        depth = 0
        with open(os.path.join(data_path, data_name), mode = 'r') as f:
            DPF_boron_normal = []
            DPF_boron_skin = []
            DPF_boron_mucosa = []
            DPF_boron_tumor = []
            DPF_boron_normal_err = []
            DPF_boron_skin_err = []
            DPF_boron_mucosa_err = []
            DPF_boron_tumor_err = []
            while True:
                if depth >= TALLIES:
                    break
                dum = f.readline()
                dum = dum.replace("\n", "")
                dum = dum.split()
                if len(dum) < 4:
                    continue
                if dum[1] == 'e-lower':
                    dum_boron_normal = 0
                    dum_boron_skin = 0
                    dum_boron_mucosa = 0
                    dum_boron_tumor = 0
                    dum_boron_normal_err = 0
                    dum_boron_skin_err = 0
                    dum_boron_mucosa_err = 0
                    dum_boron_tumor_err = 0
                    for j in range(BINS_SPE_N):
                        dum = f.readline()
                        dum = dum.replace("\n", "")
                        dum = dum.split()
                        dum_boron_normal += float(dum[2])*math.exp(bkerma(math.log((float(dum[0])*float(dum[1]))**0.5)))*CbN*CBEn*PI*5*5
                        dum_boron_skin += float(dum[2])*math.exp(bkerma(math.log((float(dum[0])*float(dum[1]))**0.5)))*CbN*CBEs*PI*5*5
                        dum_boron_mucosa += float(dum[2])*math.exp(bkerma(math.log((float(dum[0])*float(dum[1]))**0.5)))*CbN*CBEm*PI*5*5
                        dum_boron_tumor += float(dum[2])*math.exp(bkerma(math.log((float(dum[0])*float(dum[1]))**0.5)))*CbN*TNR*CBEt*PI*5*5
                        dum_boron_normal_err += (float(dum[2])*float(dum[3])*math.exp(bkerma(math.log((float(dum[0])*float(dum[1]))**0.5)))*CbN*CBEn*PI*5*5)**2
                        dum_boron_skin_err += (float(dum[2])*float(dum[3])*math.exp(bkerma(math.log((float(dum[0])*float(dum[1]))**0.5)))*CbN*CBEs*PI*5*5)**2
                        dum_boron_mucosa_err += (float(dum[2])*float(dum[3])*math.exp(bkerma(math.log((float(dum[0])*float(dum[1]))**0.5)))*CbN*CBEm*PI*5*5)**2
                        dum_boron_tumor_err += (float(dum[2])*float(dum[3])*math.exp(bkerma(math.log((float(dum[0])*float(dum[1]))**0.5)))*CbN*TNR*CBEt*PI*5*5)**2
                    DPF_boron_normal.append(dum_boron_normal)
                    DPF_boron_skin.append(dum_boron_skin)
                    DPF_boron_mucosa.append(dum_boron_mucosa)
                    DPF_boron_tumor.append(dum_boron_tumor)
                    DPF_boron_normal_err.append(dum_boron_normal_err**0.5)
                    DPF_boron_skin_err.append(dum_boron_skin_err**0.5)
                    DPF_boron_mucosa_err.append(dum_boron_mucosa_err**0.5)
                    DPF_boron_tumor_err.append(dum_boron_tumor_err**0.5)
                    depth += 1
        out_path = os.path.join(dpf_path, r"boron_normal.dpf")
        with open(out_path, mode = 'a')  as f:
            for depth in range(TALLIES):
                f.write('%10.4e\t' %(DPF_boron_normal[depth]))
            f.write('\n')
        out_path = os.path.join(dpf_path, r"boron_skin.dpf")
        with open(out_path, mode = 'a')  as f:
            for depth in range(TALLIES):
                f.write('%10.4e\t' %(DPF_boron_skin[depth]))
            f.write('\n')
        out_path = os.path.join(dpf_path, r"boron_mucosa.dpf")
        with open(out_path, mode = 'a')  as f:
            for depth in range(TALLIES):
                f.write('%10.4e\t' %(DPF_boron_mucosa[depth]))
            f.write('\n')
        out_path = os.path.join(dpf_path, r"boron_tumor.dpf")
        with open(out_path, mode = 'a')  as f:
            for depth in range(TALLIES):
                f.write('%10.4e\t' %(DPF_boron_tumor[depth]))
            f.write('\n')
        out_path = os.path.join(dpf_path, r"boron_normal_err.dpf")
        with open(out_path, mode = 'a')  as f:
            for depth in range(TALLIES):
                f.write('%10.4e\t' %(DPF_boron_normal_err[depth]))
            f.write('\n')
        out_path = os.path.join(dpf_path, r"boron_skin_err.dpf")
        with open(out_path, mode = 'a')  as f:
            for depth in range(TALLIES):
                f.write('%10.4e\t' %(DPF_boron_skin_err[depth]))
            f.write('\n')
        out_path = os.path.join(dpf_path, r"boron_mucosa_err.dpf")
        with open(out_path, mode = 'a')  as f:
            for depth in range(TALLIES):
                f.write('%10.4e\t' %(DPF_boron_mucosa_err[depth]))
            f.write('\n')
        out_path = os.path.join(dpf_path, r"boron_tumor_err.dpf")
        with open(out_path, mode = 'a')  as f:
            for depth in range(TALLIES):
                f.write('%10.4e\t' %(DPF_boron_tumor_err[depth]))
            f.write('\n')
    print(' spe2dpf (boron DPFs) successfully finished!')
    statusBar["text"] = 'Finished computing boron DPFs! (Processing time: %.4f min)' %((time.time()-start_time)/60)
# spe2dpf for primary gamma-ray DPF
def spe2dpf4gammaPriDPFFunc(inBox_workDir,
                             inBox_tabN, inBox_tabGpri, inBox_tabGsec,
                             inBox_engListN, inBox_engListG,
                             inBox_gkerma, inBox_nkerma, inBox_bkerma,
                             inBox_CbN, inBox_TNR,
                             inBox_CBEt, inBox_CBEm, inBox_CBEs, inBox_CBEn,
                             inBox_RBEgamma, inBox_RBEthermal, inBox_RBEfast,
                             inBox_BINS_MONO_N, inBox_BINS_MONO_G,
                             inBox_BINS_SPE_N, inBox_BINS_SPE_G_PRI, inBox_BINS_SPE_G_SEC,
                             inBox_EMIN_SPE_N, inBox_EMIN_SPE_G_PRI, inBox_EMIN_SPE_G_SEC,
                             inBox_EMAX_SPE_N, inBox_EMAX_SPE_G_PRI, inBox_EMAX_SPE_G_SEC,
                             inBox_TALLIES,
                             inBox_gPriDPFName):
    # check PATH
    if inBox_workDir.get() != '' and os.path.isdir(inBox_workDir.get()) == True:
        pass
    else:
        print(" ERROR: Working directory does not exist!")
        statusBar["text"] = 'ERROR: Working directory does not exist!'
        return
    if inBox_tabGpri.get() != '' and os.path.isdir(inBox_tabGpri.get()) == True:
        pass
    else:
        print(" ERROR: Directory of primary gamma-ray table data does not exist!")
        statusBar["text"] = 'ERROR: Directory of primary gamma-ray table data does not exist!'
        return
    if inBox_gkerma.get() != '' and os.path.isfile(os.path.join(r"kerma", inBox_gkerma.get())) == True:
        pass
    else:
        print(" ERROR: Gamma-ray kerma file does not exist!")
        statusBar["text"] = 'ERROR: Gamma-ray kerma file does not exist!'
        return
    if (inBox_BINS_MONO_N.get() != ''
        and inBox_BINS_MONO_G.get() != ''
        and inBox_BINS_SPE_N.get() != ''
        and inBox_BINS_SPE_G_PRI.get() != ''
        and inBox_BINS_SPE_G_SEC.get() != ''
        and inBox_EMIN_SPE_N.get() != ''
        and inBox_EMIN_SPE_G_PRI.get() != ''
        and inBox_EMIN_SPE_G_SEC.get() != ''
        and inBox_EMAX_SPE_N.get() != ''
        and inBox_EMAX_SPE_G_PRI.get() != ''
        and inBox_EMAX_SPE_G_SEC.get() != ''
        and inBox_TALLIES.get() != ''):
        pass
    else:
        print(" ERROR: Some or all of parameters in table data are missing!")
        statusBar["text"] = 'ERROR: Some or all of parameters in table data are missing!'
        return
    if (inBox_CbN.get() != ''
        and inBox_TNR.get() != ''
        and inBox_CBEt.get() != ''
        and inBox_CBEm.get() != ''
        and inBox_CBEs.get() != ''
        and inBox_CBEn.get() != ''
        and inBox_RBEgamma.get() != ''
        and inBox_RBEthermal.get() != ''
        and inBox_RBEfast.get() != ''):
        pass
    else:
        print(" ERROR: Some or all of radiobiological coefficients are missing!")
        statusBar["text"] = 'ERROR: Some or all of radiobiological coefficients are missing!'
        return
    if inBox_engListG.get() != '' and os.path.isfile(os.path.join(r"phits", inBox_engListG.get())):
        pass
    else:
        print(" ERROR: Monoenergy list of gamma-ray does not exist!")
        statusBar["text"] = 'ERROR: Monoenergy list of gamma-ray does not exist!'
        return
    if inBox_gPriDPFName.get() != '':
        pass
    else:
        print(" ERROR: Directory name to generate DPF file is missing!")
        statusBar["text"] = 'ERROR: Directory name to generate DPF file is missing!'
        return
    # main
    start_time = time.time()
    print(' #####################################################')
    print(' ########## spe2dpf (primary gamma-ray DPF) ##########')
    print(' #####################################################')
    print(' Computing & outputting...')
    # generate DPFs storaging directory
    dpf_path = os.path.join(inBox_workDir.get(), inBox_gPriDPFName.get())
    os.makedirs(dpf_path, exist_ok=True)
    # constants
    PI = PI = 3.141592653589793
    # read parameters
    CbN = float(inBox_CbN.get())
    TNR = float(inBox_TNR.get())
    CBEt = float(inBox_CBEt.get())
    CBEm = float(inBox_CBEm.get())
    CBEs = float(inBox_CBEs.get())
    CBEn = float(inBox_CBEn.get())
    RBEgamma = float(inBox_RBEgamma.get())
    RBEthermal = float(inBox_RBEthermal.get())
    RBEfast = float(inBox_RBEfast.get())
    BINS_MONO_N = int(inBox_BINS_MONO_N.get())
    BINS_MONO_G = int(inBox_BINS_MONO_G.get())
    BINS_SPE_N = int(inBox_BINS_SPE_N.get())
    BINS_SPE_G_PRI = int(inBox_BINS_SPE_G_PRI.get())
    BINS_SPE_G_SEC = int(inBox_BINS_SPE_G_SEC.get())
    EMIN_SPE_N = float(inBox_EMIN_SPE_N.get())
    EMIN_SPE_G_PRI = float(inBox_EMIN_SPE_G_PRI.get())
    EMIN_SPE_G_SEC = float(inBox_EMIN_SPE_G_SEC.get())
    EMAX_SPE_N = float(inBox_EMAX_SPE_N.get())
    EMAX_SPE_G_PRI = float(inBox_EMAX_SPE_G_PRI.get())
    EMAX_SPE_G_SEC = float(inBox_EMAX_SPE_G_SEC.get())
    TALLIES = int(inBox_TALLIES.get())
    # loading values of monoenergy (gamma-ray)
    LIST_E_mono_g = []
    LIST_E_mono_g_flt = []
    data_path = os.path.join(r"phits", inBox_engListG.get())
    with open(data_path, mode = 'r') as f:
        for i in range(BINS_MONO_G):
            dum = f.readline()
            dum = dum.replace("\n", "")
            LIST_E_mono_g.append(dum)
            LIST_E_mono_g_flt.append(float(dum))
    # loading and inter/extrapolating gamma-ray kerma coefficients
    LIST_eng_gkerma = []
    LIST_val_gkerma = []
    LIST_gkerma = []
    with open(os.path.join(r"kerma", inBox_gkerma.get()), mode = 'r') as f:
        dum = f.readline() # skipping header part
        for i in range(36):
            dum = f.readline()
            dum = dum.replace("\n", "")
            dum = dum.split()
            LIST_eng_gkerma.append(math.log(float(dum[0])))
            LIST_val_gkerma.append(math.log(float(dum[4])))
        gkerma = interpolate.interp1d(LIST_eng_gkerma, LIST_val_gkerma, kind='linear', fill_value='extrapolate') # usage: math.exp(kerma(math.log(ENG)))
    plt.plot(np.exp(np.array(LIST_eng_gkerma)), np.exp(np.array(LIST_val_gkerma)), linestyle="solid", color="green", label="ICRU soft tissue")
    plt.legend()
    plt.title('gamma-ray kerma coeff. ')
    plt.xlabel('gamma-ray energy [MeV]')
    plt.ylabel('kerma coeff. [pGy $\mathregular{cm^2}$]')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(color='r', linestyle='dotted', linewidth=1)
    plt.savefig(inBox_workDir.get() + r"\gkerma.png")
    plt.close()
    #generating primary gamma-ray DPF
    data_path = inBox_tabGpri.get()
    out_path = os.path.join(dpf_path, r"g_pri.dpf")
    with open(out_path, mode = 'w')  as f:
        f.write('#constants: borofalan concentartion in normal tissue %.1f, tumor/normal tissue ratio of borofalan concentration %.1f, RBE for gamma-ray %.1f, RBE for nitrogen dose %.1f, RBE for hydrogen dose and other neutron dose %.1f, CBE for tumor %.1f, CBE for mucosa %.1f, CBE for skin %.1f, CBE for other normal tissue %.2f\n' %(CbN, TNR, RBEgamma, RBEthermal, RBEfast, CBEt, CBEm, CBEs, CBEn))
    out_path = os.path.join(dpf_path, r"g_pri_err.dpf")
    with open(out_path, mode = 'w')  as f:
        f.write('#constants: borofalan concentartion in normal tissue %.1f, tumor/normal tissue ratio of borofalan concentration %.1f, RBE for gamma-ray %.1f, RBE for nitrogen dose %.1f, RBE for hydrogen dose and other neutron dose %.1f, CBE for tumor %.1f, CBE for mucosa %.1f, CBE for skin %.1f, CBE for other normal tissue %.2f\n' %(CbN, TNR, RBEgamma, RBEthermal, RBEfast, CBEt, CBEm, CBEs, CBEn))
    for i in range(BINS_MONO_G):
        data_name = LIST_E_mono_g[i]+'MeV-gspe_pri.dat'
        depth = 0
        with open(os.path.join(data_path, data_name), mode = 'r') as f:
            DPF_g_pri = []
            DPF_g_pri_err = []
            while True:
                if depth >= TALLIES:
                    break
                dum = f.readline()
                dum = dum.replace("\n", "")
                dum = dum.split()
                if len(dum) < 4:
                    continue
                if dum[1] == 'e-lower':
                    dum_g_pri = 0
                    dum_g_pri_err = 0
                    for k in range(BINS_SPE_G_PRI):
                        dum = f.readline()
                        dum = dum.replace("\n", "")
                        dum = dum.split()
                        dum_g_pri += float(dum[2])*math.exp(gkerma(math.log((float(dum[0])*float(dum[1]))**0.5)))*RBEgamma*PI*5*5
                        dum_g_pri_err += (float(dum[2])*float(dum[3])*math.exp(gkerma(math.log((float(dum[0])*float(dum[1]))**0.5)))*RBEgamma*PI*5*5)**2
                    DPF_g_pri.append(dum_g_pri)
                    DPF_g_pri_err.append(dum_g_pri_err**0.5)
                    depth += 1
        out_path = os.path.join(dpf_path, r"g_pri.dpf")
        with open(out_path, mode = 'a')  as f:
            for depth in range(TALLIES):
                f.write('%10.4e\t' %(DPF_g_pri[depth]))
            f.write('\n')
        out_path = os.path.join(dpf_path, r"g_pri_err.dpf")
        with open(out_path, mode = 'a')  as f:
            for depth in range(TALLIES):
                f.write('%10.4e\t' %(DPF_g_pri_err[depth]))
            f.write('\n')
    print(' spe2dpf (primary gamma-ray DPF) successfully finished!')
    statusBar["text"] = 'Finished computing primary gamma-ray DPF! (Processing time: %.4f min)' %((time.time()-start_time)/60)
# spe2dpf for secondary gamma-ray DPF
def spe2dpf4gammaSecDPFFunc(inBox_workDir,
                             inBox_tabN, inBox_tabGpri, inBox_tabGsec,
                             inBox_engListN, inBox_engListG,
                             inBox_gkerma, inBox_nkerma, inBox_bkerma,
                             inBox_CbN, inBox_TNR,
                             inBox_CBEt, inBox_CBEm, inBox_CBEs, inBox_CBEn,
                             inBox_RBEgamma, inBox_RBEthermal, inBox_RBEfast,
                             inBox_BINS_MONO_N, inBox_BINS_MONO_G,
                             inBox_BINS_SPE_N, inBox_BINS_SPE_G_PRI, inBox_BINS_SPE_G_SEC,
                             inBox_EMIN_SPE_N, inBox_EMIN_SPE_G_PRI, inBox_EMIN_SPE_G_SEC,
                             inBox_EMAX_SPE_N, inBox_EMAX_SPE_G_PRI, inBox_EMAX_SPE_G_SEC,
                             inBox_TALLIES,
                             inBox_gSecDPFName):
    # check PATH
    if inBox_workDir.get() != '' and os.path.isdir(inBox_workDir.get()) == True:
        pass
    else:
        print(" ERROR: Working directory does not exist!")
        statusBar["text"] = 'ERROR: Working directory does not exist!'
        return
    if inBox_tabGsec.get() != '' and os.path.isdir(inBox_tabGsec.get()) == True:
        pass
    else:
        print(" ERROR: Directory of secondary gamma-ray table data does not exist!")
        statusBar["text"] = 'ERROR: Directory of secondary gamma-ray table data does not exist!'
        return
    if inBox_gkerma.get() != '' and os.path.isfile(os.path.join(r"kerma", inBox_gkerma.get())) == True:
        pass
    else:
        print(" ERROR: Gamma-ray kerma file does not exist!")
        statusBar["text"] = 'ERROR: Gamma-ray kerma file does not exist!'
        return
    if (inBox_BINS_MONO_N.get() != ''
        and inBox_BINS_MONO_G.get() != ''
        and inBox_BINS_SPE_N.get() != ''
        and inBox_BINS_SPE_G_PRI.get() != ''
        and inBox_BINS_SPE_G_SEC.get() != ''
        and inBox_EMIN_SPE_N.get() != ''
        and inBox_EMIN_SPE_G_PRI.get() != ''
        and inBox_EMIN_SPE_G_SEC.get() != ''
        and inBox_EMAX_SPE_N.get() != ''
        and inBox_EMAX_SPE_G_PRI.get() != ''
        and inBox_EMAX_SPE_G_SEC.get() != ''
        and inBox_TALLIES.get() != ''):
        pass
    else:
        print(" ERROR: Some or all of parameters in table data are missing!")
        statusBar["text"] = 'ERROR: Some or all of parameters in table data are missing!'
        return
    if (inBox_CbN.get() != ''
        and inBox_TNR.get() != ''
        and inBox_CBEt.get() != ''
        and inBox_CBEm.get() != ''
        and inBox_CBEs.get() != ''
        and inBox_CBEn.get() != ''
        and inBox_RBEgamma.get() != ''
        and inBox_RBEthermal.get() != ''
        and inBox_RBEfast.get() != ''):
        pass
    else:
        print(" ERROR: Some or all of radiobiological coefficients are missing!")
        statusBar["text"] = 'ERROR: Some or all of radiobiological coefficients are missing!'
        return
    if inBox_engListN.get() != '' and os.path.isfile(os.path.join(r"phits", inBox_engListN.get())):
        pass
    else:
        print(" ERROR: Monoenergy list of neutron does not exist!")
        statusBar["text"] = 'ERROR: Monoenergy list of neutron does not exist!'
        return
    if inBox_gSecDPFName.get() != '':
        pass
    else:
        print(" ERROR: Directory name to generate DPF file is missing!")
        statusBar["text"] = 'ERROR: Directory name to generate DPF file is missing!'
        return
    # main
    start_time = time.time()
    print(' #######################################################')
    print(' ########## spe2dpf (secondary gamma-ray DPF) ##########')
    print(' #######################################################')
    print(' Computing & outputting...')
    # generate DPFs storaging directory
    dpf_path = os.path.join(inBox_workDir.get(), inBox_gSecDPFName.get())
    os.makedirs(dpf_path, exist_ok=True)
    # constants
    PI = PI = 3.141592653589793
    # read parameters
    CbN = float(inBox_CbN.get())
    TNR = float(inBox_TNR.get())
    CBEt = float(inBox_CBEt.get())
    CBEm = float(inBox_CBEm.get())
    CBEs = float(inBox_CBEs.get())
    CBEn = float(inBox_CBEn.get())
    RBEgamma = float(inBox_RBEgamma.get())
    RBEthermal = float(inBox_RBEthermal.get())
    RBEfast = float(inBox_RBEfast.get())
    BINS_MONO_N = int(inBox_BINS_MONO_N.get())
    BINS_MONO_G = int(inBox_BINS_MONO_G.get())
    BINS_SPE_N = int(inBox_BINS_SPE_N.get())
    BINS_SPE_G_PRI = int(inBox_BINS_SPE_G_PRI.get())
    BINS_SPE_G_SEC = int(inBox_BINS_SPE_G_SEC.get())
    EMIN_SPE_N = float(inBox_EMIN_SPE_N.get())
    EMIN_SPE_G_PRI = float(inBox_EMIN_SPE_G_PRI.get())
    EMIN_SPE_G_SEC = float(inBox_EMIN_SPE_G_SEC.get())
    EMAX_SPE_N = float(inBox_EMAX_SPE_N.get())
    EMAX_SPE_G_PRI = float(inBox_EMAX_SPE_G_PRI.get())
    EMAX_SPE_G_SEC = float(inBox_EMAX_SPE_G_SEC.get())
    TALLIES = int(inBox_TALLIES.get())
    # loading values of monoenergy (neutron)
    LIST_E_mono_n = []
    LIST_E_mono_n_flt = []
    data_path = os.path.join(r"phits", inBox_engListN.get())
    with open(data_path, mode = 'r') as f:
        for i in range(BINS_MONO_N):
            dum = f.readline()
            dum = dum.replace("\n", "")
            LIST_E_mono_n.append(dum)
            LIST_E_mono_n_flt.append(float(dum))
    # loading and inter/extrapolating gamma-ray kerma coefficients
    LIST_eng_gkerma = []
    LIST_val_gkerma = []
    LIST_gkerma = []
    with open(os.path.join(r"kerma", inBox_gkerma.get()), mode = 'r') as f:
        dum = f.readline() # skipping header part
        for i in range(36):
            dum = f.readline()
            dum = dum.replace("\n", "")
            dum = dum.split()
            LIST_eng_gkerma.append(math.log(float(dum[0])))
            LIST_val_gkerma.append(math.log(float(dum[4])))
        gkerma = interpolate.interp1d(LIST_eng_gkerma, LIST_val_gkerma, kind='linear', fill_value='extrapolate') # usage: math.exp(kerma(math.log(ENG)))
    plt.plot(np.exp(np.array(LIST_eng_gkerma)), np.exp(np.array(LIST_val_gkerma)), linestyle="solid", color="green", label="ICRU soft tissue")
    plt.legend()
    plt.title('gamma-ray kerma coeff. ')
    plt.xlabel('gamma-ray energy [MeV]')
    plt.ylabel('kerma coeff. [pGy $\mathregular{cm^2}$]')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(color='r', linestyle='dotted', linewidth=1)
    plt.savefig(inBox_workDir.get() + r"\gkerma.png")
    plt.close()
    #generating secondary gamma-ray DPF
    data_path = inBox_tabGsec.get()
    out_path = os.path.join(dpf_path, r"g_sec.dpf")
    with open(out_path, mode = 'w')  as f:
        f.write('#constants: borofalan concentartion in normal tissue %.1f, tumor/normal tissue ratio of borofalan concentration %.1f, RBE for gamma-ray %.1f, RBE for nitrogen dose %.1f, RBE for hydrogen dose and other neutron dose %.1f, CBE for tumor %.1f, CBE for mucosa %.1f, CBE for skin %.1f, CBE for other normal tissue %.2f\n' %(CbN, TNR, RBEgamma, RBEthermal, RBEfast, CBEt, CBEm, CBEs, CBEn))
    out_path = os.path.join(dpf_path, r"g_sec_err.dpf")
    with open(out_path, mode = 'w')  as f:
        f.write('#constants: borofalan concentartion in normal tissue %.1f, tumor/normal tissue ratio of borofalan concentration %.1f, RBE for gamma-ray %.1f, RBE for nitrogen dose %.1f, RBE for hydrogen dose and other neutron dose %.1f, CBE for tumor %.1f, CBE for mucosa %.1f, CBE for skin %.1f, CBE for other normal tissue %.2f\n' %(CbN, TNR, RBEgamma, RBEthermal, RBEfast, CBEt, CBEm, CBEs, CBEn))
    for i in range(BINS_MONO_N):
        data_name = LIST_E_mono_n[i]+'MeV-gspe_sec.dat'
        depth = 0
        with open(os.path.join(data_path, data_name), mode = 'r') as f:
            DPF_g_sec = []
            DPF_g_sec_err = []
            while True:
                if depth >= TALLIES:
                    break
                dum = f.readline()
                dum = dum.replace("\n", "")
                dum = dum.split()
                if len(dum) < 4:
                    continue
                if dum[1] == 'e-lower':
                    dum_g_sec = 0
                    dum_g_sec_err = 0
                    for k in range(BINS_SPE_G_SEC):
                        dum = f.readline()
                        dum = dum.replace("\n", "")
                        dum = dum.split()
                        dum_g_sec += float(dum[2])*math.exp(gkerma(math.log((float(dum[0])*float(dum[1]))**0.5)))*RBEgamma*PI*5*5
                        dum_g_sec_err += (float(dum[2])*float(dum[3])*math.exp(gkerma(math.log((float(dum[0])*float(dum[1]))**0.5)))*RBEgamma*PI*5*5)**2
                    DPF_g_sec.append(dum_g_sec)
                    DPF_g_sec_err.append(dum_g_sec_err**0.5)
                    depth += 1
        out_path = os.path.join(dpf_path, r"g_sec.dpf")
        with open(out_path, mode = 'a')  as f:
            for depth in range(TALLIES):
                f.write('%10.4e\t' %(DPF_g_sec[depth]))
            f.write('\n')
        out_path = os.path.join(dpf_path, r"g_sec_err.dpf")
        with open(out_path, mode = 'a')  as f:
            for depth in range(TALLIES):
                f.write('%10.4e\t' %(DPF_g_sec_err[depth]))
            f.write('\n')
    print(' spe2dpf (secondary gamma-ray DPF) successfully finished!')
    statusBar["text"] = 'Finished computing secondary gamma-ray DPF! (Processing time: %.4f min)' %((time.time()-start_time)/60)
# dpf2dist
def dpf2distFunc(inBox_workDir,
                 inBox_engListN, inBox_engListG,
                 inBox_BINS_MONO_N, inBox_BINS_MONO_G,
                 inBox_TALLIES,
                 inBox_DPFDir,
                 inBox_nspeDir, inBox_nspeFile,
                 inBox_gspeDir, inBox_gspeFile,
                 inBox_distName):
    # check PATH
    if inBox_workDir.get() != '' and os.path.isdir(inBox_workDir.get()) == True:
        pass
    else:
        print(" ERROR: Working directory does not exist!")
        statusBar["text"] = 'ERROR: Working directory does not exist!'
        return
    if inBox_nspeDir.get() != '' and os.path.isfile(os.path.join(inBox_nspeDir.get(), inBox_nspeFile.get())):
        pass
    else:
        print(" ERROR: Input neutron spectrum file does not exist!")
        statusBar["text"] = 'ERROR: Input neutron spectrum file does not exist!'
        return
    if inBox_gspeDir.get() != '' and os.path.isfile(os.path.join(inBox_gspeDir.get(), inBox_gspeFile.get())):
        pass
    else:
        print(" ERROR: Input gamma-ray spectrum file does not exist!")
        statusBar["text"] = 'ERROR: Input gamma-ray spectrum file does not exist!'
        return
    if inBox_DPFDir.get() != '' and os.path.isdir(inBox_DPFDir.get()) == True:
        pass
    else:
        print(" ERROR: DPF directory does not exist!")
        statusBar["text"] = 'ERROR: DPF directory does not exist!'
        return
    if (os.path.exists(os.path.join(inBox_DPFDir.get(), 'hyd.dpf')) == True
        and os.path.exists(os.path.join(inBox_DPFDir.get(), 'nit.dpf')) == True
        and os.path.exists(os.path.join(inBox_DPFDir.get(), 'oth.dpf')) == True
        and os.path.exists(os.path.join(inBox_DPFDir.get(), 'g_pri.dpf')) == True
        and os.path.exists(os.path.join(inBox_DPFDir.get(), 'g_sec.dpf')) == True
        and os.path.exists(os.path.join(inBox_DPFDir.get(), 'boron_normal.dpf')) == True
        and os.path.exists(os.path.join(inBox_DPFDir.get(), 'boron_skin.dpf')) == True
        and os.path.exists(os.path.join(inBox_DPFDir.get(), 'boron_mucosa.dpf')) == True
        and os.path.exists(os.path.join(inBox_DPFDir.get(), 'boron_tumor.dpf')) == True
        and os.path.exists(os.path.join(inBox_DPFDir.get(), 'hyd_err.dpf')) == True
        and os.path.exists(os.path.join(inBox_DPFDir.get(), 'nit_err.dpf')) == True
        and os.path.exists(os.path.join(inBox_DPFDir.get(), 'oth_err.dpf')) == True
        and os.path.exists(os.path.join(inBox_DPFDir.get(), 'g_pri_err.dpf')) == True
        and os.path.exists(os.path.join(inBox_DPFDir.get(), 'g_sec_err.dpf')) == True
        and os.path.exists(os.path.join(inBox_DPFDir.get(), 'boron_normal_err.dpf')) == True
        and os.path.exists(os.path.join(inBox_DPFDir.get(), 'boron_skin_err.dpf')) == True
        and os.path.exists(os.path.join(inBox_DPFDir.get(), 'boron_mucosa_err.dpf')) == True
        and os.path.exists(os.path.join(inBox_DPFDir.get(), 'boron_tumor_err.dpf')) == True):
        pass
    else:
        print(" ERROR: Some or all of DPF files does not exist!")
        statusBar["text"] = 'ERROR: Some or all of DPF files does not exist!'
        return
    if inBox_distName.get() != '':
        pass
    else:
        print(" ERROR: Output file name is missing!")
        statusBar["text"] = 'ERROR: Output file name is missing!'
        return
    # main
    start_time = time.time()
    print(' ##############################')
    print(' ########## dpf2dist ##########')
    print(' ##############################')
    print(' Computing...')
    # read parameters
    BINS_MONO_N = int(inBox_BINS_MONO_N.get())
    BINS_MONO_G = int(inBox_BINS_MONO_G.get())
    TALLIES = int(inBox_TALLIES.get())
    # loading values of monoenergy (neutron)
    LIST_E_mono_n = []
    LIST_E_mono_n_flt = []
    log_LIST_E_mono_n_flt = []
    with open(os.path.join(r"phits", inBox_engListN.get()), mode = 'r') as f:
        for i in range(BINS_MONO_N):
            dum = f.readline()
            dum = dum.replace("\n", "")
            LIST_E_mono_n.append(dum)
            LIST_E_mono_n_flt.append(float(dum))
            log_LIST_E_mono_n_flt.append(math.log(float(dum)))
    # loading values of monoenergy (gamma-ray)
    LIST_E_mono_g = []
    LIST_E_mono_g_flt = []
    with open(os.path.join(r"phits", inBox_engListG.get()), mode = 'r') as f:
        for i in range(BINS_MONO_G):
            dum = f.readline()
            dum = dum.replace("\n", "")
            LIST_E_mono_g.append(dum)
            LIST_E_mono_g_flt.append(float(dum))
    # loading DPFs and DPF errors
    # hydrogen dose DPF
    LIST_DPF_hyd = np.empty(0)
    dpf_path = os.path.join(inBox_DPFDir.get(), r"hyd.dpf")
    with open(dpf_path, mode = 'r') as f:
        dum = f.readline()
        for i in range(BINS_MONO_N):
            dum = f.readline()
            dum = dum.replace("\n", "")
            dum = dum.split()
            dum = np.array(dum)
            LIST_DPF_hyd = np.append(LIST_DPF_hyd, dum)
    LIST_DPF_hyd = LIST_DPF_hyd.reshape(BINS_MONO_N, TALLIES)
    LIST_DPF_hyd = LIST_DPF_hyd.T # LIST_DPF[i][j], where i is monoenergy and j is depth
    # hydrogen dose DPF error
    LIST_DPF_hyd_err = np.empty(0)
    dpf_path = os.path.join(inBox_DPFDir.get(), r"hyd_err.dpf")
    with open(dpf_path, mode = 'r') as f:
        dum = f.readline()
        for i in range(BINS_MONO_N):
            dum = f.readline()
            dum = dum.replace("\n", "")
            dum = dum.split()
            dum = np.array(dum)
            LIST_DPF_hyd_err = np.append(LIST_DPF_hyd_err, dum)
    LIST_DPF_hyd_err = LIST_DPF_hyd_err.reshape(BINS_MONO_N, TALLIES)
    LIST_DPF_hyd_err = LIST_DPF_hyd_err.T # LIST_DPF[i][j], where i is monoenergy and j is depth
    # nitrogen dose DPF
    LIST_DPF_nit = np.empty(0)
    dpf_path = os.path.join(inBox_DPFDir.get(), r"nit.dpf")
    with open(dpf_path, mode = 'r') as f:
        dum = f.readline()
        for i in range(BINS_MONO_N):
            dum = f.readline()
            dum = dum.replace("\n", "")
            dum = dum.split()
            dum = np.array(dum)
            LIST_DPF_nit = np.append(LIST_DPF_nit, dum)
    LIST_DPF_nit = LIST_DPF_nit.reshape(BINS_MONO_N, TALLIES)
    LIST_DPF_nit = LIST_DPF_nit.T # LIST_DPF[i][j], where i is monoenergy and j is depth
    # nitrogen dose DPF error
    LIST_DPF_nit_err = np.empty(0)
    dpf_path = os.path.join(inBox_DPFDir.get(), r"nit_err.dpf")
    with open(dpf_path, mode = 'r') as f:
        dum = f.readline()
        for i in range(BINS_MONO_N):
            dum = f.readline()
            dum = dum.replace("\n", "")
            dum = dum.split()
            dum = np.array(dum)
            LIST_DPF_nit_err = np.append(LIST_DPF_nit_err, dum)
    LIST_DPF_nit_err = LIST_DPF_nit_err.reshape(BINS_MONO_N, TALLIES)
    LIST_DPF_nit_err = LIST_DPF_nit_err.T # LIST_DPF[i][j], where i is monoenergy and j is depth
    # other neutron dose DPF
    LIST_DPF_oth = np.empty(0)
    dpf_path = os.path.join(inBox_DPFDir.get(), r"oth.dpf")
    with open(dpf_path, mode = 'r') as f:
        dum = f.readline()
        for i in range(BINS_MONO_N):
            dum = f.readline()
            dum = dum.replace("\n", "")
            dum = dum.split()
            dum = np.array(dum)
            LIST_DPF_oth = np.append(LIST_DPF_oth, dum)
    LIST_DPF_oth = LIST_DPF_oth.reshape(BINS_MONO_N, TALLIES)
    LIST_DPF_oth = LIST_DPF_oth.T # LIST_DPF[i][j], where i is monoenergy and j is depth
    # other neutron dose DPF error
    LIST_DPF_oth_err = np.empty(0)
    dpf_path = os.path.join(inBox_DPFDir.get(), r"oth_err.dpf")
    with open(dpf_path, mode = 'r') as f:
        dum = f.readline()
        for i in range(BINS_MONO_N):
            dum = f.readline()
            dum = dum.replace("\n", "")
            dum = dum.split()
            dum = np.array(dum)
            LIST_DPF_oth_err = np.append(LIST_DPF_oth_err, dum)
    LIST_DPF_oth_err = LIST_DPF_oth_err.reshape(BINS_MONO_N, TALLIES)
    LIST_DPF_oth_err = LIST_DPF_oth_err.T # LIST_DPF[i][j], where i is monoenergy and j is depth
    # primary gamma-ray dose DPF
    LIST_DPF_g_pri = np.empty(0)
    dpf_path = os.path.join(inBox_DPFDir.get(), r"g_pri.dpf")
    with open(dpf_path, mode = 'r') as f:
        dum = f.readline()
        for i in range(BINS_MONO_G):
            dum = f.readline()
            dum = dum.replace("\n", "")
            dum = dum.split()
            dum = np.array(dum)
            LIST_DPF_g_pri = np.append(LIST_DPF_g_pri, dum)
    LIST_DPF_g_pri = LIST_DPF_g_pri.reshape(BINS_MONO_G, TALLIES)
    LIST_DPF_g_pri = LIST_DPF_g_pri.T # LIST_DPF[i][j], where i is monoenergy and j is depth
    # primary gamma-ray dose DPF error
    LIST_DPF_g_pri_err = np.empty(0)
    dpf_path = os.path.join(inBox_DPFDir.get(), r"g_pri_err.dpf")
    with open(dpf_path, mode = 'r') as f:
        dum = f.readline()
        for i in range(BINS_MONO_G):
            dum = f.readline()
            dum = dum.replace("\n", "")
            dum = dum.split()
            dum = np.array(dum)
            LIST_DPF_g_pri_err = np.append(LIST_DPF_g_pri_err, dum)
    LIST_DPF_g_pri_err = LIST_DPF_g_pri_err.reshape(BINS_MONO_G, TALLIES)
    LIST_DPF_g_pri_err = LIST_DPF_g_pri_err.T # LIST_DPF[i][j], where i is monoenergy and j is depth
    # secondary gamma-ray dose DPF
    LIST_DPF_g_sec = np.empty(0)
    dpf_path = os.path.join(inBox_DPFDir.get(), r"g_sec.dpf")
    with open(dpf_path, mode = 'r') as f:
        dum = f.readline()
        for i in range(BINS_MONO_N):
            dum = f.readline()
            dum = dum.replace("\n", "")
            dum = dum.split()
            dum = np.array(dum)
            LIST_DPF_g_sec = np.append(LIST_DPF_g_sec, dum)
    LIST_DPF_g_sec = LIST_DPF_g_sec.reshape(BINS_MONO_N, TALLIES)
    LIST_DPF_g_sec = LIST_DPF_g_sec.T # LIST_DPF[i][j], where i is monoenergy and j is depth
    # secondary gamma-ray dose DPF error
    LIST_DPF_g_sec_err = np.empty(0)
    dpf_path = os.path.join(inBox_DPFDir.get(), r"g_sec_err.dpf")
    with open(dpf_path, mode = 'r') as f:
        dum = f.readline()
        for i in range(BINS_MONO_N):
            dum = f.readline()
            dum = dum.replace("\n", "")
            dum = dum.split()
            dum = np.array(dum)
            LIST_DPF_g_sec_err = np.append(LIST_DPF_g_sec_err, dum)
    LIST_DPF_g_sec_err = LIST_DPF_g_sec_err.reshape(BINS_MONO_N, TALLIES)
    LIST_DPF_g_sec_err = LIST_DPF_g_sec_err.T # LIST_DPF[i][j], where i is monoenergy and j is depth
    # other normal tissue boron dose DPF
    LIST_DPF_boron_normal = np.empty(0)
    dpf_path = os.path.join(inBox_DPFDir.get(), r"boron_normal.dpf")
    with open(dpf_path, mode = 'r') as f:
        dum = f.readline()
        for i in range(BINS_MONO_N):
            dum = f.readline()
            dum = dum.replace("\n", "")
            dum = dum.split()
            dum = np.array(dum)
            LIST_DPF_boron_normal = np.append(LIST_DPF_boron_normal, dum)
    LIST_DPF_boron_normal = LIST_DPF_boron_normal.reshape(BINS_MONO_N, TALLIES)
    LIST_DPF_boron_normal = LIST_DPF_boron_normal.T # LIST_DPF[i][j], where i is monoenergy and j is depth
    # other normal tissue boron dose DPF error
    LIST_DPF_boron_normal_err = np.empty(0)
    dpf_path = os.path.join(inBox_DPFDir.get(), r"boron_normal_err.dpf")
    with open(dpf_path, mode = 'r') as f:
        dum = f.readline()
        for i in range(BINS_MONO_N):
            dum = f.readline()
            dum = dum.replace("\n", "")
            dum = dum.split()
            dum = np.array(dum)
            LIST_DPF_boron_normal_err = np.append(LIST_DPF_boron_normal_err, dum)
    LIST_DPF_boron_normal_err = LIST_DPF_boron_normal_err.reshape(BINS_MONO_N, TALLIES)
    LIST_DPF_boron_normal_err = LIST_DPF_boron_normal_err.T # LIST_DPF[i][j], where i is monoenergy and j is depth
    # skin boron dose DPF
    LIST_DPF_boron_skin = np.empty(0)
    dpf_path = os.path.join(inBox_DPFDir.get(), r"boron_skin.dpf")
    with open(dpf_path, mode = 'r') as f:
        dum = f.readline()
        for i in range(BINS_MONO_N):
            dum = f.readline()
            dum = dum.replace("\n", "")
            dum = dum.split()
            dum = np.array(dum)
            LIST_DPF_boron_skin = np.append(LIST_DPF_boron_skin, dum)
    LIST_DPF_boron_skin = LIST_DPF_boron_skin.reshape(BINS_MONO_N, TALLIES)
    LIST_DPF_boron_skin = LIST_DPF_boron_skin.T # LIST_DPF[i][j], where i is monoenergy and j is depth
    # skin boron dose DPF error
    LIST_DPF_boron_skin_err = np.empty(0)
    dpf_path = os.path.join(inBox_DPFDir.get(), r"boron_skin_err.dpf")
    with open(dpf_path, mode = 'r') as f:
        dum = f.readline()
        for i in range(BINS_MONO_N):
            dum = f.readline()
            dum = dum.replace("\n", "")
            dum = dum.split()
            dum = np.array(dum)
            LIST_DPF_boron_skin_err = np.append(LIST_DPF_boron_skin_err, dum)
    LIST_DPF_boron_skin_err = LIST_DPF_boron_skin_err.reshape(BINS_MONO_N, TALLIES)
    LIST_DPF_boron_skin_err = LIST_DPF_boron_skin_err.T # LIST_DPF[i][j], where i is monoenergy and j is depth
    # mucosa boron dose DPF
    LIST_DPF_boron_mucosa = np.empty(0)
    dpf_path = os.path.join(inBox_DPFDir.get(), r"boron_mucosa.dpf")
    with open(dpf_path, mode = 'r') as f:
        dum = f.readline()
        for i in range(BINS_MONO_N):
            dum = f.readline()
            dum = dum.replace("\n", "")
            dum = dum.split()
            dum = np.array(dum)
            LIST_DPF_boron_mucosa = np.append(LIST_DPF_boron_mucosa, dum)
    LIST_DPF_boron_mucosa = LIST_DPF_boron_mucosa.reshape(BINS_MONO_N, TALLIES)
    LIST_DPF_boron_mucosa = LIST_DPF_boron_mucosa.T # LIST_DPF[i][j], where i is monoenergy and j is depth
    # mucosa boron dose DPF error
    LIST_DPF_boron_mucosa_err = np.empty(0)
    dpf_path = os.path.join(inBox_DPFDir.get(), r"boron_mucosa_err.dpf")
    with open(dpf_path, mode = 'r') as f:
        dum = f.readline()
        for i in range(BINS_MONO_N):
            dum = f.readline()
            dum = dum.replace("\n", "")
            dum = dum.split()
            dum = np.array(dum)
            LIST_DPF_boron_mucosa_err = np.append(LIST_DPF_boron_mucosa_err, dum)
    LIST_DPF_boron_mucosa_err = LIST_DPF_boron_mucosa_err.reshape(BINS_MONO_N, TALLIES)
    LIST_DPF_boron_mucosa_err = LIST_DPF_boron_mucosa_err.T # LIST_DPF[i][j], where i is monoenergy and j is depth
    # tumor boron dose DPF
    LIST_DPF_boron_tumor = np.empty(0)
    dpf_path = os.path.join(inBox_DPFDir.get(), r"boron_tumor.dpf")
    with open(dpf_path, mode = 'r') as f:
        dum = f.readline()
        for i in range(BINS_MONO_N):
            dum = f.readline()
            dum = dum.replace("\n", "")
            dum = dum.split()
            dum = np.array(dum)
            LIST_DPF_boron_tumor = np.append(LIST_DPF_boron_tumor, dum)
    LIST_DPF_boron_tumor = LIST_DPF_boron_tumor.reshape(BINS_MONO_N, TALLIES)
    LIST_DPF_boron_tumor = LIST_DPF_boron_tumor.T # LIST_DPF[i][j], where i is monoenergy and j is depth
    # tumor boron dose DPF error
    LIST_DPF_boron_tumor_err = np.empty(0)
    dpf_path = os.path.join(inBox_DPFDir.get(), r"boron_tumor_err.dpf")
    with open(dpf_path, mode = 'r') as f:
        dum = f.readline()
        for i in range(BINS_MONO_N):
            dum = f.readline()
            dum = dum.replace("\n", "")
            dum = dum.split()
            dum = np.array(dum)
            LIST_DPF_boron_tumor_err = np.append(LIST_DPF_boron_tumor_err, dum)
    LIST_DPF_boron_tumor_err = LIST_DPF_boron_tumor_err.reshape(BINS_MONO_N, TALLIES)
    LIST_DPF_boron_tumor_err = LIST_DPF_boron_tumor_err.T # LIST_DPF[i][j], where i is monoenergy and j is depth
    print(' Loaded DPFs!')
    # calculating neutron, secondary gamma-ray, and boron dose rates
    LIST_doserate_hyd = np.zeros(TALLIES)
    LIST_doserate_nit = np.zeros(TALLIES)
    LIST_doserate_oth = np.zeros(TALLIES)
    LIST_doserate_g_sec = np.zeros(TALLIES)
    LIST_doserate_boron_normal = np.zeros(TALLIES)
    LIST_doserate_boron_skin = np.zeros(TALLIES)
    LIST_doserate_boron_mucosa = np.zeros(TALLIES)
    LIST_doserate_boron_tumor = np.zeros(TALLIES)
    LIST_doserate_hyd_err = np.zeros(TALLIES)
    LIST_doserate_nit_err = np.zeros(TALLIES)
    LIST_doserate_oth_err = np.zeros(TALLIES)
    LIST_doserate_g_sec_err = np.zeros(TALLIES)
    LIST_doserate_boron_normal_err = np.zeros(TALLIES)
    LIST_doserate_boron_skin_err = np.zeros(TALLIES)
    LIST_doserate_boron_mucosa_err = np.zeros(TALLIES)
    LIST_doserate_boron_tumor_err = np.zeros(TALLIES)
    in_path_n = os.path.join(inBox_nspeDir.get(), inBox_nspeFile.get())
    with open(in_path_n, mode = 'r') as f:
        dum = f.readline() #skipping header part
        while True:
            dum = f.readline()
            dum = dum.replace("\n", "")
            dum = dum.split()
            if len(dum) == 0:
                break
            if len(dum) == 4:
                E_S_n = (float(dum[0])*float(dum[1]))**0.5
                F_S_n = float(dum[2])
                F_S_n_err = float(dum[2])*float(dum[3])
                for i in range(TALLIES):
                    # inter/extrapolate DPF
                    DPF_hyd = interpolate.interp1d(LIST_E_mono_n_flt, LIST_DPF_hyd[i], kind='linear', fill_value='extrapolate')
                    DPF_nit = interpolate.interp1d(LIST_E_mono_n_flt, LIST_DPF_nit[i], kind='linear', fill_value='extrapolate')
                    DPF_oth = interpolate.interp1d(LIST_E_mono_n_flt, LIST_DPF_oth[i], kind='linear', fill_value='extrapolate')
                    DPF_g_sec = interpolate.interp1d(LIST_E_mono_n_flt, LIST_DPF_g_sec[i], kind='linear', fill_value='extrapolate')
                    DPF_boron_normal = interpolate.interp1d(LIST_E_mono_n_flt, LIST_DPF_boron_normal[i], kind='linear', fill_value='extrapolate')
                    DPF_boron_skin = interpolate.interp1d(LIST_E_mono_n_flt, LIST_DPF_boron_skin[i], kind='linear', fill_value='extrapolate')
                    DPF_boron_mucosa = interpolate.interp1d(LIST_E_mono_n_flt, LIST_DPF_boron_mucosa[i], kind='linear', fill_value='extrapolate')
                    DPF_boron_tumor = interpolate.interp1d(LIST_E_mono_n_flt, LIST_DPF_boron_tumor[i], kind='linear', fill_value='extrapolate')
                    # inter/extrapolate DPF error
                    DPF_hyd_err = interpolate.interp1d(LIST_E_mono_n_flt, LIST_DPF_hyd_err[i], kind='linear', fill_value='extrapolate')
                    DPF_nit_err = interpolate.interp1d(LIST_E_mono_n_flt, LIST_DPF_nit_err[i], kind='linear', fill_value='extrapolate')
                    DPF_oth_err = interpolate.interp1d(LIST_E_mono_n_flt, LIST_DPF_oth_err[i], kind='linear', fill_value='extrapolate')
                    DPF_g_sec_err = interpolate.interp1d(LIST_E_mono_n_flt, LIST_DPF_g_sec_err[i], kind='linear', fill_value='extrapolate')
                    DPF_boron_normal_err = interpolate.interp1d(LIST_E_mono_n_flt, LIST_DPF_boron_normal_err[i], kind='linear', fill_value='extrapolate')
                    DPF_boron_skin_err = interpolate.interp1d(LIST_E_mono_n_flt, LIST_DPF_boron_skin_err[i], kind='linear', fill_value='extrapolate')
                    DPF_boron_mucosa_err = interpolate.interp1d(LIST_E_mono_n_flt, LIST_DPF_boron_mucosa_err[i], kind='linear', fill_value='extrapolate')
                    DPF_boron_tumor_err = interpolate.interp1d(LIST_E_mono_n_flt, LIST_DPF_boron_tumor_err[i], kind='linear', fill_value='extrapolate')
                    # calculating dose rates
                    LIST_doserate_hyd[i] += DPF_hyd(E_S_n)*F_S_n*3600*1e-12
                    LIST_doserate_nit[i] += DPF_nit(E_S_n)*F_S_n*3600*1e-12
                    LIST_doserate_oth[i] += DPF_oth(E_S_n)*F_S_n*3600*1e-12
                    LIST_doserate_g_sec[i] += DPF_g_sec(E_S_n)*F_S_n*3600*1e-12
                    LIST_doserate_boron_normal[i] += DPF_boron_normal(E_S_n)*F_S_n*3600*1e-12
                    LIST_doserate_boron_skin[i] += DPF_boron_skin(E_S_n)*F_S_n*3600*1e-12
                    LIST_doserate_boron_mucosa[i] += DPF_boron_mucosa(E_S_n)*F_S_n*3600*1e-12
                    LIST_doserate_boron_tumor[i] += DPF_boron_tumor(E_S_n)*F_S_n*3600*1e-12
                    # calculating dose rate errors (CAUTION: THESE ARE SQUARED VALUES!!)
                    LIST_doserate_hyd_err[i] += (DPF_hyd(E_S_n)*F_S_n_err)**2 + (DPF_hyd_err(E_S_n)*F_S_n)**2
                    LIST_doserate_nit_err[i] += (DPF_nit(E_S_n)*F_S_n_err)**2 + (DPF_nit_err(E_S_n)*F_S_n)**2
                    LIST_doserate_oth_err[i] += (DPF_oth(E_S_n)*F_S_n_err)**2 + (DPF_oth_err(E_S_n)*F_S_n)**2
                    LIST_doserate_g_sec_err[i] += (DPF_g_sec(E_S_n)*F_S_n_err)**2 + (DPF_g_sec_err(E_S_n)*F_S_n)**2
                    LIST_doserate_boron_normal_err[i] += (DPF_boron_normal(E_S_n)*F_S_n_err)**2 + (DPF_boron_normal_err(E_S_n)*F_S_n)**2
                    LIST_doserate_boron_skin_err[i] += (DPF_boron_skin(E_S_n)*F_S_n_err)**2 + (DPF_boron_skin_err(E_S_n)*F_S_n)**2
                    LIST_doserate_boron_mucosa_err[i] += (DPF_boron_mucosa(E_S_n)*F_S_n_err)**2 + (DPF_boron_mucosa_err(E_S_n)*F_S_n)**2
                    LIST_doserate_boron_tumor_err[i] += (DPF_boron_tumor(E_S_n)*F_S_n_err)**2 + (DPF_boron_tumor_err(E_S_n)*F_S_n)**2
    LIST_doserate_hyd_err = np.sqrt(LIST_doserate_hyd_err)*3600*1e-12
    LIST_doserate_nit_err = np.sqrt(LIST_doserate_nit_err)*3600*1e-12
    LIST_doserate_oth_err = np.sqrt(LIST_doserate_oth_err)*3600*1e-12
    LIST_doserate_g_sec_err = np.sqrt(LIST_doserate_g_sec_err)*3600*1e-12
    LIST_doserate_boron_normal_err = np.sqrt(LIST_doserate_boron_normal_err)*3600*1e-12
    LIST_doserate_boron_skin_err = np.sqrt(LIST_doserate_boron_skin)*3600*1e-12
    LIST_doserate_boron_mucosa_err = np.sqrt(LIST_doserate_boron_mucosa_err)*3600*1e-12
    LIST_doserate_boron_tumor_err = np.sqrt(LIST_doserate_boron_tumor_err)*3600*1e-12
    # calculating primary gamma-ray dose rate
    LIST_doserate_g_pri = np.zeros(TALLIES)
    LIST_doserate_g_pri_err = np.zeros(TALLIES)
    in_path_g = os.path.join(inBox_gspeDir.get(), inBox_gspeFile.get())
    with open(in_path_g, mode = 'r') as f:
        dum = f.readline() #skipping header part
        while True:
            dum = f.readline()
            dum = dum.replace("\n", "")
            dum = dum.split()
            if len(dum) == 0:
                break
            if len(dum) == 4:
                E_S_g = (float(dum[0])*float(dum[1]))**0.5
                F_S_g = float(dum[2])
                F_S_g_err = float(dum[2])*float(dum[3])
                for i in range(TALLIES):
                    # inter/extrapolate DPF
                    DPF_g_pri = interpolate.interp1d(LIST_E_mono_g_flt, LIST_DPF_g_pri[i], kind='linear', fill_value='extrapolate')
                    # inter/extrapolate DPF error
                    DPF_g_pri_err = interpolate.interp1d(LIST_E_mono_g_flt, LIST_DPF_g_pri_err[i], kind='linear', fill_value='extrapolate')
                    #calculating dose rates
                    LIST_doserate_g_pri[i] += DPF_g_pri(E_S_g)*F_S_g*3600*1e-12
                    #calculating dose rate errors  (CAUTION: THIS IS SQUARED VALUE!!)
                    LIST_doserate_g_pri_err[i] += (DPF_g_pri(E_S_g)*F_S_g_err)**2 + (DPF_g_pri_err(E_S_g)*F_S_g)**2
    LIST_doserate_g_pri_err = np.sqrt(LIST_doserate_g_pri_err)*3600*1e-12
    # output dose rate distributions
    # calculating summation of dose rates
    LIST_doserate_n = LIST_doserate_hyd + LIST_doserate_nit + LIST_doserate_oth
    LIST_doserate_g = LIST_doserate_g_pri + LIST_doserate_g_sec
    LIST_doserate_total_normal = LIST_doserate_n + LIST_doserate_g + LIST_doserate_boron_normal
    LIST_doserate_total_skin = LIST_doserate_n + LIST_doserate_g + LIST_doserate_boron_skin
    LIST_doserate_total_mucosa = LIST_doserate_n + LIST_doserate_g + LIST_doserate_boron_mucosa
    LIST_doserate_total_tumor = LIST_doserate_n + LIST_doserate_g + LIST_doserate_boron_tumor
    # calculating propagation of dose rate errors
    LIST_doserate_n_err = np.sqrt(LIST_doserate_hyd_err**2 + LIST_doserate_nit_err**2 + LIST_doserate_oth_err**2)
    LIST_doserate_g_err = np.sqrt(LIST_doserate_g_pri_err**2 + LIST_doserate_g_sec_err**2)
    LIST_doserate_total_normal_err = np.sqrt(LIST_doserate_n_err**2 + LIST_doserate_g_err**2 + LIST_doserate_boron_normal_err**2)
    LIST_doserate_total_skin_err = np.sqrt(LIST_doserate_n_err**2 + LIST_doserate_g_err**2 + LIST_doserate_boron_skin_err**2)
    LIST_doserate_total_mucosa_err = np.sqrt(LIST_doserate_n_err**2 + LIST_doserate_g_err**2 + LIST_doserate_boron_mucosa_err**2)
    LIST_doserate_total_tumor_err = np.sqrt(LIST_doserate_n_err**2 + LIST_doserate_g_err**2 + LIST_doserate_boron_tumor_err**2)
    print(' Calculated dose rate distributions!')
    plt.plot(range(TALLIES), LIST_doserate_hyd, linestyle="dotted", color="cyan", label="hydrogen dose")
    plt.plot(range(TALLIES), LIST_doserate_nit, linestyle="dotted", color="yellowgreen", label="nitrogen dose")
    plt.plot(range(TALLIES), LIST_doserate_oth, linestyle="dotted", color="magenta", label="other neutron dose")
    plt.plot(range(TALLIES), LIST_doserate_n, linestyle="solid", color="red", label="neutron dose")
    plt.plot(range(TALLIES), LIST_doserate_g_pri, linestyle="dotted", color="gray", label="primary gamma-ray dose")
    plt.plot(range(TALLIES), LIST_doserate_g_sec, linestyle="dotted", color="green", label="secondary gamma-ray dose")
    plt.plot(range(TALLIES), LIST_doserate_g, linestyle="solid", color="green", label="gamma-ray dose")
    plt.plot(range(TALLIES), LIST_doserate_boron_normal, linestyle="solid", color="blue", label="boron dose in other normal tissue")
    plt.plot(range(TALLIES), LIST_doserate_total_normal, linestyle="solid", color="black", label="total dose in other normal tissue")
    plt.legend()
    plt.title('RBE dose rate distribution')
    plt.xlabel('depth from phantom surface [mm]')
    plt.ylabel('RBE dose rate [$\mathregular{Gy_{RBE}}$/h]')
    plt.grid(color='r', linestyle='dotted', linewidth=1)
    out_path = os.path.join(inBox_workDir.get(), str(os.path.splitext(inBox_distName.get())[0]) + r"_rate.png")
    plt.savefig(out_path)
    plt.close()
    print(' Outputting...')
    out_path = os.path.join(inBox_workDir.get(), inBox_distName.get())
    with open(out_path, mode = 'w') as f:
        f.write('#version: %s\n' %(VERSION))
        f.write('#input neutron energy spectrum: %s\n' %(in_path_n))
        f.write('#input gamma-ray energy spectrum: %s\n' %(in_path_g))
        f.write('#DPFs: %s\n' %(inBox_DPFDir.get()))
        f.write('#column1: Depth from phantom surface [mm]\n')
        f.write('#column2: Hydrogen dose rate [Gy-eq/h]\n')
        f.write('#column3: Uncertainty of hydrogen dose rate [Gy-eq/h]\n')
        f.write('#column4: Nitrogen dose rate [Gy-eq/h]\n')
        f.write('#column5: Uncertainty of nitrogen dose rate [Gy-eq/h]\n')
        f.write('#column6: Other neutron dose rate [Gy-eq/h]\n')
        f.write('#column7: Uncertainty of other neutron dose rate [Gy-eq/h]\n')
        f.write('#column8: Neutron dose rate [Gy-eq/h]\n')
        f.write('#column9: Uncertainty of neutron dose rate [Gy-eq/h]\n')
        f.write('#column10: Primary gamma-ray dose rate [Gy/h]\n')
        f.write('#column11: Uncertainty of primary gamma-ray dose rate [Gy/h]\n')
        f.write('#column12: Secondary gamma-ray dose rate [Gy/h]\n')
        f.write('#column13: Uncertainty of secondary gamma-ray dose rate [Gy/h]\n')
        f.write('#column14: Gamma-ray dose rate [Gy/h]\n')
        f.write('#column15: Uncertainty of gamma-ray dose rate [Gy/h]\n')
        f.write('#column16: Boron dose rate in other normal tissue [Gy-eq/h]\n')
        f.write('#column17: Uncertainty of boron dose rate in other normal tissue [Gy-eq/h]\n')
        f.write('#column18: Boron dose rate in skin [Gy-eq/h]\n')
        f.write('#column19: Uncertainty of boron dose rate in skin [Gy-eq/h]\n')
        f.write('#column20: Boron dose rate in mucosa [Gy-eq/h]\n')
        f.write('#column21: Uncertainty of boron dose rate in mucosa [Gy-eq/h]\n')
        f.write('#column22: Boron dose rate in tumor [Gy-eq/h]\n')
        f.write('#column23: Uncertainty of boron dose rate in tumor [Gy-eq/h]\n')
        f.write('#column24: Total dose rate in other normal tissue [Gy-eq/h]\n')
        f.write('#column25: Uncertainty of total dose rate in other normal tissue [Gy-eq/h]\n')
        f.write('#column26: Total dose rate in skin [Gy-eq/h]\n')
        f.write('#column27: Uncertainty of total dose rate in skin [Gy-eq/h]\n')
        f.write('#column28: Total dose rate in mucosa [Gy-eq/h]\n')
        f.write('#column29: Uncertainty of total dose rate in mucosa [Gy-eq/h]\n')
        f.write('#column30: Total dose rate in tumor [Gy-eq/h]\n')
        f.write('#column31: Uncertainty of total dose rate in tumor [Gy-eq/h]\n')
        for i in range(TALLIES):
            f.write('%d\t' %(i))
            f.write('%10.4e\t' %(LIST_doserate_hyd[i]))
            f.write('%10.4e\t' %(LIST_doserate_hyd_err[i]))
            f.write('%10.4e\t' %(LIST_doserate_nit[i]))
            f.write('%10.4e\t' %(LIST_doserate_nit_err[i]))
            f.write('%10.4e\t' %(LIST_doserate_oth[i]))
            f.write('%10.4e\t' %(LIST_doserate_oth_err[i]))
            f.write('%10.4e\t' %(LIST_doserate_n[i]))
            f.write('%10.4e\t' %(LIST_doserate_n_err[i]))
            f.write('%10.4e\t' %(LIST_doserate_g_pri[i]))
            f.write('%10.4e\t' %(LIST_doserate_g_pri_err[i]))
            f.write('%10.4e\t' %(LIST_doserate_g_sec[i]))
            f.write('%10.4e\t' %(LIST_doserate_g_sec_err[i]))
            f.write('%10.4e\t' %(LIST_doserate_g[i]))
            f.write('%10.4e\t' %(LIST_doserate_g_err[i]))
            f.write('%10.4e\t' %(LIST_doserate_boron_normal[i]))
            f.write('%10.4e\t' %(LIST_doserate_boron_normal_err[i]))
            f.write('%10.4e\t' %(LIST_doserate_boron_skin[i]))
            f.write('%10.4e\t' %(LIST_doserate_boron_skin_err[i]))
            f.write('%10.4e\t' %(LIST_doserate_boron_mucosa[i]))
            f.write('%10.4e\t' %(LIST_doserate_boron_mucosa_err[i]))
            f.write('%10.4e\t' %(LIST_doserate_boron_tumor[i]))
            f.write('%10.4e\t' %(LIST_doserate_boron_tumor_err[i]))
            f.write('%10.4e\t' %(LIST_doserate_total_normal[i]))
            f.write('%10.4e\t' %(LIST_doserate_total_normal_err[i]))
            f.write('%10.4e\t' %(LIST_doserate_total_skin[i]))
            f.write('%10.4e\t' %(LIST_doserate_total_skin_err[i]))
            f.write('%10.4e\t' %(LIST_doserate_total_mucosa[i]))
            f.write('%10.4e\t' %(LIST_doserate_total_mucosa_err[i]))
            f.write('%10.4e\t' %(LIST_doserate_total_tumor[i]))
            f.write('%10.4e\n' %(LIST_doserate_total_tumor_err[i]))
    print(' dpf2dist successfully finished!')
    statusBar["text"] = 'Finished computing dose distributions! (Processing time: %.4f sec)' %(time.time()-start_time)
# dist2idx
def dist2idxFunc(inBox_workDir, 
                 inBox_TALLIES,
                 inBox_distDir, inBox_distFile,
                 inBox_LIMIT_MUCOSA, inBox_LIMIT_SKIN, inBox_LIMIT_NORMAL,
                 inBox_idxName):
    # check PATH
    if inBox_workDir.get() != '' and os.path.isdir(inBox_workDir.get()) == True:
        pass
    else:
        print(" ERROR: Working directory does not exist!")
        statusBar["text"] = 'ERROR: Working directory does not exist!'
        return
    if inBox_distDir.get() != '' and os.path.isfile(os.path.join(inBox_distDir.get(), inBox_distFile.get())) == True:
        pass
    else:
        print(" ERROR: Input distributions file does not exist!")
        statusBar["text"] = 'ERROR: Input distributions file does not exist!'
        return
    if (inBox_LIMIT_MUCOSA.get() != ''
        and inBox_LIMIT_SKIN.get() != ''
        and inBox_LIMIT_NORMAL.get() != ''):
        pass
    else:
        print(" ERROR: Some or all of dose restrictions are missing!")
        statusBar["text"] = 'ERROR: Some or all of dose restrictions are missing!'
        return
    if inBox_idxName.get() != '':
        pass
    else:
        print(" ERROR: Output file name is missing!")
        statusBar["text"] = 'ERROR: Output file name is missing!'
        return
    # main
    start_time = time.time()
    print(' ##############################')
    print(' ########## dist2idx ##########')
    print(' ##############################')
    print(' Computing...')
    # read parameters
    TALLIES = int(inBox_TALLIES.get())
    LIMIT_MUCOSA = float(inBox_LIMIT_MUCOSA.get())
    LIMIT_SKIN = float(inBox_LIMIT_SKIN.get())
    LIMIT_NORMAL = float(inBox_LIMIT_NORMAL.get())
    # read dose rate distributions
    LIST_doserate = np.empty(0, dtype = float)
    with open(os.path.join(inBox_distDir.get(), inBox_distFile.get()), mode='r') as f:
        for i in range(35): # skip header part
            dum = f.readline()
        for i in range(TALLIES):
            dum = f.readline()
            dum = dum.replace("\n", "")
            dum = dum.split()
            dum = np.array(dum)
            dum = dum.astype(float)
            LIST_doserate = np.append(LIST_doserate, dum)
    LIST_doserate = LIST_doserate.reshape(TALLIES, 31)
    LIST_doserate = LIST_doserate.T
    LIST_doserate_hyd = LIST_doserate[1]
    LIST_doserate_hyd_err = LIST_doserate[2]
    LIST_doserate_nit = LIST_doserate[3]
    LIST_doserate_nit_err = LIST_doserate[4]
    LIST_doserate_oth = LIST_doserate[5]
    LIST_doserate_oth_err = LIST_doserate[6]
    LIST_doserate_n = LIST_doserate[7]
    LIST_doserate_n_err = LIST_doserate[8]
    LIST_doserate_g_pri = LIST_doserate[9]
    LIST_doserate_g_pri_err = LIST_doserate[10]
    LIST_doserate_g_sec = LIST_doserate[11]
    LIST_doserate_g_sec_err = LIST_doserate[12]
    LIST_doserate_g = LIST_doserate[13]
    LIST_doserate_g_err = LIST_doserate[14]
    LIST_doserate_boron_normal = LIST_doserate[15]
    LIST_doserate_boron_normal_err = LIST_doserate[16]
    LIST_doserate_boron_skin = LIST_doserate[17]
    LIST_doserate_boron_skin_err = LIST_doserate[18]
    LIST_doserate_boron_mucosa = LIST_doserate[19]
    LIST_doserate_boron_mucosa_err = LIST_doserate[20]
    LIST_doserate_boron_tumor = LIST_doserate[21]
    LIST_doserate_boron_tumor_err = LIST_doserate[22]
    LIST_doserate_total_normal = LIST_doserate[23]
    LIST_doserate_total_normal_err = LIST_doserate[24]
    LIST_doserate_total_skin = LIST_doserate[25]
    LIST_doserate_total_skin_err = LIST_doserate[26]
    LIST_doserate_total_mucosa = LIST_doserate[27]
    LIST_doserate_total_mucosa_err = LIST_doserate[28]
    LIST_doserate_total_tumor = LIST_doserate[29]
    LIST_doserate_total_tumor_err = LIST_doserate[30]
    #output normalized dose distiburtions and dose indexes
    dum = [LIMIT_MUCOSA/np.max(LIST_doserate_total_mucosa), LIMIT_SKIN/np.max(LIST_doserate_total_skin), LIMIT_NORMAL/np.max(LIST_doserate_total_normal)]
    TT = min(dum)
    if dum.index(min(dum)) == 0:
        Limitation_part = 'Mucosa'
        Limitation_dose = LIMIT_MUCOSA
    elif dum.index(min(dum)) == 1:
        Limitation_part = 'Skin'
        Limitation_dose = LIMIT_SKIN
    elif dum.index(min(dum)) == 2:
        Limitation_part = 'Other normal tissue'
        Limitation_dose = LIMIT_NORMAL

#    #####################################
#    ######### from Prof. Tanaka #########
#    TT = 12.0/np.max(LIST_doserate_total_normal)
#    Limitation_part = 'Other normal tissue'
#    Limitation_dose = 12.0
#    #####################################
#    #####################################
    print(' TT was derived!')
    LIST_dose_hyd = LIST_doserate_hyd*TT
    LIST_dose_hyd_err = LIST_doserate_hyd_err*TT
    LIST_dose_nit = LIST_doserate_nit*TT
    LIST_dose_nit_err = LIST_doserate_nit_err*TT
    LIST_dose_oth = LIST_doserate_oth*TT
    LIST_dose_oth_err = LIST_doserate_oth_err*TT
    LIST_dose_n = LIST_doserate_n*TT
    LIST_dose_n_err = LIST_doserate_n_err*TT
    LIST_dose_g_pri = LIST_doserate_g_pri*TT
    LIST_dose_g_pri_err = LIST_doserate_g_pri_err*TT
    LIST_dose_g_sec = LIST_doserate_g_sec*TT
    LIST_dose_g_sec_err = LIST_doserate_g_sec_err*TT
    LIST_dose_g = LIST_doserate_g*TT
    LIST_dose_g_err = LIST_doserate_g_err*TT
    LIST_dose_boron_normal = LIST_doserate_boron_normal*TT
    LIST_dose_boron_normal_err = LIST_doserate_boron_normal_err*TT
    LIST_dose_boron_skin = LIST_doserate_boron_skin*TT
    LIST_dose_boron_skin_err = LIST_doserate_boron_skin_err*TT
    LIST_dose_boron_mucosa = LIST_doserate_boron_mucosa*TT
    LIST_dose_boron_mucosa_err = LIST_doserate_boron_mucosa_err*TT
    LIST_dose_boron_tumor = LIST_doserate_boron_tumor*TT
    LIST_dose_boron_tumor_err = LIST_doserate_boron_tumor_err*TT
    LIST_dose_total_normal = LIST_doserate_total_normal*TT
    LIST_dose_total_normal_err = LIST_doserate_total_normal_err*TT
    LIST_dose_total_skin = LIST_doserate_total_skin*TT
    LIST_dose_total_skin_err = LIST_doserate_total_skin_err*TT
    LIST_dose_total_mucosa = LIST_doserate_total_mucosa*TT
    LIST_dose_total_mucosa_err = LIST_doserate_total_mucosa_err*TT
    LIST_dose_total_tumor = LIST_doserate_total_tumor*TT
    LIST_dose_total_tumor_err = LIST_doserate_total_tumor_err*TT
    print(' Normalized dose distributions were derived!')
    plt.plot(range(TALLIES), LIST_dose_total_normal, linestyle="dotted", color="red", label="total dose in other normal tissue")
    plt.plot(range(TALLIES), LIST_dose_total_skin, linestyle="dotted", color="blue", label="total dose in skin")
    plt.plot(range(TALLIES), LIST_dose_total_mucosa, linestyle="dotted", color="green", label="total dose in mucosa")
    plt.plot(range(TALLIES), LIST_dose_total_tumor, linestyle="solid", color="black", label="total dose in tumor")
    plt.legend()
    plt.title('normalized RBE dose distribution')
    plt.xlabel('depth from phantom surface [mm]')
    plt.ylabel('normalized RBE dose [$\mathregular{Gy_{RBE}}$]')
    plt.grid(color='r', linestyle='dotted', linewidth=1)
    out_path = os.path.join(inBox_workDir.get(), str(os.path.splitext(inBox_distFile.get())[0]) + r"_norm.png")
    plt.savefig(out_path)
    plt.close()
    LIST_dose_total_tumor_extract_reverse = list(reversed(LIST_dose_total_tumor[LIST_dose_total_tumor.tolist().index(np.max(LIST_dose_total_tumor)):]))
    LIST_depth_extract_reverse = list(reversed(range(TALLIES)[LIST_dose_total_tumor.tolist().index(np.max(LIST_dose_total_tumor)):]))
    ADs = interpolate.interp1d(LIST_dose_total_tumor_extract_reverse, LIST_depth_extract_reverse)
    AD = ADs('{:.4e}'.format(max(LIST_dose_total_normal)))
    if max(LIST_dose_total_tumor) < 30:
        print(' AD30 was not available!')
        AD30 = r"N/A"
    else:
        print(' AD30 was derived!')
        AD30 = ADs('{:.4e}'.format(30))
    if max(LIST_dose_total_tumor) < 25:
        print(' AD25 was not available!')
        AD25 = r"N/A"
    else:
        print(' AD25 was derived!')
        AD25 = ADs('{:.4e}'.format(25))
    if max(LIST_dose_total_tumor) < 20:
        print(' AD20 was not available!')
        AD20 = r"N/A"
    else:
        print(' AD20 was derived!')
        AD20 = ADs('{:.4e}'.format(20))
    PTD = max(LIST_dose_total_tumor)
    print(' PTD was derived!')
    IND = 0
    for i in range(TALLIES):
        IND += LIST_dose_total_normal[i]
    print(' IND was derived!')
    print(' Outputting...')
    out_path = os.path.join(inBox_workDir.get(), inBox_idxName.get())
    with open(out_path, mode = 'w') as f:
        f.write('#version: %s\n' %(VERSION))
        f.write('#input dose rate distributions: %s\n' %(os.path.join(inBox_distDir.get(), inBox_distFile.get())))
        f.write('#limitation: %s %4.1f Gy(RBE)\n' %(Limitation_part, Limitation_dose))
        f.write('#TT [h]: %10.4e\n' %(TT))
        f.write('#AD [mm]: %.2f\n' %(AD))
        if type(AD20) is str:
            f.write('#AD20 [mm]: %s\n' %(AD20))
        else:
            f.write('#AD20 [mm]: %.2f\n' %(AD20))
        if type(AD25) is str:
            f.write('#AD25 [mm]: %s\n' %(AD25))
        else:
            f.write('#AD25 [mm]: %.2f\n' %(AD25))
        if type(AD30) is str:
            f.write('#AD30 [mm]: %s\n' %(AD30))
        else:
            f.write('#AD30 [mm]: %.2f\n' %(AD30))
        f.write('#PTD [Gy(RBE)]: %.2f\n' %(PTD))
        f.write('#IND [Gy(RBE) mm]: %.2f\n' %(IND))
        f.write('#column1: Depth from phantom surface [mm]\n')
        f.write('#column2: Hydrogen dose [Gy-eq]\n')
        f.write('#column3: Uncertainty of hydrogen dose [Gy-eq]\n')
        f.write('#column4: Nitrogen dose [Gy-eq]\n')
        f.write('#column5: Uncertainty of nitrogen dose [Gy-eq]\n')
        f.write('#column6: Other neutron dose [Gy-eq]\n')
        f.write('#column7: Uncertainty of other neutron dose [Gy-eq]\n')
        f.write('#column8: Neutron dose [Gy-eq]\n')
        f.write('#column9: Uncertainty of neutron dose [Gy-eq]\n')
        f.write('#column10: Primary gamma-ray dose [Gy]\n')
        f.write('#column11: Uncertainty of primary gamma-ray dose [Gy]\n')
        f.write('#column12: Secondary gamma-ray dose [Gy]\n')
        f.write('#column13: Uncertainty of secondary gamma-ray dose [Gy]\n')
        f.write('#column14: Gamma-ray dose [Gy]\n')
        f.write('#column15: Uncertainty of gamma-ray dose [Gy]\n')
        f.write('#column16: Boron dose in other normal tissue [Gy-eq]\n')
        f.write('#column17: Uncertainty of boron dose in other normal tissue [Gy-eq]\n')
        f.write('#column18: Boron dose in skin [Gy-eq]\n')
        f.write('#column19: Uncertainty of boron dose in skin [Gy-eq]\n')
        f.write('#column20: Boron dose in mucosa [Gy-eq]\n')
        f.write('#column21: Uncertainty of boron dose in mucosa [Gy-eq]\n')
        f.write('#column22: Boron dose in tumor [Gy-eq]\n')
        f.write('#column23: Uncertainty of boron dose in tumor [Gy-eq]\n')
        f.write('#column24: Total dose in other normal tissue [Gy-eq]\n')
        f.write('#column25: Uncertainty of total dose in other normal tissue [Gy-eq]\n')
        f.write('#column26: Total dose in skin [Gy-eq]\n')
        f.write('#column27: Uncertainty of total dose in skin [Gy-eq]\n')
        f.write('#column28: Total dose in mucosa [Gy-eq]\n')
        f.write('#column29: Uncertainty of total dose in mucosa [Gy-eq]\n')
        f.write('#column30: Total dose in tumor [Gy-eq]\n')
        f.write('#column31: Uncertainty of total dose in tumor [Gy-eq]\n')
        for i in range(TALLIES):
            f.write('%d\t' %(i))
            f.write('%10.4e\t' %(LIST_dose_hyd[i]))
            f.write('%10.4e\t' %(LIST_dose_hyd_err[i]))
            f.write('%10.4e\t' %(LIST_dose_nit[i]))
            f.write('%10.4e\t' %(LIST_dose_nit_err[i]))
            f.write('%10.4e\t' %(LIST_dose_oth[i]))
            f.write('%10.4e\t' %(LIST_dose_oth_err[i]))
            f.write('%10.4e\t' %(LIST_dose_n[i]))
            f.write('%10.4e\t' %(LIST_dose_n_err[i]))
            f.write('%10.4e\t' %(LIST_dose_g_pri[i]))
            f.write('%10.4e\t' %(LIST_dose_g_pri_err[i]))
            f.write('%10.4e\t' %(LIST_dose_g_sec[i]))
            f.write('%10.4e\t' %(LIST_dose_g_sec_err[i]))
            f.write('%10.4e\t' %(LIST_dose_g[i]))
            f.write('%10.4e\t' %(LIST_dose_g_err[i]))
            f.write('%10.4e\t' %(LIST_dose_boron_normal[i]))
            f.write('%10.4e\t' %(LIST_dose_boron_normal_err[i]))
            f.write('%10.4e\t' %(LIST_dose_boron_skin[i]))
            f.write('%10.4e\t' %(LIST_dose_boron_skin_err[i]))
            f.write('%10.4e\t' %(LIST_dose_boron_mucosa[i]))
            f.write('%10.4e\t' %(LIST_dose_boron_mucosa_err[i]))
            f.write('%10.4e\t' %(LIST_dose_boron_tumor[i]))
            f.write('%10.4e\t' %(LIST_dose_boron_tumor_err[i]))
            f.write('%10.4e\t' %(LIST_dose_total_normal[i]))
            f.write('%10.4e\t' %(LIST_dose_total_normal_err[i]))
            f.write('%10.4e\t' %(LIST_dose_total_skin[i]))
            f.write('%10.4e\t' %(LIST_dose_total_skin_err[i]))
            f.write('%10.4e\t' %(LIST_dose_total_mucosa[i]))
            f.write('%10.4e\t' %(LIST_dose_total_mucosa_err[i]))
            f.write('%10.4e\t' %(LIST_dose_total_tumor[i]))
            f.write('%10.4e\n' %(LIST_dose_total_tumor_err[i]))
    print(' dist2idx sucessfully finished!')
    statusBar["text"] = 'Finished computing dose indexes! (Processing time: %.4f sec)' %(time.time()-start_time)

#########################
########## GUI ##########
#########################

##### general #####
# generate window
root = tk.Tk()
# window size
root.geometry("1080x720")
# text in title-bar
root.title(u"SiDE4BNCT_"+VERSION)
# version
txt = tk.Label(root, text=VERSION, fg="blue")
txt.place(x=980, y=0)
# status bar
statusBar = tk.Label(root, text='Set parameters to compute dose distributions & dose indexes!', bd=1, relief=tk.SUNKEN, anchor=tk.W)
statusBar.pack(side=tk.BOTTOM, fill=tk.X)
# set working directory
txt_workDir = tk.Label(root, text="PATH to working directory:")
txt_workDir.place(x=10, y=10)
inBox_workDir = tk.Entry(width=97)
inBox_workDir.insert(0, r"./")
inBox_workDir.place(x=10, y=30)
Btn_workDir = tk.Button(root, text="search", bg="spring green", command=lambda: DirectoryPathFunc(inBox_workDir))
Btn_workDir.place(x=600, y=30)
# copy right
txt = tk.Label(root, text="© 2022 Akihisa ISHIKAWA, Nagoya University", fg="gray")
txt.place(x=410, y=670)

##### spe2dpf #####
# origin coordinate
X0_spe2dpf = 30
Y0_spe2dpf = 70
txt = tk.Label(root, text="< spe2dpf >")
txt.place(x=X0_spe2dpf, y=Y0_spe2dpf)
# set table data to generate DPFs
txt = tk.Label(root, text="[ table data ]")
txt.place(x=X0_spe2dpf, y=Y0_spe2dpf+20)
# neutron DPFs
txt = tk.Label(root, text="neutron:")
txt.place(x=X0_spe2dpf+20, y=Y0_spe2dpf+40)
inBox_tabN = tk.Entry(width=27)
inBox_tabN.insert(0, r"./phits/nspe")
inBox_tabN.place(x=X0_spe2dpf+20, y=Y0_spe2dpf+60)
Btn_tabN = tk.Button(root, text="search", bg="spring green", command=lambda: DirectoryPathFunc(inBox_tabN))
Btn_tabN.place(x=X0_spe2dpf+190, y=Y0_spe2dpf+60)
# primary gamma-ray DPF
txt = tk.Label(root, text="primary gamma-ray:")
txt.place(x=X0_spe2dpf+20, y=Y0_spe2dpf+80)
inBox_tabGpri = tk.Entry(width=27)
inBox_tabGpri.insert(0, r"./phits/gspe_pri")
inBox_tabGpri.place(x=X0_spe2dpf+20, y=Y0_spe2dpf+100)
Btn_tabGpri = tk.Button(root, text="search", bg="spring green", command=lambda: DirectoryPathFunc(inBox_tabGpri))
Btn_tabGpri.place(x=X0_spe2dpf+190, y=Y0_spe2dpf+100)
# secondary gamma-ray DPF
txt = tk.Label(root, text="secondary gamma-ray:")
txt.place(x=X0_spe2dpf+20, y=Y0_spe2dpf+120)
inBox_tabGsec = tk.Entry(width=27)
inBox_tabGsec.insert(0, r"./phits/gspe_sec")
inBox_tabGsec.place(x=X0_spe2dpf+20, y=Y0_spe2dpf+140)
Btn_tabGsec = tk.Button(root, text="search", bg="spring green", command=lambda: DirectoryPathFunc(inBox_tabGsec))
Btn_tabGsec.place(x=X0_spe2dpf+190, y=Y0_spe2dpf+140)
# kerma
txt = tk.Label(root, text="[ kerma coefficients ]")
txt.place(x=X0_spe2dpf, y=Y0_spe2dpf+180)
txt = tk.Label(root, text="gamma-ray:")
txt.place(x=X0_spe2dpf+20, y=Y0_spe2dpf+200)
inBox_gkerma = tk.Entry(width=27)
inBox_gkerma.insert(0, r"gamma.txt")
inBox_gkerma.place(x=X0_spe2dpf+20, y=Y0_spe2dpf+220)
Btn_gkerma = tk.Button(root, text="search", bg="spring green", command=lambda: FilePathFunc(r"./kerma", inBox_gkerma))
Btn_gkerma.place(x=X0_spe2dpf+190, y=Y0_spe2dpf+220)
txt = tk.Label(root, text="neutron:")
txt.place(x=X0_spe2dpf+20, y=Y0_spe2dpf+240)
inBox_nkerma = tk.Entry(width=27)
inBox_nkerma.insert(0, r"neutron.txt")
inBox_nkerma.place(x=X0_spe2dpf+20, y=Y0_spe2dpf+260)
Btn_nkerma = tk.Button(root, text="search", bg="spring green", command=lambda: FilePathFunc(r"./kerma", inBox_nkerma))
Btn_nkerma.place(x=X0_spe2dpf+190, y=Y0_spe2dpf+260)
txt = tk.Label(root, text="boron:")
txt.place(x=X0_spe2dpf+20, y=Y0_spe2dpf+280)
inBox_bkerma = tk.Entry(width=27)
inBox_bkerma.insert(0, r"boron.txt")
inBox_bkerma.place(x=X0_spe2dpf+20, y=Y0_spe2dpf+300)
Btn_bkerma = tk.Button(root, text="search", bg="spring green", command=lambda: FilePathFunc(r"./kerma", inBox_bkerma))
Btn_bkerma.place(x=X0_spe2dpf+190, y=Y0_spe2dpf+300)
# parameters in table data
txt = tk.Label(root, text="[ parameters in table data ]")
txt.place(x=X0_spe2dpf+250, y=Y0_spe2dpf+20)
txt = tk.Label(root, text="number of monoenergy:")
txt.place(x=X0_spe2dpf+270, y=Y0_spe2dpf+40)
txt = tk.Label(root, text="neutron:")
txt.place(x=X0_spe2dpf+290, y=Y0_spe2dpf+60)
inBox_BINS_MONO_N = tk.Entry(width=5)
inBox_BINS_MONO_N.insert(0, 117)
inBox_BINS_MONO_N.place(x=X0_spe2dpf+420, y=Y0_spe2dpf+60)
txt = tk.Label(root, text="gamma-ray:")
txt.place(x=X0_spe2dpf+290, y=Y0_spe2dpf+80)
inBox_BINS_MONO_G = tk.Entry(width=5)
inBox_BINS_MONO_G.insert(0, 45)
inBox_BINS_MONO_G.place(x=X0_spe2dpf+420, y=Y0_spe2dpf+80)
txt = tk.Label(root, text="number of energy bins:")
txt.place(x=X0_spe2dpf+270, y=Y0_spe2dpf+100)
txt = tk.Label(root, text="neutron:")
txt.place(x=X0_spe2dpf+290, y=Y0_spe2dpf+120)
inBox_BINS_SPE_N = tk.Entry(width=5)
inBox_BINS_SPE_N.insert(0, 126)
inBox_BINS_SPE_N.place(x=X0_spe2dpf+420, y=Y0_spe2dpf+120)
txt = tk.Label(root, text="primary gamma-ray:")
txt.place(x=X0_spe2dpf+290, y=Y0_spe2dpf+140)
inBox_BINS_SPE_G_PRI = tk.Entry(width=5)
inBox_BINS_SPE_G_PRI.insert(0, 45)
inBox_BINS_SPE_G_PRI.place(x=X0_spe2dpf+420, y=Y0_spe2dpf+140)
txt = tk.Label(root, text="secondary gamma-ray:")
txt.place(x=X0_spe2dpf+290, y=Y0_spe2dpf+160)
inBox_BINS_SPE_G_SEC = tk.Entry(width=5)
inBox_BINS_SPE_G_SEC.insert(0, 45)
inBox_BINS_SPE_G_SEC.place(x=X0_spe2dpf+420, y=Y0_spe2dpf+160)
txt = tk.Label(root, text="range of energy [MeV]:")
txt.place(x=X0_spe2dpf+270, y=Y0_spe2dpf+180)
txt = tk.Label(root, text="min")
txt.place(x=X0_spe2dpf+420, y=Y0_spe2dpf+180)
txt = tk.Label(root, text="max")
txt.place(x=X0_spe2dpf+460, y=Y0_spe2dpf+180)
txt = tk.Label(root, text="neutron:")
txt.place(x=X0_spe2dpf+290, y=Y0_spe2dpf+200)
inBox_EMIN_SPE_N = tk.Entry(width=5)
inBox_EMIN_SPE_N.insert(0, '1e-12')
inBox_EMIN_SPE_N.place(x=X0_spe2dpf+420, y=Y0_spe2dpf+200)
inBox_EMAX_SPE_N = tk.Entry(width=5)
inBox_EMAX_SPE_N.insert(0, '1e2')
inBox_EMAX_SPE_N.place(x=X0_spe2dpf+460, y=Y0_spe2dpf+200)
txt = tk.Label(root, text="primary gamma-ray:")
txt.place(x=X0_spe2dpf+290, y=Y0_spe2dpf+220)
inBox_EMIN_SPE_G_PRI = tk.Entry(width=5)
inBox_EMIN_SPE_G_PRI.insert(0, '1e-3')
inBox_EMIN_SPE_G_PRI.place(x=X0_spe2dpf+420, y=Y0_spe2dpf+220)
inBox_EMAX_SPE_G_PRI = tk.Entry(width=5)
inBox_EMAX_SPE_G_PRI.insert(0, '1e2')
inBox_EMAX_SPE_G_PRI.place(x=X0_spe2dpf+460, y=Y0_spe2dpf+220)
txt = tk.Label(root, text="secondary gamma-ray:")
txt.place(x=X0_spe2dpf+290, y=Y0_spe2dpf+240)
inBox_EMIN_SPE_G_SEC = tk.Entry(width=5)
inBox_EMIN_SPE_G_SEC.insert(0, '1e-3')
inBox_EMIN_SPE_G_SEC.place(x=X0_spe2dpf+420, y=Y0_spe2dpf+240)
inBox_EMAX_SPE_G_SEC = tk.Entry(width=5)
inBox_EMAX_SPE_G_SEC.insert(0, '1e2')
inBox_EMAX_SPE_G_SEC.place(x=X0_spe2dpf+460, y=Y0_spe2dpf+240)
txt = tk.Label(root, text="number of tallies:")
txt.place(x=X0_spe2dpf+270, y=Y0_spe2dpf+260)
inBox_TALLIES = tk.Entry(width=5)
inBox_TALLIES.insert(0, 200)
inBox_TALLIES.place(x=X0_spe2dpf+420, y=Y0_spe2dpf+260)
# radiobiological coefficients
txt = tk.Label(root, text="[ radiobiological coefficients ]")
txt.place(x=X0_spe2dpf+500, y=Y0_spe2dpf+20)
txt = tk.Label(root, text="boron concentration")
txt.place(x=X0_spe2dpf+520, y=Y0_spe2dpf+40)
txt = tk.Label(root, text="in normal tissue [ppm]:")
txt.place(x=X0_spe2dpf+520, y=Y0_spe2dpf+60)
inBox_CbN = ttk.Combobox(root, width=4, values=('25.0', '24.0', '23.0', '22.0', '21.0', '20.0', '15.0', '10.0'))
inBox_CbN.insert(0, 25.0)
inBox_CbN.place(x=X0_spe2dpf+650, y=Y0_spe2dpf+60)
txt = tk.Label(root, text="TNR:")
txt.place(x=X0_spe2dpf+520, y=Y0_spe2dpf+80)
inBox_TNR = ttk.Combobox(root, width=4, values=('3.5', '3.4', '3.3', '3.2', '3.1', '3.0'))
inBox_TNR.insert(0, 3.5)
inBox_TNR.place(x=X0_spe2dpf+650, y=Y0_spe2dpf+80)
txt = tk.Label(root, text="CBE:")
txt.place(x=X0_spe2dpf+520, y=Y0_spe2dpf+100)
txt = tk.Label(root, text="tumor:")
txt.place(x=X0_spe2dpf+540, y=Y0_spe2dpf+120)
inBox_CBEt = ttk.Combobox(root, width=4, values=('4.0'))
inBox_CBEt.insert(0, 4.0)
inBox_CBEt.place(x=X0_spe2dpf+650, y=Y0_spe2dpf+120)
txt = tk.Label(root, text="mucosa:")
txt.place(x=X0_spe2dpf+540, y=Y0_spe2dpf+140)
inBox_CBEm = ttk.Combobox(root, width=4, values=('4.9'))
inBox_CBEm.insert(0, 4.9)
inBox_CBEm.place(x=X0_spe2dpf+650, y=Y0_spe2dpf+140)
txt = tk.Label(root, text="skin:")
txt.place(x=X0_spe2dpf+540, y=Y0_spe2dpf+160)
inBox_CBEs = ttk.Combobox(root, width=4, values=('2.5'))
inBox_CBEs.insert(0, 2.5)
inBox_CBEs.place(x=X0_spe2dpf+650, y=Y0_spe2dpf+160)
txt = tk.Label(root, text="other normal tissue:")
txt.place(x=X0_spe2dpf+540, y=Y0_spe2dpf+180)
inBox_CBEn = ttk.Combobox(root, width=4, values=('1.34'))
inBox_CBEn.insert(0, 1.34)
inBox_CBEn.place(x=X0_spe2dpf+650, y=Y0_spe2dpf+180)
txt = tk.Label(root, text="RBE:")
txt.place(x=X0_spe2dpf+520, y=Y0_spe2dpf+200)
txt = tk.Label(root, text="gamma-ray:")
txt.place(x=X0_spe2dpf+540, y=Y0_spe2dpf+220)
inBox_RBEgamma = ttk.Combobox(root, width=4, values=('1.0'))
inBox_RBEgamma.insert(0, 1.0)
inBox_RBEgamma.place(x=X0_spe2dpf+650, y=Y0_spe2dpf+220)
txt = tk.Label(root, text="thermal neutron:")
txt.place(x=X0_spe2dpf+540, y=Y0_spe2dpf+240)
inBox_RBEthermal = ttk.Combobox(root, width=4, values=('2.9'))
inBox_RBEthermal.insert(0, 2.9)
inBox_RBEthermal.place(x=X0_spe2dpf+650, y=Y0_spe2dpf+240)
txt = tk.Label(root, text="fast neutron:")
txt.place(x=X0_spe2dpf+540, y=Y0_spe2dpf+260)
inBox_RBEfast = ttk.Combobox(root, width=4, values=('2.4'))
inBox_RBEfast.insert(0, 2.4)
inBox_RBEfast.place(x=X0_spe2dpf+650, y=Y0_spe2dpf+260)
# monoenergy lists
txt = tk.Label(root, text="[ monoenergy lists ]")
txt.place(x=X0_spe2dpf+750, y=Y0_spe2dpf+20)
txt = tk.Label(root, text="neutron:")
txt.place(x=X0_spe2dpf+770, y=Y0_spe2dpf+40)
inBox_engListN = tk.Entry(width=27)
inBox_engListN.insert(0, r"002_monoenergy_neutron.txt")
inBox_engListN.place(x=X0_spe2dpf+770, y=Y0_spe2dpf+60)
Btn_engListN = tk.Button(root, text="search", bg="spring green", command=lambda: FilePathFunc(r"./phits", inBox_engListN))
Btn_engListN.place(x=X0_spe2dpf+940, y=Y0_spe2dpf+60)
txt = tk.Label(root, text="gamma-ray:")
txt.place(x=X0_spe2dpf+770, y=Y0_spe2dpf+80)
inBox_engListG = tk.Entry(width=27)
inBox_engListG.insert(0, r"002_monoenergy_gamma.txt")
inBox_engListG.place(x=X0_spe2dpf+770, y=Y0_spe2dpf+100)
Btn_engListG = tk.Button(root, text="search", bg="spring green", command=lambda: FilePathFunc(r"./phits", inBox_engListG))
Btn_engListG.place(x=X0_spe2dpf+940, y=Y0_spe2dpf+100)
# set output
txt = tk.Label(root, text="[ output directory names ]")
txt.place(x=X0_spe2dpf+750, y=Y0_spe2dpf+140)
txt = tk.Label(root, text="neutron:")
txt.place(x=X0_spe2dpf+770, y=Y0_spe2dpf+160)
inBox_nDPFName = tk.Entry(width=27)
inBox_nDPFName.insert(0, r"dpf")
inBox_nDPFName.place(x=X0_spe2dpf+770, y=Y0_spe2dpf+180)
txt = tk.Label(root, text="boron:")
txt.place(x=X0_spe2dpf+770, y=Y0_spe2dpf+200)
inBox_bDPFName = tk.Entry(width=27)
inBox_bDPFName.insert(0, r"dpf")
inBox_bDPFName.place(x=X0_spe2dpf+770, y=Y0_spe2dpf+220)
txt = tk.Label(root, text="primary gamma-ray:")
txt.place(x=X0_spe2dpf+770, y=Y0_spe2dpf+240)
inBox_gPriDPFName = tk.Entry(width=27)
inBox_gPriDPFName.insert(0, r"dpf")
inBox_gPriDPFName.place(x=X0_spe2dpf+770, y=Y0_spe2dpf+260)
txt = tk.Label(root, text="secondary gamma-ray:")
txt.place(x=X0_spe2dpf+770, y=Y0_spe2dpf+280)
inBox_gSecDPFName = tk.Entry(width=27)
inBox_gSecDPFName.insert(0, r"dpf")
inBox_gSecDPFName.place(x=X0_spe2dpf+770, y=Y0_spe2dpf+300)
Btn_spe2dpfN = tk.Button(root, text="generate", bg="cornflower blue",
                         command=lambda: spe2dpf4neutronDPFsFunc(inBox_workDir,
                                                                 inBox_tabN, inBox_tabGpri, inBox_tabGsec,
                                                                 inBox_engListN, inBox_engListG,
                                                                 inBox_gkerma, inBox_nkerma, inBox_bkerma,
                                                                 inBox_CbN, inBox_TNR,
                                                                 inBox_CBEt, inBox_CBEm, inBox_CBEs, inBox_CBEn,
                                                                 inBox_RBEgamma, inBox_RBEthermal, inBox_RBEfast,
                                                                 inBox_BINS_MONO_N, inBox_BINS_MONO_G,
                                                                 inBox_BINS_SPE_N, inBox_BINS_SPE_G_PRI, inBox_BINS_SPE_G_SEC,
                                                                 inBox_EMIN_SPE_N, inBox_EMIN_SPE_G_PRI, inBox_EMIN_SPE_G_SEC,
                                                                 inBox_EMAX_SPE_N, inBox_EMAX_SPE_G_PRI, inBox_EMAX_SPE_G_SEC,
                                                                 inBox_TALLIES,
                                                                 inBox_nDPFName))
Btn_spe2dpfN.place(x=X0_spe2dpf+940, y=Y0_spe2dpf+180)
Btn_spe2dpfB = tk.Button(root, text="generate", bg="cornflower blue",
                         command=lambda: spe2dpf4boronDPFsFunc(inBox_workDir,
                                                               inBox_tabN, inBox_tabGpri, inBox_tabGsec,
                                                               inBox_engListN, inBox_engListG,
                                                               inBox_gkerma, inBox_nkerma, inBox_bkerma,
                                                               inBox_CbN, inBox_TNR,
                                                               inBox_CBEt, inBox_CBEm, inBox_CBEs, inBox_CBEn,
                                                               inBox_RBEgamma, inBox_RBEthermal, inBox_RBEfast,
                                                               inBox_BINS_MONO_N, inBox_BINS_MONO_G,
                                                               inBox_BINS_SPE_N, inBox_BINS_SPE_G_PRI, inBox_BINS_SPE_G_SEC,
                                                               inBox_EMIN_SPE_N, inBox_EMIN_SPE_G_PRI, inBox_EMIN_SPE_G_SEC,
                                                               inBox_EMAX_SPE_N, inBox_EMAX_SPE_G_PRI, inBox_EMAX_SPE_G_SEC,
                                                               inBox_TALLIES,
                                                               inBox_bDPFName))
Btn_spe2dpfB.place(x=X0_spe2dpf+940, y=Y0_spe2dpf+220)
Btn_spe2dpfGpri = tk.Button(root, text="generate", bg="cornflower blue",
                            command=lambda: spe2dpf4gammaPriDPFFunc(inBox_workDir,
                                                                    inBox_tabN, inBox_tabGpri, inBox_tabGsec,
                                                                    inBox_engListN, inBox_engListG,
                                                                    inBox_gkerma, inBox_nkerma, inBox_bkerma,
                                                                    inBox_CbN, inBox_TNR,
                                                                    inBox_CBEt, inBox_CBEm, inBox_CBEs, inBox_CBEn,
                                                                    inBox_RBEgamma, inBox_RBEthermal, inBox_RBEfast,
                                                                    inBox_BINS_MONO_N, inBox_BINS_MONO_G,
                                                                    inBox_BINS_SPE_N, inBox_BINS_SPE_G_PRI, inBox_BINS_SPE_G_SEC,
                                                                    inBox_EMIN_SPE_N, inBox_EMIN_SPE_G_PRI, inBox_EMIN_SPE_G_SEC,
                                                                    inBox_EMAX_SPE_N, inBox_EMAX_SPE_G_PRI, inBox_EMAX_SPE_G_SEC,
                                                                    inBox_TALLIES,
                                                                    inBox_gPriDPFName))
Btn_spe2dpfGpri.place(x=X0_spe2dpf+940, y=Y0_spe2dpf+260)
Btn_spe2dpfGsec = tk.Button(root, text="generate", bg="cornflower blue",
                            command=lambda: spe2dpf4gammaSecDPFFunc(inBox_workDir,
                                                                    inBox_tabN, inBox_tabGpri, inBox_tabGsec,
                                                                    inBox_engListN, inBox_engListG,
                                                                    inBox_gkerma, inBox_nkerma, inBox_bkerma,
                                                                    inBox_CbN, inBox_TNR,
                                                                    inBox_CBEt, inBox_CBEm, inBox_CBEs, inBox_CBEn,
                                                                    inBox_RBEgamma, inBox_RBEthermal, inBox_RBEfast,
                                                                    inBox_BINS_MONO_N, inBox_BINS_MONO_G,
                                                                    inBox_BINS_SPE_N, inBox_BINS_SPE_G_PRI, inBox_BINS_SPE_G_SEC,
                                                                    inBox_EMIN_SPE_N, inBox_EMIN_SPE_G_PRI, inBox_EMIN_SPE_G_SEC,
                                                                    inBox_EMAX_SPE_N, inBox_EMAX_SPE_G_PRI, inBox_EMAX_SPE_G_SEC,
                                                                    inBox_TALLIES,
                                                                    inBox_gSecDPFName))
Btn_spe2dpfGsec.place(x=X0_spe2dpf+940, y=Y0_spe2dpf+300)

##### dpf2dist #####
# origin coordinate
X0_dpf2dist = 30
Y0_dpf2dist = 410
txt = tk.Label(root, text="< dpf2dist >")
txt.place(x=X0_dpf2dist, y=Y0_dpf2dist)
# set inputs
txt = tk.Label(root, text="[ input energy spectra ]")
txt.place(x=X0_dpf2dist, y=Y0_dpf2dist+20)
# neutron spectrum
txt = tk.Label(root, text="neutron spectra:")
txt.place(x=X0_dpf2dist+20, y=Y0_dpf2dist+40)
inBox_nspeDir = tk.Entry(width=27)
inBox_nspeDir.insert(0, r"./spectra")
inBox_nspeDir.place(x=X0_dpf2dist+20, y=Y0_dpf2dist+60)
Btn_nspeDir = tk.Button(root, text="search", bg="spring green", command=lambda: DirectoryPathFunc(inBox_nspeDir))
Btn_nspeDir.place(x=X0_dpf2dist+190, y=Y0_dpf2dist+60)
txt = tk.Label(root, text="input neutron spectrum:")
txt.place(x=X0_dpf2dist+20, y=Y0_dpf2dist+80)
inBox_nspeFile = tk.Entry(width=27)
inBox_nspeFile.insert(0, r"NUANS-LiF(nat)-nspe.dat")
inBox_nspeFile.place(x=X0_dpf2dist+20, y=Y0_dpf2dist+100)
Btn_nspeFile = tk.Button(root, text="search", bg="spring green", command=lambda: FilePathFunc(inBox_nspeDir.get(), inBox_nspeFile))
Btn_nspeFile.place(x=X0_dpf2dist+190, y=Y0_dpf2dist+100)
# gamma-ray spectrum
txt = tk.Label(root, text="gamma-ray spectra:")
txt.place(x=X0_dpf2dist+20, y=Y0_dpf2dist+120)
inBox_gspeDir = tk.Entry(width=27)
inBox_gspeDir.insert(0, r"./spectra")
inBox_gspeDir.place(x=X0_dpf2dist+20, y=Y0_dpf2dist+140)
Btn_gspeDir = tk.Button(root, text="search", bg="spring green", command=lambda: DirectoryPathFunc(inBox_gspeDir))
Btn_gspeDir.place(x=X0_dpf2dist+190, y=Y0_dpf2dist+140)
txt = tk.Label(root, text="input gamma-ray spectrum:")
txt.place(x=X0_dpf2dist+20, y=Y0_dpf2dist+160)
inBox_gspeFile = tk.Entry(width=27)
inBox_gspeFile.insert(0, r"NUANS-LiF(nat)-gspe.dat")
inBox_gspeFile.place(x=X0_dpf2dist+20, y=Y0_dpf2dist+180)
Btn_gspeFile = tk.Button(root, text="search", bg="spring green", command=lambda: FilePathFunc(inBox_gspeDir.get(), inBox_gspeFile))
Btn_gspeFile.place(x=X0_dpf2dist+190, y=Y0_dpf2dist+180)
# DPFs
txt = tk.Label(root, text="[ DPFs ]")
txt.place(x=X0_dpf2dist+250, y=Y0_dpf2dist+20)
txt = tk.Label(root, text="directory of DPFs:")
txt.place(x=X0_dpf2dist+270, y=Y0_dpf2dist+40)
inBox_DPFDir = tk.Entry(width=27)
inBox_DPFDir.insert(0, r"./dpf")
inBox_DPFDir.place(x=X0_dpf2dist+270, y=Y0_dpf2dist+60)
Btn_DPFDir = tk.Button(root, text="search", bg="spring green", command=lambda: DirectoryPathFunc(inBox_DPFDir))
Btn_DPFDir.place(x=X0_dpf2dist+440, y=Y0_dpf2dist+60)
# set output
txt = tk.Label(root, text="[ ouput settings ]")
txt.place(x=X0_dpf2dist+250, y=Y0_dpf2dist+100)
txt = tk.Label(root, text="output file name:")
txt.place(x=X0_dpf2dist+270, y=Y0_dpf2dist+120)
inBox_distName = tk.Entry(width=27)
inBox_distName.insert(0, r"dist.dat")
inBox_distName.place(x=X0_dpf2dist+270, y=Y0_dpf2dist+140)
Btn_dpf2dist = tk.Button(root, text="compute", bg="yellow",
                         command=lambda: dpf2distFunc(inBox_workDir,
                                                      inBox_engListN, inBox_engListG,
                                                      inBox_BINS_MONO_N, inBox_BINS_MONO_G,
                                                      inBox_TALLIES,
                                                      inBox_DPFDir,
                                                      inBox_nspeDir, inBox_nspeFile,
                                                      inBox_gspeDir, inBox_gspeFile,
                                                      inBox_distName))
Btn_dpf2dist.place(x=X0_dpf2dist+440, y=Y0_dpf2dist+140)

##### dist2idx #####
# origin coordinate
X0_dist2idx = 530
Y0_dist2idx = 410
txt = tk.Label(root, text="< dist2idx >")
txt.place(x=X0_dist2idx, y=Y0_dist2idx)
# set input
txt = tk.Label(root, text="[ input dose rate distributions ]")
txt.place(x=X0_dist2idx, y=Y0_dist2idx+20)
txt = tk.Label(root, text="directory of distributions:")
txt.place(x=X0_dist2idx+20, y=Y0_dist2idx+40)
inBox_distDir = tk.Entry(width=27)
inBox_distDir.insert(0, r"./")
inBox_distDir.place(x=X0_dist2idx+20, y=Y0_dist2idx+60)
Btn_DPFDir = tk.Button(root, text="search", bg="spring green", command=lambda: DirectoryPathFunc(inBox_distDir))
Btn_DPFDir.place(x=X0_dist2idx+190, y=Y0_dist2idx+60)
txt = tk.Label(root, text="input distributions:")
txt.place(x=X0_dist2idx+20, y=Y0_dist2idx+80)
inBox_distFile = tk.Entry(width=27)
inBox_distFile.insert(0, r"dist.dat")
inBox_distFile.place(x=X0_dist2idx+20, y=Y0_dist2idx+100)
Btn_distFile = tk.Button(root, text="search", bg="spring green", command=lambda: FilePathFunc(inBox_distDir.get(), inBox_distFile))
Btn_distFile.place(x=X0_dist2idx+190, y=Y0_dist2idx+100)
# dose restrictions
txt = tk.Label(root, text="[ dose restrictions ]")
txt.place(x=X0_dist2idx+250, y=Y0_dist2idx+20)
txt = tk.Label(root, text=" tolerance dose [Gy-eq]:")
txt.place(x=X0_dist2idx+270, y=Y0_dist2idx+40)
txt = tk.Label(root, text="mucosa:")
txt.place(x=X0_dist2idx+290, y=Y0_dist2idx+60)
inBox_LIMIT_MUCOSA = tk.Entry(width=5)
inBox_LIMIT_MUCOSA.insert(0, 12.0)
inBox_LIMIT_MUCOSA.place(x=X0_dist2idx+420, y=Y0_dist2idx+60)
txt = tk.Label(root, text="skin:")
txt.place(x=X0_dist2idx+290, y=Y0_dist2idx+80)
inBox_LIMIT_SKIN = tk.Entry(width=5)
inBox_LIMIT_SKIN.insert(0, 15.0)
inBox_LIMIT_SKIN.place(x=X0_dist2idx+420, y=Y0_dist2idx+80)
txt = tk.Label(root, text="other normal tissue:")
txt.place(x=X0_dist2idx+290, y=Y0_dist2idx+100)
inBox_LIMIT_NORMAL = tk.Entry(width=5)
inBox_LIMIT_NORMAL.insert(0, 10.0)
inBox_LIMIT_NORMAL.place(x=X0_dist2idx+420, y=Y0_dist2idx+100)
txt = tk.Label(root, text=" prescribed dose [Gy-eq]:")
txt.place(x=X0_dist2idx+270, y=Y0_dist2idx+120)
inBox_PRESCRIBED = tk.Entry(width=5)
inBox_PRESCRIBED.insert(0, 100.0)
inBox_PRESCRIBED.place(x=X0_dist2idx+420, y=Y0_dist2idx+120)
# set output
txt = tk.Label(root, text="[ output settings ]")
txt.place(x=X0_dist2idx+250, y=Y0_dist2idx+160)
txt = tk.Label(root, text="output file name:")
txt.place(x=X0_dist2idx+270, y=Y0_dist2idx+180)
inBox_idxName = tk.Entry(width=27)
inBox_idxName.insert(0, r"idx.dat")
inBox_idxName.place(x=X0_dist2idx+270, y=Y0_dist2idx+200)
Btn_dist2idx = tk.Button(root, text="compute", bg="yellow",
                         command=lambda: dist2idxFunc(inBox_workDir, 
                                                      inBox_TALLIES,
                                                      inBox_distDir, inBox_distFile,
                                                      inBox_LIMIT_MUCOSA, inBox_LIMIT_SKIN, inBox_LIMIT_NORMAL,
                                                      inBox_idxName))
Btn_dist2idx.place(x=X0_dist2idx+440, y=Y0_dist2idx+200)

# keep window
root.mainloop()
