import os
import numpy as np
import scipy as sci
from scipy import interpolate
import matplotlib
from matplotlib import pyplot as plt
import time
import math

#VERSION
VERSION = "ver220310"

#CONSTANTS
PI = 3.141592653589793    #pi

C_NORMAL = 25.0           #BPA (generic name: borofalan) concentration in normal tisuue
RBE_GAMMA = 1.0           #relative biological effectiveness for gamma-ray dose
RBE_THERMAL = 2.9         #relative biological effectiveness for nitrogen dose
RBE_FAST = 2.4            #relative biological effectiveness for hydrogen dose and other neutron dose
CBE_TUMOR = 4.0           #compound biological effectiveness of borofalan for tumor
CBE_MUCOSA = 4.9          #compound biological effectiveness of borofalan for mucosa
CBE_SKIN = 2.5            #compound biological effectiveness of borofalan for skin
CBE_NORMAL = 1.34         #compound biological effectiveness of borofalan for other normal tissue
TNR = 3.5                 #tumor/normal tissue ratio of borofalan concentration
LIMIT_MUCOSA = 12.0       #dose limit for mucosa in unit of Gy(RBE)
LIMIT_SKIN = 15.0         #dose limit for skin in unit of Gy(RBE)
LIMIT_NORMAL = 10.0       #dose limit for other normal tissue in unit of Gy(RBE)

HEADER_N = 16             #number of header lines in ndose.dat file
HEADER_G = 4              #number of header lines in gdose.dat file
TALLIES = 200             #number of tallies in phantom

#load dose rate distributions
start_time = time.time()
LIST_doserate_hyd = []
LIST_doserate_nit = []
LIST_doserate_oth = []
LIST_doserate_n = []
LIST_doserate_g_pri = []
LIST_doserate_g_sec = []
LIST_doserate_g = []
LIST_doserate_boron_normal = []
LIST_doserate_boron_skin = []
LIST_doserate_boron_mucosa = []
LIST_doserate_boron_tumor = []
LIST_doserate_total_normal = []
LIST_doserate_total_skin = []
LIST_doserate_total_mucosa = []
LIST_doserate_total_tumor = []
ndose_path = './depVSndose25ppm-fc-2.2e5.dat'
with open(ndose_path, mode = 'r') as f:
    for i in range(HEADER_N):
        dum = f.readline() #skipping header part
    for i in range(TALLIES):
        dum = f.readline()
        dum = dum.replace("\n", "")
        dum = dum.split()
        LIST_doserate_hyd.append(float(dum[1]))
        LIST_doserate_nit.append(float(dum[2]))
        LIST_doserate_oth.append(float(dum[3]))
        LIST_doserate_n.append(float(dum[4]))
        LIST_doserate_g_sec.append(float(dum[5]))
        LIST_doserate_boron_normal.append(float(dum[6]))
        LIST_doserate_boron_skin.append(float(dum[7]))
        LIST_doserate_boron_mucosa.append(float(dum[8]))
        LIST_doserate_boron_tumor.append(float(dum[9]))
gdose_path = 'C:/Lab/001_Study/003_Prog/phantom-dose-new-ver4/002_gamma/001_main/depVSgdose_x100.dat'
with open(gdose_path, mode = 'r') as f:
    for i in range(HEADER_G):
        dum = f.readline() #skipping header part
    for i in range(TALLIES):
        dum = f.readline()
        dum = dum.replace("\n", "")
        dum = dum.split()
        LIST_doserate_g_pri.append(float(dum[1]))
LIST_doserate_g = (np.array(LIST_doserate_g_pri)+np.array(LIST_doserate_g_sec)).tolist()
LIST_doserate_total_normal = (np.array(LIST_doserate_n)+np.array(LIST_doserate_g)+np.array(LIST_doserate_boron_normal)).tolist()
LIST_doserate_total_skin = (np.array(LIST_doserate_n)+np.array(LIST_doserate_g)+np.array(LIST_doserate_boron_skin)).tolist()
LIST_doserate_total_mucosa = (np.array(LIST_doserate_n)+np.array(LIST_doserate_g)+np.array(LIST_doserate_boron_mucosa)).tolist()
LIST_doserate_total_tumor = (np.array(LIST_doserate_n)+np.array(LIST_doserate_g)+np.array(LIST_doserate_boron_tumor)).tolist()

#output normalized dose distiburtions and dose indexes
start_time = time.time()
dum = [LIMIT_MUCOSA/max(LIST_doserate_total_mucosa), LIMIT_SKIN/max(LIST_doserate_total_skin), LIMIT_NORMAL/max(LIST_doserate_total_normal)]
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

###### from Tanaka sensei ######
TT = 12.0/max(LIST_doserate_total_normal)
Limitation_part = 'Other normal tissue'
Limitation_dose = 12.0
################################

LIST_dose_hyd = (np.array(LIST_doserate_hyd)*TT).tolist()
LIST_dose_nit = (np.array(LIST_doserate_nit)*TT).tolist()
LIST_dose_oth = (np.array(LIST_doserate_oth)*TT).tolist()
LIST_dose_n = (np.array(LIST_doserate_n)*TT).tolist()
LIST_dose_g_pri = (np.array(LIST_doserate_g_pri)*TT).tolist()
LIST_dose_g_sec = (np.array(LIST_doserate_g_sec)*TT).tolist()
LIST_dose_g = (np.array(LIST_doserate_g)*TT).tolist()
LIST_dose_boron_normal = (np.array(LIST_doserate_boron_normal)*TT).tolist()
LIST_dose_boron_skin = (np.array(LIST_doserate_boron_skin)*TT).tolist()
LIST_dose_boron_mucosa = (np.array(LIST_doserate_boron_mucosa)*TT).tolist()
LIST_dose_boron_tumor = (np.array(LIST_doserate_boron_tumor)*TT).tolist()
LIST_dose_total_normal = (np.array(LIST_doserate_total_normal)*TT).tolist()
LIST_dose_total_skin = (np.array(LIST_doserate_total_skin)*TT).tolist()
LIST_dose_total_mucosa = (np.array(LIST_doserate_total_mucosa)*TT).tolist()
LIST_dose_total_tumor = (np.array(LIST_doserate_total_tumor)*TT).tolist()
plt.plot(range(TALLIES), LIST_dose_total_normal, linestyle="dotted", color="red", label="total dose in other normal tissue")
plt.plot(range(TALLIES), LIST_dose_total_skin, linestyle="dotted", color="blue", label="total dose in skin")
plt.plot(range(TALLIES), LIST_dose_total_mucosa, linestyle="dotted", color="green", label="total dose in mucosa")
plt.plot(range(TALLIES), LIST_dose_total_tumor, linestyle="solid", color="black", label="total dose in tumor")
plt.legend()
plt.title('normalized RBE/CBE-weighted dose distribution')
plt.xlabel('depth from phantom surface [mm]')
plt.ylabel('normalized RBE/CBE-weighted dose [$\mathregular{Gy_{eq}}$]')
#plt.yscale('log')
plt.grid(color='r', linestyle='dotted', linewidth=1)
plt.savefig('dist_norm.png')
print('normalized dose distributions were calculated!')
print('processing time was %f sec!' %(time.time()-start_time))
plt.close()
LIST_dose_total_tumor_extract_reverse = list(reversed(LIST_dose_total_tumor[LIST_dose_total_tumor.index(max(LIST_dose_total_tumor)):]))
LIST_depth_extract_reverse = list(reversed(range(TALLIES)[LIST_dose_total_tumor.index(max(LIST_dose_total_tumor)):]))
ADs = interpolate.interp1d(LIST_dose_total_tumor_extract_reverse, LIST_depth_extract_reverse)
AD = ADs('{:.4e}'.format(max(LIST_dose_total_normal)))
if max(LIST_dose_total_tumor) < 30:
    print('AD30 was not available!')
    AD30 = r"N/A"
else:
    print('AD30 was derived!')
    AD30 = ADs('{:.4e}'.format(30))
if max(LIST_dose_total_tumor) < 25:
    print('AD25 was not available!')
    AD25 = r"N/A"
else:
    print('AD25 was derived!')
    AD25 = ADs('{:.4e}'.format(25))
if max(LIST_dose_total_tumor) < 20:
    print('AD20 was not available!')
    AD20 = r"N/A"
else:
    print('AD20 was derived!')
    AD20 = ADs('{:.4e}'.format(20))
PTD = max(LIST_dose_total_tumor)
IND = 0
for i in range(TALLIES):
    IND += LIST_dose_total_normal[i]
out_path = r"combination contamination/idx-fc-2.2e5_gc_x100.dat"
with open(out_path, mode = 'w') as f:
    f.write('#version: %s\n' %(VERSION))
    f.write('#input dose distribution (N): %s\n' %(ndose_path))
    f.write('#input dose distribution (G): %s\n' %(gdose_path))
    f.write('#constants : borofalan concentartion in normal tissue %.1f, RBE for gamma-ray %.1f, RBE for nitrogen dose %.1f, RBE for hydrogen dose and other neutron dose %.1f, CBE for tumor %.1f, CBE for mucosa %.1f, CBE for skin %.1f, CBE for other normal tissue %.2f, tumor/normal tissue ratio of borofalan concentration %.1f, maximum mucosa dose %.1f, maximum skin dose %.1f, maximum other normal tissue dose %.1f\n' %(C_NORMAL, RBE_GAMMA, RBE_THERMAL, RBE_FAST, CBE_TUMOR, CBE_MUCOSA, CBE_SKIN, CBE_NORMAL, TNR, LIMIT_MUCOSA, LIMIT_SKIN, LIMIT_NORMAL))
    f.write('#limitation : %s %4.1f Gy-eq\n' %(Limitation_part, Limitation_dose))
    f.write('#TT [h]: %10.4e\n' %(TT))
    f.write('#AD [mm]: %5.2f\n' %(AD))
    if type(AD20) is np.ndarray:
        f.write('#AD20 [mm]: %.2f\n' %(AD20))
    else:
        f.write('#AD20 [mm]: %s\n' %(AD20))
    if type(AD25) is np.ndarray:
        f.write('#AD25 [mm]: %.2f\n' %(AD25))
    else:
        f.write('#AD25 [mm]: %s\n' %(AD25))
    if type(AD30) is np.ndarray:
        f.write('#AD30 [mm]: %.2f\n' %(AD30))
    else:
        f.write('#AD30 [mm]: %s\n' %(AD30))
    f.write('#PTD [Gy-eq]: %.2f\n' %(PTD))
    f.write('#IND [Gy-eq mm]: %.2f\n' %(IND))
    f.write('#column1: Depth from phantom surface [mm]\n')
    f.write('#column2: Hydrogen dose [Gy-eq]\n')
    f.write('#column3: Nitrogen dose [Gy-eq]\n')
    f.write('#column4: Other neutron dose [Gy-eq]\n')
    f.write('#column5: Neutron dose [Gy-eq]\n')
    f.write('#column6: Primary gamma-ray dose [Gy]\n')
    f.write('#column7: Secondary gamma-ray dose [Gy]\n')
    f.write('#column8: Gamma-ray dose [Gy]\n')
    f.write('#column9: Boron dose in other normal tissue [Gy-eq]\n')
    f.write('#column10: Boron dose in skin [Gy-eq]\n')
    f.write('#column11: Boron dose in mucosa [Gy-eq]\n')
    f.write('#column12: Boron dose in tumor [Gy-eq]\n')
    f.write('#column13: Total dose in other normal tissue [Gy-eq]\n')
    f.write('#column14: Total dose in skin [Gy-eq]\n')
    f.write('#column15: Total dose in mucosa [Gy-eq]\n')
    f.write('#column16: Total dose in tumor [Gy-eq]\n')
    for i in range(TALLIES):
        f.write('%d\t' %(i))
        f.write('%10.4e\t' %(LIST_dose_hyd[i]))
        f.write('%10.4e\t' %(LIST_dose_nit[i]))
        f.write('%10.4e\t' %(LIST_dose_oth[i]))
        f.write('%10.4e\t' %(LIST_dose_n[i]))
        f.write('%10.4e\t' %(LIST_dose_g_pri[i]))
        f.write('%10.4e\t' %(LIST_dose_g_sec[i]))
        f.write('%10.4e\t' %(LIST_dose_g[i]))
        f.write('%10.4e\t' %(LIST_dose_boron_normal[i]))
        f.write('%10.4e\t' %(LIST_dose_boron_skin[i]))
        f.write('%10.4e\t' %(LIST_dose_boron_mucosa[i]))
        f.write('%10.4e\t' %(LIST_dose_boron_tumor[i]))
        f.write('%10.4e\t' %(LIST_dose_total_normal[i]))
        f.write('%10.4e\t' %(LIST_dose_total_skin[i]))
        f.write('%10.4e\t' %(LIST_dose_total_mucosa[i]))
        f.write('%10.4e\n' %(LIST_dose_total_tumor[i]))

#ACKNOWLEDGEMENTS
acknowledgements_1 = r"The predecessor of SiDE was devised by Dr. Nakao, an ex-researcher at Nagoya University."
acknowledgements_2 = r"This work was partially supported by JSPS KAKENHI Grant Number JP20J21215."

#REFERENCES
referrences_1 = r"[1] W. H. Sweet, Early history of development of boron neutron capture therapy of tumors, Journal of Neuro-Oncology 33: 19-26, 1997."
referrences_2 = r"[2] Y. Nakagawa, K. Pooh, T. Kobayashi, T. Kageji, S. Uyama, A. Matsumura, H. Kumada, Clinical review of the Japanese experience with boron neutron capture therapy and a proposed strategy using epithermal neutron beams, Journal of Neuro-Oncology 62: 87-99, 2003."
referrences_3 = r"[3] INTERNATIONAL ATOMIC ENERGY AGENCY, Current Status of Neutron Capture Therapy, IAEA-TECDOC-1223, IAEA, Vienna (2001). https://www.iaea.org/publications/6168/current-status-of-neutron-capture-therapy"
referrences_4 = r"[4] K. Hirose, A. Konno, J. Hiratsuka, S. Yoshimoto, T. Kato, K. Ono, N. Otsuki, J. Hatazawa, H. Tanaka, K. Takayama, H. Wada, M. Suzuki, M. Sato, H. Yamaguchi, I Seto, Y. Ueki, S. Iketani, S. Imai, T. Nakamura, T. Ono, H. Endo, Y. Azami, Y. Kikuchi, M. Murakami, Y. Takai, Boron neutron capture therapy using cyclotron-based epithermal neutron source and borofalan (10B) for recurrent or locally advanced head and neck cancer (JHN002): An open-label phase II trial, Radiotherapy and Oncology 155 (2021) 182-187."
referrences_5 = r"[5] D. R. White, R. V. Griffith, I. J. Wilson, Report 46, Journal of the International Commission on Radiation Units and Measurements, Volume os24, Issue 1, 28 February 1992.　https://doi.org/10.1093/jicru/os24.1.Report46"
referrences_6 = r"[6] Y. Sakurai, T. Kobayashi, Spectrum evaluation at the filter-modified neutron irradiation field for neutron capture therapy in Kyoto University Research Reactor, Nuclear Instruments and Methods in Physics Research A 531 (2004) 585–595. https://doi.org/10.1016/j.nima.2004.05.084"

#ending messages
print('all calculation process has finished!')
print('\n-----Acknowledgements-----\n')
print(acknowledgements_1)
print(acknowledgements_2)
print('\n--------References--------\n')
print(referrences_1)
print(referrences_2)
print(referrences_3)
print(referrences_4)
print(referrences_5)
print(referrences_6)
print('\nSiDE %s' %(VERSION))
