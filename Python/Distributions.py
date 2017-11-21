# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import scipy.io
import numpy as np
import os


def save(name='', fmt='png'):
    pwd = os.getcwd()
    path = './pic/{}'.format(fmt)
    if not os.path.exists(path):
        os.mkdir(path)
    os.chdir(path)
    plt.savefig('{}.{}'.format(name, fmt), fmt='png')#, bbox_inches='tight')
    os.chdir(pwd)

labelsize = 22
legendsize = 14
matname = 'NOZ3_1_7000_OSC2_EX3_REC1_'

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
data = scipy.io.loadmat('../Matlab/MAT/' + matname + '.mat')
#data1T = scipy.io.loadmat('./MAT/1T_100_7000_distr.mat')

xr = [0, 1, 2, 3, 5, 10, 15, 20, 25, 30, 40, 50]
#xr = [0, 2, 5, 15, 25, 40, 50]

n_i_N2 = data['u_N2']
e_i_N2 = data['e_i_N2']
X = data['X']
#n_i_N2_1T = data1T['u_N2_1t']
#X_1T = data1T['X_1t']

i_N2 = range(0, 48)
#i_N2 = range(0, 34)
plt.figure()
plt.set_cmap('jet_r')
plt.gca().set_prop_cycle(plt.cycler('color', plt.get_cmap()(np.linspace(0.0, 1.0, len(xr)))))

for i in range(0, len(xr)):
    ch = str(xr[i])
    l = r'$x/{r^*} ={ }$'
    lab = l + ch
   # plt.semilogy(e_i_N2[:, 0], n_i_N2[:, i], label=lab)
    plt.semilogy(i_N2, n_i_N2[i, :], label=lab)

#for i in range(0, n_i_N2.shape[1]):
    # plt.semilogy(e_i_N2[:, 0], n_i_N2_1T[i, :], '--')
#    plt.semilogy(i_N2, n_i_N2_1T[i, :], '--')

#plt.xlabel(r'$\varepsilon_i^{\mathrm{N}_2}$, eV', fontsize=16)
plt.xlabel(r'$i$', fontsize=labelsize)
plt.ylabel(r'$n_{\mathrm{N}_2, i}/n$', rotation=0, fontsize=labelsize)
plt.title(r'$\mathrm{N}_2$', fontsize=24)
plt.gca().yaxis.set_label_coords(0.07, 1)
plt.gca().xaxis.set_label_coords(1.03, 0.05)
#plt.legend(bbox_to_anchor=(1.005, 1), loc=2, borderaxespad=0., fontsize=12)
plt.legend(loc=0, fontsize=legendsize)
plt.tick_params(labelsize=labelsize)
#plt.gca().set_xlim([min(e_i_N2), max(e_i_N2)])
plt.gca().set_xlim([min(i_N2), max(i_N2)])
plt.gca().set_ylim([n_i_N2.min(), n_i_N2.max()])
plt.locator_params(axis='y', numticks=10)

save(name='n_i_N2_' + matname, fmt='pdf')
plt.tight_layout()
plt.show()

plt.close()

n_i_O2 = data['u_O2']
i_O2 = range(0, 37)
#i_O2 = range(0, 27)
plt.figure()
plt.set_cmap('jet_r')
plt.gca().set_prop_cycle(plt.cycler('color', plt.get_cmap()(np.linspace(0.0, 1.0, len(xr)))))

for i in range(0, len(xr)):
    ch = str(xr[i])
    l = r'$x/{r^*} ={ }$'
    lab = l + ch
   # plt.semilogy(e_i_N2[:, 0], n_i_N2[:, i], label=lab)
    plt.semilogy(i_O2, n_i_O2[i, :], label=lab)

#for i in range(0, n_i_N2.shape[1]):
    # plt.semilogy(e_i_N2[:, 0], n_i_N2_1T[i, :], '--')
#    plt.semilogy(i_N2, n_i_N2_1T[i, :], '--')

#plt.xlabel(r'$\varepsilon_i^{\mathrm{N}_2}$, eV', fontsize=16)
plt.xlabel(r'$i$', fontsize=labelsize)
plt.ylabel(r'$n_{\mathrm{O}_2, i}/n$', rotation=0, fontsize=labelsize)
#plt.title(r'$\mathrm{O}_2$', fontsize=24)
plt.gca().yaxis.set_label_coords(0.08, 1)
plt.gca().xaxis.set_label_coords(1.03, 0.05)
#plt.legend(bbox_to_anchor=(1.005, 1), loc=2, borderaxespad=0., fontsize=12)
plt.legend(loc=3, fontsize=legendsize)
plt.tick_params(labelsize=labelsize)
#plt.gca().set_xlim([min(e_i_N2), max(e_i_N2)])
plt.gca().set_xlim([min(i_O2), max(i_O2)])
plt.gca().set_ylim([n_i_O2.min(), n_i_O2.max()])
plt.locator_params(axis='y', numticks=10)

save(name='n_i_O2_' + matname, fmt='pdf')
plt.tight_layout()
plt.show()

plt.close()
