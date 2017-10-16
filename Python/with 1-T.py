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
    plt.savefig('{}.{}'.format(name, fmt), fmt='png')
    os.chdir(pwd)

matname = 'NOZ1_100_7000_OSC2_EX3_REC1_'

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
data = scipy.io.loadmat('../Matlab/MAT/' + matname + '.mat')
data1T = scipy.io.loadmat('../Matlab/MAT/1T_100_7000_distr.mat')

xr = [0, 1, 2, 3, 5, 10, 15, 20, 25, 30, 40, 50]

n_i_N2 = data['u_N2']
#e_i_N2 = data['e_i_N2']
X = data['X']
n_i_N2_1T = data1T['u_N2_1t']
X_1T = data1T['X_1t']
i_N2 = range(0, 48)
plt.figure()
plt.set_cmap('jet_r')
plt.gca().set_prop_cycle(plt.cycler('color', plt.get_cmap()(np.linspace(0.0, 1.0, len(xr)))))

for i in range(0, len(xr)):
    ch = str(xr[i])
    l = r'$x/{r^*} ={ }$'
    lab = l + ch
   # plt.semilogy(e_i_N2[:, 0], n_i_N2[:, i], label=lab)
    plt.semilogy(i_N2, n_i_N2[i, :], label=lab)

for i in range(0, len(xr)):
    #plt.semilogy(e_i_N2[:, 0], n_i_N2_1T[i, :], '--')
    plt.semilogy(i_N2, n_i_N2_1T[i, :], '--')

#plt.xlabel(r'$\varepsilon_i^{\mathrm{N}_2}$, eV', fontsize=16)
plt.xlabel(r'$i$', fontsize=16)
plt.ylabel(r'$n_{\mathrm{N}_2, i}/n$', rotation=0, fontsize=16)
plt.title(r'$\mathrm{N}_2$', fontsize=24)
plt.gca().yaxis.set_label_coords(0.05, 1)
plt.gca().xaxis.set_label_coords(1.03, 0.05)
#plt.legend(bbox_to_anchor=(1.005, 1), loc=2, borderaxespad=0., fontsize=12)
plt.legend(loc=0, fontsize=12)
plt.tick_params(labelsize=14)
#plt.gca().set_xlim([min(e_i_N2), max(e_i_N2)])
plt.gca().set_xlim([min(i_N2), max(i_N2)])
#plt.gca().set_ylim([n_i_N2.min(), n_i_N2.max()])
plt.gca().set_ylim([1e-30, n_i_N2.max()])
plt.locator_params(axis='y', numticks=20)

save(name='n_i_N2_1T_' + matname, fmt='pdf')
plt.tight_layout()
plt.show()

plt.close()

n_i_O2 = data['u_O2']
n_i_O2_1T = data1T['u_O2_1t']
i_O2 = range(0, 37)
plt.figure()
plt.set_cmap('jet_r')
plt.gca().set_prop_cycle(plt.cycler('color', plt.get_cmap()(np.linspace(0.0, 1.0, len(xr)))))

for i in range(0, len(xr)):
    ch = str(xr[i])
    l = r'$x/{r^*} ={ }$'
    lab = l + ch
   # plt.semilogy(e_i_N2[:, 0], n_i_N2[:, i], label=lab)
    plt.semilogy(i_O2, n_i_O2[i, :], label=lab)

for i in range(0, len(xr)):
    # plt.semilogy(e_i_N2[:, 0], n_i_N2_1T[i, :], '--')
    plt.semilogy(i_O2, n_i_O2_1T[i, :], '--')

#plt.xlabel(r'$\varepsilon_i^{\mathrm{N}_2}$, eV', fontsize=16)
plt.xlabel(r'$i$', fontsize=16)
plt.ylabel(r'$n_{\mathrm{O}_2, i}/n$', rotation=0, fontsize=16)
plt.title(r'$\mathrm{O}_2$', fontsize=24)
plt.gca().yaxis.set_label_coords(0.05, 1)
plt.gca().xaxis.set_label_coords(1.03, 0.05)
#plt.legend(bbox_to_anchor=(1.005, 1), loc=2, borderaxespad=0., fontsize=12)
plt.legend(loc=3, fontsize=12)
plt.tick_params(labelsize=14)
#plt.gca().set_xlim([min(e_i_N2), max(e_i_N2)])
plt.gca().set_xlim([min(i_O2), max(i_O2)])
plt.gca().set_ylim([1e-15, n_i_O2.max()])
plt.locator_params(axis='y', numticks=20)

save(name='n_i_O2_1T_' + matname, fmt='pdf')
plt.tight_layout()
plt.show()

plt.close()

#### FROM X/R^*

n_i_N2_1T = data1T['n_i_N2_1t']
n_i_O2_1T = data1T['n_i_O2_1t']
n_i_N2 = data['n_i_N2']
n_i_O2 = data['n_i_O2']

ii_N2 = [0, 1, 5, 10, 20, 30, 40, 45]
plt.figure()
plt.set_cmap('rainbow')
plt.gca().set_prop_cycle(plt.cycler('color', plt.get_cmap()(np.linspace(0.0, 1.0, len(ii_N2)))))

for i in range(0, len(ii_N2)):
    plt.semilogy(X, n_i_N2[:, ii_N2[i]])

for i in range(0, len(ii_N2)):
    plt.semilogy(X_1T, n_i_N2_1T[:, ii_N2[i]], '--')

plt.xlabel(r'$x/r^*$', fontsize=16)
plt.ylabel(r'$n_{\mathrm{N}_2, i}/n$', rotation=0, fontsize=16)
plt.title(r'$\mathrm{N}_2$', fontsize=24)
plt.annotate(r'$i = 0$', xy=(45, 1), fontsize=12)
plt.annotate(r'$i = 1$', xy=(45, 0.0793), fontsize=12)
plt.annotate(r'$i = 5$', xy=(45, 0.00021), fontsize=12)
plt.annotate(r'$i = 10$', xy=(45, 1.4e-8), fontsize=12)
plt.annotate(r'$i = 20$', xy=(45, 4.6e-17), fontsize=12)
plt.annotate(r'$i = 30$', xy=(45, 7.2e-25), fontsize=12)
plt.annotate(r'$i = 40$', xy=(45, 1.7e-26), fontsize=12)
plt.annotate(r'$i = 45$', xy=(45, 4e-28), fontsize=12)
plt.gca().yaxis.set_label_coords(0.05, 1)
plt.gca().xaxis.set_label_coords(1.05, 0.05)
plt.tick_params(labelsize=14)
plt.gca().set_xlim([0, 50])
plt.gca().set_ylim([1e-30, 1])
plt.locator_params(axis='y', numticks=20)

save(name='n_i_N2_1T_x_' + matname, fmt='pdf')
plt.tight_layout()
plt.show()
plt.close()

ii_O2 = [0, 1, 5, 10, 20, 30, 32]
plt.figure()
plt.set_cmap('rainbow')
plt.gca().set_prop_cycle(plt.cycler('color', plt.get_cmap()(np.linspace(0.0, 1.0, len(ii_O2)))))

for i in range(0, len(ii_O2)):
    plt.semilogy(X, n_i_O2[:, ii_O2[i]])

for i in range(0, len(ii_O2)):
    plt.semilogy(X_1T, n_i_O2_1T[:, ii_O2[i]], '--')

plt.xlabel(r'$x/r^*$', fontsize=16)
plt.ylabel(r'$n_{\mathrm{O}_2, i}/n$', rotation=0, fontsize=16)
plt.title(r'$\mathrm{O}_2$', fontsize=24)
plt.annotate(r'$i = 0$', xy=(45, 0.11), fontsize=12)
plt.annotate(r'$i = 1$', xy=(45, 0.045), fontsize=12)
plt.annotate(r'$i = 5$', xy=(45, 0.00135), fontsize=12)
plt.annotate(r'$i = 10$', xy=(45, 2.8e-5), fontweight='bold', fontsize=12)
plt.annotate(r'$i = 20$', xy=(45, 1e-5), fontsize=12)
plt.annotate(r'$i = 30$', xy=(45, 3.7327e-6), fontsize=12)
plt.annotate(r'$i = 32$', xy=(45, 5.9902e-7), fontsize=12)
plt.gca().yaxis.set_label_coords(0.05, 1)
plt.gca().xaxis.set_label_coords(1.05, 0.05)
plt.tick_params(labelsize=14)
plt.gca().set_xlim([0, 50])
plt.gca().set_ylim([1e-10, 1])
plt.locator_params(axis='y', numticks=20)

save(name='n_i_O2_1T_x_' + matname, fmt='pdf')
plt.tight_layout()
plt.show()

plt.close()
