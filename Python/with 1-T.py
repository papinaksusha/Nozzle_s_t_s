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

matname = 'NOZ1_100_7000_OSC2_EX3_REC1_'

labelsize = 22
legendsize = 18

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

T = data['T']
v = data['v']
T_1t = data1T['T_1t']
v_1t = data1T['v_1t']


epsV = max(abs(v - v_1t)/v)*100
epsT = max(abs(T - T_1t)/T)*100
print('epsV', epsV)
print('epsT', epsT)

# # N2

col = len(xr)/2 + 1;
plt.figure()
plt.set_cmap('jet_r')
plt.gca().set_prop_cycle(plt.cycler('color', plt.get_cmap()(np.linspace(0.0, 1.0, col))))

for i in range(0, len(xr), 2):
    ch = str(xr[i])
    l = r'$x/{r^*} ={ }$'
    lab = l + ch
   # plt.semilogy(e_i_N2[:, 0], n_i_N2[:, i], label=lab)
    plt.semilogy(i_N2, n_i_N2[i, :], label=lab)

ch = str(xr[-1])
l = r'$x/{r^*} ={ }$'
lab = l + ch
# plt.semilogy(e_i_N2[:, 0], n_i_N2[:, i], label=lab)
plt.semilogy(i_N2, n_i_N2[-1, :], label=lab)

for i in range(0, len(xr), 2):
    #plt.semilogy(e_i_N2[:, 0], n_i_N2_1T[i, :], '--')
    plt.semilogy(i_N2, n_i_N2_1T[i, :], '--')

plt.semilogy(i_N2, n_i_N2_1T[-1, :], '--')
#plt.xlabel(r'$\varepsilon_i^{\mathrm{N}_2}$, eV', fontsize=16)
plt.xlabel(r'$i$', fontsize=labelsize)
plt.ylabel(r'$n_{\mathrm{N}_2, i}/n$', rotation=0, fontsize=labelsize)
#plt.title(r'$\mathrm{N}_2$', fontsize=24)
plt.gca().yaxis.set_label_coords(0.08, 1)
plt.gca().xaxis.set_label_coords(1.03, 0.05)
#plt.legend(bbox_to_anchor=(1.005, 1), loc=2, borderaxespad=0., fontsize=12)
plt.legend(loc='lower left', fontsize=legendsize)
plt.tick_params(labelsize=labelsize)
#plt.gca().set_xlim([min(e_i_N2), max(e_i_N2)])
plt.gca().set_xlim([min(i_N2), max(i_N2)])
#plt.gca().set_ylim([n_i_N2.min(), n_i_N2.max()])
plt.gca().set_ylim([1e-30, n_i_N2.max()])
plt.locator_params(axis='y', numticks=7)

save(name='n_i_N2_1T_' + matname, fmt='pdf')
plt.tight_layout()
plt.show()

plt.close()

n_i_O2 = data['u_O2']
n_i_O2_1T = data1T['u_O2_1t']
i_O2 = range(0, 37)

# O2

plt.figure()
plt.set_cmap('jet_r')
plt.gca().set_prop_cycle(plt.cycler('color', plt.get_cmap()(np.linspace(0.0, 1.0, col))))

for i in range(0, len(xr), 2):
    ch = str(xr[i])
    l = r'$x/{r^*} ={ }$'
    lab = l + ch
   # plt.semilogy(e_i_N2[:, 0], n_i_N2[:, i], label=lab)
    plt.semilogy(i_O2, n_i_O2[i, :], label=lab)

ch = str(xr[-1])
l = r'$x/{r^*} ={ }$'
lab = l + ch
# plt.semilogy(e_i_N2[:, 0], n_i_N2[:, i], label=lab)
plt.semilogy(i_O2, n_i_O2[-1, :], label=lab)

for i in range(0, len(xr), 2):
    # plt.semilogy(e_i_N2[:, 0], n_i_N2_1T[i, :], '--')
    plt.semilogy(i_O2, n_i_O2_1T[i, :], '--')

plt.semilogy(i_O2, n_i_O2_1T[-1, :], '--')
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
plt.gca().set_ylim([1e-15, n_i_O2.max()])
plt.locator_params(axis='y', numticks=7)

save(name='n_i_O2_1T_' + matname, fmt='pdf')
plt.tight_layout()
plt.show()

plt.close()

#### FROM X/R^*

# N2

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

plt.xlabel(r'$x/r^*$', fontsize=labelsize)
plt.ylabel(r'$n_{\mathrm{N}_2, i}/n$', rotation=0, fontsize=labelsize)
#plt.title(r'$\mathrm{N}_2$', fontsize=24)
plt.annotate(r'$i = 0$', xy=(43, 1), fontsize=legendsize)
plt.annotate(r'$i = 1$', xy=(43, 0.0793), fontsize=legendsize)
plt.annotate(r'$i = 5$', xy=(43, 0.00021), fontsize=legendsize)
plt.annotate(r'$i = 10$', xy=(43, 1.4e-8), fontsize=legendsize)
plt.annotate(r'$i = 20$', xy=(43, 4.6e-17), fontsize=legendsize)
plt.annotate(r'$i = 30$', xy=(43, 8.88e-23), fontsize=legendsize)
plt.annotate(r'$i = 40$', xy=(43, 7.2e-25), fontsize=legendsize)
plt.annotate(r'$i = 45$', xy=(43, 6.9e-27), fontsize=legendsize)
plt.gca().yaxis.set_label_coords(0.08, 1)
plt.gca().xaxis.set_label_coords(1.07, 0.07)
plt.tick_params(labelsize=labelsize)
plt.gca().set_xlim([0, 50])
plt.gca().set_ylim([1e-30, 1])
plt.locator_params(axis='y', numticks=7)

save(name='n_i_N2_1T_x_' + matname, fmt='pdf')
plt.tight_layout()
plt.show()
plt.close()

ii_O2 = [0, 1, 5, 10, 20, 30, 32]

# O2
plt.figure()
plt.set_cmap('rainbow')
plt.gca().set_prop_cycle(plt.cycler('color', plt.get_cmap()(np.linspace(0.0, 1.0, len(ii_O2)))))

for i in range(0, len(ii_O2)):
    plt.semilogy(X, n_i_O2[:, ii_O2[i]])

for i in range(0, len(ii_O2)):
    plt.semilogy(X_1T, n_i_O2_1T[:, ii_O2[i]], '--')

plt.xlabel(r'$x/r^*$', fontsize=labelsize )
plt.ylabel(r'$n_{\mathrm{O}_2, i}/n$', rotation=0, fontsize=labelsize)
#plt.title(r'$\mathrm{O}_2$', fontsize=24)
plt.annotate(r'$i = 0$', xy=(43, 0.11), fontsize=legendsize)
plt.annotate(r'$i = 1$', xy=(43, 0.045), fontsize=legendsize)
plt.annotate(r'$i = 5$', xy=(43, 0.00135), fontsize=legendsize)
plt.annotate(r'$i = 10$', xy=(43, 2.8e-5), fontsize=legendsize)
plt.annotate(r'$i = 20$', xy=(43, 1e-5), fontsize=legendsize)
plt.annotate(r'$i = 30$', xy=(43, 2.7327e-6), fontsize=legendsize)
plt.annotate(r'$i = 32$', xy=(43, 3.9902e-7), fontsize=legendsize)
plt.gca().yaxis.set_label_coords(0.08, 1)
plt.gca().xaxis.set_label_coords(1.07, 0.07)
plt.tick_params(labelsize=labelsize)
plt.gca().set_xlim([0, 50])
plt.gca().set_ylim([1e-10, 1])
plt.locator_params(axis='y', numticks=7)

save(name='n_i_O2_1T_x_' + matname, fmt='pdf')
plt.tight_layout()
plt.show()

plt.close()
