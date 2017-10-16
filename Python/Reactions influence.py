# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import scipy.io
import numpy as np
import os

plt.rc('font', **{'family': 'serif'})
plt.rc('text', usetex=True)
plt.rc('text.latex', unicode=True)
plt.rc('text.latex', preamble=r"\usepackage[utf8]{inputenc}")
plt.rc('text.latex', preamble=r"\usepackage[russian]{babel}")

def save(name='', fmt='png'):
    pwd = os.getcwd()
    path = './pic/{}'.format(fmt)
    if not os.path.exists(path):
        os.mkdir(path)
    os.chdir(path)
    plt.savefig('{}.{}'.format(name, fmt), fmt='png')
    os.chdir(pwd)

matname = 'NOZ1_1_7000_OSC2_EX3_REC1_'
matname_ex = 'NOZ1_1_7000_OSC2_EX5_REC1_'
matname_rec = 'NOZ1_100_7000_OSC2_EX3_REC0_'

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
data = scipy.io.loadmat('../Matlab/MAT/' + matname + '.mat')
data_ex = scipy.io.loadmat('../Matlab/MAT/' + matname_ex + '.mat')
data_rec = scipy.io.loadmat('../Matlab/MAT/' + matname_rec + '.mat')
dataN2 = scipy.io.loadmat('../Matlab/MAT/NOZ1_sp1_1_100_7000_OSC2_.mat')
dataO2 = scipy.io.loadmat('../Matlab/MAT/NOZ1_sp2_1_100_7000_OSC2_.mat')

X = data['X']
X_ex = data_ex['X']
X_rec = data_rec['X']
n_N2 = data['n_N2']
n_N2_ex = data_ex['n_N2']
n_N2_rec = data_rec['n_N2']
n_O2 = data['n_O2']
n_O2_ex = data_ex['n_O2']
n_O2_rec = data_rec['n_O2']
n_NO = data['n_NO']
n_NO_ex = data_ex['n_NO']
n_N = data['n_N']
n_N_ex = data_ex['n_N']
n_O = data['n_O']
n_O_ex = data_ex['n_O']

# all species

plt.figure()
plt.semilogy(X, n_N2, label=r'$\mathrm{N_2}$')
plt.semilogy(X, n_O2, label=r'$\mathrm{O_2}$')
plt.semilogy(X, n_NO, label=r'$\mathrm{NO}$')
plt.semilogy(X, n_N, label=r'$\mathrm{N}$')
plt.semilogy(X, n_O, label=r'$\mathrm{O}$')
plt.xlabel(r'$x/r^*$', fontsize=16)
plt.ylabel(r'$n_{c}/n$', rotation=0, fontsize=16)
plt.gca().yaxis.set_label_coords(0.05, 1)
plt.gca().xaxis.set_label_coords(1.05, 0.05)
#plt.legend(loc=0, fontsize=22)
plt.annotate(r'$\mathrm{N_2}$', xy=(4.5, 0.71), fontsize=22)
plt.annotate(r'$\mathrm{O_2}$', xy=(4.5, 0.21), fontsize=22)
plt.annotate(r'$\mathrm{NO}$', xy=(4.5, 0.0078), fontsize=22)
plt.annotate(r'$\mathrm{N}$', xy=(4.5, 0.0019), fontsize=22)
plt.annotate(r'$\mathrm{O}$', xy=(4.5, 0.0004), fontsize=22)
plt.tick_params(labelsize=14)
plt.gca().set_xlim(0, 5)
#plt.gca().set_ylim([n_i_N2.min(), n_i_N2.max()])
plt.locator_params(axis='y', numticks=20)

save(name='all_species_' + matname, fmt='pdf')
plt.tight_layout()
plt.show()
plt.close()

# all species without exchange

plt.figure()
plt.semilogy(X_ex, n_N2_ex, label=r'$\mathrm{N_2}$')
plt.semilogy(X_ex, n_O2_ex, label=r'$\mathrm{O_2}$')
plt.semilogy(X_ex, n_NO_ex, label=r'$\mathrm{NO}$')
plt.semilogy(X_ex, n_N_ex, label=r'$\mathrm{N}$')
plt.semilogy(X_ex, n_O_ex, label=r'$\mathrm{O}$')
plt.xlabel(r'$x/r^*$', fontsize=16)
plt.ylabel(r'$n_{c}/n$', rotation=0, fontsize=16)
plt.gca().yaxis.set_label_coords(0.05, 1)
plt.gca().xaxis.set_label_coords(1.05, 0.05)
#plt.legend(loc=0, fontsize=22)
plt.annotate(r'$\mathrm{N_2}$', xy=(4.5, 0.82), fontsize=22)
plt.annotate(r'$\mathrm{O_2}$', xy=(4.5, 0.12), fontsize=22)
plt.annotate(r'$\mathrm{NO}$', xy=(4.5, 0.0072), fontsize=22)
plt.annotate(r'$\mathrm{N}$', xy=(4.5, 1.5e-5), fontsize=22)
plt.annotate(r'$\mathrm{O}$', xy=(4.5, 6.4e-9), fontsize=22)
plt.tick_params(labelsize=14)
plt.gca().set_xlim(0, 5)
#plt.gca().set_ylim([n_i_N2.min(), n_i_N2.max()])
plt.locator_params(axis='y', numticks=20)

save(name='all_species_' + matname_ex, fmt='pdf')
plt.tight_layout()
plt.show()
plt.close()

# N2

matname = 'NOZ1_100_7000_OSC2_EX3_REC1_'
matname_ex = 'NOZ1_100_7000_OSC2_EX5_REC1_'
data = scipy.io.loadmat('../Matlab/MAT/' + matname + '.mat')
data_ex = scipy.io.loadmat('../Matlab/MAT/' + matname_ex + '.mat')
n_N2 = data['n_N2']
n_O2 = data['n_O2']
n_N2_ex = data_ex['n_N2']
n_O2_ex = data_ex['n_O2']

plt.figure()
plt.plot(X, n_N2, 'k', label=u'все реакции')
plt.plot(X_ex, n_N2_ex, 'b', label=u'без обмена')
plt.plot(X_rec, n_N2_rec, 'r',  label=u'без рекомбинации')
plt.xlabel(r'$x/r^*$', fontsize=16)
plt.ylabel(r'$n_{\mathrm{N}_2}/n$', rotation=0, fontsize=16)
plt.gca().yaxis.set_label_coords(0.05, 1)
plt.gca().xaxis.set_label_coords(1.05, 0.05)
plt.legend(loc='upper right', bbox_to_anchor=(1, 0.5), fontsize=18)
plt.tick_params(labelsize=14)
plt.gca().set_xlim(0, 5)
#plt.gca().set_ylim([n_i_N2.min(), n_i_N2.max()])
#plt.locator_params(axis='y', numticks=20)

save(name='N2_react', fmt='pdf')
plt.tight_layout()
plt.show()
plt.close()

# O2

plt.figure()
plt.plot(X, n_O2, 'k', label=u'все реакции')
plt.plot(X_ex, n_O2_ex, 'b', label=u'без обмена')
plt.plot(X_rec, n_O2_rec, 'r',  label=u'без рекомбинации')
plt.xlabel(r'$x/r^*$', fontsize=16)
plt.ylabel(r'$n_{\mathrm{O}_2}/n$', rotation=0, fontsize=16)
plt.gca().yaxis.set_label_coords(0.05, 1)
plt.gca().xaxis.set_label_coords(1.05, 0.05)
plt.legend(loc='upper right', bbox_to_anchor=(1, 0.5), fontsize=18)
plt.tick_params(labelsize=14)
plt.gca().set_xlim(0, 5)
#plt.gca().set_ylim([n_i_N2.min(), n_i_N2.max()])
#plt.locator_params(axis='y', numticks=20)

save(name='O2_react', fmt='pdf')
plt.tight_layout()
plt.show()
plt.close()

# T

T = data['T']
T_ex = data_ex['T']
T_rec = data_rec['T']
T_N2 = dataN2['T']
T_O2 = dataO2['T']
X_N2 = dataN2['X']
X_O2 = dataO2['X']

plt.figure()
plt.plot(X, T, 'k', label=u'все реакции')
plt.plot(X_ex, T_ex, 'y--', label=u'без обмена')
plt.plot(X_rec, T_rec, 'r',  label=u'без рекомбинации')
plt.plot(X_N2, T_N2, 'g',  label=r'$\mathrm{N_2/N}$')
plt.plot(X_O2, T_O2, 'b',  label=r'$\mathrm{O_2/O}$')
plt.xlabel(r'$x/r^*$', fontsize=16)
plt.ylabel(r'$T$, K', rotation=0, fontsize=16)
plt.gca().yaxis.set_label_coords(0.05, 1)
plt.gca().xaxis.set_label_coords(1.05, 0.05)
plt.legend(loc=0, fontsize=18)
plt.tick_params(labelsize=14)
plt.gca().set_xlim(0, 50)
#plt.gca().set_ylim([n_i_N2.min(), n_i_N2.max()])
plt.locator_params(axis='y', numticks=20)

save(name='T_react', fmt='pdf')
plt.tight_layout()
plt.show()
plt.close()

plt.figure()
plt.plot(X, T, 'k', label=u'все реакции')
plt.plot(X_ex, T_ex, 'y--', label=u'без обмена')
plt.plot(X_rec, T_rec, 'r',  label=u'без рекомбинации')

plt.xlabel(r'$x/r^*$', fontsize=16)
plt.ylabel(r'$T$, K', rotation=0, fontsize=16)
plt.gca().yaxis.set_label_coords(0.05, 1)
plt.gca().xaxis.set_label_coords(1.05, 0.05)
plt.legend(loc=0, fontsize=18)
plt.tick_params(labelsize=14)
plt.gca().set_xlim(0, 0.01)
plt.gca().set_ylim([6200, 7000])
plt.locator_params(axis='y', numticks=20)

save(name='T_react3', fmt='pdf')
plt.tight_layout()
plt.show()
plt.close()
