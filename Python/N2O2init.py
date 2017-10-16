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

matname1_7000 = 'NOZ1_1_7000_OSC2_EX3_REC1_'
matname10_7000 = 'NOZ1_10_7000_OSC2_EX3_REC1_'
matname100_7000 = 'NOZ1_100_7000_OSC2_EX3_REC1_'

plt.rc('font', **{'family':'serif'})
plt.rc('text', usetex=True)
plt.rc('text.latex', unicode=True)
plt.rc('text.latex', preamble=r"\usepackage[utf8]{inputenc}")
plt.rc('text.latex', preamble=r"\usepackage[russian]{babel}")

data1_7000 = scipy.io.loadmat('../Matlab/MAT/' + matname1_7000 + '.mat')
data10_7000 = scipy.io.loadmat('../Matlab/MAT/' + matname10_7000 + '.mat')
data100_7000 = scipy.io.loadmat('../Matlab/MAT/' + matname100_7000 + '.mat')

n_N217 = data1_7000['n_N2']
n_N2107 = data10_7000['n_N2']
n_N21007 = data100_7000['n_N2']

n_O217 = data1_7000['n_O2']
n_O2107 = data10_7000['n_O2']
n_O21007 = data100_7000['n_O2']

X17 = data1_7000['X']
X107 = data10_7000['X']
X1007 = data100_7000['X']

plt.figure()

plt.plot(X17, n_N217, 'k', label=r'$p^* = 1$ atm')
plt.plot(X107, n_N2107, 'b', label=r'$p^* = 10$ atm')
plt.plot(X1007, n_N21007, 'r', label=r'$p^* = 100$ atm')

plt.title(r'$\mathrm{N}_2$', fontsize=24)
plt.xlabel(r'$x/r^*$', fontsize=16)
plt.ylabel(r'$n_{\mathrm{N}_2}/n$', rotation=0, fontsize=16)
plt.gca().yaxis.set_label_coords(0.07, 1)
plt.gca().xaxis.set_label_coords(1.05, 0.05)
plt.legend(loc='upper right', bbox_to_anchor=(1, 0.9), fontsize=18)
plt.tick_params(labelsize=14)
plt.gca().set_xlim([0, 5])
plt.gca().set_ylim([n_N21007.min(), n_N21007.max()])

save(name='N2_init', fmt='pdf')
plt.tight_layout()
plt.show()
plt.close()


plt.figure()

plt.plot(X17, n_O217, 'k', label=r'$p^* = 1$ atm')
plt.plot(X107, n_O2107, 'b', label=r'$p^* = 10$ atm')
plt.plot(X1007, n_O21007, 'r', label=r'$p^* = 100$ atm')

plt.title(r'$\mathrm{O}_2$', fontsize=24)
plt.xlabel(r'$x/r^*$', fontsize=16)
plt.ylabel(r'$n_{\mathrm{O}_2}/n$', rotation=0, fontsize=16)
plt.gca().yaxis.set_label_coords(0.07, 1)
plt.gca().xaxis.set_label_coords(1.05, 0.05)
plt.legend(loc='upper right', bbox_to_anchor=(1, 0.9), fontsize=18)
plt.tick_params(labelsize=14)
plt.gca().set_xlim([0, 5])
plt.gca().set_ylim([n_O21007.min(), n_O21007.max()])


save(name='O2_init', fmt='pdf')
plt.tight_layout()
plt.show()

plt.close()
