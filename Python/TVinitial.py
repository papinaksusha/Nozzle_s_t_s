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
matname1_5000 = 'NOZ1_1_5000_OSC2_EX3_REC1_'
matname100_7000 = 'NOZ1_100_7000_OSC2_EX3_REC1_'
matname100_5000 = 'NOZ1_100_5000_OSC2_EX3_REC1_'

plt.rc('font',**{'family':'serif'})
plt.rc('text', usetex=True)
plt.rc('text.latex',unicode=True)
plt.rc('text.latex',preamble=r"\usepackage[utf8]{inputenc}")
plt.rc('text.latex',preamble=r"\usepackage[russian]{babel}")

data1_7000 = scipy.io.loadmat('../Matlab/MAT/' + matname1_7000 + '.mat')
data1_5000 = scipy.io.loadmat('../Matlab/MAT/' + matname1_5000 + '.mat')
data100_7000 = scipy.io.loadmat('../Matlab/MAT/' + matname100_7000 + '.mat')
data100_5000 = scipy.io.loadmat('../Matlab/MAT/' + matname100_5000 + '.mat')


V17 = data1_7000['v']
V1007 = data100_7000['v']
V15 = data1_5000['v']
V1005 = data100_5000['v']

T17 = data1_7000['T']
T1007 = data100_7000['T']
T15 = data1_5000['T']
T1005 = data100_5000['T']

X17 = data1_7000['X']
X1007 = data100_7000['X']
X15 = data1_5000['X']
X1005 = data100_5000['X']

plt.figure()

plt.plot(X17, T17, 'r', label=r'$T^* = 7000$ K, $p^* = 1$ atm')
plt.plot(X1007, T1007, 'r--', label=r'$T^* = 7000$ K, $p^* = 100$ atm')
plt.plot(X15, T15, 'b', label=r'$T^* = 5000$ K, $p^* = 1$ atm$')
plt.plot(X1005, T1005, 'b--', label=r'$T^* = 5000$ K, $p^* = 100$ atm')

plt.xlabel(r'$x/r^*$', fontsize=16)
plt.ylabel(r'$T$, K', rotation=0, fontsize=16)
plt.gca().yaxis.set_label_coords(0.05, 1)
plt.gca().xaxis.set_label_coords(1.05, 0.05)
plt.legend(loc=1, fontsize=18)
plt.tick_params(labelsize=14)
plt.gca().set_xlim([min(X17), max(X17)])
plt.gca().set_ylim([0, 7000])


save(name='T_init', fmt='pdf')
plt.tight_layout()
plt.show()

plt.close()

plt.figure()

plt.plot(X17, V17, 'r', label=r'$T^* = 7000$ K, $p^* = 1$ atm')
plt.plot(X1007, V1007, 'r--', label=r'$T^* = 7000$ K, $p^* = 100$ atm')
plt.plot(X15, V15, 'b', label=r'$T^* = 5000$ K, $p^* = 1$ atm$')
plt.plot(X1005, V1005, 'b--', label=r'$T^* = 5000$ K, $p^* = 100$ atm')

plt.xlabel(r'$x/r^*$', fontsize=16)
plt.ylabel(r'$v$, ' + u'м/с', rotation=0, fontsize=16)
plt.gca().yaxis.set_label_coords(0.05, 1)
plt.gca().xaxis.set_label_coords(1.05, 0.05)
plt.legend(loc=0, fontsize=18)
plt.tick_params(labelsize=14)
plt.gca().set_xlim([min(X17), max(X17)])
plt.gca().set_ylim([1500, 4500])


save(name='v_init', fmt='pdf')
plt.tight_layout()
plt.show()

plt.close()
