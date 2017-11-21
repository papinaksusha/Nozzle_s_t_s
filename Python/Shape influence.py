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
legendsize = 22
matname1 = 'NOZ1_1_7000_OSC2_EX3_REC1_'
matname2 = 'NOZ2_1_7000_OSC2_EX3_REC1_'
matname3 = 'NOZ3_1_7000_OSC2_EX3_REC1_'

plt.rc('font', **{'family':'serif'})
plt.rc('text', usetex=True)
plt.rc('text.latex', unicode=True)
plt.rc('text.latex', preamble=r"\usepackage[utf8]{inputenc}")
plt.rc('text.latex', preamble=r"\usepackage[russian]{babel}")

data1 = scipy.io.loadmat('../Matlab/MAT/' + matname1 + '.mat')
data2 = scipy.io.loadmat('../Matlab/MAT/' + matname2 + '.mat')
data3 = scipy.io.loadmat('../Matlab/MAT/' + matname3 + '.mat')

n_N21 = data1['n_N2']
n_N22 = data2['n_N2']
n_N23 = data3['n_N2']

n_O21 = data1['n_O2']
n_O22 = data2['n_O2']
n_O23 = data3['n_O2']

T1 = data1['T']
T2 = data2['T']
T3 = data3['T']

X1 = data1['X']
X2 = data2['X']
X3 = data3['X']

# T

plt.figure()

plt.plot(X1, T1, 'k', label=u'коническое')
plt.plot(X2, T2, 'b', label=u'гиперболическое')
plt.plot(X3, T3, 'r', label=u'F4')

plt.xlabel(r'$x/r^*$', fontsize=labelsize)
plt.ylabel(r'$T$, K', rotation=0, fontsize=labelsize)
plt.gca().yaxis.set_label_coords(0.05, 1)
plt.gca().xaxis.set_label_coords(1.07, 0.07)
#plt.legend(loc='upper right', bbox_to_anchor=(1, 0.9), fontsize=18)
plt.legend(loc='upper right', fontsize=legendsize)
plt.tick_params(labelsize=labelsize)
plt.gca().set_xlim([0, 50])
plt.gca().set_ylim([0, 7000])

save(name='T_shapes', fmt='pdf')
plt.tight_layout()
plt.show()
plt.close()

plt.figure()

plt.plot(X1, n_N21, 'k', label=u'коническое')
plt.plot(X2, n_N22, 'b', label=u'гиперболическое')
plt.plot(X3, n_N23, 'r', label=u'F4')

plt.title(r'$\mathrm{N}_2$', fontsize=24)
plt.xlabel(r'$x/r^*$', fontsize=labelsize)
plt.ylabel(r'$n_{\mathrm{N}_2}/n$', rotation=0, fontsize=labelsize)
plt.gca().yaxis.set_label_coords(0.07, 1)
plt.gca().xaxis.set_label_coords(1.07, 0.07)
plt.legend(loc='upper right', bbox_to_anchor=(1, 0.9), fontsize=legendsize)
plt.tick_params(labelsize=labelsize)
plt.gca().set_xlim([0, 5])
plt.locator_params(axis='y', numticks=10)
#plt.gca().set_ylim([n_O21007.min(), n_O21007.max()])

save(name='N2_shapes', fmt='pdf')
plt.tight_layout()
plt.show()

plt.close()

plt.figure()

plt.plot(X1, n_O21, 'k', label=u'коническое')
plt.plot(X2, n_O22, 'b', label=u'гиперболическое')
plt.plot(X3, n_O23, 'r', label=u'F4')

plt.title(r'$\mathrm{O}_2$', fontsize=24)
plt.xlabel(r'$x/r^*$', fontsize=labelsize)
plt.ylabel(r'$n_{\mathrm{O}_2}/n$', rotation=0, fontsize=labelsize)
plt.gca().yaxis.set_label_coords(0.07, 1)
plt.gca().xaxis.set_label_coords(1.07, 0.07)
plt.legend(loc='upper right', bbox_to_anchor=(1, 0.9), fontsize=legendsize)
plt.tick_params(labelsize=labelsize)
plt.locator_params(axis='y', numticks=10)
plt.gca().set_xlim([0, 5])
#plt.gca().set_ylim([n_O21007.min(), n_O21007.max()])

save(name='O2_shapes', fmt='pdf')
plt.tight_layout()
plt.show()

plt.close()
