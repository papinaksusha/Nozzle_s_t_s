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
matname1 = 'NOZ1_100_7000_OSC2_EX1_REC1_'
matname2 = 'NOZ1_100_7000_OSC2_EX2_REC1_'
matname3 = 'NOZ1_100_7000_OSC2_EX3_REC1_'

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

data1 = scipy.io.loadmat('../Matlab/MAT/' + matname1 + '.mat')
data2 = scipy.io.loadmat('../Matlab/MAT/' + matname2 + '.mat')
data3 = scipy.io.loadmat('../Matlab/MAT/' + matname3 + '.mat')

n_i_N21 = data1['u_N2']
n_i_O21 = data1['u_O2']
nn_i_N21 = data1['n_i_N2']
nn_i_O21 = data1['n_i_O2']
X1 = data1['X']
T1 = data1['T']
n_i_N22 = data2['u_N2']
n_i_O22 = data2['u_O2']
nn_i_N22 = data2['n_i_N2']
nn_i_O22 = data2['n_i_O2']
X2 = data2['X']
T2 = data2['T']
n_i_N23 = data3['u_N2']
n_i_O23 = data3['u_O2']
nn_i_N23 = data3['n_i_N2']
nn_i_O23 = data3['n_i_O2']
X3 = data3['X']
T3 = data3['T']

i_N2 = range(0, 48)
i_O2 = range(0, 37)

#N2_i

plt.figure()

plt.semilogy(X1, nn_i_N21[:, 0], 'k', label=u'(1)')
plt.semilogy(X1, nn_i_N21[:, 20], 'k--')
plt.semilogy(X1, nn_i_N21[:, 40], 'k:')

plt.semilogy(X2, nn_i_N22[:, 0], 'b', label=u'(2)')
plt.semilogy(X2, nn_i_N22[:, 20], 'b--')
plt.semilogy(X2, nn_i_N22[:, 40], 'b:')

plt.semilogy(X3, nn_i_N23[:, 0], 'r', label=u'(3)')
plt.semilogy(X3, nn_i_N23[:, 20], 'r--')
plt.semilogy(X3, nn_i_N23[:, 40], 'r:')

#plt.title(r'$\mathrm{N}_2$', fontsize=24)
plt.xlabel(r'$x/r^*$', fontsize=labelsize)
plt.ylabel(r'$n_{\mathrm{N}_2}/n$', rotation=0, fontsize=labelsize)
plt.gca().yaxis.set_label_coords(0.08, 1)
plt.gca().xaxis.set_label_coords(1.07, 0.07)
plt.legend(loc='lower left' , fontsize=legendsize) # bbox_to_anchor=(1, 0.9)
plt.tick_params(labelsize=labelsize)
plt.gca().set_xlim([0, 50])
plt.gca().set_ylim(1e-25, 1)# ni_N217[:, 0].max()])

save(name='N2_i_ex', fmt='pdf')
plt.tight_layout()
plt.show()
plt.close()

#O2_i

plt.figure()

plt.semilogy(X1, nn_i_O21[:, 0], 'k', label=u'(1)')
plt.semilogy(X1, nn_i_O21[:, 15], 'k--')
plt.semilogy(X1, nn_i_O21[:, 30], 'k:')

plt.semilogy(X2, nn_i_O22[:, 0], 'b', label=u'(2)')
plt.semilogy(X2, nn_i_O22[:, 15], 'b--')
plt.semilogy(X2, nn_i_O22[:, 30], 'b:')

plt.semilogy(X3, nn_i_O23[:, 0], 'r', label=u'(3)')
plt.semilogy(X3, nn_i_O23[:, 15], 'r--')
plt.semilogy(X3, nn_i_O23[:, 30], 'r:')

#plt.title(r'$\mathrm{O}_2$', fontsize=24)
plt.xlabel(r'$x/r^*$', fontsize=labelsize)
plt.ylabel(r'$n_{\mathrm{N}_2}/n$', rotation=0, fontsize=labelsize)
plt.gca().yaxis.set_label_coords(0.08, 1)
plt.gca().xaxis.set_label_coords(1.07, 0.07)
plt.legend(loc='lower left' , fontsize=legendsize) # bbox_to_anchor=(1, 0.9)
plt.tick_params(labelsize=labelsize)
plt.gca().set_xlim([0, 50])
plt.gca().set_ylim(1e-7, 1)# ni_N217[:, 0].max()])

save(name='O2_i_ex', fmt='pdf')
plt.tight_layout()
plt.show()
plt.close()

# N2_i x=50

plt.figure()

plt.semilogy(i_N2, n_i_N21[-1, :], 'k', label=u'(1)')
plt.semilogy(i_N2, n_i_N22[-1, :], 'b', label=u'(2)')
plt.semilogy(i_N2, n_i_N23[-1, :], 'r', label=u'(3)')

#plt.title(r'$\mathrm{N}_2$', fontsize=24)
plt.xlabel(r'$i$', fontsize=labelsize)
plt.ylabel(r'$n_{\mathrm{N}_2,i}/n$', rotation=0, fontsize=labelsize)
plt.gca().yaxis.set_label_coords(0.08, 1)
plt.gca().xaxis.set_label_coords(1.05, 0.07)
plt.legend(loc='lower left' , fontsize=legendsize) # bbox_to_anchor=(1, 0.9)
plt.tick_params(labelsize=labelsize)
plt.gca().set_xlim([0, max(i_N2)])
#plt.gca().set_ylim(1e-14, 1)# ni_N217[:, 0].max()])

save(name='N2_i_ex_50', fmt='pdf')
plt.tight_layout()
plt.show()
plt.close()

# O2_i_ex_50
plt.figure()

plt.semilogy(i_O2, n_i_O21[-1, :], 'k', label=u'(1)')
plt.semilogy(i_O2, n_i_O22[-1, :], 'b', label=u'(2)')
plt.semilogy(i_O2, n_i_O23[-1, :], 'r', label=u'(3)')

#plt.title(r'$\mathrm{O}_2$', fontsize=24)
plt.xlabel(r'$i$', fontsize=labelsize)
plt.ylabel(r'$n_{\mathrm{O}_2,i}/n$', rotation=0, fontsize=labelsize)
plt.gca().yaxis.set_label_coords(0.08, 1)
plt.gca().xaxis.set_label_coords(1.05, 0.07)
plt.legend(loc='lower left' , fontsize=legendsize) # bbox_to_anchor=(1, 0.9)
plt.tick_params(labelsize=labelsize)
plt.gca().set_xlim([0, max(i_O2)])
#plt.gca().set_ylim(1e-14, 1)# ni_N217[:, 0].max()])

save(name='O2_i_ex_50', fmt='pdf')
plt.tight_layout()
plt.show()
plt.close()

# T

plt.figure()

plt.plot(X1, T1, 'k', label=u'(1)')
plt.plot(X2, T2, 'b', label=u'(2)')
plt.plot(X3, T3, 'r', label=u'$(3)')

Teps12 = max(abs(T1 - T2)/T1)*100
Teps13 = max(abs(T1 - T3)/T1)*100
Teps32 = max(abs(T3 - T2)/T3)*100

print('Teps12 = ', Teps12, 'Teps13 = ', Teps13, 'Teps32 = ', Teps32)

plt.xlabel(r'$x/r^*$', fontsize=labelsize)
plt.ylabel(r'$T$, K', rotation=0, fontsize=labelsize)
plt.gca().yaxis.set_label_coords(0.07, 1)
plt.gca().xaxis.set_label_coords(1.07, 0.07)
plt.legend(loc=1, fontsize=legendsize)
plt.tick_params(labelsize=labelsize)
plt.gca().set_xlim([min(X1), 5])
plt.gca().set_ylim([0, 7000])

save(name='T_ex', fmt='pdf')
plt.tight_layout()
plt.show()

plt.close()

