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

plt.rc('font', **{'family': 'serif'})
plt.rc('text', usetex=True)
plt.rc('text.latex', unicode=True)
plt.rc('text.latex', preamble=r"\usepackage[utf8]{inputenc}")
plt.rc('text.latex', preamble=r"\usepackage[russian]{babel}")

x = np.arange(0, 50, 0.01)

# conical

r_cr1 = 1e-3
alpha1 = 0.117 * np.pi
r1 = r_cr1 + x * r_cr1 * np.tan(alpha1)
S1 = np.pi * r1 ** 2 / 2 / r_cr1 ** 2

# hyperbolic

r_cr2 = 3e-3
alpha2 = 1 / 18 * np.pi
r2 = r_cr2 * (1 + (x * r_cr2) ** 2 * (np.tan(alpha2) / r_cr2) ** 2) ** 0.5
S2 = np.pi * r2 ** 2 / 2 / r_cr2 ** 2

# F4

a = 0.3599
bb = 0.2277
cc = 0.1884
d = 0.0184
e = 0.1447
r = a - bb - d / e
r_cr3 = a - bb - d / e
xx = x * r_cr3
r3 = a - bb * np.exp(-cc * xx ** 2) - d / (xx ** 2 + e)
S3 = np.pi * r3 ** 2 / 2 / r_cr3 ** 2

plt.figure()

plt.semilogy(x, S1, 'k', label=u'коническое')
plt.semilogy(x, S2, 'b', label=u'гиперболическое')
plt.semilogy(x, S3, 'r', label=u'F4')

plt.xlabel(r'$x/r^*$', fontsize=labelsize)
plt.ylabel(r'$S(x)/r^*^2$', rotation=0, fontsize=labelsize)
plt.gca().yaxis.set_label_coords(0.12, 1.01)
plt.gca().xaxis.set_label_coords(1.07, 0.07)
plt.legend(loc=0, fontsize=legendsize)
plt.tick_params(labelsize=labelsize)
plt.gca().set_xlim([min(x), max(x)])
plt.gca().set_ylim([0, S1.max()])

save(name='squares', fmt='pdf')
plt.tight_layout()
plt.show()

plt.close()

plt.figure()

plt.plot(x, r1, 'k', label=u'коническое')
plt.plot(x, r2, 'b', label=u'гиперболическое')
plt.plot(x, r3, 'r', label=u'F4')

plt.xlabel(r'$x/r^*$', fontsize=labelsize)
plt.ylabel(r'$r(x),$' + u' м', rotation=0, fontsize=labelsize)
plt.gca().yaxis.set_label_coords(0.1, 1.01)
plt.gca().xaxis.set_label_coords(1.07, 0.07)
plt.legend(loc=0, fontsize=legendsize )
plt.tick_params(labelsize=labelsize)
plt.gca().set_xlim([min(x), max(x)])
plt.gca().set_ylim([0, r3.max()])

save(name='profiles', fmt='pdf')
plt.tight_layout()
plt.show()

plt.close()
