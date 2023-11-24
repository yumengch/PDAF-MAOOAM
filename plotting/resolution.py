import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

#  n_exp           9 nt        2000 n_sample        1000
y = [7.7422326912817974E-002,
     4.7561770811221454E-002,
     9.5612101521581269E-003,
     6.8200650626671071E-004,
     1.0978467893127872E-005,
     4.3202805836482769E-008,
     1.0362144475051586E-009,
     1.0362145993658874E-009,
     1.0362145993862703E-009]

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
fig = plt.figure(1, figsize=(12, 6))
fig.clf()
ax = fig.add_subplot(111)
ax.plot(range(1, 10), np.log(y), 'k', linewidth=2, marker='+', markersize=8)
f = lambda x, pos: r'{} $\times$ {}'.format(int(2**x + 1), int(2**x + 1))
ax.xaxis.set_major_locator(mticker.MultipleLocator(1))
ax.xaxis.set_major_formatter(f)
ax.set_xlabel('No. of grid points')
ax.set_ylabel('Maximum error from transformation')
fig.savefig('trans_error.eps', bbox_inches='tight')