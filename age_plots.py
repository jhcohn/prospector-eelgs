import matplotlib.pyplot as plt

mass_frac = [[], []]
t_h = [[], []]
sigup = [[], []]
sigdown = [[], []]

name = 'out.txt'
name_l = 'out/mass_noelg/7065_uds_noelg_mass.txt'

names = [name, name_l]
for i in range(len(names)):
    idx = 0
    if i % 2 == 1:
        idx = 1
    with open(names[i], 'r') as f:
        for line in f:
            if line[0] == '#':
                pass
            else:
                cols = line.split()
                print(cols)
                mass_frac[idx].append(float(cols[0]))
                t_h[idx].append(float(cols[1]))
                sigup[idx].append(float(cols[2]))
                sigdown[idx].append(float(cols[3]))

err = [[sigdown[0], sigup[0]], [sigdown[1], sigup[1]]]
print(t_h)
print(sigup)
print(sigdown)
print(mass_frac)

# Figure stuff
fig = plt.figure()
ax1 = plt.subplot(1, 2, 1)
ax2 = plt.subplot(1, 2, 2)

ax1.set_yscale('log')  # y is currently time
ax1.set_xlim(0.0, 1.0)
ax1.set_ylim(1e-4, 10)
ax2.set_yscale('log')  # y is currently time
ax2.set_xlim(0.0, 1.0)
ax2.set_ylim(1e-4, 10)

# Labels etc.
fsz = 20  # fontsize
fig.text(0.5, 0.066, r'Fraction of Stellar Mass Formed', ha='center', va='center', fontsize=fsz)
# ax1.set_xlabel(r'Fraction of Stellar Mass Formed')
ax1.set_ylabel(r'Lookback time [Gyr]', fontsize=fsz)
ax1.text(0.8, 3, 'EELGs', fontsize=fsz)  # use if uvj_in
ax2.text(0.8, 3, 'LBGs', fontsize=fsz)  # use if uvj_in

ax1.axvline(x=0.5, linestyle='--', color='k')
ax1.errorbar(mass_frac[0], t_h[0], yerr=err[0])

ax2.axvline(x=0.5, linestyle='--', color='k')
ax2.errorbar(mass_frac[1], t_h[1], yerr=err[1])
plt.show()
