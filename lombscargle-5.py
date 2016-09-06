import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')

from pylab import *
import pylab as P

from astropy.stats import LombScargle

from astropy import table

# read table
#data = table.Table.read('/home/ksilva/Desktop/Var-python/rx0154/fotometria/rx0154_vmag.dat',format='ascii.fixed_width_no_header') #This reads a table without header
data = table.Table.read('/home/ksilva/Desktop/Var-python/rx0154_vflux.dat',format='ascii') #This reads a table without header
print data

t   = data['col1'] 
mag =data['col2']
dmag =data['col3']

# compute LS
#freq, PLS = LombScargle(t, mag).autopower(minimum_frequency=1 / 0.2,maximum_frequency=1 / 0.001, method='chi2')
freq, PLS = LombScargle(t, mag).autopower(minimum_frequency=1 / 0.105,maximum_frequency=1 / 0.01, method='chi2')
best_freq = freq[np.argmax(PLS)]
phase = (t * best_freq) % 1

# compute the best-fit model
phase_fit = np.linspace(0, 1)
mag_fit = LombScargle(t, mag, dmag).model(t=phase_fit / best_freq,
                                          frequency=best_freq)

# set up the figure & axes for plotting
fig, ax = plt.subplots(1, 3, figsize=(12, 5))
fig.suptitle('Lomb-Scargle Periodogram (period=' +str(1/best_freq)+ ' days)')
fig.subplots_adjust(bottom=0.12, left=0.07, right=0.95)
	#inset = fig.add_axes([0.78, 0.56, 0.15, 0.3])

	# plot the raw data
	
ax[0].errorbar(t, mag, dmag, fmt='ok', elinewidth=1.5, capsize=0)
#ax[0].invert_yaxis()
ax[0].set(
		#xlim=(0, 50)
xlabel='HJD (days)',ylabel='Observed Flux (Relative)')

	# plot the periodogram
ax[1].plot(1. / freq, PLS)
ax[1].set(xlabel='period (days)',ylabel='Lomb-Scargle Power', ylim=(0, 0.8))

	# plot the phased data & model in the inset
	#inset.errorbar(phase, mag, dmag, fmt='.k', capsize=0)
	#inset.plot(phase_fit, mag_fit)
	#inset.invert_yaxis()
	#inset.set_xlabel('phase')
	#inset.set_ylabel('mag')

	# plot the periodogram
ax[2].errorbar(phase, mag, dmag, fmt='.k', capsize=0)
ax[2].errorbar(1+phase, mag, dmag, fmt='.k', capsize=0)
ax[2].plot(phase_fit, mag_fit, color='red')
ax[2].plot(1+phase_fit, mag_fit, color='red')
#plt.gca().invert_yaxis()
ax[2].set(xlabel='phase')#,
	          #xlim=(0.2, 2),
	          #ylim=(0, 1));
#lg=P.legend(loc=1,prop={'size':7})
#lg.draw_frame(False)
#plt.show()
PLOTFILENAME = '/home/ksilva/Desktop/Var-python/RX0154.jpg' #accepts pdf,png
P.savefig(PLOTFILENAME)
#print
#print 'Plot saved on: ', PLOTFILENAME


