import numpy as np
import matplotlib.pyplot as plt
f1='./compare_amp/obs_amp_1675_1725.lst'
f2='./compare_amp/pre_amp_1675_1725.lst'

inArr1=np.loadtxt(f1)
inArr2=np.loadtxt(f2)
amp1=inArr1[:,3]
az1=inArr1[:,0]
amp2=inArr2[:,3]
az2=inArr2[:,0]
amp1=amp1/amp1.max()
amp2=amp2/amp2.max()
# ax=plt.subplot(111)
line1, =plt.plot(az1, amp1,'ob', markersize=15)
line2, =plt.plot(az2, amp2,'or', markersize=15)
plt.legend([line1, line2], ['3d', 'ray'])
plt.show()