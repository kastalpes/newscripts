from GL import *
import numpy.fft as fpack
class thing():
    def __init__(self,data):
        self.data=data
        self._kk = []
        l = len(self.data.shape)
        for i,dim in enumerate(self.data.shape):
            sl = i*(1,)+(dim,)+(l-i-1)*(1,)
            k = fpack.fftfreq(dim)
            k.resize(sl)
            self._kk.append(k)
        self.kradius = np.sqrt(np.sum((k**2 for k in self._kk)))
        self.shortind = (np.array([a.size for a in self._kk])).argmin()
        self.bins = self._kk[self.shortind].flatten()
        self.nx = int(self.bins.size/2)
        self.dk = self.bins[1] - self.bins[0]
box = np.zeros([16]*3)
t1 = thing(box)
