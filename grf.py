import numpy as np
import argparse
import csv

p = argparse.ArgumentParser()
p.add_argument('-a',"--alpha",dest="alpha",type=float,default=-0.5, help="Exponent of power spectrum. Default: -0.5")
p.add_argument('-n',"--norm",dest="amp",type=float,default=0.3, help="Amplitude of fluctuations. Default: 0.3")
p.add_argument('-r',"--resolution",dest="res",type=int,default=100, help="Number of grid points; total is res**(ndim). Default: 100")
p.add_argument('-d',"--ndim",dest="dim",type=int,default=3, help="Number of Cartesian dimensions. Default: 3")
p.add_argument('-b',"--nbase",dest="base",type=float,default=1.0, help="Zero fluctuation value for GRF. Default: 1.0")
p.add_argument('-f',"--filename",dest="fn",type=str,default="grf.txt", help="Name of file to save to. Default: \"grf.txt\"")
args = p.parse_args()

alpha = args.alpha
dim = args.dim
res = args.res
amp = args.amp
base = args.base
fn = args.fn

def fftIndgen(n):
    a = range(0, n/2+1)
    b = range(1, n/2)
    b.reverse()
    b = [-i for i in b]
    return a + b

def gaussian_random_field(Pk = lambda k : k**-3.0, size = 100):
    def Pk3(kx, ky, kz):
        if kx == 0 and ky == 0 and kz == 0:
            return 0.0
        return np.sqrt(Pk(np.sqrt(kx**2 + ky**2 + kz**2)))
    noise = np.fft.fftn(np.random.normal(size = (size, size, size)))
    amplitude = np.zeros((size,size, size))
    for i, kx in enumerate(fftIndgen(size)):
        for j, ky in enumerate(fftIndgen(size)):
            for k, kz in enumerate(fftIndgen(size)):
                amplitude[i, j, k] = Pk3(kx, ky, kz)
    return np.fft.ifftn(noise * amplitude)

out = gaussian_random_field(Pk = lambda k: (k)**alpha, size=res)
phi = amp*out.real/((np.max(out.real) - np.min(out.real))/2.0) + base
phi_output = phi.flatten().astype(str)

fi = open(fn,"w")
for i in range(len(phi_output)):
    fi.write(phi_output[i] + "\n")
fi.close()
