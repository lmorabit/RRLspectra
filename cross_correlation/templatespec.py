import matplotlib.pyplot as plt
import numpy as np


def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return idx

def find_array(myArray, refArray):
    out1 = np.empty(myArray.size, dtype=int)
    for i, value in enumerate(myArray):
        index = find_nearest(refArray, value)
        out1[i] = index
    return out1

def templatespec(fmin, fmax, c_freq, vwidth):
    speedoflight = 2.99792458

    cfreqs = [i for i in c_freq if i >= fmin and i <=fmax]
    cfreqs = np.asarray(cfreqs[::-1])
    gwidth = cfreqs * vwidth / speedoflight
    gsig = gwidth / 2.0 / np.sqrt(2.0* np.log(2.0))

    nabs = (fmax -fmin) * 1e5
    cont_abscissa = np.arange(nabs)/1e5 + fmin
    cont_values = np.zeros(cont_abscissa.size)

    print('cont_abscissa',cont_abscissa.size)

    cfreqs_index = find_array(cfreqs,cont_abscissa)

    print('cfreqs_index',cfreqs_index)

    nvals = 10000
    xvals = np.arange(nvals) / 1e5 - nvals/ 2 / 1e5
    maxval = 1e0

    for ii in range(0,len(cfreqs_index)-1):

        gs = gsig[ii]
        tmpg = maxval*np.exp(-(xvals-0.0)**2/gs**2)
        cont_values[cfreqs_index[ii] - int((nvals / 2)):cfreqs_index[ii] + int((nvals / 2))] = cont_values[cfreqs_index[ii] - int((nvals / 2)):cfreqs_index[ii] + int((nvals / 2))] - tmpg

    tempspec = np.zeros((2, cont_abscissa.size))
    tempspec[0, :] = cont_abscissa
    tempspec[1, :] = cont_values

    return tempspec

# def main():
#     fl = open('LLINE/RRL_CI_alpha.txt', 'r')
#     freq = []
#     for line in fl:
#         if line.startswith('#'):
#             continue
#         data = line.split(' ')
#         freq.append(float(data[-1].strip()))
#     res =templatespec(200,300,freq,20)
#
#     plt.plot(res[0],res[1])
#     plt.show()
# if __name__ == "__main__":
#     main()