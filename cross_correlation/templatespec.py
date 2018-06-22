import numpy as np

def find_nearest(array, value):
    """
    Finds the index of an array which corresponds to the closest requested value 
    :param array: array of values
    :param value: desired value
    :return: index of desired value within the array
    """
    idx = (np.abs(array - value)).argmin()
    return idx

def find_array(myArray, refArray):
    """
    Finds nearest positions in a reference array, given an input array of values
    :param myArray: input array of values 
    :param refArray: reference array within which the nearest values of myArray are found
    :return: array of index values 
    """
    out1 = np.empty(myArray.size, dtype=int)
    for i, value in enumerate(myArray):
        index = find_nearest(refArray, value)
        out1[i] = index
    return out1

def templatespec(fmin, fmax, c_freq, vwidth):
    """
    function for creating template spectra
    :param fmin: minimum frequency (MHz)
    :param fmax: maximum frequency (MHz)
    :param c_freq: array of frequency values (MHz)
    :param vwidth: frequency width
    :return: template array [0]:frequency [1]: spectrum
    """
    speedoflight = 2.99792458

    cfreqs = [nu for nu in c_freq if nu >= fmin and nu <=fmax]
    cfreqs = np.asarray(cfreqs[::-1])
    gwidth = cfreqs * vwidth / speedoflight
    gsig = gwidth / 2.0 / np.sqrt(2.0 * np.log(2.0))

    nabs = (fmax - fmin) * 1e5
    cont_abscissa = np.arange(nabs)/1e5 + fmin
    cont_values = np.zeros(cont_abscissa.size)
    cfreqs_index = find_array(cfreqs, cont_abscissa)

    nvals = 10000
    xvals = np.arange(nvals) / 1e5 - nvals / 2 / 1e5
    maxval = 1e0

    for ii in range(0,len(cfreqs_index)-1):
        gs = gsig[ii]
        tmpg = maxval*np.exp(-(xvals-0.0)**2/gs**2)
        cont_values[cfreqs_index[ii] - int((nvals / 2)):cfreqs_index[ii] + int((nvals / 2))] = \
            cont_values[cfreqs_index[ii] - int((nvals / 2)):cfreqs_index[ii] + int((nvals / 2))] - tmpg

    tempspec = np.zeros((2, cont_abscissa.size))
    tempspec[0, :] = cont_abscissa
    tempspec[1, :] = cont_values

    return tempspec
