import os, os.path
import numpy
import asciitable
datadir= '../data'
def read_valentini():
    data1= asciitable.read(os.path.join(datadir,'arcs_a.dat'),
                           readme="../data/ReadMe",
                           Reader=asciitable.cds.Cds,
                           guess=False,
                    fill_values=[('', '-999')])
    data2= asciitable.read(os.path.join(datadir,'arcs_b.dat'),
                           readme="../data/ReadMe",
                           Reader=asciitable.cds.Cds,
                           guess=False,
                    fill_values=[('', '-999')])
    #Remove duplicates
    indxarray1= numpy.zeros(len(data1),dtype='bool')+True
    hipnums= []
    for ii in range(len(data1)):
        if data1['HIP'][ii] in hipnums:
            indxarray1[ii]= False
            continue
        hipnums.append(data1['HIP'][ii])
    data1= data1[indxarray1]
    indxarray2= numpy.zeros(len(data2),dtype='bool')+True
    hipnums= []
    for ii in range(len(data2)):
        if data2['HIP'][ii] in hipnums:
            indxarray2[ii]= False
            continue
        hipnums.append(data2['HIP'][ii])
    data2= data2[indxarray2]
    return data2
