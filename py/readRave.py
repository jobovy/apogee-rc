#Module to read the RAVE DR4 data
import os, os.path
import asciitable
_DATAFILE= os.path.join(os.getenv('DATADIR'),'rave','ravedr4.dat')
_DATAREADME= os.path.join(os.getenv('DATADIR'),'rave','ReadMe')
def readRave():
    data= asciitable.read(_DATAFILE,
                          readme=_DATAREADME,
                          Reader=asciitable.cds.Cds,
                          guess=False,
                          fill_values=[('', '-999')])
    return data
