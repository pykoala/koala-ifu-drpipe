"""
Module for reading KOALA-related files
"""

from glob import glob

from pykoala.instruments.koala_ifu import koala_rss
from pykoala import vprint

def list_koala_fits_files_in_folder(path, verbose = True, use2=True, use3=False, ignore_offsets=True, 
                              skyflat_names=None, ignore_list=None, return_list=False):  
    """
    This task just reads the fits files in a folder and prints the KOALA-related fits files.
    
    The files are shown organised by name (in header). Exposition times are also given.
    
    Option of returning the list if return_list == True.

    Parameters
    ----------
    path : string
        path to data
    verbose : Boolean, optional
        Print the list. The default is True.
    use2 : Boolean, optional
        If True, the SECOND word of the name in fit files will be used
    use3 : Boolean, optional
        If True, the THIRD word of the name in fit files will be used
    ignore_offsets : Boolean, optional
        If True it will show all the fits files with the same name but different offsets together. The default is True.
    skyflat_names : list of strings, optional
        List with the names of the skyflat in fits file
    ignore_list : list of strings, optional
        List of words to ignore. The default is None.
    return_list : Boolean, optional
        Return the list. The default is False.

    Raises
    ------
    NameError
        DESCRIPTION.

    Returns
    -------
    if return_list == True, it returns the list of files as: 
        list_of_objetos, list_of_files, list_of_exptimes, date, grating

    """
    
    #FIXME This task needs to be checked, it was used by Ángel old automatic script
    #      but it is still useful for easy printing what the user has in a particular folder.
    
    list_of_objetos=[]
    list_of_files=[]
    list_of_exptimes=[]
    if skyflat_names is None: skyflat_names = ["skyflat", "SKYFLAT", "SkyFlat"]

    if ignore_list is None:  ignore_list = ["a", "b", "c", "d", "e", "f", "p", "pos", "Pos",
                                             "A", "B", "C", "D", "E", "F", "P", "POS",
                                             "p1", "p2","p3","p4","p5","p6",
                                             "P1", "P2","P3","P4","P5","P6",
                                             "pos1", "pos2","pos3","pos4","pos5","pos6",
                                             "Pos1", "Pos2","Pos3","Pos4","Pos5","Pos6",
                                             "POS1", "POS2","POS3","POS4","POS5","POS6"] 
        
    if verbose: print("\n> Listing 2dFdr fits files in folder",path,":\n")
    
    if path[-1] != "/" : path=path+"/"
    date_ = ''  # TODO: This must be filled somehow. It is just for printing data
    files_ = glob.glob(path + '*.fits')
    if len(files_) == 0: raise NameError('No files found within folder '+path)
 
    # Ignore fits products from 2dFdr, darks, flats, arcs...
    files=[]
    for fitsName in sorted(files_):
        include_this = True
        if fitsName[-8:] == "tlm.fits" : include_this = False
        if fitsName[-7:] == "im.fits" : include_this = False
        if fitsName[-7:] == "ex.fits" : include_this = False  
        
        if include_this: 
            hdulist = fits.open(fitsName)
            try:
                object_class = hdulist[0].header['NDFCLASS']
                if object_class ==  "MFOBJECT": 
                    files.append(fitsName) 
                #else: print(object_class, fitsName)
            except Exception:
                pass

    
    for fitsName in sorted(files):
                
        check_file = True
        if fitsName[-8:] != "red.fits" : 
            check_file = False
        if fitsName[0:8] == "combined" and check_file == False: 
            check_file = True
        for skyflat_name in skyflat_names:
            if skyflat_name in fitsName : check_file = True
        
        hdulist = fits.open(fitsName)   # it was pyfits

        object_fits = hdulist[0].header['OBJECT'].split(" ")
    
        if object_fits[0] in ["HD", "NGC", "IC"] or use2:
            try:
                if not ignore_offsets:
                    object_fits[0]=object_fits[0]+object_fits[1]
                elif object_fits[1] not in ignore_list:
                    object_fits[0]=object_fits[0]+object_fits[1]
            except Exception:
                pass
        if use3:
            try:
                if not ignore_offsets:
                    object_fits[0]=object_fits[0]+object_fits[2]
                elif object_fits[2] not in ignore_list:
                    object_fits[0]=object_fits[0]+object_fits[2]
            except Exception:
                pass
            
        try:
            exptime = hdulist[0].header['EXPOSED']
        except Exception:
            check_file = False 

        grating = hdulist[0].header['GRATID']
        try:
            date_ = hdulist[0].header['UTDATE']
        except Exception:
            check_file = False 
        hdulist.close()
        
        if check_file:
            found=False
            for i in range(len(list_of_objetos)):          
                if list_of_objetos[i] == object_fits[0]:
                    found=True
                    list_of_files[i].append(fitsName)               
                    list_of_exptimes[i].append(exptime)
            if not found:
                list_of_objetos.append(object_fits[0])
                list_of_files.append([fitsName])
                list_of_exptimes.append([exptime])
             
    date =date_[0:4]+date_[5:7]+date_[8:10]   
             
    if verbose:
        for i in range(len(list_of_objetos)):
            for j in range(len(list_of_files[i])):
                if j == 0: 
                    vprint("  {:15s}  {}          {:.1f} s".format(list_of_objetos[i], list_of_files[i][0], list_of_exptimes[i][0]))
                else:
                    vprint("                   {}          {:.1f} s".format(list_of_files[i][j], list_of_exptimes[i][j]))
                        
        vprint("\n  They were obtained on {} using the grating {}".format(date, grating))

    if return_list: return list_of_objetos, list_of_files, list_of_exptimes, date, grating
