



def process_koala_rss(filename=None, path=None,
                      rss = None,  rss_object_name = None,
                      save_rss_to_fits_file = None,
                      rss_clean=False,
                      # Calibration of the night
                      calibration_night = None,
                      # MASK
                      apply_mask = False,
                      mask = None,   # This can be from a file or a mask
                      make_zeros_in_mask = False, #plot_mask=False,  # Mask if given
                      #valid_wave_min=0, valid_wave_max=0,  # These two are not needed if Mask is given
                      apply_throughput=False,    # ----------------- THROUGHPUT (T)
                      throughput = None,
                      #throughput_2D=[], throughput_2D_file="", throughput_2D_wavecor=False,
                      correct_ccd_defects=False,  # ----------------- CORRECT NANs in CCD (C)
                      #remove_5577=False, kernel_correct_ccd_defects=51, fibre_p=-1, plot_suspicious_fibres=False,
                      
                      fix_wavelengths=False,     # ----------------- FIX WAVELENGTH SHIFTS   (W)
                      wavelength_shift_correction = None, 
                      sky_lines_for_wavelength_shifts=None, 
                      sky_lines_file_for_wavelength_shifts = None, 
                      n_sky_lines_for_wavelength_shifts  = 3,
                      maxima_sigma_for_wavelength_shifts = 2.5, 
                      maxima_offset_for_wavelength_shifts = 1.5,
                      median_fibres_for_wavelength_shifts = 7, 
                      index_fit_for_wavelength_shifts = 2, 
                      kernel_fit_for_wavelength_shifts= None, 
                      clip_fit_for_wavelength_shifts =0.4,
                      fibres_to_plot_for_wavelength_shifts=None,
                      show_fibres_for_wavelength_shifts=None,
                      median_offset_per_skyline_weight = 1.,     # 1 is the BLUE line (median offset per skyline), 0 is the GREEN line (median of solutions), anything between [0,1] is a combination.
                      show_skylines_for_wavelength_shifts = None,
                      plot_wavelength_shift_correction_solution = None,
                      
                      correct_for_extinction=False, # ----------------- EXTINCTION  (X)
                      
                      apply_telluric_correction = False,    # ----------------- TELLURIC  (U)
                      telluric_correction=None, 
                      width_for_telluric_correction = 30,
                      clean_5577 = False,
                      
                      sky_method=None,        # ----------------- SKY (S)
                      skycorrection = None,
                      n_sky=None, 
                      sky_wave_min=None, sky_wave_max=None,
                      sky_fibres=None,  # do_sky=True
                      sky_spectrum = None,  #sky_spectrum_file=None      ### These should be together
                      bright_emission_lines_to_substract_in_sky = None,
                      list_of_skylines_to_fit_near_bright_emission_lines = None,
                      list_of_skylines_to_fit = None,
                      fix_edges = None,
                      fix_edges_wavelength_continuum = None,
                      fix_edges_index_fit=None, 
                      fix_edges_kernel_fit=None,
                      
                      scale_sky = None,
                      #sky_rss=[0], scale_sky_rss=0, scale_sky_1D=0.,
                      #maxima_sigma=3.,
                
                      #sky_lines_file=None, exclude_wlm=[[0, 0]], emission_line_file = None,
                      is_sky=False, win_sky=0, #auto_scale_sky=False, ranges_with_emission_lines=[0], cut_red_end=0,
                      
                      correct_negative_sky=False,           # ----------------- NEGATIVE SKY  (N)
                      min_percentile_for_negative_sky = 5,
                      kernel_for_negative_sky=21,
                      order_fit_for_negative_sky=7,  
                      clip_fit_for_negative_sky = 0.8,
                      individual_check_for_negative_sky=True, # NOT IMPLEMENTED YET #TODO IF NEEDED
                      use_fit_for_negative_sky=True,
                      check_only_sky_fibres = False,
                      force_sky_fibres_to_zero=False,
                      show_fibres_for_negative_sky = None,
                      plot_rss_map_for_negative_sky = False, 
                      
                      id_el=False,     # ----------------- ID emission lines   (I)
                      brightest_line=None, # "Ha",
                      brightest_line_wavelength=None,
                      brightest_fibres_to_combine = None,    #high_fibres=20, 
                      #lowest_fibres_to_combine = None, low_fibres=10,  using n_sky if needed
                      
                      #clean_sky_residuals=False, # ----------------- CLEAN SKY RESIDUALS   (R)
                      big_telluric_residua_correction = False ,
                      max_dispersion_for_big_telluric = 1.4,
                      min_value_per_wave_for_fitting_big_telluric = 50, 
                      telluric_residua_at_6860_correction = False,
                      max_dispersion_for_6860 = 1.4,
                      min_value_per_wave_for_for_6860 = None, 
                      continuum_model_after_sky_correction = None,
                      fibres_to_fix=None,
                      #features_to_fix=[], sky_fibres_for_residuals=[],
                      #remove_negative_median_values=False,
                      
                      correct_extreme_negatives=False,    # ----------- EXTREME NEGATIVES   (R)
                      percentile_min_for_extreme_negatives=0.05,
                      clean_cosmics=False,                # ----------- CLEAN COSMICS     (R)
                      width_bl=20., kernel_median_cosmics=5, cosmic_higher_than=100., extra_factor=1.,
                      max_number_of_cosmics_per_fibre=12, 
                      only_plot_cosmics_cleaned = False,
                      
                      print_summary=False, 
                      plot_final_rss=None,   # None: if a correction is done, it will plot it at the end 
                      plot_final_rss_title = None,   # log= True, gamma = 0.,fig_size=12,
                      **kwargs):    # verbose, plot, warnings should be there     
    """
    This is the most important task, as it calls the rests.
    
    It keeps almost all the parameters that are used in the different tasks. Their description are there, but it will be included here ASAP.
    
    """         
    if rss_clean:                        # Just read file if rss_clean = True
        apply_mask = False
        apply_throughput = False                             # T
        correct_ccd_defects = False                          # C
        fix_wavelengths = False                              # W
        correct_for_extinction = False                       # X
        apply_telluric_correction = False                    # T
        clean_5577 = False                                   # R
        sky_method = None                                    # S
        correct_negative_sky = False                         # N
        id_el = False                                        # E
        big_telluric_residua_correction = False              # R01
        telluric_residua_at_6860_correction = False          # R02
        clean_cosmics = False                                # R04
        correct_extreme_negatives = False                    # R08
        # plot_final_rss = plot
        plot = False
        #verbose = False 
    elif calibration_night is not None:
        if throughput is None and calibration_night.throughput is not None: throughput = calibration_night.throughput
        if wavelength_shift_correction is None and calibration_night.wavelength_shift_correction is not None: wavelength_shift_correction = calibration_night.wavelength_shift_correction     
        if telluric_correction is None and calibration_night.telluric_correction is not None: telluric_correction=calibration_night.telluric_correction
        #if flux_calibration is not None
                   
    verbose = kwargs.get('verbose', False)
    #warnings = kwargs.get('warnings', verbose)
    plot =  kwargs.get('plot', False)
    #plot_all = kwargs.get('plot_all', False)
    
    if plot is False:
        only_plot_cosmics_cleaned = False
        plot_rss_map_for_negative_sky = False

    # Reading the file or the object
    if filename is not None:
        rss = koalaRSS(filename, path = path,
                       rss_object_name = rss_object_name, **kwargs)
    elif rss is None:
        raise RuntimeError("  No rss provided !!!!")
    else:  # rss is an object
        rss=copy.deepcopy(rss)
        if rss_object_name is not None:
            rss.koala.info["rss_object_name"]=rss_object_name
    
    # Get name of original file in case we need it for saving rss at the end
    if filename is None: filename = rss.koala.info["path_to_file"]   

    # Check the number of corrections to be applied
    if (apply_throughput == False and correct_ccd_defects == False and fix_wavelengths == False
        and correct_for_extinction == False and apply_telluric_correction == False and clean_5577 == False
        and sky_method == None and correct_negative_sky == False and id_el == False
        and big_telluric_residua_correction == False and telluric_residua_at_6860_correction == False #and clean_sky_residuals == False    
        and clean_cosmics == False and correct_extreme_negatives == False #and fix_edges == False # and remove_negative_median_values == False
        and is_sky == False and apply_mask == False):
        # If nothing is selected to do, we assume that the RSS file is CLEAN
        rss_clean = True
        #plot_final_rss = plot   
        #plot = False
        #verbose = False
    elif verbose:
        if filename is not None: 
            print("\n> Processing file {} as requested... ".format(filename))
        else:
            filename = rss.koala.info["path_to_file"]
            
            if rss_clean is False:
                if rss_object_name is None:
                    print("\n> Processing rss object as requested... ")
                else:
                    print("\n> Processing rss object {} as requested... ".format(rss_object_name))
                if calibration_night is not None: print("  Calibration of the night provided!")
 
    # Check wavelength range to guess brightest emission line 
    if brightest_line is None:
        if rss.wavelength[0] < 6562.82 and 6562.82 < rss.wavelength[-1]: brightest_line="Ha"
        if rss.wavelength[0] < 5006.84 and 5006.84 < rss.wavelength[-1]: brightest_line="[OIII]"
    if brightest_line_wavelength is None:
        brightest_line_wavelength=quick_find_brightest_line(rss, brightest_fibres_to_combine=brightest_fibres_to_combine, lowest_fibres_to_combine=n_sky)
              
    # Corrections:
    corrections_done = []
    
    if apply_throughput:   # -------------------------------------------------------------------  (T)
        throughput_corr = ThroughputCorrection(throughput=throughput, **kwargs)
        rss = throughput_corr.apply(rss) 
        corrections_done.append("apply_throughput") 
    
    if correct_ccd_defects:  # -----------------------------------------------------------------  (C)
        rss = clean_nan(rss, **kwargs)
        corrections_done.append("correct_ccd_defects")

    if fix_wavelengths:    # -------------------------------------------------------------------  (W)
        # Find correction if not provided
        if wavelength_shift_correction is None:
            wavelength_shift_correction = WavelengthShiftCorrection.wavelength_shift_using_skylines(rss, 
                                                                                                    sky_lines =sky_lines_for_wavelength_shifts,
                                                                                                    sky_lines_file = sky_lines_file_for_wavelength_shifts,
                                                                                                    n_sky_lines = n_sky_lines_for_wavelength_shifts,
                                                                                                    valid_wave_min = rss.koala.info["valid_wave_min"],
                                                                                                    valid_wave_max = rss.koala.info["valid_wave_max"],
                                                                                                    maxima_sigma=maxima_sigma_for_wavelength_shifts, 
                                                                                                    maxima_offset=maxima_offset_for_wavelength_shifts,
                                                                                                    median_fibres = median_fibres_for_wavelength_shifts,
                                                                                                    index_fit = index_fit_for_wavelength_shifts, 
                                                                                                    kernel_fit = kernel_fit_for_wavelength_shifts, 
                                                                                                    clip_fit =clip_fit_for_wavelength_shifts,
                                                                                                    fibres_to_plot=fibres_to_plot_for_wavelength_shifts,
                                                                                                    show_fibres=show_fibres_for_wavelength_shifts,
                                                                                                    **kwargs) #, plot=True, verbose =True)
        elif plot and plot_wavelength_shift_correction_solution is None: 
            plot_wavelength_shift_correction_solution = True
        else: plot_wavelength_shift_correction_solution = False
     
        # Apply solution
        rss = wavelength_shift_correction.apply(rss, 
                                                wavelength_shift_correction=wavelength_shift_correction,
                                                median_offset_per_skyline_weight = median_offset_per_skyline_weight,
                                                show_fibres_for_wavelength_shifts=show_fibres_for_wavelength_shifts,
                                                show_skylines_for_wavelength_shifts = show_skylines_for_wavelength_shifts,
                                                plot_wavelength_shift_correction_solution = plot_wavelength_shift_correction_solution,
                                                **kwargs)   # verbose = True, plot=False)

        corrections_done.append("wavelength_shift_correction")
    
    if correct_for_extinction: # ---------------------------------------------------------------  (X)
        atm_ext_corr = AtmosphericExtCorrection(verbose=verbose)
        rss = atm_ext_corr.apply(rss)
        corrections_done.append("extinction_correction")
        
    if apply_telluric_correction: # ------------------------------------------------------------  (U)
        if telluric_correction is None:
            telluric_correction = TelluricCorrection(rss, verbose=verbose)
            _, fig = telluric_correction.telluric_from_model(plot=plot, width=width_for_telluric_correction)
        elif str(type(telluric_correction)) == "<class 'str'>":    # It is a file
            if path is not None: telluric_correction = os.path.join(path,telluric_correction)
            if verbose: print(" - Reading telluric correction from file",telluric_correction)
            telluric_correction = TelluricCorrection(telluric_correction_file = telluric_correction)
        
            
        rss = telluric_correction.apply(rss, verbose=verbose)
        corrections_done.append("telluric_correction")
        
    if rss.wavelength[0] < 5577 and rss.wavelength[-1] > 5577 and clean_5577: # ----------------  (S)
        rss = clean_skyline(rss, skyline = 5577, **kwargs)
        corrections_done.append("clean_5577")
    
    if sky_method is not None:  # --------------------------------------------------------------  (S)
    
        #TODO: Perhaps this section should be a separated task...
    
        if verbose: print("> Correcting sky using the",sky_method,"method.")

        if skycorrection is None:  
            if sky_spectrum is not None:   # Valid for 1D, 2D
                if verbose: print("  Sky spectrum provided ...")
                skymodel = SkyModel(wavelength=rss.wavelength, intensity = sky_spectrum, variance=np.sqrt(sky_spectrum))  
            else:  # sky_method in ["self", "selffit"]: DEFAULT if sky_method is not ["2D", "1D", "1Dfit"]
                if n_sky is None:
                    if verbose: print("  Using Pablo's method for obtaining self sky spectrum ...")
                    skymodel = SkyFromObject(rss, bckgr_estimator='mad', source_mask_nsigma=5, remove_cont=False)   #TODO
                else:
                    skymodel = SkyFrom_n_sky(rss, n_sky, sky_fibres=sky_fibres, 
                                              sky_wave_min=sky_wave_min, sky_wave_max=sky_wave_max, 
                                              bright_emission_lines_to_substract_in_sky = bright_emission_lines_to_substract_in_sky,
                                              list_of_skylines_to_fit_near_bright_emission_lines = list_of_skylines_to_fit_near_bright_emission_lines,
                                              fix_edges = fix_edges,
                                              fix_edges_wavelength_continuum = fix_edges_wavelength_continuum,
                                              fix_edges_index_fit=fix_edges_index_fit, 
                                              fix_edges_kernel_fit=fix_edges_kernel_fit,
                                              **kwargs)
                    sky_fibres = skymodel.sky_fibres
                    
            #TODO Perhaps it is best to do here the identification of emission lines to know where they are....
            
            if sky_method in ["1Dfit", "selffit"]:   
                sky_spectrum=skymodel.intensity
                sky2D, gaussian_model = model_sky_fitting_gaussians(rss, sky_spectrum, 
                                                                    list_of_skylines_to_fit = list_of_skylines_to_fit, 
                                                                    fibre_list= sky_fibres,
                                                                    plot_continuum = False, **kwargs)
                skymodel = SkyModel(wavelength=rss.wavelength, intensity = sky2D, variance=np.sqrt(sky2D)) 
                skymodel.gaussian_model = gaussian_model
            skycorrection = SkySubsCorrection(skymodel)
        
        rss, _ = skycorrection.apply(rss, verbose=verbose)
        rss.skymodel = skycorrection.skymodel
        if sky_fibres is not None: rss.skymodel.sky_fibres = sky_fibres
        corrections_done.append("sky_correction")
        
    # Correct negative sky  # ------------------------------------------------------------------  (N)
    if is_sky is False and correct_negative_sky is True:   #TODO: This has to be a correction applied to data container
    
        rss= correcting_negative_sky(rss, 
                                      min_percentile_for_negative_sky = min_percentile_for_negative_sky,
                                      individual_check_for_negative_sky = individual_check_for_negative_sky,
                                      kernel_for_negative_sky = kernel_for_negative_sky,
                                      order_fit_for_negative_sky = order_fit_for_negative_sky,
                                      clip_fit_for_negative_sky = clip_fit_for_negative_sky,
                                      use_fit_for_negative_sky = use_fit_for_negative_sky,
                                      check_only_sky_fibres = check_only_sky_fibres,
                                      force_sky_fibres_to_zero=force_sky_fibres_to_zero,
                                      sky_fibres= sky_fibres,
                                      show_fibres = show_fibres_for_negative_sky,
                                      plot_rss_map_for_negative_sky = plot_rss_map_for_negative_sky,
                                      **kwargs)  
        corrections_done.append("negative_sky_correction")
    
    
    # Identify emission lines     # ------------------------------------------------------------  (E)  
    if id_el:
        find_emission_lines_in_koala(rss, brightest_fibres_to_combine = brightest_fibres_to_combine, 
                                      brightest_line = brightest_line,  **kwargs)
        # Update brightest_line_wavelength
        brightest_line_wavelength = rss.koala.info["brightest_line_wavelength"]          
        corrections_done.append("emission_line_identification")
        
    # Clean telluric residua, extreme negatives & cosmics    # ---------------------------------  (R)  
    if big_telluric_residua_correction or telluric_residua_at_6860_correction or correct_extreme_negatives or clean_cosmics:
        # Get continuum image if any of these have been requested
        if continuum_model_after_sky_correction is None:
            if rss.koala.continuum_model_after_sky_correction is not None:
                continuum_model_after_sky_correction = rss.koala.continuum_model_after_sky_correction
            else:
                continuum_model_after_sky_correction = rss_continuum_image(rss, **kwargs)
                rss.koala.continuum_model_after_sky_correction = continuum_model_after_sky_correction
    
    if big_telluric_residua_correction:               #TODO: This has to be a correction applied to data container
        rss = clean_telluric_residuals(rss,                                    
                                        fibre_list = fibres_to_fix,       
                                        continuum = continuum_model_after_sky_correction,
                                        max_dispersion = max_dispersion_for_big_telluric,
                                        min_value_per_wave_for_fitting = min_value_per_wave_for_fitting_big_telluric, 
                                        **kwargs)   
        corrections_done.append("big_telluric_residua_correction")

    if telluric_residua_at_6860_correction:          #TODO: This has to be a correction applied to data container
        rss =     clean_telluric_residuals (rss, continuum=continuum_model_after_sky_correction,    
                                            max_dispersion = max_dispersion_for_6860,
                                            min_value_per_wave_for_fitting = min_value_per_wave_for_for_6860, #10, #50,
                                            interval_to_clean = [6850,6876],
                                            lines_to_fit=[6857,6867],
                                            sigman = [[3.5,3.5]],
                                            max_sigma = [5,5],
                                            max_wave_disp = [15,15],
                                            #fibre_list=use_list,
                                            #plot_fibre_list=use_list,
                                            #plot_individual_comparison 
                                            **kwargs)        
        corrections_done.append("telluric_residua_at_6860_correction")
    
    if correct_extreme_negatives:        #TODO: This has to be a correction applied to data container
        rss = clean_extreme_negatives(rss,                                       
                                      fibre_list=fibres_to_fix, 
                                      percentile_min=percentile_min_for_extreme_negatives,  
                                      continuum = continuum_model_after_sky_correction,
                                      **kwargs)
        corrections_done.append("correct_extreme_negatives")

    # Clean cosmics    (R)
    if clean_cosmics:                    #TODO: This has to be a correction applied to data container
        rss=kill_cosmics(rss,                                                  
                          brightest_line_wavelength, 
                          width_bl=width_bl, 
                          kernel_median_cosmics=kernel_median_cosmics,
                          cosmic_higher_than=cosmic_higher_than, extra_factor=extra_factor,
                          max_number_of_cosmics_per_fibre=max_number_of_cosmics_per_fibre,
                          continuum = continuum_model_after_sky_correction,
                          fibre_list=fibres_to_fix, 
                          only_plot_cosmics_cleaned = only_plot_cosmics_cleaned, **kwargs)  #plot_cosmic_image=plot, plot_RSS_images=plot, verbose=verbose)   
        corrections_done.append("clean_cosmics")
    
    # Apply mask for edges if corrections applied or requested
    if apply_telluric_correction or sky_method is not None or correct_negative_sky or clean_cosmics or correct_extreme_negatives or apply_mask:
        rss = apply_mask_to_rss(rss, mask=mask, make_zeros=make_zeros_in_mask, verbose=verbose) 
        
        
    # Add corrections_done and history:    
    if rss.koala.corrections_done is None:
        rss.koala.corrections_done = corrections_done
    else:
        rss.koala.corrections_done.append(corrections_done)
    
    if rss.koala.info["history"] is None: rss.koala.info["history"] = []  # TODO: needs to do proper history 
    for item in corrections_done: 
        #print(item)
        rss.koala.info["history"].append(item)

    # Summary:
    
    if plot_final_rss is None and len(corrections_done) > 0 and plot is not False : plot_final_rss= True
    
    if len(corrections_done) > 0: 
        if plot_final_rss:
            if plot_final_rss_title is None:
                plot_final_rss_title = rss.info['name'] + " - RSS image - "+str(len(corrections_done))+" corrections applied"
            rss_image(rss, title=plot_final_rss_title, **kwargs)
    
        if verbose or print_summary:    
            if filename is not None:
                print("\n> Summary of processing rss file", '"' + filename + '"', ":")
            elif rss_object_name is not None:
                print("\n> Summary of processing rss object", '"' + rss_object_name + '"', ":")
            else: 
                print("\n> Summary of processing this rss:")
            print('  Name of the observation = "{}",   Name of this Python RSS object = "{}".'.format(rss.info['name'],rss.koala.info['rss_object_name']))
            if len(corrections_done) > 0:
                print(f"  Corrections applied: {len(corrections_done)} in total:")
                for correction in corrections_done: print(f"  - {correction}")
            if rss.koala.info['rss_object_name'] is not None:
                print(f"\n  All applied corrections are stored in {rss.koala.info['rss_object_name']}.intensity !")
    
    # if save_rss_to_fits_file is not None:
    #     if save_rss_to_fits_file == "auto": # These two options, "auto" and "clean", should go in task save_rss_to_fits_file
    #         save_rss_to_fits_file = name_keys(filename, path= path, 
    #                                           apply_throughput = apply_throughput,                                       # T
    #                                           correct_ccd_defects = correct_ccd_defects,                                 # C
    #                                           fix_wavelengths = fix_wavelengths,                                         # W        
    #                                           correct_for_extinction = correct_for_extinction,                           # X
    #                                           apply_telluric_correction = apply_telluric_correction,                     # T
    #                                           sky_method = sky_method,                                                   # S
    #                                           correct_negative_sky = correct_negative_sky,                               # N
    #                                           id_el = id_el,                                                             # E
    #                                           big_telluric_residua_correction = big_telluric_residua_correction,         # R01
    #                                           telluric_residua_at_6860_correction = telluric_residua_at_6860_correction, # R02
    #                                           clean_cosmics = clean_cosmics,                                             # R04
    #                                           correct_extreme_negatives = correct_extreme_negatives)                     # R08 
        
    #     if save_rss_to_fits_file == "clean": save_rss_to_fits_file = filename[:-5]+"_clean.fits"
            
    #     koala_rss_to_fits(rss, fits_file=save_rss_to_fits_file, path = path, verbose=verbose)
    
    return rss



def process_n_koala_rss_files(filename_list = None,
                              path = None,
                              rss_list = None,
                              rss_object_name_list = None,
                              save_rss_to_fits_file_list = None,
                              # more things to add, e.g., sky.. #TODO
                              **kwargs):
    """
    This task process several koala rss files.
    
    As **kwargs, any parameter in process_koala_rss().
    
    Note that this task will use the same parameters for ALL rss files (we will add more options, e.g., adding different sky spectra to used in each file, soon)

    Parameters
    ----------
    filename_list : list of strings, optional
        List with the names of the fits files to process
    path : string or list of strings, optional
        Path to the fits files. 
        If only 1 string is given, it assumes all files listed in filename_list are in the same folder
    rss_list : list of objects, optional
        lift of rss objects
    rss_object_name_list : list of strings, optional
        list with the names of the rss objects
    save_rss_to_fits_file_list : list of strings, optional
        if provided, the processed rss files will be saved in these files
        it will use path if provided

    Raises
    ------
    RuntimeError
        If not filename_list or rss_list is provided.

    Returns
    -------
    processed_rss_files : list of rss objects
        list with the processed rss objects.

    """
    verbose = kwargs.get('verbose', False)
    
    number_rss_files = None
    if filename_list is not None:  number_rss_files = len(filename_list)
    if rss_list is not None:  number_rss_files = len(rss_list)
    if number_rss_files is None:
        raise RuntimeError("NO filename_list or rss_list provided!!!") 
            
    if rss_list is None: rss_list = [None] * number_rss_files
    if filename_list is None: filename_list = [None] * number_rss_files
    if rss_object_name_list is None: rss_object_name_list = [None] * number_rss_files
    if path is None: 
        path = [None] * number_rss_files
    elif np.isscalar(path):
        path = [path] * number_rss_files
    if save_rss_to_fits_file_list is None: 
        save_rss_to_fits_file_list = [None] * number_rss_files
    elif np.isscalar(save_rss_to_fits_file_list):
        save_rss_to_fits_file_list = [save_rss_to_fits_file_list] * number_rss_files
        
    processed_rss_files = []
    
    for i in range(number_rss_files):
        if verbose: print("\n> Processing rss {} of {} :    ---------------------------------------------".format(i+1,number_rss_files))

        _rss_ = process_koala_rss(filename=filename_list[i], 
                                  path = path[i], 
                                  rss = rss_list[i],
                                  rss_object_name=rss_object_name_list[i], 
                                  save_rss_to_fits_file = save_rss_to_fits_file_list[i],
                                  # more things to add #TODO
                                  **kwargs)
        processed_rss_files.append(_rss_)   

    return processed_rss_files
