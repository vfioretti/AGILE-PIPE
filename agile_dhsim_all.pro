; agile_dhsim_all.pro - Description
; ---------------------------------------------------------------------------------
; Compacting the BoGEMMS AGILE processed files to a unique file:
; - Tracker
; - AC
; - Calorimeter
; ---------------------------------------------------------------------------------
; Output:
; - all files are created in the same directory of the input files
; ---------> ASCII files
; - G4.RAW.KALMAN.AGILE<version>.<phys>List.<sim_type>.<strip>.<point>.<n_in>ph.<energy>MeV.<theta>.<phi>.all.dat
; - G4.RAW.GENERAL.AGILE<version>.<phys>List.<sim_type>.<strip>.<point>.<n_in>ph.<energy>MeV.<theta>.<phi>.all.dat
; - G4.DIGI.KALMAN.AGILE<version>.<phys>List.<sim_type>.<strip>.<point>.<n_in>ph.<energy>MeV.<theta>.<phi>.all.dat
; - G4.DIGI.GENERAL.AGILE<version>.<phys>List.<sim_type>.<strip>.<point>.<n_in>ph.<energy>MeV.<theta>.<phi>.all.dat
; - G4_GAMS_XPLANE_AGILE<version>_<phys>List_<sim_type>_<strip>_<point>_<n_in>ph_<energy>MeV_<theta>_<phi>.all.dat
; - G4_GAMS_YPLANE_AGILE<version>_<phys>List_<sim_type>_<strip>_<point>_<n_in>ph_<energy>MeV_<theta>_<phi>.all.dat
; - G4_GAMS_CAL_AGILE<version>_<phys>List_<sim_type>_<strip>_<point>_<n_in>ph_<energy>MeV_<theta>_<phi>.all.dat
; - G4_GAMS_AC_AGILE<version>_<phys>List_<sim_type>_<strip>_<point>_<n_in>ph_<energy>MeV_<theta>_<phi>.all.dat
; - stripx.dat
; - stripy.dat
; - ac.dat
; - calo.dat
; ---------> FITS files
; - G4.RAW.AGILE<version>.<phys>List.<sim_type>.<strip>.<point>.<n_in>ph.<energy>MeV.<theta>.<phi>.all.fits
; - L0.AGILE<version>.<phys>List.<sim_type>.<strip>.<point>.<n_in>ph.<energy>MeV.<theta>.<phi>.all.fits
; - L0.5.DIGI.AGILE<version>.<phys>List.<sim_type>.<strip>.<point>.<n_in>ph.<energy>MeV.<theta>.<phi>.all.fits
; - G4.AGILE<version>.AC.<phys>List.<sim_type>.<strip>.<point>.<n_in>ph.<energy>MeV.<theta>.<phi>.all.fits
; - G4.AGILE<version>.CAL.<phys>List.<sim_type>.<strip>.<point>.<n_in>ph.<energy>MeV.<theta>.<phi>.all.fits 
; ---------------------------------------------------------------------------------
; copyright            : (C) 2013 Valentina Fioretti
; email                : fioretti@iasfbo.inaf.it
; ----------------------------------------------
; Usage:
; agile_dhsim_all
; ---------------------------------------------------------------------------------


pro agile_dhsim_all, $
		agile_version, $
		out_type, $
		sim_type, $
		py_list, $
		N_in, $
		n_files, $
		part_type, $
		ene_range, $
		ene_min, $
		ene_max, $
		theta_type, $
		phi_type, $
		source_g, $
		isStrip, $
		repli, $
		cal_flag, $
		ac_flag


; Variables initialization
;N_in = 0UL            ;--> Number of emitted photons
;n_fits = 0           ;--> Number of FITS files produced by the simulation

;agile_version = ''
;out_type = 0
;sim_type = 0
;py_list = 0
;ene_range = 0
;ene_type = 0
;ene_min = 0
;ene_max = 0
;theta_type = 0
;phi_type = 0
;source_g = 0
;cal_flag = 0
;ac_flag = 0

;read, agile_version, PROMPT='% - Enter AGILE release (e.g. V1.4):'
;read, out_type, PROMPT='% - Enter simulation output type (0: FULL, 1: only DHSim):'
;read, sim_type, PROMPT='% - Enter simulation type [0 = Mono, 1 = Chen, 2: Vela, 3: Crab, 4: G400, 5 = SS]:'
;read, py_list, PROMPT='% - Enter the Physics List [0 = QGSP_BERT_EMV, 100 = ARGO, 300 = FERMI, 400 = ASTROMEV]:'
;read, N_in, PROMPT='% - Enter the number of emitted photons:'
;read, n_fits, PROMPT='% - Enter number of FITS files:'
;read, ene_range, PROMPT='% - Enter energy distribution [0 = mono, 1 = range]:'

if (ene_range EQ 0) then begin
  ene_dis = 'MONO'
  ene_type = ene_min
  if (ene_type GE 1) then ene_type = strtrim(string(long(ene_type)),1)
  if (ene_type LT 1) then ene_type = STRMID(STRTRIM(STRING(ene_type),1),0,5)
  if (size(ene_min, /TYPE) NE 2) then begin
      nstring = strlen(ene_type)
      ene_type_notzero = ene_type
      flag = 1 
      for ichar_reverse=0, nstring-1 do begin
        ichar = (nstring-1) - ichar_reverse
        if ((strmid(ene_type, ichar, 1) EQ '0') or (strmid(ene_type, ichar, 1) EQ '.')) then begin
            if (flag EQ 1) then ene_type_notzero = STRMID(ene_type_notzero, 0, ichar)
        endif else begin
            flag = 0
        endelse
      endfor
      ene_type = ene_type_notzero
  endif
endif
if (ene_range EQ 1) then begin
    ene_dis = 'POW'
    
    ene_min_string = strtrim(string(ene_min),1)
    if (size(ene_min, /TYPE) NE 2) then begin
      nstring = strlen(ene_min_string)
      ene_min_string_notzero = ene_min_string
      flag = 1 
      for ichar_reverse=0, nstring-1 do begin
        ichar = (nstring-1) - ichar_reverse
        if ((strmid(ene_min_string, ichar, 1) EQ '0') or (strmid(ene_min_string, ichar, 1) EQ '.')) then begin
            if (flag EQ 1) then ene_min_string_notzero = STRMID(ene_min_string_notzero, 0, ichar)
        endif else begin
            flag = 0
        endelse
      endfor
      ene_min_string = ene_min_string_notzero
    endif
    
    ene_max_string = strtrim(string(ene_max),1)
    if (size(ene_max, /TYPE) NE 2) then begin
      nstring = strlen(ene_max_string)
      ene_max_string_notzero = ene_max_string
      flag = 1
      for ichar_reverse=0, nstring-1 do begin
        ichar = (nstring-1) - ichar_reverse
        if ((strmid(ene_max_string, ichar, 1) EQ '0') or (strmid(ene_max_string, ichar, 1) EQ '.')) then begin
          if (flag EQ 1) then ene_max_string_notzero = STRMID(ene_max_string_notzero, 0, ichar)
        endif else begin
          flag = 0
        endelse
      endfor
      ene_max_string = ene_max_string_notzero
    endif
    
    ene_type = ene_min_string+'.'+ene_max_string
endif
if (ene_range EQ 2) then begin
  ene_dis = 'EXP'
    ene_min_string = strtrim(string(ene_min),1)
    if (size(ene_min_string, /TYPE) NE 2) then begin
      nstring = strlen(ene_min_string)
      ene_min_string_notzero = ene_min_string
      flag = 1 
      for ichar_reverse=0, nstring-1 do begin
        ichar = (nstring-1) - ichar_reverse
        if ((strmid(ene_min_string, ichar, 1) EQ '0') or (strmid(ene_min_string, ichar, 1) EQ '.')) then begin
            if (flag EQ 1) then ene_min_string_notzero = STRMID(ene_min_string_notzero, 0, ichar)
        endif else begin
            flag = 0
        endelse
      endfor
      ene_min_string = ene_min_string_notzero
    endif
    
    ene_max_string = strtrim(string(ene_max),1)
    if (size(ene_max_string, /TYPE) NE 2) then begin
      nstring = strlen(ene_max_string)
      ene_max_string_notzero = ene_max_string
      flag = 1
      for ichar_reverse=0, nstring-1 do begin
        ichar = (nstring-1) - ichar_reverse
        if ((strmid(ene_max_string, ichar, 1) EQ '0') or (strmid(ene_max_string, ichar, 1) EQ '.')) then begin
          if (flag EQ 1) then ene_max_string_notzero = STRMID(ene_max_string_notzero, 0, ichar)
        endif else begin
          flag = 0
        endelse
      endfor
      ene_max_string = ene_max_string_notzero
    endif
  ene_type = ene_min_string+'.'+ene_max_string
endif
if (ene_range EQ 3) then begin
  ene_dis = 'LIN'
    ene_min_string = strtrim(string(ene_min),1)
    if (size(ene_min_string, /TYPE) NE 2) then begin
      nstring = strlen(ene_min_string)
      ene_min_string_notzero = ene_min_string
      flag = 1 
      for ichar_reverse=0, nstring-1 do begin
        ichar = (nstring-1) - ichar_reverse
        if ((strmid(ene_min_string, ichar, 1) EQ '0') or (strmid(ene_min_string, ichar, 1) EQ '.')) then begin
            if (flag EQ 1) then ene_min_string_notzero = STRMID(ene_min_string_notzero, 0, ichar)
        endif else begin
            flag = 0
        endelse
      endfor
      ene_min_string = ene_min_string_notzero
    endif
    
    ene_max_string = strtrim(string(ene_max),1)
    if (size(ene_max_string, /TYPE) NE 2) then begin
      nstring = strlen(ene_max_string)
      ene_max_string_notzero = ene_max_string
      flag = 1
      for ichar_reverse=0, nstring-1 do begin
        ichar = (nstring-1) - ichar_reverse
        if ((strmid(ene_max_string, ichar, 1) EQ '0') or (strmid(ene_max_string, ichar, 1) EQ '.')) then begin
          if (flag EQ 1) then ene_max_string_notzero = STRMID(ene_max_string_notzero, 0, ichar)
        endif else begin
          flag = 0
        endelse
      endfor
      ene_max_string = ene_max_string_notzero
    endif
  ene_type = ene_min_string+'.'+ene_max_string
endif

;read, theta_type, PROMPT='% - Enter theta:'
;read, phi_type, PROMPT='% - Enter phi:'
;read, source_g, PROMPT='% - Enter source geometry [0 = Point, 1 = Plane]:'

if (py_list EQ 0) then begin
   py_dir = 'QGSP_BERT_EMV'
   py_name = 'QGSP_BERT_EMV'
endif
if (py_list EQ 100) then begin
   py_dir = '100List'
   py_name = '100List'
endif
if (py_list EQ 300) then begin
   py_dir = '300List'
   py_name = '300List'
endif
if (py_list EQ 400) then begin
   py_dir = 'ASTROMEV'
   py_name = 'ASTROMEV'
endif

if (sim_type EQ 0) then begin
   sim_name = 'MONO'
endif
if (sim_type EQ 1) then begin
   sim_name = 'CHEN'
endif
if (sim_type EQ 2) then begin
   sim_name = 'VELA'
endif
if (sim_type EQ 3) then begin
   sim_name = 'CRAB'
endif
if (sim_type EQ 4) then begin
   sim_name = 'G400'
endif
if (sim_type EQ 5) then begin
   sim_name = 'SS'
endif
if (sim_type EQ 6) then begin
   sim_name = 'RANGE'
endif

if (source_g EQ 0) then begin
 sdir = '/Point'
 sname = 'Point'
endif
if (source_g EQ 1) then begin
 sdir = '/Plane'
 sname = 'Plane'
endif

;read, isStrip, PROMPT='% - Strip activated?:'
;read, repli, PROMPT='% - Strips replicated?:'

if (isStrip EQ 0) then stripDir = 'NoStrip/'
if ((isStrip EQ 1) AND (repli EQ 0)) then stripDir = 'StripNoRepli/'
if ((isStrip EQ 1) AND (repli EQ 1)) then stripDir = 'StripRepli/'

if (isStrip EQ 0) then stripname = 'NOSTRIP'
if ((isStrip EQ 1) AND (repli EQ 0)) then stripname = 'STRIP'
if ((isStrip EQ 1) AND (repli EQ 1)) then stripname = 'STRIP.REPLI'

;read, cal_flag, PROMPT='% - Cal present?:'
;read, ac_flag, PROMPT='% - AC present?:'

; Reading the FITS files

; G4.RAW.AGILE<version>.<phys>List.<strip>.<point>.<n_in><part_type>.<energy>MeV.<theta>.<phi>.all.fits
rawData_event_id = -1l
rawData_Si_id = -1l
rawData_tray_id =  -1l
rawData_plane_id =  -1l
rawData_Strip_id =  -1l
rawData_energy_dep =  -1.
rawData_ent_x =  -1.
rawData_ent_y =  -1.
rawData_ent_z =  -1.
rawData_exit_x =  -1.
rawData_exit_y =  -1.
rawData_exit_z =  -1.

; G4.FIRST.AGILE<version>.<phys>List.<strip>.<point>.<n_in><part_type>.<energy>MeV.<theta>.<phi>.all.fits
firstData_event_id = -1l
firstData_Si_id = -1l
firstData_tray_id =  -1l
firstData_plane_id =  -1l
firstData_Strip_id =  -1l
firstData_energy_dep =  -1.
firstData_ent_x =  -1.
firstData_ent_y =  -1.
firstData_ent_z =  -1.
firstData_exit_x =  -1.
firstData_exit_y =  -1.
firstData_exit_z =  -1.

; L0.AGILE<version>.<phys>List.<strip>.<point>.<n_in><part_type>.<energy>MeV.<theta>.<phi>.all.fits
L0TRACKERGLOBAL_Glob_event_id_test = -1l
L0TRACKERGLOBAL_Glob_vol_id_test = -1l
L0TRACKERGLOBAL_Glob_moth_id_test = -1l
L0TRACKERGLOBAL_Glob_tray_id_test = -1l
L0TRACKERGLOBAL_Glob_plane_id_test = -1l
L0TRACKERGLOBAL_Glob_Si_id_test = -1l
L0TRACKERGLOBAL_Glob_Strip_id_test = -1l
L0TRACKERGLOBAL_Glob_Strip_type_test = -1l
L0TRACKERGLOBAL_Glob_pos_test = -1.
L0TRACKERGLOBAL_Glob_zpos_test = -1.
L0TRACKERGLOBAL_Glob_energy_dep_test = -1.

; L0.5.DIGI.AGILE<version>.<phys>List.<strip>.<point>.<n_in><part_type>.<energy>MeV.<theta>.<phi>.all.fits
L05TRACKERGLOBAL_Glob_event_id_acap = -1l
L05TRACKERGLOBAL_Glob_vol_id_acap = -1l
L05TRACKERGLOBAL_Glob_moth_id_acap = -1l
L05TRACKERGLOBAL_Glob_tray_id_acap = -1l
L05TRACKERGLOBAL_Glob_plane_id_acap = -1l
L05TRACKERGLOBAL_Glob_Si_id_acap = -1l
L05TRACKERGLOBAL_Glob_Strip_id_acap = -1l
L05TRACKERGLOBAL_Glob_Strip_type_acap = -1l
L05TRACKERGLOBAL_Glob_pos_acap = -1.
L05TRACKERGLOBAL_Glob_zpos_acap = -1.
L05TRACKERGLOBAL_Glob_energy_dep_acap = -1.

; G4.CAL.AGILE<version>.<phys>List.<strip>.<point>.<n_in><part_type>.<energy>MeV.<theta>.<phi>.all.fits 
calInput_event_id_tot_cal = -1l
calInput_bar_plane_tot = -1l
calInput_bar_id_tot = -1l
calInput_ene_a_tot = -1.
calInput_ene_b_tot = -1.

; G4.AC.AGILE<version>.<phys>List.<strip>.<point>.<n_in><part_type>.<energy>MeV.<theta>.<phi>.all.fits
acInput_event_id_tot_ac = -1l
acInput_AC_panel = ''
acInput_AC_subpanel = -1l
acInput_energy_dep_tot_ac = -1.

; - G4.RAW.KALMAN.AGILE<version>.<phys>List.<strip>.<point>.<n_in><part_type>.<energy>MeV.<theta>.<phi>.all.dat
; ASCII Columns:
; - c1 = event ID
; - c2 = Silicon layer ID
; - c3 = x/y pos [cm]
; - c4 = z pos [cm]
; - c5 = plane ID
; - c6 = strip ID
; - c7 = energy dep [keV]

raw_kalman_event_id = -1l
raw_kalman_Si_id = -1l
raw_kalman_pos = -1.
raw_kalman_zpos = -1.
raw_kalman_plane_id = -1l
raw_kalman_strip_id = -1l
raw_kalman_edep = -1.

; - G4.RAW.GENERAL.AGILE<version>.<phys>List.<strip>.<point>.<n_in><part_type>.<energy>MeV.<theta>.<phi>.all.dat
; ASCII Columns:
; - c1 = event ID
; - c2 = Silicon layer ID
; - c3 = x/y pos [cm]
; - c4 = z pos [cm]
; - c5 = tray ID
; - c6 = plane ID
; - c7 = strip ID 
; - c8 = energy dep [keV]    

raw_general_event_id = -1l
raw_general_Si_id = -1l
raw_general_pos = -1.
raw_general_zpos = -1.
raw_general_tray_id = -1l
raw_general_plane_id = -1l
raw_general_strip_id = -1l
raw_general_edep = -1.    

; - G4.DIGI.KALMAN.AGILE<version>.<phys>List.<strip>.<point>.<n_in><part_type>.<energy>MeV.<theta>.<phi>.all.dat
; ASCII Columns:
; - c1 = event ID
; - c2 = Silicon layer ID
; - c3 = x/y pos [cm]
; - c4 = z pos [cm]
; - c5 = plane ID
; - c6 = strip ID
; - c7 = energy dep [keV]

digi_kalman_event_id = -1l
digi_kalman_Si_id = -1l
digi_kalman_pos = -1.
digi_kalman_zpos = -1.
digi_kalman_plane_id = -1l
digi_kalman_strip_id = -1l
digi_kalman_edep = -1.

; - G4.DIGI.GENERAL.AGILE<version>.<phys>List.<strip>.<point>.<n_in><part_type>.<energy>MeV.<theta>.<phi>.all.dat
; ASCII Columns:
; - c1 = event ID
; - c2 = Silicon layer ID
; - c3 = x/y pos [cm]
; - c4 = z pos [cm]
; - c5 = tray ID
; - c6 = plane ID
; - c7 = strip ID 
; - c8 = energy dep [keV]    

digi_general_event_id = -1l
digi_general_Si_id = -1l
digi_general_pos = -1.
digi_general_zpos = -1.
digi_general_tray_id = -1l
digi_general_plane_id = -1l
digi_general_strip_id = -1l
digi_general_edep = -1.    

; - /G4.AA.STRIP.AGILE'+agile_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+strmid(strtrim(string(N_in),1),0,10)+'<part_type>.'+ene_type+'MeV.'+strmid(strtrim(string(theta_type),1),0,10)+'.'+strmid(strtrim(string(phi_type),1),0,10)+'.'+strtrim(string(ifile),1)+'.dat'  ; ASCII Columns:
; - c1 = event ID
; - c2 = theta input
; - c3 = phi input
; - c4 = energy input
; - c5 = plane ID
; - c6 = Pos Z
; - c7 = X/Y flag (X = 0, Y = 1)
; - c8 = Strip ID
; - c9 = strip position (reference system center at the Silicon layer center)
; - c10 = energy deposition (keV)

aa_strip_event_id = -1l
aa_strip_theta_in = -1l
aa_strip_phi_in = -1l
aa_strip_ene_in = -1l
aa_strip_plane_id = -1l
aa_strip_zpos = -1.
aa_strip_si_id = -1l
aa_strip_strip_id = -1l
aa_strip_pos = -1.
aa_strip_edep = -1.
aa_strip_pair = -1

; - /G4.AA.CLUSTER.AGILE'+agile_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+strmid(strtrim(string(N_in),1),0,10)+'<part_type>.'+ene_type+'MeV.'+strmid(strtrim(string(theta_type),1),0,10)+'.'+strmid(strtrim(string(phi_type),1),0,10)+'.'+strtrim(string(ifile),1)+'.dat'  ; ASCII Columns:
; - c1 = event ID
; - c2 = theta input
; - c3 = phi input
; - c4 = energy input
; - c5 = plane ID
; - c6 = Pos Z
; - c7 = X/Y flag (X = 0, Y = 1)
; - c8 = Cluster position (reference system center at the Silicon layer center)
; - c9 = energy deposition (keV)
; - c10 = number of strips composing the cluster
; - c11 = pair flag

aa_kalman_event_id = -1l
aa_kalman_theta_in = -1l
aa_kalman_phi_in = -1l
aa_kalman_ene_in = -1l
aa_kalman_plane_id = -1l
aa_kalman_zpos = -1.
aa_kalman_si_id = -1l
aa_kalman_pos = -1.
aa_kalman_edep = -1.
aa_kalman_strip_number = -1.
aa_kalman_pair = -1l

if (out_type EQ 0) then begin
	; - G4_AIV_XPLANE_AGILE<version>_<phys>List_<strip>_<point>_<n_in><part_type>_<energy>MeV_<theta>_<phi>.all.dat
	; ASCII Columns:
	; - c1 = event ID
	; - c2 = plane ID
	; - c3 = readout strip ID
	; - c4 = -999
	; - c5 = -999
	; - c6 = energy dep in keV
	; - c7 = -999    event_start = -1
	
	aiv_xplane_event_id = -1l
	aiv_xplane_plane_id = -1l
	aiv_xplane_readstreap_id = -1l
	aiv_xplane_c4 = -1l
	aiv_xplane_c5 = -1l
	aiv_xplane_edep = -1.
	aiv_xplane_c7 = -1l

	; - G4_AIV_YPLANE_AGILE<version>_<phys>List_<strip>_<point>_<n_in><part_type>_<energy>MeV_<theta>_<phi>.all.dat
	; ASCII Columns:
	; - c1 = event ID
	; - c2 = plane ID
	; - c3 = readout strip ID
	; - c4 = -999
	; - c5 = -999
	; - c6 = energy dep in keV
	; - c7 = -999
	
	aiv_yplane_event_id = -1l
	aiv_yplane_plane_id = -1l
	aiv_yplane_readstreap_id = -1l
	aiv_yplane_c4 = -1l
	aiv_yplane_c5 = -1l
	aiv_yplane_edep = -1.
	aiv_yplane_c7 = -1l
endif

; - G4_GAMS_XPLANE_AGILE<version>_<phys>List_<strip>_<point>_<n_in><part_type>_<energy>MeV_<theta>_<phi>.all.dat
; ASCII Columns:
; - c1 = event ID
; - c2 = plane ID
; - c3 = readout strip ID
; - c4 = -999
; - c5 = -999
; - c6 = energy dep in MIP
; - c7 = -999    event_start = -1
    
gams_xplane_event_id = -1l
gams_xplane_plane_id = -1l
gams_xplane_readstreap_id = -1l
gams_xplane_c4 = -1l
gams_xplane_c5 = -1l
gams_xplane_edep = -1.
gams_xplane_c7 = -1l

; - G4_GAMS_YPLANE_AGILE<version>_<phys>List_<strip>_<point>_<n_in><part_type>_<energy>MeV_<theta>_<phi>.all.dat
; ASCII Columns:
; - c1 = event ID
; - c2 = plane ID
; - c3 = readout strip ID
; - c4 = -999
; - c5 = -999
; - c6 = energy dep in MIP
; - c7 = -999
    
gams_yplane_event_id = -1l
gams_yplane_plane_id = -1l
gams_yplane_readstreap_id = -1l
gams_yplane_c4 = -1l
gams_yplane_c5 = -1l
gams_yplane_edep = -1.
gams_yplane_c7 = -1l


; - G4_GAMS_CAL_AGILE<version>_<phys>List_<strip>_<point>_<n_in><part_type>_<energy>MeV_<theta>_<phi>.all.dat
; ASCII Columns:
; - c1 = event ID
; - c2 = bar plane 
; - c3 = bar id
; - c4 = 0
; - c5 = energy A
; - c6 = energy B

gams_cal_event_id = -1l
gams_cal_bar_plane = -1l
gams_cal_bar_id= -1l
gams_cal_c4 = -1l
gams_cal_ene_a = -1.
gams_cal_ene_b = -1.

; - G4_GAMS_AC_AGILE<version>_<phys>List_<strip>_<point>_<n_in><part_type>_<energy>MeV_<theta>_<phi>.all.dat
; ASCII Columns:
; - c1 = event ID
; - c2 = AC panel
; - c3 = AC subpanel
; - c4 = energy deposit

gams_ac_event_id = -1l
gams_ac_panel = ''
gams_ac_subpanel = -1l
gams_ac_edep = -1.

filepath = './AGILE'+agile_version+sdir+'/theta'+strtrim(string(theta_type),1)+'/'+stripDir+py_dir+'/'+sim_name+'/'+ene_type+'MeV/'+strtrim(string(N_in),1)+part_type+'/'
print, 'LEVEL0 file path: ', filepath

for ifile=0, n_files-1 do begin

 	
  if (out_type EQ 0) then begin 
    filenamedat_raw_kalman = filepath+'G4.RAW.KALMAN.AGILE'+agile_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+strmid(strtrim(string(N_in),1),0,10)+part_type+'.'+ene_type+'MeV.'+strmid(strtrim(string(theta_type),1),0,10)+'.'+strmid(strtrim(string(phi_type),1),0,10)+'.'+strtrim(string(ifile),1)+'.dat'   
    readcol, filenamedat_raw_kalman, c1, c2, c3, c4, c5, c6, c7, format='(l,l,d,d,l,l,d)'
    raw_kalman_event_id = [raw_kalman_event_id, c1]
    raw_kalman_Si_id = [raw_kalman_Si_id, c1]
    raw_kalman_pos = [raw_kalman_pos, c2]
    raw_kalman_zpos = [raw_kalman_zpos, c3]
    raw_kalman_plane_id = [raw_kalman_plane_id, c4]
    raw_kalman_strip_id = [raw_kalman_strip_id, c5]
    raw_kalman_edep = [raw_kalman_edep, c6]  
    
    filenamedat_raw_general = filepath+'G4.RAW.GENERAL.AGILE'+agile_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+strmid(strtrim(string(N_in),1),0,10)+part_type+'.'+ene_type+'MeV.'+strmid(strtrim(string(theta_type),1),0,10)+'.'+strmid(strtrim(string(phi_type),1),0,10)+'.'+strtrim(string(ifile),1)+'.dat'   
    readcol, filenamedat_raw_general, c1, c2, c3, c4, c5, c6, c7, c8, format='(l,l,d,d,l,l,l,d)'
    raw_general_event_id = [raw_general_event_id, c1]
    raw_general_Si_id = [raw_general_Si_id, c2]
    raw_general_pos = [raw_general_pos, c3]
    raw_general_zpos = [raw_general_zpos, c4]
    raw_general_tray_id = [raw_general_tray_id, c5]
    raw_general_plane_id = [raw_general_plane_id, c6]
    raw_general_strip_id = [raw_general_strip_id, c7]
    raw_general_edep = [raw_general_edep, c8]  

    filenamedat_digi_kalman = filepath+'G4.DIGI.KALMAN.AGILE'+agile_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+strmid(strtrim(string(N_in),1),0,10)+part_type+'.'+ene_type+'MeV.'+strmid(strtrim(string(theta_type),1),0,10)+'.'+strmid(strtrim(string(phi_type),1),0,10)+'.'+strtrim(string(ifile),1)+'.dat'   
    readcol, filenamedat_digi_kalman, c1, c2, c3, c4, c5, c6, c7, format='(l,l,d,d,l,l,d)'
    digi_kalman_event_id = [digi_kalman_event_id, c1]
    digi_kalman_Si_id = [digi_kalman_Si_id, c1]
    digi_kalman_pos = [digi_kalman_pos, c2]
    digi_kalman_zpos = [digi_kalman_zpos, c3]
    digi_kalman_plane_id = [digi_kalman_plane_id, c4]
    digi_kalman_strip_id = [digi_kalman_strip_id, c5]
    digi_kalman_edep = [digi_kalman_edep, c6]  
    
    filenamedat_digi_general = filepath+'G4.DIGI.GENERAL.AGILE'+agile_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+strmid(strtrim(string(N_in),1),0,10)+part_type+'.'+ene_type+'MeV.'+strmid(strtrim(string(theta_type),1),0,10)+'.'+strmid(strtrim(string(phi_type),1),0,10)+'.'+strtrim(string(ifile),1)+'.dat'   
    readcol, filenamedat_digi_general, c1, c2, c3, c4, c5, c6, c7, c8, format='(l,l,d,d,l,l,l,d)'
    digi_general_event_id = [digi_general_event_id, c1]
    digi_general_Si_id = [digi_general_Si_id, c2]
    digi_general_pos = [digi_general_pos, c3]
    digi_general_zpos = [digi_general_zpos, c4]
    digi_general_tray_id = [digi_general_tray_id, c5]
    digi_general_plane_id = [digi_general_plane_id, c6]
    digi_general_strip_id = [digi_general_strip_id, c7]
    digi_general_edep = [digi_general_edep, c8]
    
    filenamedat_aiv_xplane = filepath+'G4_AIV_XPLANE_AGILE'+agile_version+'_'+py_name+'_'+sim_name+'_'+stripname+'_'+sname+'_'+strmid(strtrim(string(N_in),1),0,10)+part_type+'_'+ene_type+'MeV_'+strmid(strtrim(string(theta_type),1),0,10)+'_'+strmid(strtrim(string(phi_type),1),0,10)+'.'+strtrim(string(ifile),1)+'.dat'   
    readcol, filenamedat_aiv_xplane, c1, c2, c3, c4, c5, c6, c7, format='(l,l,l,l,l,d,l)' 
    aiv_xplane_event_id = [aiv_xplane_event_id, c1]
    aiv_xplane_plane_id = [aiv_xplane_plane_id, c2]
    aiv_xplane_readstreap_id = [aiv_xplane_readstreap_id, c3]
    aiv_xplane_c4 = [aiv_xplane_c4, c4]
    aiv_xplane_c5 = [aiv_xplane_c5, c5]
    aiv_xplane_edep = [aiv_xplane_edep, c6]
    aiv_xplane_c7 = [aiv_xplane_c7, c7]
    
    filenamedat_aiv_yplane = filepath+'G4_AIV_YPLANE_AGILE'+agile_version+'_'+py_name+'_'+sim_name+'_'+stripname+'_'+sname+'_'+strmid(strtrim(string(N_in),1),0,10)+part_type+'_'+ene_type+'MeV_'+strmid(strtrim(string(theta_type),1),0,10)+'_'+strmid(strtrim(string(phi_type),1),0,10)+'.'+strtrim(string(ifile),1)+'.dat'   
    readcol, filenamedat_aiv_yplane, c1, c2, c3, c4, c5, c6, c7, format='(l,l,l,l,l,d,l)' 
    aiv_yplane_event_id = [aiv_yplane_event_id, c1]
    aiv_yplane_plane_id = [aiv_yplane_plane_id, c2]
    aiv_yplane_readstreap_id = [aiv_yplane_readstreap_id, c3]
    aiv_yplane_c4 = [aiv_yplane_c4, c4]
    aiv_yplane_c5 = [aiv_yplane_c5, c5]
    aiv_yplane_edep = [aiv_yplane_edep, c6]
    aiv_yplane_c7 = [aiv_yplane_c7, c7]    
      
  endif

    filenamedat_aa_strip = filepath+'/G4.AA.STRIP.AGILE'+agile_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+strmid(strtrim(string(N_in),1),0,10)+part_type+'.'+ene_type+'MeV.'+strmid(strtrim(string(theta_type),1),0,10)+'.'+strmid(strtrim(string(phi_type),1),0,10)+'.'+strtrim(string(ifile),1)+'.dat'
    readcol, filenamedat_aa_strip, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, format='(l,l,l,l,l,f,l,l,f,f,l)'
    print, filenamedat_aa_strip
   
     aa_strip_event_id = [aa_strip_event_id, c1]
     aa_strip_theta_in = [aa_strip_theta_in, c2]
     aa_strip_phi_in = [aa_strip_phi_in, c3]
     aa_strip_ene_in = [aa_strip_ene_in, c4]
     aa_strip_plane_id = [aa_strip_plane_id, c5]
     aa_strip_zpos = [aa_strip_zpos, c6]
     aa_strip_si_id = [aa_strip_si_id, c7]
     aa_strip_strip_id = [aa_strip_strip_id, c8]
     aa_strip_pos = [aa_strip_pos, c9]
     aa_strip_edep = [aa_strip_edep, c10]
     aa_strip_pair = [aa_strip_pair, c11]
     
    filenamedat_aa_kalman = filepath+'/G4.AA.CLUSTER.AGILE'+agile_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+strmid(strtrim(string(N_in),1),0,10)+part_type+'.'+ene_type+'MeV.'+strmid(strtrim(string(theta_type),1),0,10)+'.'+strmid(strtrim(string(phi_type),1),0,10)+'.'+strtrim(string(ifile),1)+'.dat'
    readcol, filenamedat_aa_kalman, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, format='(l,l,l,l,l,f,l,f,f,l,l)'
    print, filenamedat_aa_kalman   
     aa_kalman_event_id = [aa_kalman_event_id, c1]
     aa_kalman_theta_in = [aa_kalman_theta_in, c2]
     aa_kalman_phi_in = [aa_kalman_phi_in, c3]
     aa_kalman_ene_in = [aa_kalman_ene_in, c4]
     aa_kalman_plane_id = [aa_kalman_plane_id, c5]
     aa_kalman_zpos = [aa_kalman_zpos, c6]
     aa_kalman_si_id = [aa_kalman_si_id, c7]
     aa_kalman_pos = [aa_kalman_pos, c8]
     aa_kalman_edep = [aa_kalman_edep, c9]
     aa_kalman_strip_number = [aa_kalman_strip_number, c10]
     aa_kalman_pair = [aa_kalman_pair, c11]

    
    filenamedat_gams_xplane = filepath+'G4_GAMS_XPLANE_AGILE'+agile_version+'_'+py_name+'_'+sim_name+'_'+stripname+'_'+sname+'_'+strmid(strtrim(string(N_in),1),0,10)+part_type+'_'+ene_type+'MeV_'+strmid(strtrim(string(theta_type),1),0,10)+'_'+strmid(strtrim(string(phi_type),1),0,10)+'.'+strtrim(string(ifile),1)+'.dat'   
    readcol, filenamedat_gams_xplane, c1, c2, c3, c4, c5, c6, c7, format='(l,l,l,l,l,d,l)' 
    gams_xplane_event_id = [gams_xplane_event_id, c1]
    gams_xplane_plane_id = [gams_xplane_plane_id, c2]
    gams_xplane_readstreap_id = [gams_xplane_readstreap_id, c3]
    gams_xplane_c4 = [gams_xplane_c4, c4]
    gams_xplane_c5 = [gams_xplane_c5, c5]
    gams_xplane_edep = [gams_xplane_edep, c6]
    gams_xplane_c7 = [gams_xplane_c7, c7]
    
    filenamedat_gams_yplane = filepath+'G4_GAMS_YPLANE_AGILE'+agile_version+'_'+py_name+'_'+sim_name+'_'+stripname+'_'+sname+'_'+strmid(strtrim(string(N_in),1),0,10)+part_type+'_'+ene_type+'MeV_'+strmid(strtrim(string(theta_type),1),0,10)+'_'+strmid(strtrim(string(phi_type),1),0,10)+'.'+strtrim(string(ifile),1)+'.dat'   
    readcol, filenamedat_gams_yplane, c1, c2, c3, c4, c5, c6, c7, format='(l,l,l,l,l,d,l)' 
    gams_yplane_event_id = [gams_yplane_event_id, c1]
    gams_yplane_plane_id = [gams_yplane_plane_id, c2]
    gams_yplane_readstreap_id = [gams_yplane_readstreap_id, c3]
    gams_yplane_c4 = [gams_yplane_c4, c4]
    gams_yplane_c5 = [gams_yplane_c5, c5]
    gams_yplane_edep = [gams_yplane_edep, c6]
    gams_yplane_c7 = [gams_yplane_c7, c7]
    
    if (cal_flag) then begin
		filenamedat_gams_cal = filepath+'G4_GAMS_CAL_AGILE'+agile_version+'_'+py_name+'_'+sim_name+'_'+stripname+'_'+sname+'_'+strmid(strtrim(string(N_in),1),0,10)+part_type+'_'+ene_type+'MeV_'+strmid(strtrim(string(theta_type),1),0,10)+'_'+strmid(strtrim(string(phi_type),1),0,10)+'.'+strtrim(string(ifile),1)+'.dat'   
		readcol, filenamedat_gams_cal, c1, c2, c3, c4, c5, c6, format='(l,l,l,l,d,d)' 
		gams_cal_event_id = [gams_cal_event_id, c1]
		gams_cal_bar_plane = [gams_cal_bar_plane, c2]
		gams_cal_bar_id= [gams_cal_bar_id, c3]
		gams_cal_c4 = [gams_cal_c4, c4]
		gams_cal_ene_a = [gams_cal_ene_a, c5]
		gams_cal_ene_b = [gams_cal_ene_b, c6]
	endif
	if (ac_flag) then begin
		filenamedat_gams_ac = filepath+'G4_GAMS_AC_AGILE'+agile_version+'_'+py_name+'_'+sim_name+'_'+stripname+'_'+sname+'_'+strmid(strtrim(string(N_in),1),0,10)+part_type+'_'+ene_type+'MeV_'+strmid(strtrim(string(theta_type),1),0,10)+'_'+strmid(strtrim(string(phi_type),1),0,10)+'.'+strtrim(string(ifile),1)+'.dat'   
		readcol, filenamedat_gams_ac, c1, c2, c3, c4, format='(l,a,l,d)' 
		gams_ac_event_id = [gams_ac_event_id, c1]
		gams_ac_panel = [gams_ac_panel, c2]
		gams_ac_subpanel = [gams_ac_subpanel, c3]
		gams_ac_edep = [gams_ac_edep, c4]
    endif

   if (out_type EQ 0) then begin 
    filenamefits_raw = filepath+'G4.RAW.AGILE'+agile_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+strmid(strtrim(string(N_in),1),0,10)+part_type+'.'+ene_type+'MeV.'+strmid(strtrim(string(theta_type),1),0,10)+'.'+strmid(strtrim(string(phi_type),1),0,10)+'.'+strtrim(string(ifile),1)+'.fits'   
    print, filenamefits_raw
    struct_raw = mrdfits(filenamefits_raw,$ 
                     1, $
                     structyp = 'raw', $
                     /unsigned)

    rawData_event_id = [rawData_event_id, struct_raw.EVT_ID]
    rawData_Si_id = [rawData_Si_id, struct_raw.TRK_FLAG]
    rawData_tray_id = [rawData_tray_id, struct_raw.TRAY_ID]
    rawData_plane_id = [rawData_plane_id, struct_raw.PLANE_ID]
    rawData_Strip_id = [rawData_Strip_id, struct_raw.STRIP_ID]
    rawData_energy_dep = [rawData_energy_dep, struct_raw.E_DEP]
    rawData_ent_x = [rawData_ent_x, struct_raw.X_ENT]
    rawData_ent_y = [rawData_ent_y, struct_raw.Y_ENT]
    rawData_ent_z = [rawData_ent_z, struct_raw.Z_ENT]
    rawData_exit_x = [rawData_exit_x, struct_raw.X_EXIT]
    rawData_exit_y = [rawData_exit_y, struct_raw.Y_EXIT]
    rawData_exit_z = [rawData_exit_z, struct_raw.Z_EXIT]

    filenamefits_first = filepath+'G4.FIRST.AGILE'+agile_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+strmid(strtrim(string(N_in),1),0,10)+part_type+'.'+ene_type+'MeV.'+strmid(strtrim(string(theta_type),1),0,10)+'.'+strmid(strtrim(string(phi_type),1),0,10)+'.'+strtrim(string(ifile),1)+'.fits'
    print, filenamefits_first
    struct_first = mrdfits(filenamefits_first,$
      1, $
      structyp = 'first', $
      /unsigned)

    firstData_event_id = [firstData_event_id, struct_first.EVT_ID]
    firstData_Si_id = [firstData_Si_id, struct_first.TRK_FLAG]
    firstData_tray_id = [firstData_tray_id, struct_first.TRAY_ID]
    firstData_plane_id = [firstData_plane_id, struct_first.PLANE_ID]
    firstData_Strip_id = [firstData_Strip_id, struct_first.STRIP_ID]
    firstData_energy_dep = [firstData_energy_dep, struct_first.E_DEP]
    firstData_ent_x = [firstData_ent_x, struct_first.X_ENT]
    firstData_ent_y = [firstData_ent_y, struct_first.Y_ENT]
    firstData_ent_z = [firstData_ent_z, struct_first.Z_ENT]
    firstData_exit_x = [firstData_exit_x, struct_first.X_EXIT]
    firstData_exit_y = [firstData_exit_y, struct_first.Y_EXIT]
    firstData_exit_z = [firstData_exit_z, struct_first.Z_EXIT]
    
    filenamefits_l0 = filepath+'L0.AGILE'+agile_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+STRMID(STRTRIM(STRING(N_IN),1),0,10)+part_type+'.'+ene_type+'MeV.'+STRMID(STRTRIM(STRING(THETA_TYPE),1),0,10)+'.'+STRMID(STRTRIM(STRING(PHI_TYPE),1),0,10)+'.'+strtrim(string(ifile),1)+'.fits'  
    struct_l0 = mrdfits(filenamefits_l0,$ 
                     1, $
                     structyp = 'l0', $
                     /unsigned)

    L0TRACKERGLOBAL_Glob_event_id_test = [L0TRACKERGLOBAL_Glob_event_id_test, struct_l0.EVT_ID]
    L0TRACKERGLOBAL_Glob_vol_id_test = [L0TRACKERGLOBAL_Glob_vol_id_test, struct_l0.VOLUME_ID]
    L0TRACKERGLOBAL_Glob_moth_id_test = [L0TRACKERGLOBAL_Glob_moth_id_test, struct_l0.MOTHER_ID]
    L0TRACKERGLOBAL_Glob_tray_id_test = [L0TRACKERGLOBAL_Glob_tray_id_test, struct_l0.TRAY_ID]
    L0TRACKERGLOBAL_Glob_plane_id_test = [L0TRACKERGLOBAL_Glob_plane_id_test, struct_l0.PLANE_ID]
    L0TRACKERGLOBAL_Glob_Si_id_test = [L0TRACKERGLOBAL_Glob_Si_id_test, struct_l0.TRK_FLAG]    
    L0TRACKERGLOBAL_Glob_Strip_id_test = [L0TRACKERGLOBAL_Glob_Strip_id_test, struct_l0.STRIP_ID]
    L0TRACKERGLOBAL_Glob_Strip_type_test = [L0TRACKERGLOBAL_Glob_Strip_type_test, struct_l0.STRIP_TYPE]
    L0TRACKERGLOBAL_Glob_pos_test = [L0TRACKERGLOBAL_Glob_pos_test, struct_l0.POS]
    L0TRACKERGLOBAL_Glob_zpos_test = [L0TRACKERGLOBAL_Glob_zpos_test, struct_l0.ZPOS]
    L0TRACKERGLOBAL_Glob_energy_dep_test = [L0TRACKERGLOBAL_Glob_energy_dep_test, struct_l0.E_DEP]

    filenamefits_l05 = filepath+'L0.5.DIGI.AGILE'+agile_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+STRMID(STRTRIM(STRING(N_IN),1),0,10)+part_type+'.'+ene_type+'MeV.'+STRMID(STRTRIM(STRING(THETA_TYPE),1),0,10)+'.'+STRMID(STRTRIM(STRING(PHI_TYPE),1),0,10)+'.'+strtrim(string(ifile),1)+'.fits'   
    struct_l05 = mrdfits(filenamefits_l05,$ 
                     1, $
                     structyp = 'l05', $
                     /unsigned)

    L05TRACKERGLOBAL_Glob_event_id_acap = [L05TRACKERGLOBAL_Glob_event_id_acap, struct_l05.EVT_ID]
    L05TRACKERGLOBAL_Glob_vol_id_acap = [L05TRACKERGLOBAL_Glob_vol_id_acap, struct_l05.VOLUME_ID]
    L05TRACKERGLOBAL_Glob_moth_id_acap = [L05TRACKERGLOBAL_Glob_moth_id_acap, struct_l05.MOTHER_ID]
    L05TRACKERGLOBAL_Glob_tray_id_acap = [L05TRACKERGLOBAL_Glob_tray_id_acap, struct_l05.TRAY_ID]
    L05TRACKERGLOBAL_Glob_plane_id_acap = [L05TRACKERGLOBAL_Glob_plane_id_acap, struct_l05.PLANE_ID]
    L05TRACKERGLOBAL_Glob_Si_id_acap = [L05TRACKERGLOBAL_Glob_Si_id_acap, struct_l05.TRK_FLAG]    
    L05TRACKERGLOBAL_Glob_Strip_id_acap = [L05TRACKERGLOBAL_Glob_Strip_id_acap, struct_l05.STRIP_ID]
    L05TRACKERGLOBAL_Glob_Strip_type_acap = [L05TRACKERGLOBAL_Glob_Strip_type_acap, struct_l05.STRIP_TYPE]
    L05TRACKERGLOBAL_Glob_pos_acap = [L05TRACKERGLOBAL_Glob_pos_acap, struct_l05.POS]
    L05TRACKERGLOBAL_Glob_zpos_acap = [L05TRACKERGLOBAL_Glob_zpos_acap, struct_l05.ZPOS]
    L05TRACKERGLOBAL_Glob_energy_dep_acap = [L05TRACKERGLOBAL_Glob_energy_dep_acap, struct_l05.E_DEP]

	if (cal_flag) then begin
		filenamefits_cal = filepath+'G4.CAL.AGILE'+agile_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+strmid(strtrim(string(N_in),1),0,10)+part_type+'.'+ene_type+'MeV.'+strmid(strtrim(string(theta_type),1),0,10)+'.'+strmid(strtrim(string(phi_type),1),0,10)+'.'+strtrim(string(ifile),1)+'.fits'   
		print, filenamefits_cal
		struct_cal = mrdfits(filenamefits_cal,$ 
						 1, $
						 structyp = 'cal', $
						 /unsigned)

		calInput_event_id_tot_cal = [calInput_event_id_tot_cal, struct_cal.EVT_ID]
		calInput_bar_plane_tot = [calInput_bar_plane_tot, struct_cal.BAR_PLANE]
		calInput_bar_id_tot = [calInput_bar_id_tot, struct_cal.BAR_ID]
		calInput_ene_a_tot = [calInput_ene_a_tot, struct_cal.ENERGY_A]
		calInput_ene_b_tot = [calInput_ene_b_tot, struct_cal.ENERGY_B]
	endif
	if (ac_flag) then begin
		filenamefits_ac = filepath+'G4.AC.AGILE'+agile_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+strmid(strtrim(string(N_in),1),0,10)+part_type+'.'+ene_type+'MeV.'+strmid(strtrim(string(theta_type),1),0,10)+'.'+strmid(strtrim(string(phi_type),1),0,10)+'.'+strtrim(string(ifile),1)+'.fits'   
		struct_ac = mrdfits(filenamefits_ac,$ 
						 1, $
						 structyp = 'ac', $
						 /unsigned)

		acInput_event_id_tot_ac = [acInput_event_id_tot_ac, struct_ac.EVT_ID]
		acInput_AC_panel = [acInput_AC_panel, struct_ac.AC_PANEL]
		acInput_AC_subpanel = [acInput_AC_subpanel, struct_ac.AC_SUBPANEL]
		acInput_energy_dep_tot_ac = [acInput_energy_dep_tot_ac, struct_ac.E_DEP]
    endif
 endif
endfor

; -----> ASCII files

if (out_type EQ 0) then begin
	raw_kalman_event_id = raw_kalman_event_id[1:*]
	raw_kalman_Si_id = raw_kalman_Si_id[1:*]
	raw_kalman_pos = raw_kalman_pos[1:*]
	raw_kalman_zpos = raw_kalman_zpos[1:*]
	raw_kalman_plane_id = raw_kalman_plane_id[1:*]
	raw_kalman_strip_id = raw_kalman_strip_id[1:*]
	raw_kalman_edep = raw_kalman_edep[1:*]    

	openw,lun,filepath+'G4.RAW.KALMAN.AGILE'+agile_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+strmid(strtrim(string(N_in),1),0,10)+part_type+'.'+ene_type+'MeV.'+strmid(strtrim(string(theta_type),1),0,10)+'.'+strmid(strtrim(string(phi_type),1),0,10)+'.all.dat',/get_lun
	; ASCII Columns:
	; - c1 = Silicon layer ID
	; - c2 = x/y pos [cm]
	; - c3 = z pos [cm]
	; - c4 = plane ID
	; - c5 = strip ID
	; - c6 = energy dep [keV]

	event_start = -1
	j=0l
	while (1) do begin
		where_event_eq = where(raw_kalman_event_id EQ raw_kalman_event_id(j))
		raw_kalman_Si_id_temp = raw_kalman_Si_id(where_event_eq)
		raw_kalman_pos_temp = raw_kalman_pos(where_event_eq)
		raw_general_zpos_temp  = raw_general_zpos(where_event_eq)
		raw_kalman_plane_id_temp  = raw_kalman_plane_id(where_event_eq)
		raw_kalman_strip_id_temp = raw_kalman_strip_id(where_event_eq)
		raw_kalman_edep_temp = raw_kalman_edep(where_event_eq)

		printf, lun, '; Event:', raw_kalman_event_id(j)
		printf, lun, '; ', theta_type, phi_type, ene_type   
	
		for r=0l, n_elements(raw_kalman_Si_id_temp)-1 do begin
			  printf, lun, raw_kalman_Si_id_temp(r), raw_kalman_pos_temp(r), raw_general_zpos_temp(r), raw_kalman_plane_id_temp(r), raw_kalman_strip_id_temp(r), raw_kalman_edep_temp(r), format='(I6,2x,F10.5,2x,F10.5,2x,I5,2x,I5,2x,F10.5)'
		endfor
		N_event_eq = n_elements(where_event_eq)
		if where_event_eq(N_event_eq-1) LT (n_elements(raw_kalman_event_id)-1) then begin
		  j = where_event_eq(N_event_eq-1)+1
		endif else break
	endwhile

	Free_lun, lun

	raw_general_event_id = raw_general_event_id[1:*]
	raw_general_Si_id = raw_general_Si_id[1:*]
	raw_general_pos = raw_general_pos[1:*]
	raw_general_zpos = raw_general_zpos[1:*]
	raw_general_tray_id = raw_general_tray_id[1:*]
	raw_general_plane_id = raw_general_plane_id[1:*]
	raw_general_strip_id = raw_general_strip_id[1:*]
	raw_general_edep = raw_general_edep[1:*]    

	openw,lun,filepath+'G4.RAW.GENERAL.AGILE'+agile_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+strmid(strtrim(string(N_in),1),0,10)+part_type+'.'+ene_type+'MeV.'+strmid(strtrim(string(theta_type),1),0,10)+'.'+strmid(strtrim(string(phi_type),1),0,10)+'.all.dat',/get_lun
	; ASCII Columns:
	; - c1 = event ID
	; - c2 = Silicon layer ID
	; - c3 = x/y pos [cm]
	; - c4 = z pos [cm]
	; - c5 = tray ID
	; - c6 = plane ID
	; - c7 = strip ID 
	; - c8 = energy dep [keV]    

	event_start = -1
	j=0l
	while (1) do begin
		where_event_eq = where(raw_general_event_id EQ raw_general_event_id(j))
		raw_general_Si_id_temp = raw_general_Si_id(where_event_eq)
		raw_general_pos_temp = raw_general_pos(where_event_eq)
		raw_general_zpos_temp  = raw_general_zpos(where_event_eq)
		raw_general_tray_id_temp  = raw_general_tray_id(where_event_eq)
		raw_general_plane_id_temp = raw_general_plane_id(where_event_eq)    
		raw_general_strip_id_temp = raw_general_strip_id(where_event_eq)
		raw_general_edep_temp = raw_general_edep(where_event_eq)

		for r=0l, n_elements(raw_general_Si_id_temp)-1 do begin
			 printf, lun, raw_general_event_id(j), raw_general_Si_id_temp(r), raw_general_pos_temp(r), raw_general_zpos_temp(r), raw_general_tray_id_temp(r), raw_general_plane_id_temp(r), raw_general_strip_id_temp(r), raw_general_edep_temp(r), format='(I6,2x,I5,2x,F10.5,2x,F10.5,2x,I5,2x,I5,2x,I5,2x,F10.5)'
		endfor
	
		N_event_eq = n_elements(where_event_eq)
		if where_event_eq(N_event_eq-1) LT (n_elements(raw_general_event_id)-1) then begin
		  j = where_event_eq(N_event_eq-1)+1
		endif else break
	endwhile

	Free_lun, lun


	digi_kalman_event_id = digi_kalman_event_id[1:*]
	digi_kalman_Si_id = digi_kalman_Si_id[1:*]
	digi_kalman_pos = digi_kalman_pos[1:*]
	digi_kalman_zpos = digi_kalman_zpos[1:*]
	digi_kalman_plane_id = digi_kalman_plane_id[1:*]
	digi_kalman_strip_id = digi_kalman_strip_id[1:*]
	digi_kalman_edep = digi_kalman_edep[1:*]    

	openw,lun,filepath+'G4.DIGI.KALMAN.AGILE'+agile_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+strmid(strtrim(string(N_in),1),0,10)+part_type+'.'+ene_type+'MeV.'+strmid(strtrim(string(theta_type),1),0,10)+'.'+strmid(strtrim(string(phi_type),1),0,10)+'.all.dat',/get_lun
	; ASCII Columns:
	; - c1 = Silicon layer ID
	; - c2 = x/y pos [cm]
	; - c3 = z pos [cm]
	; - c4 = plane ID
	; - c5 = strip ID
	; - c6 = energy dep [keV]

	event_start = -1
	j=0l
	while (1) do begin
		where_event_eq = where(digi_kalman_event_id EQ digi_kalman_event_id(j))
		digi_kalman_Si_id_temp = digi_kalman_Si_id(where_event_eq)
		digi_kalman_pos_temp = digi_kalman_pos(where_event_eq)
		digi_general_zpos_temp  = digi_general_zpos(where_event_eq)
		digi_kalman_plane_id_temp  = digi_kalman_plane_id(where_event_eq)
		digi_kalman_strip_id_temp = digi_kalman_strip_id(where_event_eq)
		digi_kalman_edep_temp = digi_kalman_edep(where_event_eq)

		printf, lun, '; Event:', digi_kalman_event_id(j)
		printf, lun, '; ', theta_type, phi_type, ene_type   
	
		for r=0l, n_elements(digi_kalman_Si_id_temp)-1 do begin
			  printf, lun, digi_kalman_Si_id_temp(r), digi_kalman_pos_temp(r), digi_general_zpos_temp(r), digi_kalman_plane_id_temp(r), digi_kalman_strip_id_temp(r), digi_kalman_edep_temp(r), format='(I6,2x,F10.5,2x,F10.5,2x,I5,2x,I5,2x,F10.5)'
		endfor
		N_event_eq = n_elements(where_event_eq)
		if where_event_eq(N_event_eq-1) LT (n_elements(digi_kalman_event_id)-1) then begin
		  j = where_event_eq(N_event_eq-1)+1
		endif else break
	endwhile

	Free_lun, lun

	digi_general_event_id = digi_general_event_id[1:*]
	digi_general_Si_id = digi_general_Si_id[1:*]
	digi_general_pos = digi_general_pos[1:*]
	digi_general_zpos = digi_general_zpos[1:*]
	digi_general_tray_id = digi_general_tray_id[1:*]
	digi_general_plane_id = digi_general_plane_id[1:*]
	digi_general_strip_id = digi_general_strip_id[1:*]
	digi_general_edep = digi_general_edep[1:*]    

	openw,lun,filepath+'G4.DIGI.GENERAL.AGILE'+agile_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+strmid(strtrim(string(N_in),1),0,10)+part_type+'.'+ene_type+'MeV.'+strmid(strtrim(string(theta_type),1),0,10)+'.'+strmid(strtrim(string(phi_type),1),0,10)+'.all.dat',/get_lun
	; ASCII Columns:
	; - c1 = event ID
	; - c2 = Silicon layer ID
	; - c3 = x/y pos [cm]
	; - c4 = z pos [cm]
	; - c5 = tray ID
	; - c6 = plane ID
	; - c7 = strip ID 
	; - c8 = energy dep [keV]    

	event_start = -1
	j=0l
	while (1) do begin
		where_event_eq = where(digi_general_event_id EQ digi_general_event_id(j))
		digi_general_Si_id_temp = digi_general_Si_id(where_event_eq)
		digi_general_pos_temp = digi_general_pos(where_event_eq)
		digi_general_zpos_temp  = digi_general_zpos(where_event_eq)
		digi_general_tray_id_temp  = digi_general_tray_id(where_event_eq)
		digi_general_plane_id_temp = digi_general_plane_id(where_event_eq)    
		digi_general_strip_id_temp = digi_general_strip_id(where_event_eq)
		digi_general_edep_temp = digi_general_edep(where_event_eq)

		for r=0l, n_elements(digi_general_Si_id_temp)-1 do begin
			 printf, lun, digi_general_event_id(j), digi_general_Si_id_temp(r), digi_general_pos_temp(r), digi_general_zpos_temp(r), digi_general_tray_id_temp(r), digi_general_plane_id_temp(r), digi_general_strip_id_temp(r), digi_general_edep_temp(r), format='(I6,2x,I5,2x,F10.5,2x,F10.5,2x,I5,2x,I5,2x,I5,2x,F10.5)'
		endfor
	
		N_event_eq = n_elements(where_event_eq)
		if where_event_eq(N_event_eq-1) LT (n_elements(digi_general_event_id)-1) then begin
		  j = where_event_eq(N_event_eq-1)+1
		endif else break
	endwhile

	Free_lun, lun
	
	aiv_xplane_event_id = aiv_xplane_event_id[1:*]
	aiv_xplane_plane_id = aiv_xplane_plane_id[1:*]
	aiv_xplane_readstreap_id = aiv_xplane_readstreap_id[1:*]
	aiv_xplane_c4 = aiv_xplane_c4[1:*]
	aiv_xplane_c5 = aiv_xplane_c5[1:*]
	aiv_xplane_edep = aiv_xplane_edep[1:*]
	aiv_xplane_c7 = aiv_xplane_c7[1:*]

	openw,lun,filepath+'G4_AIV_XPLANE_AGILE'+agile_version+'_'+py_name+'_'+sim_name+'_'+stripname+'_'+sname+'_'+strmid(strtrim(string(N_in),1),0,10)+part_type+'_'+ene_type+'MeV_'+strmid(strtrim(string(theta_type),1),0,10)+'_'+strmid(strtrim(string(phi_type),1),0,10)+'.all.dat',/get_lun
	; ASCII Columns:
	; - c1 = event ID
	; - c2 = plane ID
	; - c3 = readout strip ID
	; - c4 = -999
	; - c5 = -999
	; - c6 = energy dep in keV
	; - c7 = -999    event_start = -1
	j=0l
	while (1) do begin
		where_event_eq = where(aiv_xplane_event_id EQ aiv_xplane_event_id(j))
		aiv_xplane_plane_id_temp = aiv_xplane_plane_id(where_event_eq)
		aiv_xplane_readstreap_id_temp = aiv_xplane_readstreap_id(where_event_eq)
		aiv_xplane_c4_temp  = aiv_xplane_c4(where_event_eq)
		aiv_xplane_c5_temp  = aiv_xplane_c5(where_event_eq)
		aiv_xplane_edep_temp = aiv_xplane_edep(where_event_eq)    
		aiv_xplane_c7_temp = aiv_xplane_c7(where_event_eq)    
	
		for r=0l, n_elements(aiv_xplane_plane_id_temp)-1 do begin
			printf, lun, aiv_xplane_event_id(j), aiv_xplane_plane_id_temp(r), aiv_xplane_readstreap_id_temp(r), aiv_xplane_c4_temp(r), aiv_xplane_c5_temp(r), aiv_xplane_edep_temp(r), aiv_xplane_c7_temp(r), format='(I6,I5,I5,I5,I5,F10.5,I5)'        
		endfor

		N_event_eq = n_elements(where_event_eq)
		if where_event_eq(N_event_eq-1) LT (n_elements(aiv_xplane_event_id)-1) then begin
		  j = where_event_eq(N_event_eq-1)+1
		endif else break
	endwhile

	Free_lun, lun

	aiv_yplane_event_id = aiv_yplane_event_id[1:*]
	aiv_yplane_plane_id = aiv_yplane_plane_id[1:*]
	aiv_yplane_readstreap_id = aiv_yplane_readstreap_id[1:*]
	aiv_yplane_c4 = aiv_yplane_c4[1:*]
	aiv_yplane_c5 = aiv_yplane_c5[1:*]
	aiv_yplane_edep = aiv_yplane_edep[1:*]
	aiv_yplane_c7 = aiv_yplane_c7[1:*]

	openw,lun,filepath+'G4_AIV_YPLANE_AGILE'+agile_version+'_'+py_name+'_'+sim_name+'_'+stripname+'_'+sname+'_'+strmid(strtrim(string(N_in),1),0,10)+part_type+'_'+ene_type+'MeV_'+strmid(strtrim(string(theta_type),1),0,10)+'_'+strmid(strtrim(string(phi_type),1),0,10)+'.all.dat',/get_lun
	; ASCII Columns:
	; - c1 = event ID
	; - c2 = plane ID
	; - c3 = readout strip ID
	; - c4 = -999
	; - c5 = -999
	; - c6 = energy dep in keV
	; - c7 = -999    event_start = -1
	j=0l
	while (1) do begin
		where_event_eq = where(aiv_yplane_event_id EQ aiv_yplane_event_id(j))
		aiv_yplane_plane_id_temp = aiv_yplane_plane_id(where_event_eq)
		aiv_yplane_readstreap_id_temp = aiv_yplane_readstreap_id(where_event_eq)
		aiv_yplane_c4_temp  = aiv_yplane_c4(where_event_eq)
		aiv_yplane_c5_temp  = aiv_yplane_c5(where_event_eq)
		aiv_yplane_edep_temp = aiv_yplane_edep(where_event_eq)    
		aiv_yplane_c7_temp = aiv_yplane_c7(where_event_eq)    
	
		for r=0l, n_elements(aiv_yplane_plane_id_temp)-1 do begin
			printf, lun, aiv_yplane_event_id(j), aiv_yplane_plane_id_temp(r), aiv_yplane_readstreap_id_temp(r), aiv_yplane_c4_temp(r), aiv_yplane_c5_temp(r), aiv_yplane_edep_temp(r), aiv_yplane_c7_temp(r), format='(I6,I5,I5,I5,I5,F10.5,I5)'        
		endfor

		N_event_eq = n_elements(where_event_eq)
		if where_event_eq(N_event_eq-1) LT (n_elements(aiv_yplane_event_id)-1) then begin
		  j = where_event_eq(N_event_eq-1)+1
		endif else break
	endwhile

	Free_lun, lun
endif

aa_strip_event_id = aa_strip_event_id[1:*]
aa_strip_theta_in = aa_strip_theta_in[1:*]
aa_strip_phi_in = aa_strip_phi_in[1:*]
aa_strip_ene_in = aa_strip_ene_in[1:*]
aa_strip_plane_id = aa_strip_plane_id[1:*]
aa_strip_zpos = aa_strip_zpos[1:*]
aa_strip_strip_id = aa_strip_strip_id[1:*]
aa_strip_si_id = aa_strip_si_id[1:*]
aa_strip_pos = aa_strip_pos[1:*]
aa_strip_edep = aa_strip_edep[1:*]
aa_strip_pair = aa_strip_pair[1:*]

aa_kalman_event_id = aa_kalman_event_id[1:*]
aa_kalman_theta_in = aa_kalman_theta_in[1:*]
aa_kalman_phi_in = aa_kalman_phi_in[1:*]
aa_kalman_ene_in = aa_kalman_ene_in[1:*]
aa_kalman_plane_id = aa_kalman_plane_id[1:*]
aa_kalman_zpos = aa_kalman_zpos[1:*]
aa_kalman_si_id = aa_kalman_si_id[1:*]
aa_kalman_pos = aa_kalman_pos[1:*]
aa_kalman_edep = aa_kalman_edep[1:*]
aa_kalman_strip_number = aa_kalman_strip_number[1:*]
aa_kalman_pair = aa_kalman_pair[1:*]

; filenamedat_aa_strip = filepath+sim_tag+'_STRIP_'+strmid(strtrim(string(N_in),1),0,10)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+strmid(strtrim(string(theta_type),1),0,10)+'_'+strmid(strtrim(string(phi_type),1),0,10)+'.'+pol_string+strtrim(string(ifile),1)+'.dat'

openw,lun,filepath+'G4.AA.STRIP.AGILE'+agile_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+strmid(strtrim(string(N_in),1),0,10)+part_type+'.'+ene_type+'MeV.'+strmid(strtrim(string(theta_type),1),0,10)+'.'+strmid(strtrim(string(phi_type),1),0,10)+'.all.dat',/get_lun
; ASCII Columns:
; - c1 = event ID
; - c2 = theta input
; - c3 = phi input
; - c4 = energy input
; - c5 = plane ID
; - c6 = Pos Z
; - c7 = X/Y flag (X = 0, Y = 1)
; - c8 = Strip ID
; - c9 = strip position (reference system center at the Silicon layer center)
; - c10 = energy deposition (keV)
; - c11 = Pair flag (1 = pair, 0 = no pair)

for r=0l, n_elements(aa_strip_event_id)-1 do begin
 printf, lun, aa_strip_event_id(r), aa_strip_theta_in(r), aa_strip_phi_in(r), aa_strip_ene_in(r), aa_strip_plane_id(r), aa_strip_zpos(r), aa_strip_si_id(r), aa_strip_strip_id(r), aa_strip_pos(r), aa_strip_edep(r), aa_strip_pair(r), format='(I7,I7,I7,I7,I7,F20.7,I7,I7,F20.7,F20.7,I7)'
endfor

Free_lun, lun

; filenamedat_aa_kalman = filepath+sim_tag+'_CLUSTER'+strmid(strtrim(string(N_in),1),0,10)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+strmid(strtrim(string(theta_type),1),0,10)+'_'+strmid(strtrim(string(phi_type),1),0,10)+'.'+pol_string+strtrim(string(ifile),1)+'.dat'

openw,lun,filepath+'G4.AA.CLUSTER.AGILE'+agile_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+strmid(strtrim(string(N_in),1),0,10)+part_type+'.'+ene_type+'MeV.'+strmid(strtrim(string(theta_type),1),0,10)+'.'+strmid(strtrim(string(phi_type),1),0,10)+'.all.dat',/get_lun
; ASCII Columns:
; - c1 = event ID
; - c2 = theta input
; - c3 = phi input
; - c4 = energy input
; - c5 = plane ID
; - c6 = Pos Z
; - c7 = X/Y flag (X = 0, Y = 1)
; - c8 = Cluster position (reference system center at the Silicon layer center)
; - c9 = energy deposition (keV)
; - c10 = number of strips composing the cluster
; - c11 = Pair flag (1 = pair, 0 = no pair)

for r=0l, n_elements(aa_kalman_event_id)-1 do begin
 printf, lun, aa_kalman_event_id(r), aa_kalman_theta_in(r), aa_kalman_phi_in(r), aa_kalman_ene_in(r), aa_kalman_plane_id(r), aa_kalman_zpos(r), aa_kalman_si_id(r), aa_kalman_pos(r), aa_kalman_edep(r), aa_kalman_strip_number(r), aa_kalman_pair(r), format='(I7,I7,I7,I7,I7,F20.7,I7,F20.7,F20.7,I7,I7)'
endfor

Free_lun, lun

   
gams_xplane_event_id = gams_xplane_event_id[1:*]
gams_xplane_plane_id = gams_xplane_plane_id[1:*]
gams_xplane_readstreap_id = gams_xplane_readstreap_id[1:*]
gams_xplane_c4 = gams_xplane_c4[1:*]
gams_xplane_c5 = gams_xplane_c5[1:*]
gams_xplane_edep = gams_xplane_edep[1:*]
gams_xplane_c7 = gams_xplane_c7[1:*]

openw,lun,filepath+'G4_GAMS_XPLANE_AGILE'+agile_version+'_'+py_name+'_'+sim_name+'_'+stripname+'_'+sname+'_'+strmid(strtrim(string(N_in),1),0,10)+part_type+'_'+ene_type+'MeV_'+strmid(strtrim(string(theta_type),1),0,10)+'_'+strmid(strtrim(string(phi_type),1),0,10)+'.all.dat',/get_lun
; ASCII Columns:
; - c1 = event ID
; - c2 = plane ID
; - c3 = readout strip ID
; - c4 = -999
; - c5 = -999
; - c6 = energy dep in MIP
; - c7 = -999    event_start = -1
j=0l
while (1) do begin
    where_event_eq = where(gams_xplane_event_id EQ gams_xplane_event_id(j))
    gams_xplane_plane_id_temp = gams_xplane_plane_id(where_event_eq)
    gams_xplane_readstreap_id_temp = gams_xplane_readstreap_id(where_event_eq)
    gams_xplane_c4_temp  = gams_xplane_c4(where_event_eq)
    gams_xplane_c5_temp  = gams_xplane_c5(where_event_eq)
    gams_xplane_edep_temp = gams_xplane_edep(where_event_eq)    
    gams_xplane_c7_temp = gams_xplane_c7(where_event_eq)    
    
    for r=0l, n_elements(gams_xplane_plane_id_temp)-1 do begin
        printf, lun, gams_xplane_event_id(j), gams_xplane_plane_id_temp(r), gams_xplane_readstreap_id_temp(r), gams_xplane_c4_temp(r), gams_xplane_c5_temp(r), gams_xplane_edep_temp(r), gams_xplane_c7_temp(r), format='(I6,I5,I5,I5,I5,F10.5,I5)'        
    endfor

    N_event_eq = n_elements(where_event_eq)
    if where_event_eq(N_event_eq-1) LT (n_elements(gams_xplane_event_id)-1) then begin
      j = where_event_eq(N_event_eq-1)+1
    endif else break
endwhile

Free_lun, lun

spawn,'cp '+filepath+'G4_GAMS_XPLANE_AGILE'+agile_version+'_'+py_name+'_'+sim_name+'_'+stripname+'_'+sname+'_'+strmid(strtrim(string(N_in),1),0,10)+part_type+'_'+ene_type+'MeV_'+strmid(strtrim(string(theta_type),1),0,10)+'_'+strmid(strtrim(string(phi_type),1),0,10)+'.all.dat '+filepath+'stripx.dat'

gams_yplane_event_id = gams_yplane_event_id[1:*]
gams_yplane_plane_id = gams_yplane_plane_id[1:*]
gams_yplane_readstreap_id = gams_yplane_readstreap_id[1:*]
gams_yplane_c4 = gams_yplane_c4[1:*]
gams_yplane_c5 = gams_yplane_c5[1:*]
gams_yplane_edep = gams_yplane_edep[1:*]
gams_yplane_c7 = gams_yplane_c7[1:*]

openw,lun,filepath+'G4_GAMS_YPLANE_AGILE'+agile_version+'_'+py_name+'_'+sim_name+'_'+stripname+'_'+sname+'_'+strmid(strtrim(string(N_in),1),0,10)+part_type+'_'+ene_type+'MeV_'+strmid(strtrim(string(theta_type),1),0,10)+'_'+strmid(strtrim(string(phi_type),1),0,10)+'.all.dat',/get_lun
; ASCII Columns:
; - c1 = event ID
; - c2 = plane ID
; - c3 = readout strip ID
; - c4 = -999
; - c5 = -999
; - c6 = energy dep in MIP
; - c7 = -999    event_start = -1
j=0l
while (1) do begin
    where_event_eq = where(gams_yplane_event_id EQ gams_yplane_event_id(j))
    gams_yplane_plane_id_temp = gams_yplane_plane_id(where_event_eq)
    gams_yplane_readstreap_id_temp = gams_yplane_readstreap_id(where_event_eq)
    gams_yplane_c4_temp  = gams_yplane_c4(where_event_eq)
    gams_yplane_c5_temp  = gams_yplane_c5(where_event_eq)
    gams_yplane_edep_temp = gams_yplane_edep(where_event_eq)    
    gams_yplane_c7_temp = gams_yplane_c7(where_event_eq)    
    
    for r=0l, n_elements(gams_yplane_plane_id_temp)-1 do begin
        printf, lun, gams_yplane_event_id(j), gams_yplane_plane_id_temp(r), gams_yplane_readstreap_id_temp(r), gams_yplane_c4_temp(r), gams_yplane_c5_temp(r), gams_yplane_edep_temp(r), gams_yplane_c7_temp(r), format='(I6,I5,I5,I5,I5,F10.5,I5)'        
    endfor

    N_event_eq = n_elements(where_event_eq)
    if where_event_eq(N_event_eq-1) LT (n_elements(gams_yplane_event_id)-1) then begin
      j = where_event_eq(N_event_eq-1)+1
    endif else break
endwhile

Free_lun, lun
spawn,'cp '+filepath+'G4_GAMS_YPLANE_AGILE'+agile_version+'_'+py_name+'_'+sim_name+'_'+stripname+'_'+sname+'_'+strmid(strtrim(string(N_in),1),0,10)+part_type+'_'+ene_type+'MeV_'+strmid(strtrim(string(theta_type),1),0,10)+'_'+strmid(strtrim(string(phi_type),1),0,10)+'.all.dat '+filepath+'stripy.dat'

if (cal_flag) then begin
	gams_cal_event_id = gams_cal_event_id[1:*]
	gams_cal_bar_plane = gams_cal_bar_plane[1:*]
	gams_cal_bar_id = gams_cal_bar_id[1:*]
	gams_cal_c4 = gams_cal_c4[1:*]
	gams_cal_ene_a = gams_cal_ene_a[1:*]
	gams_cal_ene_b = gams_cal_ene_b[1:*]

	
	openw,lun,filepath+'G4_GAMS_CAL_AGILE'+agile_version+'_'+py_name+'_'+sim_name+'_'+stripname+'_'+sname+'_'+strmid(strtrim(string(N_in),1),0,10)+part_type+'_'+ene_type+'MeV_'+strmid(strtrim(string(theta_type),1),0,10)+'_'+strmid(strtrim(string(phi_type),1),0,10)+'.all.dat',/get_lun
	; ASCII Columns:
	; - c1 = event ID
	; - c2 = bar plane 
	; - c3 = bar_id
	; - c4 = 0
	; - c5 = energy A
	; - c6 = energy B

	j=0l
	while (1) do begin
		printf, lun, gams_cal_event_id(j), gams_cal_bar_plane(j), gams_cal_bar_id(j), gams_cal_c4(j), gams_cal_ene_a(j),gams_cal_ene_b(j)

		if (j LT (n_elements(gams_cal_event_id)-1)) then begin
		  j = j+1
		endif else break
	endwhile

	Free_lun, lun
	spawn,'cp '+filepath+'G4_GAMS_CAL_AGILE'+agile_version+'_'+py_name+'_'+sim_name+'_'+stripname+'_'+sname+'_'+strmid(strtrim(string(N_in),1),0,10)+part_type+'_'+ene_type+'MeV_'+strmid(strtrim(string(theta_type),1),0,10)+'_'+strmid(strtrim(string(phi_type),1),0,10)+'.all.dat '+filepath+'calo.dat'
endif
if (ac_flag) then begin
	gams_ac_event_id = gams_ac_event_id[1:*]
	gams_ac_panel = gams_ac_panel[1:*]
	gams_ac_subpanel = gams_ac_subpanel[1:*]
	gams_ac_edep = gams_ac_edep[1:*]

	openw,lun,filepath+'G4_GAMS_AC_AGILE'+agile_version+'_'+py_name+'_'+sim_name+'_'+stripname+'_'+sname+'_'+strmid(strtrim(string(N_in),1),0,10)+part_type+'_'+ene_type+'MeV_'+strmid(strtrim(string(theta_type),1),0,10)+'_'+strmid(strtrim(string(phi_type),1),0,10)+'.all.dat',/get_lun
	; ASCII Columns:
	; - c1 = event ID
	; - c2 = AC panel
	; - c3 = AC subpanel
	; - c4 = energy deposit

	j=0l
	while (1) do begin
		printf, lun, gams_ac_event_id(j),gams_ac_panel(j), gams_ac_subpanel(j), gams_ac_edep(j), format='(I,2x,A,2x,I,2x,F20.15)'

		if (j LT (n_elements(gams_ac_event_id)-1)) then begin
		  j = j+1
		endif else break
	endwhile

	Free_lun, lun
	spawn,'cp '+filepath+'G4_GAMS_AC_AGILE'+agile_version+'_'+py_name+'_'+sim_name+'_'+stripname+'_'+sname+'_'+strmid(strtrim(string(N_in),1),0,10)+part_type+'_'+ene_type+'MeV_'+strmid(strtrim(string(theta_type),1),0,10)+'_'+strmid(strtrim(string(phi_type),1),0,10)+'.all.dat '+filepath+'ac.dat'
endif

if (out_type EQ 0) then begin

	; -----> FITS files
	; G4.RAW.AGILE<version>.<phys>List.<strip>.<point>.<n_in>ph.<energy>MeV.<theta>.<phi>.all.fits
	rawData_event_id = rawData_event_id[1:*]
	rawData_Si_id = rawData_Si_id[1:*]
	rawData_tray_id =  rawData_tray_id[1:*]
	rawData_plane_id =  rawData_plane_id[1:*]
	rawData_Strip_id =  rawData_Strip_id[1:*]
	rawData_energy_dep =  rawData_energy_dep[1:*]
	rawData_ent_x =  rawData_ent_x[1:*]
	rawData_ent_y =  rawData_ent_y[1:*]
	rawData_ent_z =  rawData_ent_z[1:*]
	rawData_exit_x =  rawData_exit_x[1:*]
	rawData_exit_y =  rawData_exit_y[1:*]
	rawData_exit_z = rawData_exit_z[1:*]

	; G4.FIRST.AGILE<version>.<phys>List.<strip>.<point>.<n_in>ph.<energy>MeV.<theta>.<phi>.all.fits
	firstData_event_id = firstData_event_id[1:*]
	firstData_Si_id = firstData_Si_id[1:*]
	firstData_tray_id =  firstData_tray_id[1:*]
	firstData_plane_id =  firstData_plane_id[1:*]
	firstData_Strip_id =  firstData_Strip_id[1:*]
	firstData_energy_dep =  firstData_energy_dep[1:*]
	firstData_ent_x =  firstData_ent_x[1:*]
	firstData_ent_y =  firstData_ent_y[1:*]
	firstData_ent_z =  firstData_ent_z[1:*]
	firstData_exit_x =  firstData_exit_x[1:*]
	firstData_exit_y =  firstData_exit_y[1:*]
	firstData_exit_z = firstData_exit_z[1:*]
	
	; L0.AGILE<version>.<phys>List.<strip>.<point>.<n_in>ph.<energy>MeV.<theta>.<phi>.all.fits
	L0TRACKERGLOBAL_Glob_event_id_test = L0TRACKERGLOBAL_Glob_event_id_test[1:*]
	L0TRACKERGLOBAL_Glob_vol_id_test = L0TRACKERGLOBAL_Glob_vol_id_test[1:*]
	L0TRACKERGLOBAL_Glob_moth_id_test = L0TRACKERGLOBAL_Glob_moth_id_test[1:*]
	L0TRACKERGLOBAL_Glob_tray_id_test = L0TRACKERGLOBAL_Glob_tray_id_test[1:*]
	L0TRACKERGLOBAL_Glob_plane_id_test = L0TRACKERGLOBAL_Glob_plane_id_test[1:*]
	L0TRACKERGLOBAL_Glob_Si_id_test = L0TRACKERGLOBAL_Glob_Si_id_test[1:*]
	L0TRACKERGLOBAL_Glob_Strip_id_test = L0TRACKERGLOBAL_Glob_Strip_id_test[1:*]
	L0TRACKERGLOBAL_Glob_Strip_type_test = L0TRACKERGLOBAL_Glob_Strip_type_test[1:*]
	L0TRACKERGLOBAL_Glob_pos_test = L0TRACKERGLOBAL_Glob_pos_test[1:*]
	L0TRACKERGLOBAL_Glob_zpos_test = L0TRACKERGLOBAL_Glob_zpos_test[1:*]
	L0TRACKERGLOBAL_Glob_energy_dep_test = L0TRACKERGLOBAL_Glob_energy_dep_test[1:*]

	; L0.5.DIGI.AGILE<version>.<phys>List.<strip>.<point>.<n_in>ph.<energy>MeV.<theta>.<phi>.all.fits
	L05TRACKERGLOBAL_Glob_event_id_acap = L05TRACKERGLOBAL_Glob_event_id_acap[1:*]
	L05TRACKERGLOBAL_Glob_vol_id_acap = L05TRACKERGLOBAL_Glob_vol_id_acap[1:*]
	L05TRACKERGLOBAL_Glob_moth_id_acap = L05TRACKERGLOBAL_Glob_moth_id_acap[1:*]
	L05TRACKERGLOBAL_Glob_tray_id_acap = L05TRACKERGLOBAL_Glob_tray_id_acap[1:*]
	L05TRACKERGLOBAL_Glob_plane_id_acap = L05TRACKERGLOBAL_Glob_plane_id_acap[1:*]
	L05TRACKERGLOBAL_Glob_Si_id_acap = L05TRACKERGLOBAL_Glob_Si_id_acap[1:*]
	L05TRACKERGLOBAL_Glob_Strip_id_acap = L05TRACKERGLOBAL_Glob_Strip_id_acap[1:*]
	L05TRACKERGLOBAL_Glob_Strip_type_acap = L05TRACKERGLOBAL_Glob_Strip_type_acap[1:*]
	L05TRACKERGLOBAL_Glob_pos_acap = L05TRACKERGLOBAL_Glob_pos_acap[1:*]
	L05TRACKERGLOBAL_Glob_zpos_acap = L05TRACKERGLOBAL_Glob_zpos_acap[1:*]
	L05TRACKERGLOBAL_Glob_energy_dep_acap = L05TRACKERGLOBAL_Glob_energy_dep_acap[1:*]

	if (cal_flag) then begin
		; G4.CAL.AGILE<version>.<phys>List.<strip>.<point>.<n_in>ph.<energy>MeV.<theta>.<phi>.all.fits 
		calInput_event_id_tot_cal = calInput_event_id_tot_cal[1:*]
		calInput_bar_plane_tot = calInput_bar_plane_tot[1:*]
		calInput_bar_id_tot = calInput_bar_id_tot[1:*]
		calInput_ene_a_tot = calInput_ene_a_tot[1:*]
		calInput_ene_b_tot = calInput_ene_b_tot[1:*]
	endif
	if (ac_flag) then begin
		; G4.AC.AGILE<version>.<phys>List.<strip>.<point>.<n_in>ph.<energy>MeV.<theta>.<phi>.all.fits
		acInput_event_id_tot_ac = acInput_event_id_tot_ac[1:*]
		acInput_AC_panel = acInput_AC_panel[1:*]
		acInput_AC_subpanel = acInput_AC_subpanel[1:*]
		acInput_energy_dep_tot_ac = acInput_energy_dep_tot_ac[1:*]
	endif

	CREATE_STRUCT, rawData, 'rawData', ['EVT_ID', 'TRK_FLAG', 'TRAY_ID', 'PLANE_ID', 'STRIP_ID', 'E_DEP', 'X_ENT', 'Y_ENT', 'Z_ENT', 'X_EXIT', 'Y_EXIT', 'Z_EXIT'], $
	'I,I,I,I,I,F20.5,F20.5,F20.5,F20.5,F20.5,F20.5,F20.5', DIMEN = n_elements(rawData_event_id)
	rawData.EVT_ID = rawData_event_id
	rawData.TRK_FLAG = rawData_Si_id
	rawData.TRAY_ID = rawData_tray_id
	rawData.PLANE_ID = rawData_plane_id
	rawData.STRIP_ID = rawData_Strip_id
	rawData.E_DEP = rawData_energy_dep
	rawData.X_ENT = rawData_ent_x
	rawData.Y_ENT = rawData_ent_y
	rawData.Z_ENT = rawData_ent_z
	rawData.X_EXIT = rawData_exit_x
	rawData.Y_EXIT = rawData_exit_y
	rawData.Z_EXIT = rawData_exit_z


	hdr_rawData = ['COMMENT  AGILE '+agile_version+' Geant4 simulation', $
				   'N_in     = '+strtrim(string(N_in),1), $
				   'Energy     = '+ene_type, $
				   'Theta     = '+strtrim(string(theta_type),1), $
				   'Phi     = '+strtrim(string(phi_type),1), $
				   'Position unit = cm', $
				   'Energy unit = keV']

	MWRFITS, rawData, filepath+'G4.RAW.AGILE'+agile_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+strmid(strtrim(string(N_in),1),0,10)+part_type+'.'+ene_type+'MeV.'+strmid(strtrim(string(theta_type),1),0,10)+'.'+strmid(strtrim(string(phi_type),1),0,10)+'.all.fits', hdr_rawData, /create

	CREATE_STRUCT, firstData, 'firstData', ['EVT_ID', 'TRK_FLAG', 'TRAY_ID', 'PLANE_ID', 'STRIP_ID', 'E_DEP', 'X_ENT', 'Y_ENT', 'Z_ENT', 'X_EXIT', 'Y_EXIT', 'Z_EXIT'], $
	  'I,I,I,I,I,F20.5,F20.5,F20.5,F20.5,F20.5,F20.5,F20.5', DIMEN = n_elements(firstData_event_id)
	firstData.EVT_ID = firstData_event_id
	firstData.TRK_FLAG = firstData_Si_id
	firstData.TRAY_ID = firstData_tray_id
	firstData.PLANE_ID = firstData_plane_id
	firstData.STRIP_ID = firstData_Strip_id
	firstData.E_DEP = firstData_energy_dep
	firstData.X_ENT = firstData_ent_x
	firstData.Y_ENT = firstData_ent_y
	firstData.Z_ENT = firstData_ent_z
	firstData.X_EXIT = firstData_exit_x
	firstData.Y_EXIT = firstData_exit_y
	firstData.Z_EXIT = firstData_exit_z


	hdr_firstData = ['COMMENT  AGILE '+agile_version+' Geant4 simulation', $
	  'N_in     = '+strtrim(string(N_in),1), $
	  'Energy     = '+ene_type, $
	  'Theta     = '+strtrim(string(theta_type),1), $
	  'Phi     = '+strtrim(string(phi_type),1), $
	  'Position unit = cm', $
	  'Energy unit = keV']

	MWRFITS, firstData, filepath+'G4.FIRST.AGILE'+agile_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+strmid(strtrim(string(N_in),1),0,10)+part_type+'.'+ene_type+'MeV.'+strmid(strtrim(string(theta_type),1),0,10)+'.'+strmid(strtrim(string(phi_type),1),0,10)+'.all.fits', hdr_firstData, /create

	CREATE_STRUCT, L0TRACKERGLOBAL, 'GLOBALTRACKERL0', ['EVT_ID', 'VOLUME_ID', 'MOTHER_ID', 'TRAY_ID', 'PLANE_ID','TRK_FLAG', 'STRIP_ID', 'STRIP_TYPE', 'POS', 'ZPOS','E_DEP'], 'I,J,J,I,I,I,J,J,F20.5,F20.5,F20.5', DIMEN = N_ELEMENTS(L0TRACKERGLOBAL_Glob_event_id_test)
	L0TRACKERGLOBAL.EVT_ID = L0TRACKERGLOBAL_Glob_event_id_test
	L0TRACKERGLOBAL.VOLUME_ID = L0TRACKERGLOBAL_Glob_vol_id_test
	L0TRACKERGLOBAL.MOTHER_ID = L0TRACKERGLOBAL_Glob_moth_id_test
	L0TRACKERGLOBAL.TRAY_ID = L0TRACKERGLOBAL_Glob_tray_id_test
	L0TRACKERGLOBAL.PLANE_ID = L0TRACKERGLOBAL_Glob_plane_id_test
	L0TRACKERGLOBAL.TRK_FLAG = L0TRACKERGLOBAL_Glob_Si_id_test
	L0TRACKERGLOBAL.STRIP_ID = L0TRACKERGLOBAL_Glob_Strip_id_test
	L0TRACKERGLOBAL.STRIP_TYPE = L0TRACKERGLOBAL_Glob_Strip_type_test
	L0TRACKERGLOBAL.POS = L0TRACKERGLOBAL_Glob_pos_test
	L0TRACKERGLOBAL.ZPOS = L0TRACKERGLOBAL_Glob_zpos_test
	L0TRACKERGLOBAL.E_DEP = L0TRACKERGLOBAL_Glob_energy_dep_test

	HDR_L0GLOBAL = ['Creator          = Valentina Fioretti', $
			  'BoGEMMS release  = AGILE '+agile_version, $
			  'N_IN             = '+STRTRIM(STRING(N_IN),1)+'   /Number of simulated particles', $
			  'ENERGY           = '+ene_type+'   /Simulated input energy', $
			  'THETA            = '+STRTRIM(STRING(THETA_TYPE),1)+'   /Simulated input theta angle', $
			  'PHI              = '+STRTRIM(STRING(PHI_TYPE),1)+'   /Simulated input phi angle', $
			  'ENERGY UNIT      = KEV']


	MWRFITS, L0TRACKERGLOBAL, filepath+'L0.AGILE'+agile_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+STRMID(STRTRIM(STRING(N_IN),1),0,10)+part_type+'.'+ene_type+'MeV.'+STRMID(STRTRIM(STRING(THETA_TYPE),1),0,10)+'.'+STRMID(STRTRIM(STRING(PHI_TYPE),1),0,10)+'.all.fits', HDR_L0GLOBAL, /CREATE


	CREATE_STRUCT, L05TRACKERGLOBAL, 'GLOBALTRACKERL05', ['EVT_ID', 'VOLUME_ID', 'MOTHER_ID', 'TRAY_ID','PLANE_ID','TRK_FLAG', 'STRIP_ID', 'STRIP_TYPE', 'POS', 'ZPOS','E_DEP'], 'J,J,J,I,I,I,J,J,F20.5,F20.5,F20.5', DIMEN = N_ELEMENTS(L05TRACKERGLOBAL_Glob_event_id_acap)
	L05TRACKERGLOBAL.EVT_ID = L05TRACKERGLOBAL_Glob_event_id_acap
	L05TRACKERGLOBAL.VOLUME_ID = L05TRACKERGLOBAL_Glob_vol_id_acap
	L05TRACKERGLOBAL.MOTHER_ID = L05TRACKERGLOBAL_Glob_moth_id_acap
	L05TRACKERGLOBAL.TRAY_ID = L05TRACKERGLOBAL_Glob_tray_id_acap
	L05TRACKERGLOBAL.PLANE_ID = L05TRACKERGLOBAL_Glob_plane_id_acap
	L05TRACKERGLOBAL.TRK_FLAG = L05TRACKERGLOBAL_Glob_Si_id_acap
	L05TRACKERGLOBAL.STRIP_ID = L05TRACKERGLOBAL_Glob_Strip_id_acap
	L05TRACKERGLOBAL.STRIP_TYPE = L05TRACKERGLOBAL_Glob_Strip_type_acap
	L05TRACKERGLOBAL.POS = L05TRACKERGLOBAL_Glob_pos_acap
	L05TRACKERGLOBAL.ZPOS = L05TRACKERGLOBAL_Glob_zpos_acap
	L05TRACKERGLOBAL.E_DEP = L05TRACKERGLOBAL_Glob_energy_dep_acap

	HDR_L05GLOBAL = ['Creator          = Valentina Fioretti', $
			  'BoGEMMS release  = AGILE '+agile_version, $
			  'N_IN             = '+STRTRIM(STRING(N_IN),1)+'   /Number of simulated particles', $
			  'ENERGY           = '+ene_type+'   /Simulated input energy', $
			  'THETA            = '+STRTRIM(STRING(THETA_TYPE),1)+'   /Simulated input theta angle', $
			  'PHI              = '+STRTRIM(STRING(PHI_TYPE),1)+'   /Simulated input phi angle', $
			  'ENERGY UNIT      = KEV']


	MWRFITS, L05TRACKERGLOBAL, filepath+'L0.5.DIGI.AGILE'+agile_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+STRMID(STRTRIM(STRING(N_IN),1),0,10)+part_type+'.'+ene_type+'MeV.'+STRMID(STRTRIM(STRING(THETA_TYPE),1),0,10)+'.'+STRMID(STRTRIM(STRING(PHI_TYPE),1),0,10)+'.all.fits', HDR_L05GLOBAL, /CREATE

	if (cal_flag) then begin
		CREATE_STRUCT, calInput, 'input_cal_dhsim', ['EVT_ID', 'BAR_PLANE', 'BAR_ID', 'ENERGY_A', 'ENERGY_B'], $
		'I,I,I,F20.15,F20.15', DIMEN = n_elements(calInput_event_id_tot_cal)
		calInput.EVT_ID = calInput_event_id_tot_cal
		calInput.BAR_PLANE = calInput_bar_plane_tot 
		calInput.BAR_ID = calInput_bar_id_tot
		calInput.ENERGY_A = calInput_ene_a_tot
		calInput.ENERGY_B = calInput_ene_b_tot


		hdr_calInput = ['COMMENT  AGILE V2.0 Geant4 simulation', $
					   'N_in     = '+strtrim(string(N_in),1), $
					   'Energy     = '+ene_type, $
					   'Theta     = '+strtrim(string(theta_type),1), $
					   'Phi     = '+strtrim(string(phi_type),1), $
					   'Energy unit = GeV']

		MWRFITS, calInput, filepath+'G4.CAL.AGILE'+agile_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+strmid(strtrim(string(N_in),1),0,10)+part_type+'.'+ene_type+'MeV.'+strmid(strtrim(string(theta_type),1),0,10)+'.'+strmid(strtrim(string(phi_type),1),0,10)+'.all.fits', hdr_calInput, /create
	endif
	if (ac_flag) then begin

		CREATE_STRUCT, acInput, 'input_ac_dhsim', ['EVT_ID', 'AC_PANEL', 'AC_SUBPANEL', 'E_DEP'], $
		'I,A,I,F20.15', DIMEN = n_elements(acInput_event_id_tot_ac)
		acInput.EVT_ID = acInput_event_id_tot_ac
		acInput.AC_PANEL = acInput_AC_panel
		acInput.AC_SUBPANEL = acInput_AC_subpanel
		acInput.E_DEP = acInput_energy_dep_tot_ac


		hdr_acInput = ['COMMENT  AGILE V2.0 Geant4 simulation', $
					   'N_in     = '+strtrim(string(N_in),1), $
					   'Energy     = '+ene_type, $
					   'Theta     = '+strtrim(string(theta_type),1), $
					   'Phi     = '+strtrim(string(phi_type),1), $
					   'Energy unit = GeV']

		MWRFITS, acInput, filepath+'G4.AC.AGILE'+agile_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+strmid(strtrim(string(N_in),1),0,10)+part_type+'.'+ene_type+'MeV.'+strmid(strtrim(string(theta_type),1),0,10)+'.'+strmid(strtrim(string(phi_type),1),0,10)+'.all.fits', hdr_acInput, /create
	endif
endif


end
