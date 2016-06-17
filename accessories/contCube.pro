print,' Initializing'

;---------------------------------------------------------------
;  Define some parameters
;---------------------------------------------------------------
;

x_angle = 0.
y_angle = 0.
z_angle = 0.

proj_axy = 1   ;; x=1,y=2,z=3 ! attention, proj_abs est TJS sur Z ...

dist = 1.5 ;kpc Diab
cube_size =1.24;

;color_table = 41
;loadct,color_table,file='barbaras_colors.tbl' 

color_table = 1
loadct,color_table


size_im = 50l;; size of the initial images to be read

temp_lim = 500.

vel_0 = 23. ; km/s  parameter for the velocity field Bassgen
vel_1 = 50. ;km/s
n_vel =  41  ; nb pixels for the velocity spectra
velo_max =  80. ;; velocity maximum for the velocity spectra

velo_coeff =  [1.,1.,1.] 

v_turb =  2. ;km/s

atomic_mass = [0,1,1,4,4,14,16,16,16,16,16,32,32]

;liste_num_line_prof =  [2,3,4,5,6,7,8,9,10]
liste_num_line_prof =  [5]

;ap_diam =  3./1.3 ;; arcsec
;ap_centers = [[0,14],[7,14],[14,14],[0,7],[7,7],[14,7],[0.,0.],[7,0],[14,0]]/2.4/1.3
ap_diam =  3./1.9 ;; arcsec
ap_centers = [[0,14],[7,14],[14,14],[0,7],[7,7],[14,7],[0.,0.],[7,0],[14,0]]/2.4/1.9
rot_aper = [0.,0.,0.]
multi_aper =  [0,3,3]
;
;---------------------------------------------------------------
;  Define what to do
;---------------------------------------------------------------
;
do_ps = 1

do_read = 1 & if n_elements(rai1) lt 1 then do_read = 1
do_read_abund = 0
use_old =  0
raia_file = 1
do_mask = 0
do_x8 =  1
do_divise =  0
print,'works'
do_slicer =  1
do_axy_median = 0
do_turn = 0
do_velocity = 0 & if n_elements(velocity) lt 1 then do_velocity = 1
velo_symetric = 0
do_profil =  0
do_plot_profils = 0
do_im_unred = 1
do_temp_r = 0
do_image = 1
do_plot_aper = 0
do_diag = 0

double =  1
inter = 0

if do_ps then set_plot,'ps'
if do_ps then device,/color
loadct,color_table

correc_intens = fltarr(12) + 1.
;correc_intens[1] = 0.742        ; OIII 5007
;correc_intens[7] = 0.746        ; NII 6583
;correc_intens[8] = 0.759        ; OI 6300
;
;-----------------------------------------------------------
;-----------------------------------------------------------
;	Readind the files
;-----------------------------------------------------------
;-----------------------------------------------------------
;
print,' Let s go...'

if do_read then begin
   print,' reading files'
   read_names,liste_name,liste_abond,liste_par
   liste_name =  liste_name[12*(raia_file-1):12*raia_file-1]
  
   rai1 =  0
   rai1 = read_raias('contCube.out',$
                     size_im,size_im,size_im, viewangles,$
                     do_x8=do_x8,/contcube,double=double,do_divise = do_divise)

   ; take off the star

endif ;;do_read

;
;
;
if do_x8 then size_im =  2 * size_im -1
;
;
print, 'here'
if do_x8 then mult_by_2=size_im/2 else mult_by_2=0
if do_read then for i = 1, 1 do $
  rai1[mult_by_2,*,*].im0 = rai1[mult_by_2,*,*].im0 *2
if do_read then for i = 1, 1 do $
  rai1[*,mult_by_2,*].im0 = rai1[*,mult_by_2,*].im0 *2
if do_read then for i = 1, 1 do $
  rai1[*,*,mult_by_2].im0 = rai1[*,*,mult_by_2].im0 *2

; take the star off
rai1[mult_by_2,mult_by_2,mult_by_2].im0=0.


;
;-----------------------------------------------------------
;	Defining size of arrays and images
;-----------------------------------------------------------
;
print,' defining some variables'
size_tab_r = n_tags(rai1) - 1 

if double then im = dblarr(size_im,size_im,size_tab_r,/nozero) $
 else im = fltarr(size_im,size_im,size_tab_r,/nozero)

nb_im = size_tab_r


size_prof =  size_im
;
;-----------------------------------------------------------
;	Rotating the datas
;-----------------------------------------------------------
;
rai1_a=rai1
if do_turn then begin
   rai1_a = turn_3d(temporary(rai1_a),x_angle,y_angle,z_angle, $
                    /res,/cons,/verb,/int)
endif ;;do_turn
;
;-----------------------------------------------------------
;	Image unredened
;-----------------------------------------------------------
;
if do_im_unred then begin
   print,' Performing the images'
   im_unred = im
   
   for i_im = 0,nb_im-1 do im_unred[*,*,i_im] = total(rai1_a.(i_im+1),proj_axy)

   print,'im_unred done'
endif ;;do_im_unred
;
;-----------------------------------------------------------
; Intensite / Hbeta
;-----------------------------------------------------------
;
intens_hb =  fltarr(12)
ibeta =  total(rai1_a.(1))
for i_im = 0,nb_im-1 do intens_hb[i_im] = $
 total(rai1_a.(i_im+1))/ibeta*correc_intens[i_im]
for i_im = 0,nb_im-1 do print,liste_name[i_im],intens_hb[i_im]
;
;-----------------------------------------------------------
; Images 
;-----------------------------------------------------------
;
if do_image then begin

x = ((findgen(size_im)-(size_im-1)/2.)/(size_im-1)*2.)* cube_size 
y = x
im_hbeta = (im_unred[*,*,0]) 

!p.multi =  0;[0,3,4]
for i_im = 0,nb_im-1 do begin
   loadct,color_table

   tit = (liste_name.id)[i_im]+' '+strcompress((liste_name.lam)[i_im],/re)
   tvim,im_unred[*,*,i_im],title=tit,/sc,inter=inter
   if do_plot_aper and n_elements(n_aper) ne 0 then $
    for i_aper= 0,n_aper-1 do contour,/over,apertu[*,*,i_aper], $
    levels=[0.5],thick=i_aper+1
   
endfor ;;i_im
endif
;
;-----------------------------------------------------------
; Misc
;-----------------------------------------------------------
;

flux_tot_beta =  total(rai1_a.im0[0]) ;/4./!pi/(3.086d18*dist*1e3)^2

print,'flux_tot',flux_tot_beta;,'  masse_tot_Hioni',masse_tot_Hioni


print,'terminado'
if do_ps then device,/close
end ;;main





