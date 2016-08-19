
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

;---------------------------------------------------------------
;  Define what to do
;---------------------------------------------------------------
;
do_ps = 1

do_read = 1 & if n_elements(rai1) lt 1 then do_read = 1
do_read_abund = 1
raia_file = 1
do_x8 =  1
do_divise =  0
print,'works'
do_slicer =  1
do_turn = 1
do_im_unred = 1
do_image = 1
do_plot_aper = 0
do_diag = 0

double =  0
inter = 0

if do_ps then set_plot,'ps'
if do_ps then device,/color
loadct,color_table,file='barbaras_colors.tbl' 

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
   rai1 = read_raias('plot.out',$
                     size_im,size_im,size_im, $
                     do_x8=do_x8,/plot,double=double,do_divise = do_divise)
   
   black= 0
   black = read_raias('grid0.out',$
                      size_im,size_im,size_im, $
                      do_x8=do_x8,/black,double=double,do_divise = do_divise)
   rai1.r1 = rai1.r1*(abs(black.bla-1))

   par =  0
   par = read_raias('grid1.out',$
                    size_im,size_im,size_im, $
                    do_x8=do_x8,/param,double=double,do_divise = do_divise)

   abond = 0
   if do_read_abund then abond = read_raias('grid2.out',$
                    size_im,size_im,size_im, $
                    do_x8=do_x8,/abond,double=double,do_divise = do_divise)


endif ;;do_read

;
;
;
if do_x8 then size_im =  2 * size_im -1
if do_divise then size_im =  size_im / 2
;
;
print, 'here'
if do_x8 then mult_by_2=size_im/2 else mult_by_2=0
if do_read then for i = 1, 1 do $
  rai1[mult_by_2,*,*].(i) = rai1[mult_by_2,*,*].(i) *2
if do_read then for i = 1, 1 do $
  rai1[*,mult_by_2,*].(i) = rai1[*,mult_by_2,*].(i) *2
if do_read then for i = 1, 1 do $
  rai1[*,*,mult_by_2].(i) = rai1[*,*,mult_by_2].(i) *2

print, '1'

par(where(par.(0) eq 0.)).(0)=1.

par.(0) = alog(par.(0))

;-----------------------------------------------------------
; Recording datas in sli.dat to use slicer3
;-----------------------------------------------------------
;
if do_slicer then begin 
   create_3d,(liste_name.id)[0]+strcompress((liste_name.lam)[0]),$
    rai1.(1),'sli.dat'
;   for i =  2,11 do create_3d,(liste_name.id)[i-1]+ $
;    strcompress((liste_name.lam)[i-1]),$
;    rai1.(i),'sli.dat',/append
   for i =  0,n_elements(liste_par)-1 do $
     create_3d,liste_par[i],par.(i),'sli.dat',/append
   for i = 0, n_tags(abond)-1 do  $
     create_3d,(tag_names(abond))[i],abond.(i),'sli.dat',/append
endif ;;do_slicer
print, '2'
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
secparpix =  conv_arc(cube_size,dist)/(size_im/2)

size_prof =  size_im


par_a =  par
rai1_a =  rai1

;
;-----------------------------------------------------------
;	Rotating the datas
;-----------------------------------------------------------
;

if do_turn then begin
   rai1_a = turn_3d(temporary(rai1_a),x_angle,y_angle,z_angle, $
                    /res,/cons,/verb,/int)
   par_a = turn_3d(temporary(par_a),x_angle,y_angle,z_angle, $
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
xarc = conv_arc(x,dist)
yarc = conv_arc(y,dist)

im_hbeta = (im_unred[*,*,0]) 

!p.multi =  0;[0,3,4]
for i_im = 0,nb_im-1 do begin
   loadct,color_table,file='barbaras_colors.tbl' ,/silent

   tit = (liste_name.id)[i_im]+' '+strcompress((liste_name.lam)[i_im],/re)

   tvim,im_unred[*,*,i_im],title=tit,/sc,inter=inter

;,colors=!p.color*indgen(20)/19./15.,range=[0,max(im_unred[*,*,i_im])]
;   confill,im_unred[*,*,i_im]/max(im_unred[*,*,i_im]),title=tit,/aspect
;   tvim,im_unred[*,*,i_im]/(im_hbeta>1e20),title=tit,/sc,inter=int,$
;    range=[0,3*intens_hb[i_im]]

endfor ;;i_im
endif
;
;-----------------------------------------------------------
; Misc
;-----------------------------------------------------------
;

flux_tot_beta =  total(rai1_a.r1) ;/4./!pi/(3.086d18*dist*1e3)^2
if do_divise then flux_tot_beta =  flux_tot_beta * 8.

print,'flux_tot_beta',flux_tot_beta;,'  masse_tot_Hioni',masse_tot_Hioni


print,'tutto fatto'
if do_ps then device,/close
end ;;main





