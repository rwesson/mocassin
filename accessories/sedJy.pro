;
; calculates emerging SED in Jy at earth, plots the results and writes 
; them out to an ascii file
; 
; Barbara Ercolano, December 2004
; 

function sm,tab
coeff=5
n_interm=n_elements(reform(tab[*,0,0]))*coeff
return,smooth(rebin(reform(tab[*,0,0]),n_interm),coeff)
end

defsysv,'!col',EXISTS = col_exists
if not col_exists then $
defsysv,'!col',{black     : 0, $
                white     : 1, $
                red       : 2, $
                green     : 3, $
                blue      : 4, $
                yellow    : 5, $
                magenta   : 6, $
                cyan      : 7, $
                orange    : 8, $
                l_green   : 9, $
                violet    : 10,$
                turquoise : 11,$
                l_blue    : 12,$
                l_red     : 13 $
               }
;        0  1    2   3   4   5   6   7   8   9   10  11  12  13
tvlct,[  0,255, 255,  0,  0,255,255,  0,255,125,125,  0,  0,255],$ ;red
      [  0,255,   0,255,  0,255,  0,255,125,255,  0,255,125,  0],$ ;green
      [  0,255,   0,  0,255,  0,255,255,  0,  0,255,125,255,125]   ;blue


rstar=1500.
pc=3.086e18
ryd2hz = 3.2898e15

mocassinfile1='SED.out'
outputfile='SED.dat'
r_r1=1000.
vp=0
nbins=101

data1=replicate({ryd:fltarr(1),ang:fltarr(1),flux:fltarr(vp+1)},nbins)

openr,lun,mocassinfile1,/get_lun
nothing=' '
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf,lun,data1
close,lun 
free_lun,lun

set_plot, 'ps'
device, /color

const = 1.e36/pc^2.

data1.flux = data1.flux*const*8.*!pi/(4.*!pi*rstar^2.)

plot,data1.ang, data1.flux[0] ,psym=0,xtitle='lambda [um]', ytitle='lambda*Flambda [erg cm^-2]',/xsty,/xlog,/ylog

plot,data1.ang, 1.e23*data1.flux[0]/(data1.ryd*ryd2hz) ,psym=-1,xtitle='lambda [um]', ytitle='Fnu [Jy]',xrange=[5,500],yrange=[0.1,100],/xsty,/xlog,/ylog

plot,data1.ang, 1.e23*data1.flux[0]/(data1.ryd*ryd2hz) ,psym=0,xtitle='lambda [um]', ytitle='Fnu [Jy]',/xsty

openw,lun,outputfile,/get_lun
printf, lun,'    lambda [um]       Fnu [Jy] '
for i =0, nbins-1 do begin
printf, lun, data1[i].ang, 1.e23*data1[i].flux[0]/(data1[i].ryd*ryd2hz)
endfor
close,lun
free_lun,lun

device, /close
end