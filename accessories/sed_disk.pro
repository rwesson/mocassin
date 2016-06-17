; calculates emerging SED from a disk model at a distance equal to the
; stellar radius.
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

rstar=6.9599e10
rout=3.27e17
t0001='SED000.1.out'
t0010='SED000.1.out'
t0100='SED000.1.out'
t1000='SED000.1.out'

r_r1=1000.
vp=2


;benchfile='bench_sp.b-0.1'
;openr,lun,benchfile,/get_lun
;bench=replicate({flux:dblarr(3,61)},180)
;readf,lun,bench
;close,lun
;free_lun,lun

data0001=replicate({ryd:fltarr(1),ang:fltarr(1),flux:fltarr(vp+1)},101)
openr,lun,t0001,/get_lun
nothing=' '
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf,lun,data0001
close,lun 
free_lun,lun

data0010=replicate({ryd:fltarr(1),ang:fltarr(1),flux:fltarr(vp+1)},101)
openr,lun,t0010,/get_lun
nothing=' '
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf,lun,data0010
close,lun 
free_lun,lun

data0100=replicate({ryd:fltarr(1),ang:fltarr(1),flux:fltarr(vp+1)},101)
openr,lun,t0100,/get_lun
nothing=' '
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf,lun,data0100
close,lun 
free_lun,lun

data1000=replicate({ryd:fltarr(1),ang:fltarr(1),flux:fltarr(vp+1)},101)
openr,lun,t1000,/get_lun
nothing=' '
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf,lun,data1000
close,lun 
free_lun,lun

set_plot, 'ps'
device, /color

data0001.flux[*] = data0001.flux[*]*2.*1.e36*1.e-7*1.e4/(rstar^2.)
data0010.flux[*] = data0010.flux[*]*2.*1.e36*1.e-7*1.e4/(rstar^2.)
data0100.flux[*] = data0100.flux[*]*2.*1.e36*1.e-7*1.e4/(rstar^2.)
data1000.flux[*] = data1000.flux[*]*2.*1.e36*1.e-7*1.e4/(rstar^2.)

plot,data0001.ang, data0001.flux[0] ,psym=0,xtitle='!4k!X [!4l!Xm]', ytitle='!4k!X F!D!4k!X!N [W m!U-2!N]',title='average direction', xrange=[0.1,10000.],yrange=[1.e3,1.e8],/xsty,/xlog,/ylog
oplot,data0010.ang, data0010.flux[0], linestyle=1, col=!col.red
oplot,data0100.ang, data0100.flux[0], linestyle=2, col=!col.blue
oplot,data1000.ang, data1000.flux[0], linestyle=3, col=!col.green

;plot,data0001.ang, data0001.flux[1] ,psym=0,xtitle='lambda [um]', ytitle='lambda*Flambda [W m^-2]',title='i=12.5deg',xrange=[0.1,10000.],yrange=[1.e2,1.e8],/xsty,/xlog,/ylog
;oplot,data0010.ang, data0010.flux[1], linestyle=1, col=!col.red
;oplot,data0100.ang, data0100.flux[1], linestyle=2, col=!col.blue
;oplot,data1000.ang, data1000.flux[1], linestyle=3, col=!col.green
plot,data0001.ang, data0001.flux[1]*cos(12.5*!pi/180.) ,psym=0,xtitle='!4k!X [!4l!Xm]', ytitle='!4k!X F!D!4k!X!N [W m!U-2!N]',title='i=12.5!9%!X',xrange=[0.1,10000.],yrange=[1.e3,1.e8],/xsty,/xlog,/ylog
oplot,data0010.ang, data0010.flux[1]*cos(12.5*!pi/180.), linestyle=1, col=!col.red
oplot,data0100.ang, data0100.flux[1]*cos(12.5*!pi/180.), linestyle=2, col=!col.blue
oplot,data1000.ang, data1000.flux[1]*cos(12.5*!pi/180.), linestyle=3, col=!col.green



;plot,data0001.ang, data0001.flux[2] ,psym=0,xtitle='lambda [um]', ytitle='lambda*Flambda [W m^-2]',title='i=77.5deg',xrange=[0.1,10000.],yrange=[1.e2,1.e8],/xsty,/xlog,/ylog
;oplot,data0010.ang, data0010.flux[2], linestyle=1, col=!col.red
;oplot,data0100.ang, data0100.flux[2], linestyle=2, col=!col.blue
;oplot,data1000.ang, data1000.flux[2], linestyle=3, col=!col.green
plot,data0001.ang, data0001.flux[2]*cos(77.5*!pi/180.) ,psym=0,xtitle='!4k!X [!4l!Xm]', ytitle='!4k!X F!D!4k!X!N [W m!U-2!N]',title='i=77.5!9%!X',xrange=[0.1,10000.],yrange=[1.e3,1.e8],/xsty,/xlog,/ylog
oplot,data0010.ang, data0010.flux[2]*cos(77.5*!pi/180.), linestyle=1, col=!col.red
oplot,data0100.ang, data0100.flux[2]*cos(77.5*!pi/180.), linestyle=2, col=!col.blue
oplot,data1000.ang, data1000.flux[2]*cos(77.5*!pi/180.), linestyle=3, col=!col.green


device, /close
end
