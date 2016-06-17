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

rout=2.18e17
dustyfile='p0tau1.stb'
mocassinfile='disk0.1.SED'
Lstar = 1. 
r_r1=1000.

data=fltarr(3,61)
openr,lun,mocassinfile,/get_lun
nothing=' '
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf,lun,data
close,lun 
free_lun,lun

bench=fltarr(8,105)
openr,lun,dustyfile,/get_lun
nothing=' '
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf, lun, bench
close,lun 
free_lun,lun

set_plot, 'ps'
device, /color

data(2,*) = data(2,*)*1.e36*8./Lstar

plot,data(0,*), data(2,*) ,psym=-1,xtitle='nu [Ryd]', ytitle='relative flux',title='tau(10um)=1',xrange=[1.e-4,0.4],/xsty,/xlog

plot,data(0,*), data(2,*),xrange=[1.e-4,0.4],yrange=[1.e-6,1],psym=-1,xtitle='nu [Ryd]', ytitle='relative flux',title='p=0 tau=1',/xsty,/xlog,/ylog

plot,data(1,*), data(2,*) ,psym=-1,xtitle='lambda [um]', ytitle='relative flux',title='p=0, tau=1',yrange=[0.,1.],/xsty
oplot,bench(0,*), bench(1,*),col=!col.red,psym=-3

plot,data(1,*), data(2,*),psym=-1,xtitle='lambda [um]', ytitle='relative flux',xrange=[0.1,10000.],yrange=[1.e28,1.e34],title='tau(1um)=1',/xsty,/xlog,/ylog
oplot,bench(0,*), bench(1,*),col=!col.red,psym=-3

ones=fltarr(105)
for i = 0,104 do begin
ones(i)=1.
endfor

;plot,data(1,*), data(2,*)/bench(1,*),psym=-1,xrange=[0.12,120.],yrange=[0,5.],xtitle=' lambda [um]', ytitle='departure',title='bench1 : departure',/xsty
;yrange=[0,2],
;oplot, data(1,*),ones,col=!col.red,psym=0
;oplot, data(1,*),data(2,*)/data_save(2,*),col=!col.green,psym=0


device, /close
end
