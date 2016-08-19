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
dustyfile1='p0tau10.stb'
dustyfile2='p0tau100.stb'
mocassinfile='p0tau1.SED'
mocassinfile1='p0tau10.SED'
mocassinfile2='p0tau100.SED'
Lstar = 38.26e36 
r_r1=1000.

data=fltarr(3,101)
openr,lun,mocassinfile,/get_lun
nothing=' '
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf,lun,data
close,lun 
free_lun,lun

data1=fltarr(3,101)
openr,lun,mocassinfile1,/get_lun
nothing=' '
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf,lun,data1
close,lun 
free_lun,lun

data2=fltarr(3,101)
openr,lun,mocassinfile2,/get_lun
nothing=' '
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf,lun,data2
close,lun 
free_lun,lun

bench=fltarr(8,101)
openr,lun,dustyfile,/get_lun
nothing=' '
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf, lun, bench
close,lun 
free_lun,lun

bench1=fltarr(8,101)
openr,lun,dustyfile1,/get_lun
nothing=' '
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf, lun, bench1
close,lun 
free_lun,lun

bench2=fltarr(8,101)
openr,lun,dustyfile2,/get_lun
nothing=' '
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf, lun, bench2
close,lun 
free_lun,lun

set_plot, 'ps'
device, /color

data(2,*) = data(2,*)*1.e36*8./Lstar
data1(2,*) = data1(2,*)*1.e36*8./Lstar
data2(2,*) = data2(2,*)*1.e36*8./Lstar

plot,data(1,*), data(2,*)*!pi ,xrange=[0.12,1000],yrange=[1.e-4,1],psym=0,xtitle='lambda [um]', ytitle='relative flux',title='tau(1um)=1, 10, 100',/xsty,/xlog,/ylog
oplot,bench(0,*), bench(1,*),col=!col.red,psym=0,linestyle=2

oplot,data1(1,*),data1(2,*)*!pi,psym=0,col=!col.black,thick=4
oplot,bench1(0,*), bench1(1,*),col=!col.red,psym=0,linestyle=2,thick=4

oplot,data2(1,*),data2(2,*)*!pi,col=!col.black,psym=0,thick=8
oplot,bench2(0,*), bench2(1,*),col=!col.red,psym=0,linestyle=2,thick=8

device, /close
end
