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
nx=30
rout=1.49e16
rin=2.18e17/1000.
dustyfile='p0tau1.rtb'
mocassinfile='disk0.1.dustGrid'
axisfile='disk0.1.grid0'
r_r1=1000.

dustprop=replicate({ndust:fltarr(1),t:fltarr(1)},nx,nx,nx)
openr,lun,mocassinfile,/get_lun
readf,lun,dustprop
close,lun 
free_lun,lun
openr,lun,axisfile,/get_lun
nuthing=' '
readf, lun, nothing
x=fltarr(nx)
readf,lun,x      
close,lun 
free_lun,lun
x = x/rout
r=fltarr(nx,nx,nx)

bench=fltarr(7,18)
openr,lun,dustyfile,/get_lun
nothing=' '
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf, lun,bench
close,lun
free_lun,lun

set_plot, 'ps'
device, /color
for i = 0,nx-1 do for j = 0,nx-1 do for k =0,nx-1 do $
   r[i,j,k]=sqrt(x[i]^2+x[j]^2+x[k]^2) 
;plot,r,dustprop.t[0],psym=1,xrange=[rin/rout, 1.],yrange=[0., 800.],xtitle='r/R',$
;   ytitle='T_d.',title='silicate a=0.16um',/xsty
plot,x,dustprop[0,*,0].t[0],psym=1,xrange=[rin/rout, 1.],yrange=[0., 400.],xtitle='r/R',$
   ytitle='T_d.',title='silicate a=0.16um',/xsty
oplot,x,dustprop[*,0,0].t[0],col=!col.red
oplot,x,dustprop[0,0,*].t[0],col=!col.green
;oplot,bench(0,*)/r_r1,bench(5,*),col=!col.red,psym=-1


;oplot,r,dustprop.c[4],col=!col.cyan,psym=1

;xyouts,2.d17,0.9,'He2+',col=!col.violet
device, /close

end
