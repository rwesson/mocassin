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


;********** MODIFY BELOW ***************
mocassinfile='SED.out'   ; file name
n = 1                    ; inclination angles
D=1000.                  ; distance in pc
;********** END MODIFICATION ***********

data=fltarr(3+n,215)

deltanu=fltarr(1,215)
openr,lun,mocassinfile,/get_lun
nothing=' '
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf, lun, nothing
readf,lun,data
close,lun 
free_lun,lun


set_plot, 'ps'
device, file='SEDnew.ps',/color

plot,data(1,*), data(2,*)/D^2 ,psym=0,xrange=[.3,350.],/xs,/ys,/xlog,/ylog, xtitle= "Wavelength [um]", ytitle = 'Flux Density [Jy]', symsize = 1.5, thick=5, yrange = [max(data(2,*)/D^2)/1.e4, 2.*max(data(2,*)/D^2)]
if (n gt 0) then begin
   oplot, data(1,*), data(3,*)/D^2, color = !col.red, line = 2, thick = 5
endif
if (n gt 1) then begin
   oplot, data(1,*), data(4,*)/D^2, color = !col.green, line = 3, thick = 5
endif

device, /close
end
