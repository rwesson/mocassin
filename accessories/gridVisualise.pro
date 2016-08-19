;
; this program creates 2D projections (written out to 
; postscript files) and 3D images (to be visulaized 
; with IDL's slicer3), of a data cube. 
; needs tvim.pro and prepare_slicer3.pro
;
; Barbara Ercolano 2004
;
nx=16
ny=nx
nz=nx
densityfile='p0tau1.ndust'

set_plot,'ps'
device,/color
loadct,24
data=replicate({x:0.,y:0.,z:0.,hden:0.},nx,ny,nz)
openr, lun,densityfile,/get_lun
readf, lun, data
close,lun
free_lun,lun
maxvalue = max(data.hden)
minvalue = maxvalue/1000.
index = where( data.hden lt minvalue)
data[index].hden = 0.
tvim,data[*,*,0].hden ; projection along x
tvim,data[*,0,*].hden ; projection along y
tvim,data[0,*,*].hden ; projection along z
device,/close
set_plot, 'X'
prepare_slicer3,'density',data.hden,'slicer.dat'
slicer3

end
