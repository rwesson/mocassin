;
; this program plots the dust tempratures obtained 
; for a simulation with 20 grain sizes and 2 grain 
; species. 
; this was a multichemistry model and the results 
; are also plotted by sector.
; needs the following input files : dustGrid.out, 
; grainsize file, grainspecies file, grid0.out. 
;
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

set_plot, 'ps'
device, /color
nx = 0
ny = 0
nz = 0
gridfile='grid0.out'
openr,lun,gridfile,/get_lun
ncells=0
motherP=0
rout=0.
readf, lun, motherP
readf, lun, nx,ny,nz,ncells,motherP,rout
convergence=replicate({n:intarr(3)},nx,ny,nz)
x=fltarr(nx)
y=fltarr(ny)
z=fltarr(nz)
readf,lun,x      
readf,lun,y      
readf,lun,z      
readf,lun,convergence
close,lun 
free_lun,lun

sizefile='primary_grainsizes.dat'
nsizes=0
openr,lun,sizefile,/get_lun
readf,lun,nsizes
dustsizes=replicate({index:intarr(1),radius:fltarr(1),weight:fltarr(1)},nsizes)
readf,lun,dustsizes
close,lun 
free_lun,lun

speciesfile='primary_grainspecies.dat'
nspecies=0
openr,lun,speciesfile,/get_lun
readf,lun,nspecies
dustabundances=replicate({abund:fltarr(1)},nspecies)
readf,lun,dustabundances
close,lun
free_lun,lun

mocassinfile='dustGrid.out'
dustprop=replicate({n:fltarr(1),t:fltarr(nspecies+1,nsizes+1)},nx,ny,nz)
openr,lun,mocassinfile,/get_lun
readf,lun,dustprop
close,lun 
free_lun,lun

r=fltarr(nx,ny,nz)
celltot = nx*ny*nz
r1 = fltarr(celltot)
r2 = fltarr(celltot)
r3 = fltarr(celltot)

temp=replicate({t:fltarr(nsizes)},nx,ny,nz)
t1=replicate({t:fltarr(nsizes)},celltot)
t2=replicate({t:fltarr(nsizes)},celltot)
t3=replicate({t:fltarr(nsizes)},celltot)

ii=1
jj=1
kk=1
for i = 0,nx-1 do for j = 0,nx-1 do for k =0,nx-1 do begin
   r[i,j,k]=sqrt(x[i]^2+x[j]^2+x[k]^2) 
   if (sqrt(r[i,j,k]^2-x[i]^2)/r[i,j,k] gt 0.707107 and convergence[i,j,k].n[1] eq 1) then r1[ii] = r[i,j,k]
   if (sqrt(r[i,j,k]^2-x[i]^2)/r[i,j,k] lt 0.707107 and convergence[i,j,k].n[1] eq 1) then r2[jj] = r[i,j,k]
   if (convergence[i,j,k].n[1] eq 1) then r3[kk] = r[i,j,k]

   for l = 0,nsizes-1 do begin
       for m = 0,nspecies-1 do begin
           temp[i,j,k].t[l]=temp[i,j,k].t[l]+dustprop[i,j,k].t[m+1,l+1]*dustabundances[m].abund
           if (sqrt(r[i,j,k]^2-x[i]^2)/r[i,j,k] gt 0.707107 and convergence[i,j,k].n[1] eq 1) then t1[ii] = temp[i,j,k]
           if (sqrt(r[i,j,k]^2-x[i]^2)/r[i,j,k] lt 0.707107and convergence[i,j,k].n[1] eq 1) then t2[jj] = temp[i,j,k]
           if (convergence[i,j,k].n[1] eq 1) then t3[kk] = temp[i,j,k]
       endfor
   endfor

   if (sqrt(r[i,j,k]^2-x[i]^2)/r[i,j,k] gt 0.707107 and convergence[i,j,k].n[1] eq 1) then ii=ii+1
   if (sqrt(r[i,j,k]^2-x[i]^2)/r[i,j,k] lt 0.707107 and convergence[i,j,k].n[1] eq 1) then jj=jj+1
   if (convergence[i,j,k].n[1] eq 1) then kk=kk+1

endfor

plot,r,temp.t[0],psym=1,xrange=[0., rout],yrange=[0., 250.],xtitle='r [cm]',$
   ytitle='T_d.',title='',/xsty
oplot,r,temp.t[8],psym=1, col=!col.red
oplot,r,temp.t[19],psym=1, col=!col.green

plot,r3,t3.t[0],psym=1,xrange=[0., rout],yrange=[0., 250.],xtitle='r [cm]',$
   ytitle='T_d.',title='',/xsty
oplot,r3,t3.t[8],psym=1, col=!col.red
oplot,r3,t3.t[19],psym=1, col=!col.green

plot,r1,t1.t[0],psym=1,xrange=[0., rout],yrange=[0., 250.],xtitle='r [cm]',$
   ytitle='T_d.',title='',/xsty
oplot,r1,t1.t[8],psym=1, col=!col.red
oplot,r1,t1.t[19],psym=1, col=!col.green

plot,r2,t2.t[0],psym=1,xrange=[0., rout],yrange=[0., 250.],xtitle='r [cm]',$
   ytitle='T_d.',title='',/xsty
oplot,r2,t2.t[8],psym=1, col=!col.red
oplot,r2,t2.t[19],psym=1, col=!col.green

plot,r[*,0,0],temp[*,0,0].t[0],psym=-1,xrange=[0., rout],yrange=[0., 250.],xtitle='r [cm]',$
   ytitle='T_d.',title='',/xsty
oplot,r[*,0,0],temp[*,0,0].t[8],psym=-1, col=!col.red
oplot,r[*,0,0],temp[*,0,0].t[19],psym=-1, col=!col.green

plot,r[0,0,*],temp[0,0,*].t[0],psym=-1,xrange=[0., rout],yrange=[0., 250.],xtitle='r [cm]',$
   ytitle='T_d.',title='',/xsty
oplot,r[0,0,*],temp[0,0,*].t[8],psym=-1, col=!col.red
oplot,r[0,0,*],temp[0,0,*].t[19],psym=-1, col=!col.green

device, /close

end
