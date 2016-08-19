 pro color_key,pos=pos,ysize=ysize,range=range,inc=inc,labels=labels,$
       charsize=charsize,barwidth=barwidth,step=step,title=title,$
       colors=colors,nbotclr=nbotclr,clevels=clevels,slevels=slevels, $
       xoff=xoff,rct=rct
;+
;ROUTINE:     COLOR_KEY
;
;PURPOSE:     draws a numbered color scale.  
;             NOTE: COLOR_KEY is intended to replace PUT_COLOR_SCALE
;
;USEAGE:     
;            COLOR_KEY,pos=pos,ysize=ysize,range=range,inc=inc,$
;                    charsize=charsize,barwidth=barwidth,step=step,$
;                    labels=labels,title=title,colors=colors,nbotclr=nbotclr
;
;PURPOSE:    Draws a numbered color scale
;
;INPUT:      no required input
;
;OPTIONAL KEYWORD INPUT:
;
;  pos          2 element vector, x and y device coordinates of lower 
;               left hand corner of color bar.  If POS is not set the
;               color bar will be placed one character width to the 
;               right of the lower right corner of the last defined 
;               data window
;
;  xoff         x offset of color bar from right edge of plot frame in
;               units of default character with.  (default=1) xoff is
;               disabled when pos is set.
;
;  range        array which contains full range of physical values,
;               The number scale limits are computed fron min(range) and
;               max(range).  (default=0-255)
;
;  inc          increment step of the number scale in physical units
;               If INC is not set the number scale increment is set 
;               automatically
;
;  step         if set to one, the color table is modified with STEP_CT.
;               The color scale is stepped at the number scale intervals. 
;
;  colors       a vector of discreet colors selected from the original
;               color table.  In this case the original color table is
;               not modified and the number of elements of COLORS
;               is used to find the number scale increments. This
;               option should be used with devices that do not accept
;               user specified color tables (such as the Apple Laser
;               Writers).  BEWARE: when the COLORS option is used with
;               TVIM's SCALE option care must be taken to ensure that
;               COLOR_KEY's number scale matches the quantity values.
;               Unlike the STEP option, this is not automatic (see
;               example below)
;
;  nbotclr      number of reserved colors at the bottom of the color
;               table.  these colors will not appear in the color key.
;
;  charsize     character size on number scale
;
;  ysize        vertical size of color bar in device units. 
;               if ysize is negative, the abs of ysize is interpreted
;               as the normalized z coordinate of the top of the bar.
;
;  barwidth     width of color bar (default=1)
;
;  labels       a vector of strings used to label the color key levels
;               if not set the default color key number labels are
;               used.  If the number of elements in LABELS is the same
;               as the number of elements in COLORS then the labels
;               will appear in the middle of a given color band
;               instead of at the boundary between colors.
;
;  title        color scale title drawn to the right of the color key
;
;  slevels      selected levels at which to mark vertical scale. If LABELS
;               is also set, the number of elements in LABELS should be
;               equal to either nsl, or nsl-1 where nsl is the number of
;               elements in SLEVELS.  The SLEVELS input does not work when
;               COLORS is also specified.
;
;  rct          if set, reverse direction of grey scale increase.  This
;               keyword can be used on black & white printers to allow
;               large field values to be represented by dark colors.
;               
;               
;
;KEYWORD OUTPUT
;
;  clevels      physical values at the color key tic marks
;
;  SIDE EFFECTS:
;               if STEP is set COLOR_KEY calls STEP_CT to discreetize
;               the color scale.  To return to original color table
;               enter STEP_CT,/OFF
;
;EXAMPLE:
;
;; on devices which allow user defined color tables:
;
;      loadct,0
;      TVIM,dist(200),/scale
;      COLOR_KEY,range=dist(200),inc=20,/step ; change the color scale
;      COLOR_KEY,range=dist(200),inc=40,/step ; change it again
;
;; on devices with a single hard-wired color table:
;
;      d=dist(200)+20.
;      inc=step_size(d)            ; finds good step size
;      dmin=fix(min(d)/inc)*inc
;      dmax=fix(max(d)/inc)*inc
;      nclrs=(dmax-dmin)/inc
;      colors=!p.color*indgen(nclrs)/(nclrs-1)
;      TVIM,colors((d-dmin)/inc),/scale,/c_map
;      color_key,range=[dmin,dmax],colors=colors
;
;
;; draw one big color key which applies to all tvim images on frame
; 
;      !p.multi=[0,2,2]
;      loadct,5
;      im1=randata(128,128,s=3)
;      im2=rotate(im1,5)
;      im3=rotate(im1,7)
;      im4=rotate(im1,2)
;      range=[min(im1,max=mx),mx]
;      !x.omargin=[0,20]
;      tvim,im1,range=range
;      tvim,im2,range=range
;      tvim,im3,range=range
;      tvim,im4,range=range
;      color_key,range=range,xoff=10,ysize=-.97
;
;; logarithmic scales  
;
;      im=randata(128,128)
;      im=abs(im) & im=3*(im/max(im) > .001)
;      tvim,im,/rmarg ; here we assume im is the log10 of some plot quantity.
;
;      ndec=3       ; number of decades
;
;      slevels=[1.] ; starting value
;      for i=0,ndec-1 do slevels=[slevels,(findgen(9)+2)*10.^i]
;      slevels=alog10(slevels)
;      labels=strarr(n_elements(slevels))
;      labels(indgen(ndec+1)*9)=string(f='("10!a",i1,"!b")',indgen(ndec+1))
;      color_key,range=im,slevels=slevels,labels=labels
;
;
; AUTHOR:       Paul Ricchiazzi    jan93 
;               Earth Space Research Group, UCSB
;
; REVISIONS:
; 20sep93: put in LABELS option
; 21sep93: use pure white blankout color for postscript
; 22sep93: put in COLORS option
; 09sep96: put in SLEVELS option
; 
;-
;
xunit=!d.x_ch_size
if n_elements(xoff) eq 0 then xoff=1
if keyword_set(charsize) eq 0 then charsize=1
if keyword_set(pos) then begin
  x1=pos(0)
  y1=pos(1)
endif else begin
  px=!x.window*!d.x_vsize
  py=!y.window*!d.y_vsize
  x1=px(1)+xunit*xoff
  y1=py(0)
  ys=py(1)-py(0)
endelse

max_color=!d.n_colors-1
if keyword_set(ysize) eq 0 then begin
  if keyword_set(ys) then ysize=ys else ysize=max_color
endif else begin
  if ysize lt 0. then ysize=(-ysize*!d.y_vsize-y1) > 2
endelse

nlabels=n_elements(labels)

if keyword_set(range) then begin
  amin=min(range)
  amax=max(range)
endif else begin
  amin=0.
  amax=255.
endelse
;
s0=float(amin)
s1=float(amax)
ncolors=n_elements(colors)
if ncolors gt 1 then begin
  incc=float(amax-amin)/ncolors
endif else begin
  if keyword_set(inc) eq 0  then begin
    rng=alog10(s1-s0)
    if rng lt 0. then pt=fix(alog10(s1-s0)-.5) else pt=fix(alog10(s1-s0)+.5)
    incc=10.^pt
    tst=[.05,.1,.2,.5,1.,2.,5.,10]
    ii=where((s1-s0)/(incc*tst) le 16,n_count)
    if n_count gt 0 then incc=incc*tst(ii(0)) else incc=1
  endif else begin
    incc=inc
  endelse
  s0=fix(s0/incc)*incc     & if s0 lt amin then s0=s0+incc
  s1=fix(s1/incc)*incc     & if s1 gt amax then s1=s1-incc
endelse

;
frmt='(e9.2)'
nzs=fix(alog10(incc*1.01))
if nzs lt 0 and nzs gt -4 then begin
  frmt='(f8.'+string(form='(i1)',-nzs+1)+')'  ; used on scale
endif
if nzs ge 0 and nzs le 3 then frmt='(f8.1)'
mg=6
if nlabels eq 0 then begin
  smax=string(amax,form=frmt)
  smax=strcompress(smax,/remove_all)
  smin=string(amin,form=frmt)
  smin=strcompress(smin,/remove_all)
  lablen=strlen(smax) > strlen(smin)
endif else begin
  lablen=max(strlen(labels))
endelse

;
if keyword_set(barwidth) eq 0 then barwidth=1.
;
dx=4*xunit*barwidth         ; width of color bar
x2=x1+dx
x3=x2+.5*xunit
mg=.5*xunit                 ; black out margin
dy=ysize                    ; height of color bar
y2=y1+dy
bw=dx+2*mg+charsize*lablen*!d.x_ch_size
bh=dy+2*mg+charsize*!d.y_ch_size


if ncolors eq 0 then begin
  if keyword_set(nbotclr) eq 0 then nbotclr=0 
  if keyword_set(step) then step_ct,[amin,amax],incc,nbotclr=nbotclr
;  clrbar=bytscl(indgen(y2-y1),top=max_color)
  clrbar=nbotclr+(max_color-nbotclr)*findgen(y2-y1)/(y2-y1-1)
endif else begin
  clrbar=colors((ncolors*findgen(y2-y1)/(y2-y1-1)) < (ncolors-1))
endelse

if keyword_set(rct) then clrbar=reverse(clrbar)
;


if !d.name eq 'X' then begin
  polyfill,x1-mg+[0,bw,bw,0,0],y1-mg+[0,0,bh,bh,0],/device,$
              color=!p.background
  tv,(replicate(1b,dx) # clrbar),x1,y1,/device
endif else begin
  tvlct,rr,gg,bb,/get
  tvlct,255,255,255,255
  polyfill,x1-mg+[0,bw,bw,0,0],y1-mg+[0,0,bh,bh,0],/device,color=255
  tvlct,rr,gg,bb
  tv,(replicate(1b,2) # clrbar),x1,y1,xsize=dx,ysize=dy
endelse
;
boxx=[x1,x2,x2,x1,x1]
boxy=[y1,y1,y2,y2,y1]
plots,boxx,boxy,/device
denom=amax-amin
;
nslev=n_elements(slevels)
if nslev ne 0 then begin
  clevels=slevels
  nval=nslev-1
endif else begin
  nval=fix((s1-s0)/incc+.1)
  clevels=s0+findgen(nval+1)*incc
endelse

x4=0.

for ival=0,nval do begin
  val=clevels(ival)
  ss=(val-amin)/denom
  if ss ge 0 and ss le 1 then begin
    yval=y1+(y2-y1)*ss
    plots,[x1,x2],[yval,yval],/device
    if nlabels eq 0 or nlabels gt nval then begin
      if nlabels ge nval then begin
        sval=labels(ival) 
      endif else begin
        sval=string(val,form=frmt)
        sval=strcompress(sval,/remove_all)
      endelse
      xyouts,x3,yval,sval,width=w,/device,charsize=charsize
    endif else begin
      if ival lt nlabels then begin
        sval=labels(ival) 
;        ylab=yval+.5*incc*(y2-y1)/denom
        ylab=y1+(y2-y1)*(.5*(clevels(ival)+clevels(ival+1))-amin)/denom
        nlines=n_elements(str_sep(sval,'!c'))
        ylab=ylab+!d.y_ch_size*charsize*.5*(nlines-1)
        xyouts,x3,ylab,sval,width=w,/device,charsize=charsize
      endif
    endelse
    x4=x4 > w
  endif
endfor

if keyword_set(title) then begin
  xtitle=x3+x4*!d.x_vsize+1.5*xunit
  ytitle=y1+.5*ysize
  xyouts,xtitle,ytitle,/device,align=.5,orientation=-90,title,$
                charsize=charsize
endif
;
end
