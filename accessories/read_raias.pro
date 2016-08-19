function read_raias,filename,dim1,dim2,dim3,va,do_x8=do_x8,contcube=contcube,param=param,help=help, $
                    radio=radio,abond=abond,plot=plot,double=double,$ 
                    do_divise=do_divise

;+
;
; function read_raias,filename,dim1,dim2,dim3,param=param,help=help,
; radio=radio,abond=abond
;
;
;-



if keyword_set(help) then begin
  print,'function read_raias,filename,dim1,dim2,dim3,contcube=contcube,param=param,help=help,radio=radio,abond=abond, plot=plot' 
  return,0
endif

if keyword_set(double) then begin
case 1 of
   (keyword_set(contcube)) : $	
         res = {num:0l,xx:0.d,yy:0.d,zz:0.d,im0:fltarr(va+1)}
   (keyword_set(param)) : $
    res = {num:0l,t:0.d,de:0.d,dh:0.d}
   (keyword_set(plot)) : $
    res = {num:0l,r1:0.d}
   (keyword_set(radio)) : $ 
    res = {num:0l,r1:0.d,r2:0.d,r3:0.d,r4:0.d}
   (keyword_set(abond)) : $ 
    res = {h0:0.d,h1:0.d,he0:0.d,he1:0.d,he2:0.d,c0:0.d,c1:0.d,c2:0.d,c3:0.d,c4:0.d,c5:0.d,c6:0.d,$
	n0:0.d,n1:0.d,n2:0.d,n3:0.d,n4:0.d,n5:0.d,n6:0.d,o0:0.d,o1:0.d,o2:0.d,o3:0.d,o4:0.d,o5:0.d,o6:0.d,$
	ne0:0.d,ne1:0.d,ne2:0.d,ne3:0.d,ne4:0.d,ne5:0.d,ne6:0.d, $
mg0:0.d,mg1:0.d,mg2:0.d,mg3:0.d,mg4:0.d,mg5:0.d,mg6:0.d, $
s0:0.d,s1:0.d,s2:0.d,s3:0.d,s4:0.d,s5:0.d,s6:0.d, $
cl0:0.d,cl1:0.d,cl2:0.d,cl3:0.d,cl4:0.d,cl5:0.d,cl6:0.d, $
ar0:0.d,ar1:0.d,ar2:0.d,ar3:0.d,ar4:0.d,ar5:0.d,ar6:0.d, $
fe0:0.d,fe1:0.d,fe2:0.d,fe3:0.d,fe4:0.d,fe5:0.d,fe6:0.d}
   (keyword_set(plot)) : $ 
    res = {num:0l,r1:0.d}
   else : res = {num:0l,r1:0.d,r2:0.d,r3:0.d,r4:0.d,r5:0.d,r6:0.d,$
                 r7:0.d,r8:0.d,r9:0.d,r10:0.d,r11:0.d,r12:0.d}
endcase
endif else begin
   case 1 of
      (keyword_set(contcube)) : $
         res = {num:0l,xx:0.d,yy:0.d,zz:0.d,im0:fltarray(va+1)}
      (keyword_set(param)) : $
       res ={num:0l,t:0.,de:0.,dh:0.}
      (keyword_set(plot)) : $
       res ={num:0l,r1:0.}
      (keyword_set(radio)) : $ 
       res = {num:0l,r1:0.,r2:0.,r3:0.,r4:0.}
      (keyword_set(abond)) : $ 
       res ={h0:0.,h1:0.,he0:0.,he1:0.,he2:0.,c0:0.,c1:0.,c2:0.,c3:0.,c4:0.,c5:0.,c6:0.,n0:0.,n1:0.,n2:0.,$
	n3:0.,n4:0.,n5:0.,n6:0.,o0:0.,o1:0.,o2:0.,o3:0.,o4:0.,o5:0.,o6:0.,ne0:0.,ne1:0.,ne2:0.,$
	ne3:0.,ne4:0.,ne5:0.,ne6:0., $
mg0:0.d,mg1:0.d,mg2:0.d,mg3:0.d,mg4:0.d,mg5:0.d,mg6:0.d, $
s0:0.d,s1:0.d,s2:0.d,s3:0.d,s4:0.d,s5:0.d,s6:0.d, $
cl0:0.d,cl1:0.d,cl2:0.d,cl3:0.d,cl4:0.d,cl5:0.d,cl6:0.d, $
ar0:0.d,ar1:0.d,ar2:0.d,ar3:0.d,ar4:0.d,ar5:0.d,ar6:0.d, $
fe0:0.d,fe1:0.d,fe2:0.d,fe3:0.d,fe4:0.d,fe5:0.d,fe6:0.d}
   endcase
endelse

comment = ''
res_tab = 0
if n_params() eq 3 then res_tab = replicate(res,dim1,dim2) $
else res_tab = replicate(res,dim1,dim2,dim3)

openr,lun,filename,/get_lun

;readf,lun,comment
readf,lun,res_tab

close,lun
free_lun,lun

print,filename,' read.'

resultat =  0
if keyword_set(do_divise) then begin
   res_tab2 =  0
   res_tab2 =  replicate(res,dim1/2,dim2/2,dim3/2)
   for i= 0,n_tags(res_tab)-1 do $
    res_tab2.(i) = rebin(res_tab.(i),dim1/2,dim2/2,dim3/2)
   if keyword_set(do_x8) then resultat =  cube_x8(res,res_tab2) else resultat =  res_tab2
endif else $
 if keyword_set(do_x8) then resultat =  cube_x8(res,res_tab) $
else resultat =  res_tab



return,resultat
end


