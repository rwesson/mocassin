pro create_3d,name,data,filename,append=append
   
;   
;   Write in the FILENAME file in a format readable by slicer3
;   name: name of the datas 
;   data : 3D datas   
;   /append: will append to the FILENAME
;   WARNING: If append, all the 3D cube must have the same dimension.   
;
   
if not keyword_set(append) then begin
	spawn,'rm '+filename
	spawn,'touch '+filename
endif

openu,lun,filename,/get_lun,append=append

WRITEU, lun, SIZE(data)
WRITEU, lun, STRLEN(name)
WRITEU, lun, BYTE(name)
WRITEU, lun, data

close,lun
free_lun,lun

print,name,' in ',filename

end
