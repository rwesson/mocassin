pro prepare_slicer3,name,data,filename,append=append

;+
; Purpose: will save 3D arrays in a file with the
; apropriate format to be read by slicer3
;
; Calling sequence: prepare_slicer3,name,data,filename,[append=append]
;
; name : name of the datas (string)
; data : 3D array to save
; filename : name of the file to be read by slicer3 (string)
; append : if not set, will erase 'filename' before saving
;
; Morisset, IAG/USP, 1997 morisset@iagusp.usp.br
;
;-

if not keyword_set(append) then spawn,'rm '+filename

openu,lun,filename,/get_lun,/append

WRITEU, lun, SIZE(data)
WRITEU, lun, STRLEN(name)
WRITEU, lun, BYTE(name)
WRITEU, lun, data

close,lun
free_lun,lun

end