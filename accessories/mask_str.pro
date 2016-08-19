function mask_str,tab,mask,value,help=help

if keyword_set(help) then begin
	print,'function mask_str,tab,mask,[value],help=help'
	return,0
endif

tab_out = tab

if n_params() eq 2 then value = 0.

if n_elements(mask) ne 0 then begin
	if (size(tab))(n_elements(size(tab))-2) eq 8 then begin
	 names = tag_names(tab)
	 print,n_elements(mask),' intensities to cancel'
	 for i = 1,n_tags(tab)-1 do begin
;		tab_tmp = tab.(i)
;		tab_tmp(mask) = value
;		tab_out.(i) = tab_tmp
		tab_out[mask].(i) = value
		print,'changing intensity in ',names(i)
	 endfor
	endif else tab_out(mask) = value
endif
return,tab_out
end
