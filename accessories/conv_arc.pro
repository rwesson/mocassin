function conv_arc,dist_proj,dist,help=help

if keyword_set(help) then begin
	print,'resul = conv_arc,dist_proj,dist,help=help'
	print,'dist_proj en cm, dist en kpc, resul en arcsec'
	return,0
endif


dist_tab = (dist_proj*0.+dist)*3.086e21

return,3600.*180./3.1415927*atan(dist_proj,dist_tab)
end
