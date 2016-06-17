function cube_x8,stru,tab_in
   
;   
;  Will make an 8 times greater cube assuming 3-axis symetry
;  stru: model of the structure array of tab_in   
;
   
   s = size(tab_in)
   
   if s(1) ne s(2) or s(2) ne s(3) or s(1) ne s(3) then return,-1
   
   siz = 2 * s(1) - 1
      
   if (size(tab_in))(n_elements(size(tab_in))-2) eq 8 then begin ; tab_in is a structure
      tab_out =  0
      tab_out =  replicate(stru,siz,siz,siz)
      tagnames =  tag_names(tab_in)
      for i = 0,n_tags(tab_in)-1 do begin
         tab_out.(i) = cube_x8(stru,tab_in.(i))
         print,'cube_x8 of: ',tagnames(i)
      endfor
      return,tab_out
   endif else begin
       
      tab_out = dblarr(siz,siz,siz)
      
      tab_out(0:siz/2,0:siz/2,0:siz/2) = reverse(reverse(reverse(tab_in,3),2),1)
      tab_out(0:siz/2,0:siz/2,siz/2:siz-1) = reverse(reverse(tab_in,2),1)
      tab_out(0:siz/2,siz/2:siz-1,0:siz/2) = reverse(reverse(tab_in,3),1)
      tab_out(0:siz/2,siz/2:siz-1,siz/2:siz-1) = reverse(tab_in,1)
      tab_out(siz/2:siz-1,0:siz/2,0:siz/2) = reverse(reverse(tab_in,2),3)
      tab_out(siz/2:siz-1,0:siz/2,siz/2:siz-1) = reverse(tab_in,2)
      tab_out(siz/2:siz-1,siz/2:siz-1,0:siz/2) = reverse(tab_in,3)
      tab_out(siz/2:siz-1,siz/2:siz-1,siz/2:siz-1) = tab_in
      
      return,tab_out
   endelse
end
