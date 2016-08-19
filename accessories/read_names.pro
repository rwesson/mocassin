pro read_names,liste_name,liste_abond,liste_par
   
   name_tmp = {name,id:' ',lam:0.}
   num = 12
   liste_name= replicate(name_tmp,num)


   liste_abond = [ $
       'HI   ','HeI  ','HeII ','CII  ','CIII ','CIV  ','NI   ','NII  ',$
       'NIII ','NIV  ','OI   ','OII  ','OIII ','OIV  ','NeIII','NeIV ',$
       'NeV  ','SII  ','SIII ','SIV  ','ArIII','ArIV ','ArV  ']

   liste_par = ['t','de','dh']

   
end
