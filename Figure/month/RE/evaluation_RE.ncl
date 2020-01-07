; evaluating the reconstucted monthly data=================
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

;================set the parameters=====================

diri1="/Users/chang/Desktop/Small_work/permafrost/Feng/data/"
name1="GLDAS_FI_Feng.nc"
diri2="/Users/chang/Desktop/Small_work/permafrost/Oliver/data/"
name2 = "GLDAS_FI_Oliver.nc"
diri3="/Users/chang/Desktop/Small_work/permafrost/New_Method/monthly/data/Frost_index_data/"
name3 = "GLDAS_FI_New_ensemble.nc"
diri4="/Users/chang/Desktop/Small_work/permafrost/GLDAS/data/"
name4="GLDAS_FI_ori.nc"

Year1=1901
Year2=2011
Year=ispan(Year1,Year2,1)


out_data_ddt=new((/3,5/),"float")
out_data_ddf=new((/3,5/),"float")
out_data_fi=new((/3,5/),"float")

;=================read the data, because the value of 2012 is missing, we only use data from 1901-2011==========================

  f1  = addfile(diri1+name1,"r")
  lat=f1->lat
  lon=f1->lon
  year_1=Year(:110)
  fi_Feng = f1->fi(:110,:,:)
  ddt_Feng = f1->ddt(:110,:,:)
  ddf_Feng = f1->ddf(:110,:,:)
  printVarSummary(fi_Feng)

  f2  = addfile(diri2+name2,"r")
  fi_Oliver= f2->fi(:110,:,:)
  ddt_Oliver= f2->ddt(:110,:,:)
  ddf_Oliver= f2->ddf(:110,:,:)
  printVarSummary(fi_Oliver)

  f3  = addfile(diri3+name3,"r")
  fi_New= f3->fi(:110,:,:)
  ddt_New= f3->ddt(:110,:,:)
  ddf_New= f3->ddf(:110,:,:)
  printVarSummary(fi_New)

  f4  = addfile(diri4+name4,"r")
  fi_GLDAS= f4->fi(:110,:,:)
  ddt_GLDAS= f4->ddt(:110,:,:)
  ddf_GLDAS= f4->ddf(:110,:,:)
  printVarSummary(fi_GLDAS)


 a = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
  lsdata = a->LSMASK
  lsm = landsea_mask(lsdata,lat,lon)

  fi_GLDAS  = mask(fi_GLDAS ,lsm.eq.0,False)  
  fi_GLDAS  = mask(fi_GLDAS ,lsm.eq.2,False)  
  fi_GLDAS  = mask(fi_GLDAS ,lsm.eq.4,False)  
 
ddt_GLDAS = where(ismissing(fi_GLDAS),fi_GLDAS@_FillValue,ddt_GLDAS )
ddf_GLDAS = where(ismissing(fi_GLDAS),fi_GLDAS@_FillValue,ddf_GLDAS )

fi_New = where(ismissing(fi_GLDAS),fi_GLDAS@_FillValue,fi_New )
ddt_New = where(ismissing(fi_GLDAS),fi_GLDAS@_FillValue,ddt_New )
ddf_New = where(ismissing(fi_GLDAS),fi_GLDAS@_FillValue,ddf_New )

fi_Oliver = where(ismissing(fi_GLDAS),fi_GLDAS@_FillValue,fi_Oliver )
ddt_Oliver = where(ismissing(fi_GLDAS),fi_GLDAS@_FillValue,ddt_Oliver )
ddf_Oliver = where(ismissing(fi_GLDAS),fi_GLDAS@_FillValue,ddf_Oliver )

fi_Feng = where(ismissing(fi_GLDAS),fi_GLDAS@_FillValue,fi_Feng )
ddt_Feng = where(ismissing(fi_GLDAS),fi_GLDAS@_FillValue,ddt_Feng )
ddf_Feng = where(ismissing(fi_GLDAS),fi_GLDAS@_FillValue,ddf_Feng )


;================= using the distribution of permafrost to mask fi ddt and ddf==================

fi_GLDAS_ave=dim_avg_n_Wrap(fi_GLDAS,0)
fi_GLDAS_ave=where(fi_GLDAS_ave.eq.0,fi_GLDAS_ave@_FillValue,fi_GLDAS_ave)
;fi_GLDAS_ave=where(fi_GLDAS_ave.le.0.495,fi_GLDAS_ave@_FillValue,fi_GLDAS_ave)
fi_GLDAS_ave_test=conform_dims(dimsizes(fi_GLDAS),fi_GLDAS_ave,(/1,2/))

;ddt_GLDAS=where(ismissing(fi_GLDAS_ave_test),ddt_GLDAS@_FillValue,ddt_GLDAS)
;ddf_GLDAS=where(ismissing(fi_GLDAS_ave_test),ddf_GLDAS@_FillValue,ddf_GLDAS)

;fi_New=where(ismissing(fi_GLDAS_ave_test),fi_New@_FillValue,fi_New)
;ddt_New=where(ismissing(fi_GLDAS_ave_test),ddt_New@_FillValue,ddt_New)
;ddf_New=where(ismissing(fi_GLDAS_ave_test),ddf_New@_FillValue,ddf_New)

;fi_Oliver=where(ismissing(fi_GLDAS_ave_test),fi_Oliver@_FillValue,fi_Oliver)
;ddt_Oliver=where(ismissing(fi_GLDAS_ave_test),ddt_Oliver@_FillValue,ddt_Oliver)
;ddf_Oliver=where(ismissing(fi_GLDAS_ave_test),ddf_Oliver@_FillValue,ddf_Oliver)

;fi_Feng=where(ismissing(fi_GLDAS_ave_test),fi_Feng@_FillValue,fi_Feng)
;ddt_Feng=where(ismissing(fi_GLDAS_ave_test),ddt_Feng@_FillValue,ddt_Feng)
;ddf_Feng=where(ismissing(fi_GLDAS_ave_test),ddf_Feng@_FillValue,ddf_Feng)


;================calculate the relative error========================================================


re_Feng=dim_sum_n_Wrap(abs(fi_Feng-fi_GLDAS),0)
RE_fi_Feng=re_Feng/(dimsizes(Year)*fi_GLDAS_ave)*100
copy_VarCoords(fi_GLDAS_ave,RE_fi_Feng) 
re_Oliver=dim_sum_n_Wrap(abs(fi_Oliver-fi_GLDAS),0)
RE_fi_Oliver=re_Oliver/(dimsizes(Year)*fi_GLDAS_ave)*100
copy_VarCoords(fi_GLDAS_ave,RE_fi_Oliver) 
re_New=dim_sum_n_Wrap(abs(fi_New-fi_GLDAS),0)
RE_fi_New=re_New/(dimsizes(Year)*fi_GLDAS_ave)*100
copy_VarCoords(fi_GLDAS_ave,RE_fi_New) 


ddt_GLDAS_ave=dim_avg_n_Wrap(ddt_GLDAS,0)
ddt_GLDAS_ave=where(ddt_GLDAS_ave.eq.0,fi_GLDAS_ave@_FillValue,ddt_GLDAS_ave)
re_Feng=dim_sum_n_Wrap(abs(ddt_Feng-ddt_GLDAS),0)
RE_ddt_Feng=re_Feng/(dimsizes(Year)*ddt_GLDAS_ave)*100
copy_VarCoords(ddt_GLDAS_ave,RE_ddt_Feng) 
re_Oliver=dim_sum_n_Wrap(abs(ddt_Oliver-ddt_GLDAS),0)

RE_ddt_Oliver=re_Oliver/(dimsizes(Year)*ddt_GLDAS_ave)*100
copy_VarCoords(ddt_GLDAS_ave,RE_ddt_Oliver) 
re_New=dim_sum_n_Wrap(abs(ddt_New-ddt_GLDAS),0)
RE_ddt_New=re_New/(dimsizes(Year)*ddt_GLDAS_ave)*100
copy_VarCoords(ddt_GLDAS_ave,RE_ddt_New) 


ddf_GLDAS_ave=dim_avg_n_Wrap(ddf_GLDAS,0)
ddf_GLDAS_ave=where(ddf_GLDAS_ave.eq.0,fi_GLDAS_ave@_FillValue,ddf_GLDAS_ave)
re_Feng=dim_sum_n_Wrap(abs(ddf_Feng-ddf_GLDAS),0)
RE_ddf_Feng=re_Feng/(dimsizes(Year)*ddf_GLDAS_ave)*100
copy_VarCoords(ddf_GLDAS_ave,RE_ddf_Feng) 
re_Oliver=dim_sum_n_Wrap(abs(ddf_Oliver-ddf_GLDAS),0)
RE_ddf_Oliver=re_Oliver/(dimsizes(Year)*ddf_GLDAS_ave)*100
copy_VarCoords(ddf_GLDAS_ave,RE_ddf_Oliver) 
re_New=dim_sum_n_Wrap(abs(ddf_New-ddf_GLDAS),0)
RE_ddf_New=re_New/(dimsizes(Year)*ddf_GLDAS_ave)*100
copy_VarCoords(ddf_GLDAS_ave,RE_ddf_New) 




plot_re=new(9,"graphic")
wks = gsn_open_wks("pdf" ,"Measures_RE")      ; send graphics to PNG file
  gsn_define_colormap(wks,"amwg_blueyellowred")

  res                       = True     ; plot mods desired
  res@cnFillOn              = True     ; turn on color fill
  res@cnLinesOn             = False    ; turn of contour lines
  plot=new(10,graphic)  
  res@gsnDraw = False 
  res@gsnFrame = False 
  res@mpCenterLonF=-45
  res@gsnPolar   = "NH"  
  res@mpMinLatF            = 25   ; range to zoom in on
  res@mpFillOn             = False              
  res@cnFillOn              = True     ; turn on color fill
  res@cnLinesOn             = False    ; turn of contour lines
  ;res@cnFillMode            = "RasterFill"    
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnFillColors        =(/2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17/)
  res@lbLabelBarOn        = False
  res@gsnLeftStringFontHeightF = 0.025
  res@cnLevels = (/0.1,0.25,0.5,0.75,1.,5.,10.,20.,30.,40.,45.,50.,60.,70.,80./)
  res@gsnLeftString ="(a) FI (Frauenfeld)"
  plot_re(0)=gsn_csm_contour_map_polar(wks,RE_fi_Oliver,res)
  res@gsnLeftString ="(b) FI (Vincent)"
  plot_re(1)=gsn_csm_contour_map_polar(wks,RE_fi_Feng,res)
  res@gsnLeftString ="(c) FI (Monte Carlo monthly)"
  plot_re(2)=gsn_csm_contour_map_polar(wks,RE_fi_New,res)
  res@gsnLeftString ="(d) DDT (Frauenfeld)"
  plot_re(3)=gsn_csm_contour_map_polar(wks,RE_ddt_Oliver,res)
  res@gsnLeftString =  "(e) DDT (Vincent)"
  plot_re(4)=gsn_csm_contour_map_polar(wks,RE_ddt_Feng,res)
  res@gsnLeftString ="(f) DDT (Monte Carlo monthly)"
  plot_re(5)=gsn_csm_contour_map_polar(wks,RE_ddt_New,res)
  res@gsnLeftString ="(g) DDF (Frauenfeld)"
  plot_re(6)=gsn_csm_contour_map_polar(wks,RE_ddf_Oliver,res)
  res@gsnLeftString = "(h) DDF (Vincent)"
  plot_re(7)=gsn_csm_contour_map_polar(wks,RE_ddf_Feng,res)
  res@gsnLeftString ="(i) DDF (Monte Carlo monthly)"
  plot_re(8)=gsn_csm_contour_map_polar(wks,RE_ddf_New,res)



  resP                     = True         ; modify the panel plot
  resP@gsnFrame            = False
  resP@gsnPanelTop   = 0.85   
  resP@gsnPanelBottom   = 0.15
  resP@gsnPanelYWhiteSpacePercent = 5
  resP@gsnPanelXWhiteSpacePercent = 5
  resP@lbLabelFontHeightF  = 0.007  
  resP@gsnPanelLabelBar    = True       ; common label bar
  resP@lbOrientation       = "horizontal" ; vertical label bar
  resP@pmLabelBarWidthF    = 0.6      ; default is shorter
  resP@pmLabelBarHeightF   = 0.055        ; default is taller
  resP@lbLabelFontHeightF  = 0.008       ; make labels smaller ( default=0.02 )
  gsn_panel(wks,plot_re,(/3,3/),resP) 
  delete(resP)

;===============output data=======================

  rad    = 4.0*atan(1.0)/180.0
  re     = 6371220.0
  rr     = re*rad
  clat   = cos(lat*rad)

 RE_fi_Feng_avg = wgt_areaave(RE_fi_Feng, clat, 1.0, 0)
 RE_ddt_Feng_avg = wgt_areaave(RE_ddt_Feng, clat, 1.0, 0)
 RE_ddf_Feng_avg = wgt_areaave(RE_ddf_Feng, clat, 1.0, 0)

 RE_fi_Oliver_avg = wgt_areaave(RE_fi_Oliver, clat, 1.0, 0)
 RE_ddt_Oliver_avg = wgt_areaave(RE_ddt_Oliver, clat, 1.0, 0)
 RE_ddf_Oliver_avg = wgt_areaave(RE_ddf_Oliver, clat, 1.0, 0)

 RE_fi_New_avg = wgt_areaave(RE_fi_New, clat, 1.0, 0)
 RE_ddt_New_avg = wgt_areaave(RE_ddt_New, clat, 1.0, 0)
 RE_ddf_New_avg = wgt_areaave(RE_ddf_New, clat, 1.0, 0)

out_data=new((/3,3/),"float")

out_data(0,0)=RE_ddt_Oliver_avg
out_data(1,0)=RE_ddt_Feng_avg
out_data(2,0)=RE_ddt_New_avg

out_data(0,1)=RE_ddf_Oliver_avg
out_data(1,1)=RE_ddf_Feng_avg
out_data(2,1)=RE_ddf_New_avg

out_data(0,2)=RE_fi_Oliver_avg
out_data(1,2)=RE_fi_Feng_avg
out_data(2,2)=RE_fi_New_avg

 opt = True
 opt@title  = "        DDT        DDF       FI    "
 opt@fout   = "RE_25N.txt"
 write_matrix (out_data, "3f12.3", opt)






fi_GLDAS_ave=where(fi_GLDAS_ave.le.0.495,fi_GLDAS_ave@_FillValue,fi_GLDAS_ave)

RE_fi_New=where(ismissing(fi_GLDAS_ave),fi_New@_FillValue,RE_fi_New)
RE_ddt_New=where(ismissing(fi_GLDAS_ave),ddt_New@_FillValue,RE_ddt_New)
RE_ddf_New=where(ismissing(fi_GLDAS_ave),ddf_New@_FillValue,RE_ddf_New)

RE_fi_Oliver=where(ismissing(fi_GLDAS_ave),fi_Oliver@_FillValue,RE_fi_Oliver)
RE_ddt_Oliver=where(ismissing(fi_GLDAS_ave),ddt_Oliver@_FillValue,RE_ddt_Oliver)
RE_ddf_Oliver=where(ismissing(fi_GLDAS_ave),ddf_Oliver@_FillValue,RE_ddf_Oliver)

RE_fi_Feng=where(ismissing(fi_GLDAS_ave),fi_Feng@_FillValue,RE_fi_Feng)
RE_ddt_Feng=where(ismissing(fi_GLDAS_ave),ddt_Feng@_FillValue,RE_ddt_Feng)
RE_ddf_Feng=where(ismissing(fi_GLDAS_ave),ddf_Feng@_FillValue,RE_ddf_Feng)


  rad    = 4.0*atan(1.0)/180.0
  re     = 6371220.0
  rr     = re*rad
  clat   = cos(lat*rad)

 RE_fi_Feng_avg = wgt_areaave(RE_fi_Feng, clat, 1.0, 0)
 RE_ddt_Feng_avg = wgt_areaave(RE_ddt_Feng, clat, 1.0, 0)
 RE_ddf_Feng_avg = wgt_areaave(RE_ddf_Feng, clat, 1.0, 0)

 RE_fi_Oliver_avg = wgt_areaave(RE_fi_Oliver, clat, 1.0, 0)
 RE_ddt_Oliver_avg = wgt_areaave(RE_ddt_Oliver, clat, 1.0, 0)
 RE_ddf_Oliver_avg = wgt_areaave(RE_ddf_Oliver, clat, 1.0, 0)

 RE_fi_New_avg = wgt_areaave(RE_fi_New, clat, 1.0, 0)
 RE_ddt_New_avg = wgt_areaave(RE_ddt_New, clat, 1.0, 0)
 RE_ddf_New_avg = wgt_areaave(RE_ddf_New, clat, 1.0, 0)

out_data=new((/3,3/),"float")

out_data(0,0)=RE_ddt_Oliver_avg
out_data(1,0)=RE_ddt_Feng_avg
out_data(2,0)=RE_ddt_New_avg

out_data(0,1)=RE_ddf_Oliver_avg
out_data(1,1)=RE_ddf_Feng_avg
out_data(2,1)=RE_ddf_New_avg

out_data(0,2)=RE_fi_Oliver_avg
out_data(1,2)=RE_fi_Feng_avg
out_data(2,2)=RE_fi_New_avg

 opt = True
 opt@title  = "        DDT        DDF       FI    "
 opt@fout   = "RE_permafrost.txt"
 write_matrix (out_data, "3f12.3", opt)

exit


RE_ddf_New=where(RE_ddf_New.ge.100,RE_ddf_New@_FillValue,RE_ddf_New)
RE_ddt_New=where(RE_ddt_New.ge.100,RE_ddt_New@_FillValue,RE_ddt_New)
RE_fi_New=where(RE_fi_New.ge.100,RE_fi_New@_FillValue,RE_fi_New)

RE_ddf_Feng=where(RE_ddf_Feng.ge.100,RE_ddf_Feng@_FillValue,RE_ddf_Feng)
RE_ddt_Feng=where(RE_ddt_Feng.ge.100,RE_ddt_Feng@_FillValue,RE_ddt_Feng)
RE_fi_Feng=where(RE_fi_Feng.ge.100,RE_fi_Feng@_FillValue,RE_fi_Feng)

RE_ddf_Oliver=where(RE_ddf_Oliver.ge.100,RE_ddf_Oliver@_FillValue,RE_ddf_Oliver)
RE_ddt_Oliver=where(RE_ddt_Oliver.ge.100,RE_ddt_Oliver@_FillValue,RE_ddt_Oliver)
RE_fi_Oliver=where(RE_fi_Oliver.ge.100,RE_fi_Oliver@_FillValue,RE_fi_Oliver)

