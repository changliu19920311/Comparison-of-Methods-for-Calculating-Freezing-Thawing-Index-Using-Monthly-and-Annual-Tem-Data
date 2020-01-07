
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;================set the parameters=====================


name1="GLDAS_FI_New_ensemble.nc"
diri3="/Users/chang/Desktop/Small_work/permafrost/New_Method/annual/data/Frost_index_data/"
diri4="//Users/chang/Desktop/Small_work/permafrost/GLDAS/data/"
name2="GLDAS_FI_ori.nc"

Year1=1901
Year2=2011
Year=ispan(Year1,Year2,1)

;=================read the reconstructed data ==========================


  f3  = addfile(diri3+name1,"r")
  lat=f3->lat
  lon=f3->lon
  fi_New= f3->fi
  ddt_New= f3->ddt
  ddf_New= f3->ddf
  printVarSummary(fi_New)

  f4  = addfile(diri4+name2,"r")
  fi_GLDAS= f4->fi
  ddt_GLDAS= f4->ddt
  ddf_GLDAS= f4->ddf
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


;================calculate the relative error========================================================

fi_GLDAS_ave=dim_avg_n_Wrap(fi_GLDAS,0)
fi_GLDAS_ave=where(fi_GLDAS_ave.eq.0,fi_GLDAS_ave@_FillValue,fi_GLDAS_ave)
re_New=dim_sum_n_Wrap(abs(fi_New-fi_GLDAS),0)
RE_fi_New=re_New/(dimsizes(Year)*fi_GLDAS_ave)*100
copy_VarCoords(fi_GLDAS_ave,RE_fi_New) 


ddt_GLDAS_ave=dim_avg_n_Wrap(ddt_GLDAS,0)
ddt_GLDAS_ave=where(ddt_GLDAS_ave.eq.0,fi_GLDAS_ave@_FillValue,ddt_GLDAS_ave)
re_New=dim_sum_n_Wrap(abs(ddt_New-ddt_GLDAS),0)
RE_ddt_New=re_New/(dimsizes(Year)*ddt_GLDAS_ave)*100
copy_VarCoords(ddt_GLDAS_ave,RE_ddt_New) 


ddf_GLDAS_ave=dim_avg_n_Wrap(ddf_GLDAS,0)
ddf_GLDAS_ave=where(ddf_GLDAS_ave.eq.0,fi_GLDAS_ave@_FillValue,ddf_GLDAS_ave)
re_New=dim_sum_n_Wrap(abs(ddf_New-ddf_GLDAS),0)
RE_ddf_New=re_New/(dimsizes(Year)*ddf_GLDAS_ave)*100
copy_VarCoords(ddf_GLDAS_ave,RE_ddf_New) 
;====================calculate Root mean square error===============================================

rmse_fi_New = dim_rmsd_n(fi_New,fi_GLDAS,0)*100
copy_VarCoords(fi_GLDAS_ave,rmse_fi_New) 
rmse_ddt_New = dim_rmsd_n(ddt_New,ddt_GLDAS,0)*100
copy_VarCoords(fi_GLDAS_ave,rmse_ddt_New) 
rmse_ddf_New = dim_rmsd_n(ddf_New,ddf_GLDAS,0)*100
copy_VarCoords(fi_GLDAS_ave,rmse_ddf_New) 


;==============================correlation=====================================================

    n=dimsizes(Year)
    df=n-2
   

    rr_fi_New   = escorc_n(fi_New,fi_GLDAS,0,0) 
    rr_fi_New = where(rr_fi_New.eq.1,0.999,rr_fi_New)
    t    = rr_fi_New*sqrt((n-2)/(1-rr_fi_New^2))      
    p_rr_fi_New = student_t(t, df)
    delete(t)
    copy_VarCoords(fi_GLDAS(0,:,:),rr_fi_New) 
    copy_VarCoords(fi_GLDAS(0,:,:),p_rr_fi_New) 

    rr_ddt_New   = escorc_n(ddt_New,ddt_GLDAS,0,0) 
    rr_ddt_New = where(rr_ddt_New.eq.1,0.999,rr_ddt_New)
    t    = rr_ddt_New*sqrt((n-2)/(1-rr_ddt_New^2))      
    p_rr_ddt_New = student_t(t, df)
    delete(t)
    copy_VarCoords(ddt_GLDAS(0,:,:),rr_ddt_New) 
    copy_VarCoords(ddt_GLDAS(0,:,:),p_rr_ddt_New) 

    rr_ddf_New   = escorc_n(ddf_New,ddf_GLDAS,0,0) 
    rr_ddf_New = where(rr_ddf_New.eq.1,0.999,rr_ddf_New)
    t    = rr_ddf_New*sqrt((n-2)/(1-rr_ddf_New^2))      
    p_rr_ddf_New = student_t(t, df)
    delete(t)
    copy_VarCoords(ddf_GLDAS(0,:,:),rr_ddf_New) 
    copy_VarCoords(ddf_GLDAS(0,:,:),p_rr_ddf_New) 
;==============================Difference=====================================================


ddf_GLDAS_ave=dim_avg_n_Wrap(ddf_GLDAS,0)
ddt_GLDAS_ave=dim_avg_n_Wrap(ddt_GLDAS,0)
fi_GLDAS_ave=dim_avg_n_Wrap(fi_GLDAS,0)

ddf_New_ave=dim_avg_n_Wrap(ddf_New,0)
ddt_New_ave=dim_avg_n_Wrap(ddt_New,0)
fi_New_ave=dim_avg_n_Wrap(fi_New,0)


dif_fi_New = fi_New_ave-fi_GLDAS_ave
dif_ddf_New = ddf_New_ave-ddf_GLDAS_ave
dif_ddt_New = ddt_New_ave-ddt_GLDAS_ave

copy_VarCoords(ddf_GLDAS_ave,dif_fi_New) 
copy_VarCoords(ddf_GLDAS_ave,dif_ddf_New) 
copy_VarCoords(ddf_GLDAS_ave,dif_ddt_New)

;==============================plot======================================================================

fo=addfile("merit.nc","c")
fo->dif_fi_New=dif_fi_New
fo->dif_ddf_New=dif_ddf_New
fo->dif_ddt_New=dif_ddt_New

fo->RE_ddt_New=RE_ddt_New
fo->RE_ddf_New=RE_ddf_New
fo->RE_fi_New=RE_fi_New

fo->rr_ddf_New=rr_ddf_New
fo->rr_ddt_New=rr_ddt_New
fo->rr_fi_New=rr_fi_New

plot_re=new(3,"graphic")
plot_rmse=new(3,"graphic")
plot_rr=new(3,"graphic")
plot_dif=new(3,"graphic")


wks = gsn_open_wks("pdf" ,"Measures_annual")      ; send graphics to PNG file
  gsn_define_colormap(wks,"amwg_blueyellowred")

  res                       = True     ; plot mods desired
  res@cnFillOn              = True     ; turn on color fill
  res@cnLinesOn             = False    ; turn of contour lines
  res@gsnDraw = False 
  res@gsnFrame = False 
  res@mpCenterLonF=-45
  res@gsnPolar   = "NH"  
  res@mpMinLatF            = 25     ; range to zoom in on
  res@mpFillOn             = False              
  res@cnFillOn              = True     ; turn on color fill
  res@cnLinesOn             = False    ; turn of contour lines
  ;res@cnFillMode            = "RasterFill"    
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnFillColors        =(/2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17/)
  res@lbLabelBarOn        = False
  res@gsnLeftStringFontHeightF = 0.025
  res@cnLevels = (/0.1,0.25,0.5,0.75,1.,5.,10.,20.,30.,40.,45.,50.,60.,70.,80./)
  res@gsnLeftString ="(c) FI (RE)"
  plot_re(2)=gsn_csm_contour_map_polar(wks,RE_fi_New,res)
  res@gsnLeftString ="(b) DDT (RE)"
  plot_re(1)=gsn_csm_contour_map_polar(wks,RE_ddt_New,res)
  res@gsnLeftString ="(a) DDF (RE)"
  plot_re(0)=gsn_csm_contour_map_polar(wks,RE_ddf_New,res)
  resP                     = True         ; modify the panel plot
  resP@gsnFrame            = False
  resP@gsnPanelTop   = 0.95   
  resP@gsnPanelBottom   = 0.65
  resP@gsnPanelYWhiteSpacePercent = 5
  resP@gsnPanelXWhiteSpacePercent = 5
  resP@lbLabelFontHeightF  = 0.007  
  resP@gsnPanelLabelBar    = True       ; common label bar
  resP@lbOrientation       = "horizontal" ; vertical label bar
  resP@pmLabelBarWidthF    = 0.6      ; default is shorter
  resP@pmLabelBarHeightF   = 0.055        ; default is taller
  resP@lbLabelFontHeightF  = 0.008       ; make labels smaller ( default=0.02 )
  gsn_panel(wks,plot_re,(/1,3/),resP) 
  delete(resP)


  res@cnLevels = (/0.1,0.2,0.5,1.,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0/)
  plot_rmse(2)=gsn_csm_contour_map_polar(wks,rmse_fi_New,res)
  res@cnLevels = (/500.,1000.,1500.,2000.,2500.,3000.,3500.,4000.,4500.,5000.,5500.,6000.,6500.,7000.,7500./)
  plot_rmse(1)=gsn_csm_contour_map_polar(wks,rmse_ddt_New,res)
  res@cnLevels = (/500.,1000.,1500.,2000.,2500.,3000.,3500.,4000.,4500.,5000.,5500.,6000.,6500.,7000.,7500./)
  plot_rmse(0)=gsn_csm_contour_map_polar(wks,rmse_ddf_New,res)



  res@cnLevels := (/0.1,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.9/)
  res@gsnLeftString ="(f) FI (Cor)"
  plot_rr(2)=gsn_csm_contour_map_polar(wks,rr_fi_New,res)
  res@gsnLeftString ="(e) DDT (Cor)"
  plot_rr(1)=gsn_csm_contour_map_polar(wks,rr_ddt_New,res)
  res@gsnLeftString ="(d) DDF (Cor)"
  plot_rr(0)=gsn_csm_contour_map_polar(wks,rr_ddf_New,res)

  resP                     = True         ; modify the panel plot
  resP@gsnFrame            = False
  resP@gsnPanelTop   = 0.65   
  resP@gsnPanelBottom   = 0.35
  resP@gsnPanelYWhiteSpacePercent = 5
  resP@gsnPanelXWhiteSpacePercent = 5
  resP@lbLabelFontHeightF  = 0.007  
  resP@gsnPanelLabelBar    = True       ; common label bar
  resP@lbOrientation       = "horizontal" ; vertical label bar
  resP@pmLabelBarWidthF    = 0.6      ; default is shorter
  resP@pmLabelBarHeightF   = 0.055        ; default is taller
  resP@lbLabelFontHeightF  = 0.008       ; make labels smaller ( default=0.02 )
  gsn_panel(wks,plot_rr,(/1,3/),resP) 
  delete(resP)

  res@cnLevels :=(/-0.07,-0.06,-0.05,-0.04,-0.03,-0.02,-0.01,0,0.01,0.02,0.03,0.04,0.05,0.06,0.07/)
  res@gsnLeftString ="(i) FI (MBE)"
  plot_dif(2)=gsn_csm_contour_map_polar(wks,dif_fi_New,res)
  res@cnLevels := (/-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70/)
  res@gsnLeftString ="(h) DDT (MBE)"
  plot_dif(1)=gsn_csm_contour_map_polar(wks,dif_ddt_New,res)
  res@gsnLeftString ="(g) DDF (MBE)"
  plot_dif(0)=gsn_csm_contour_map_polar(wks,dif_ddf_New,res)


  resP                     = True         ; modify the panel plot
  resP@gsnFrame            = False
  resP@gsnPanelTop   = 0.35   
  resP@gsnPanelBottom   = 0.05
  resP@gsnPanelYWhiteSpacePercent = 5
  resP@gsnPanelXWhiteSpacePercent = 5
  resP@lbLabelFontHeightF  = 0.007  
  resP@gsnPanelLabelBar    = True       ; common label bar
  resP@lbOrientation       = "horizontal" ; vertical label bar
  resP@pmLabelBarWidthF    = 0.6      ; default is shorter
  resP@pmLabelBarHeightF   = 0.055        ; default is taller
  resP@lbLabelFontHeightF  = 0.008       ; make labels smaller ( default=0.02 )
  gsn_panel(wks,plot_dif,(/1,3/),resP) 
  delete(resP)

  nboxes  = 16
  labels  =  (/"-0.07","-0.06","-0.05","-0.04","-0.03","-0.02","-0.01","0","0.01","0.02","0.03","0.04","0.05","0.06","0.07"/)

  lbar_x1 = 0.20
  lbar_y1 = 0.06

  lbres                      = True
  lbres@vpHeightF            = 0.055
  lbres@vpWidthF             = 0.6

  lbres@lbAutoManage         = False
  lbres@lbOrientation        = "horizontal"
  lbres@lbLabelPosition      = "Bottom"
  lbres@lbMonoFillPattern         = True

  lbres@lbLabelsOn           = True
  lbres@lbPerimOn            = False
  lbres@lbTitleOn            = False

  lbres@lbLabelAlignment     = "InteriorEdges"
  lbres@lbBoxLinesOn         = True    ; True is the default
  lbres@lbFillColors         =  (/2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17/)
  lbres@lbLabelFontHeightF   = 0.008

  lbres@lbLabelJust          = "CenterLeft"
  lbid1 = gsn_create_labelbar_ndc(wks,nboxes,labels,lbar_x1,lbar_y1,lbres)
  draw(lbid1)


  txres               = True                     ; text mods desired
  txres@txFontHeightF = 0.01                     ; font smaller. default big
  gsn_text_ndc(wks,"DDT/DDF",0.16,0.09,txres) 
  gsn_text_ndc(wks,"FI",0.19,0.04,txres) 

;===============output data=======================

  rad    = 4.0*atan(1.0)/180.0
  re     = 6371220.0
  rr     = re*rad
  clat   = cos(lat*rad)

 dif_fi_New_avg = wgt_areaave(dif_fi_New, clat, 1.0, 0)
 dif_ddt_New_avg = wgt_areaave(dif_ddt_New, clat, 1.0, 0)
 dif_ddf_New_avg = wgt_areaave(dif_ddf_New, clat, 1.0, 0)

 rr_fi_New_avg = wgt_areaave(rr_fi_New, clat, 1.0, 0)
 rr_ddt_New_avg = wgt_areaave(rr_ddt_New, clat, 1.0, 0)
 rr_ddf_New_avg = wgt_areaave(rr_ddf_New, clat, 1.0, 0)

 re_fi_New_avg = wgt_areaave(RE_fi_New, clat, 1.0, 0)
 re_ddt_New_avg = wgt_areaave(RE_ddt_New, clat, 1.0, 0)
 re_ddf_New_avg = wgt_areaave(RE_ddf_New, clat, 1.0, 0)

out_data=new((/3,3/),"float")

out_data(0,0)=dif_ddt_New_avg
out_data(0,1)=dif_ddf_New_avg
out_data(0,2)=dif_fi_New_avg

out_data(1,0)=rr_ddt_New_avg
out_data(1,1)=rr_ddf_New_avg
out_data(1,2)=rr_fi_New_avg

out_data(2,0)=re_ddt_New_avg
out_data(2,1)=re_ddf_New_avg
out_data(2,2)=re_fi_New_avg


 opt = True
 opt@title  = "        DDT        DDF       FI    "
 opt@fout   = "Annual_25N.txt"
 write_matrix (out_data, "3f12.3", opt)


diri4="/Users/chang/Desktop/Small_work/permafrost/GLDAS/data/"
name4="GLDAS_FI_ori.nc"
 f4  = addfile(diri4+name4,"r")
  fi_GLDAS= f4->fi(:110,:,:)
a = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
  lsdata = a->LSMASK
  lsm = landsea_mask(lsdata,lat,lon)

  fi_GLDAS  = mask(fi_GLDAS ,lsm.eq.0,False)  
  fi_GLDAS  = mask(fi_GLDAS ,lsm.eq.2,False)  
  fi_GLDAS  = mask(fi_GLDAS ,lsm.eq.4,False)  
fi_GLDAS_ave=dim_avg_n_Wrap(fi_GLDAS,0)
fi_GLDAS_ave=where(fi_GLDAS_ave.le.0.495,fi_GLDAS_ave@_FillValue,fi_GLDAS_ave)

dif_fi_New=where(ismissing(fi_GLDAS_ave),fi_New@_FillValue,dif_fi_New)
dif_ddt_New=where(ismissing(fi_GLDAS_ave),ddt_New@_FillValue,dif_ddt_New)
dif_ddf_New=where(ismissing(fi_GLDAS_ave),ddf_New@_FillValue,dif_ddf_New)

rr_fi_New=where(ismissing(fi_GLDAS_ave),fi_New@_FillValue,rr_fi_New)
rr_ddt_New=where(ismissing(fi_GLDAS_ave),ddt_New@_FillValue,rr_ddt_New)
rr_ddf_New=where(ismissing(fi_GLDAS_ave),ddf_New@_FillValue,rr_ddf_New)

RE_fi_New=where(ismissing(fi_GLDAS_ave),fi_New@_FillValue,RE_fi_New)
RE_ddt_New=where(ismissing(fi_GLDAS_ave),ddt_New@_FillValue,RE_ddt_New)
RE_ddf_New=where(ismissing(fi_GLDAS_ave),ddf_New@_FillValue,RE_ddf_New)


;===============output data=======================

  rad    = 4.0*atan(1.0)/180.0
  re     = 6371220.0
  rr     = re*rad
  clat   = cos(lat*rad)

 dif_fi_New_avg = wgt_areaave(dif_fi_New, clat, 1.0, 0)
 dif_ddt_New_avg = wgt_areaave(dif_ddt_New, clat, 1.0, 0)
 dif_ddf_New_avg = wgt_areaave(dif_ddf_New, clat, 1.0, 0)

 rr_fi_New_avg = wgt_areaave(rr_fi_New, clat, 1.0, 0)
 rr_ddt_New_avg = wgt_areaave(rr_ddt_New, clat, 1.0, 0)
 rr_ddf_New_avg = wgt_areaave(rr_ddf_New, clat, 1.0, 0)

 re_fi_New_avg = wgt_areaave(RE_fi_New, clat, 1.0, 0)
 re_ddt_New_avg = wgt_areaave(RE_ddt_New, clat, 1.0, 0)
 re_ddf_New_avg = wgt_areaave(RE_ddf_New, clat, 1.0, 0)

out_data=new((/3,3/),"float")

out_data(0,0)=dif_ddt_New_avg
out_data(0,1)=dif_ddf_New_avg
out_data(0,2)=dif_fi_New_avg

out_data(1,0)=rr_ddt_New_avg
out_data(1,1)=rr_ddf_New_avg
out_data(1,2)=rr_fi_New_avg

out_data(2,0)=re_ddt_New_avg
out_data(2,1)=re_ddf_New_avg
out_data(2,2)=re_fi_New_avg

 opt = True
 opt@title  = "        DDT        DDF       FI    "
 opt@fout   = "Annual_permafrost.txt"
 write_matrix (out_data, "3f12.3", opt)


exit






