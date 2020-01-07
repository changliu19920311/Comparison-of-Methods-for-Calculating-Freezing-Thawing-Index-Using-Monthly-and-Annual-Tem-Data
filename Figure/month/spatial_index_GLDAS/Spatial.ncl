; evaluating the reconstucted monthly data=================
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;================set the parameters=====================



Year1=1901
Year2=2011
Year=ispan(Year1,Year2,1)




;=================read the reconstructed data ==========================

  f1  = addfile("../../../GLDAS/GLDAS_index.nc","r")
  lat=f1->lat
  lon=f1->lon
  fi= f1->fi(:110,:,:)
  ddt= f1->ddt(:110,:,:)
  ddf= f1->ddf(:110,:,:)

 a = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
  lsdata = a->LSMASK
  lsm = landsea_mask(lsdata,lat,lon)

  fi  = mask(fi ,lsm.eq.0,False)  
  fi  = mask(fi ,lsm.eq.2,False)  
  fi  = mask(fi ,lsm.eq.4,False)  


ilat=ind(lat.le.25)
fi(:,ilat,:)= fi@_FillValue
ddt = where(ismissing(fi),fi@_FillValue,ddt )
ddf = where(ismissing(fi),fi@_FillValue,ddf )
delete(ilat)

 rad    = 4.0*atan(1.0)/180.0
  re     = 6371220.0
  rr     = re*rad
  clat   = cos(lat*rad)



  fi_ave_area  = wgt_areaave(fi, clat, 1.0, 0)
  ddf_ave_area  = wgt_areaave(ddf, clat, 1.0, 0)
  ddt_ave_area  = wgt_areaave(ddt, clat, 1.0, 0)


fi_ave = dim_avg_n_Wrap(fi, 0)
ddt_ave = dim_avg_n_Wrap(ddt, 0)
ddf_ave = dim_avg_n_Wrap(ddf, 0)




  wks = gsn_open_wks("pdf" ,"Spatial_distribution")      ; send graphics to PNG file
  gsn_define_colormap(wks,"amwg_blueyellowred")
  plot=new(4,graphic)
  res = True     
  res@gsnDraw = False
  res@gsnFrame = False 
  res@mpCenterLonF=-45
  res@gsnPolar   = "NH"  
  ;res@cnFillMode = "RasterFill"
  res@gsnLeftString   = " "  
  res@mpMinLatF            = 25     ; range to zoom in on
  res@mpFillOn             = False              
  res@cnFillOn              = True     ; turn on color fill
  res@cnLinesOn             = False    ; turn of contour lines  
  res@lbLabelBarOn = False
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels = (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9/)
  res@cnFillColors        =(/5,6,7,8,9,10,11,12,13,14/)
  res@lbLabelStride         = 1 
  res@gsnPolarLabelFontHeightF = .005 
res@gsnLeftStringFontHeightF=0.008
 res@vpXF             = 0.08                  ; position and sizes
 res@vpYF             = 0.35                  ; for XY plot
 res@vpWidthF         = 0.2
 res@vpHeightF        = 0.2
  res@gsnLeftString  = "(c) FI" 
  plot_FI= gsn_csm_contour_map_polar(wks,fi_ave,res)
delete(res@cnFillColors)
  res@cnLevels := (/100,300,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000/) 

  res@gsnLeftString  = "(a) DDT"
 res@vpXF             = 0.08                  ; position and sizes
 res@vpYF             = 0.95                  ; for XY plot
 res@vpWidthF         = 0.2
 res@vpHeightF        = 0.2
  res@cnFillColors        =(/2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17/)
  plot_ddt=gsn_csm_contour_map_polar(wks,ddt_ave,res)
  res@cnLevels := (/100,300,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000/)  
  res@cnFillColors        =(/2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17/)
  res@gsnLeftString  = "(b) DDF" 
 res@vpXF             = 0.08                  ; position and sizes
 res@vpYF             = 0.65                 ; for XY plot
 res@vpWidthF         = 0.2
 res@vpHeightF        = 0.2
  plot_ddf=gsn_csm_contour_map_polar(wks,ddf_ave,res)
 ; res@cnFillColors :=(/-1,6,4,2/)
 ; res@cnLevels := (/0.495,0.570,0.625/)  
  ;res@gsnLeftString  = "(d) Permafrost zones" 
  ;res@lbLabelStrings    = (/"No permafrost","Sporadic","Discontinuous","Continuous"/)
  ;res@lbLabelAlignment = "BoxCenters"
  ;plot_FI_2=gsn_csm_contour_map_polar(wks,fi_ave,res)

delete(res)

res                  = True                     ; plot mods desired
res@gsnDraw=False
res@gsnFrame = False
res@vpHeightF          = 0.15
res@vpWidthF=  0.4
res@tmYROn            = False
res@tmXTOn            = False
res@xyMarkerSizeF    = 0.006
res@xyDashPatterns  = 0
res@xyMarkLineModes = "Lines"
res@xyMarkers      =   16 
res@trXMaxF   =  2015
res@xyLineColors="black"
res@trXMinF   =  1900
res@tmXBMode     =  "Manual"  
res@tmXBLabelDeltaF =-1
res@tmYLLabelDeltaF =-1
res@tmXBMinorOn         = False
res@xyLineThicknesses  =  1
res@xyDashPattern     =  "solid"
res@tmYLMajorOutwardLengthF = -0.0001
res@tmXBMajorOutwardLengthF =-0.0001
res@tmYLMinorOn         = False
res@tmBorderThicknessF  = 0.005
res@tmYLMajorThicknessF = 0.2
res@tmXBMajorThicknessF = 0.2
res@tmYLMajorLengthF  = 0.005
res@tmXBMajorLengthF  = 0.005

res@tiXAxisString = "Year"; xaxis string
res@tiYAxisString = "DDF"               ; yaxis string
res@gsnLeftString = ""
 res@vpXF             = 0.4                  ; position and sizes
 res@vpYF             = 0.6  
  res@gsnLeftString  = "(d)"
 plot_ddt_0 = gsn_csm_xy (wks,Year,ddf_ave_area,res)
res@tiYAxisString = "DDT"               ; yaxis string
 res@vpXF             = 0.4                  ; position and sizes
 res@vpYF             = 0.9 
  res@gsnLeftString  = "(e)"
 plot_ddf_0 = gsn_csm_xy (wks,Year,ddt_ave_area,res)
res@tiYAxisString = "FI"               ; yaxis string
 res@vpXF             = 0.4                  ; position and sizes
 res@vpYF             = 0.3  
  res@gsnLeftString  = "(f)"
 plot_fi_0 = gsn_csm_xy (wks,Year,fi_ave_area,res)

; psres = True                                                               
; maximize_output(wks,psres)
;frame(wks)

  res1P                     = True         ; modify the panel plot
  res1P@gsnFrame         = False
  res1P@gsnPanelTop   = 0.87  
  res1P@gsnPanelYWhiteSpacePercent=5
  res1P@gsnPanelXWhiteSpacePercent=5
  res1P@gsnPanelBottom    = 0.5
  res1P@gsnPanelRight =0.89
  res1P@gsnPanelLeft =0.11
  res1P@gsnPanelLabelBar    = True       ; common label bar
  res1P@pmLabelBarWidthF    = 0.5      ; default is shorter
  res1P@pmLabelBarHeightF   = 0.05        ; default is taller
  res1P@lbLabelFontHeightF  = 0.01       ; make labels smaller ( default=0.02 )
  gsn_panel(wks,(/plot_ddt,plot_ddf/),(/1,2/),res1P) 





  res2P                     = True         ; modify the panel plot
  res2P@gsnFrame         = False
  res2P@gsnPanelTop   = 0.5
  res2P@gsnPanelBottom    = 0.13
  res2P@gsnPanelRight =0.89
  res2P@gsnPanelLeft =0.11
  res2P@gsnPanelYWhiteSpacePercent=5
  res2P@gsnPanelXWhiteSpacePercent=5
  res2P@gsnPanelLabelBar    = True       ; common label bar
  res2P@pmLabelBarWidthF    = 0.5      ; default is shorter
  res2P@pmLabelBarHeightF   = 0.05        ; default is taller
  res2P@lbLabelFontHeightF  = 0.01       ; make labels smaller ( default=0.02 )
  gsn_panel(wks,(/plot_FI/),(/1,1/),res2P) 

  ;gsn_panel(wks,(/plot_ddt,plot_ddt_0,plot_ddf,plot_ddf_0,plot_FI,plot_fi_0 /),(/3,2/),resP) 
  
  
exit


