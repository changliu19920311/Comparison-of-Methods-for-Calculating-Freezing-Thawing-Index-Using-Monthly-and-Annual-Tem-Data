; evaluating the reconstucted monthly data=================
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;================set the parameters=====================

Year1=1901
Year2=2011
Year=ispan(Year1,Year2,1)

;===========GLDAS daily=====================================
name1="GLDAS_FI_ori.nc"
diri1="/Users/chang/Desktop/Small_work/permafrost/GLDAS/data/"

;===========New method (monthly)================
name2 = "GLDAS_FI_New_ensemble.nc"
diri2="/Users/chang/Desktop/Small_work/permafrost/New_Method/monthly/data/Frost_index_data/"

;===========New method (annual)================
name3="GLDAS_FI_New_ensemble.nc"
diri3="/Users/chang/Desktop/Small_work/permafrost/New_Method/annual/data/Frost_index_data/"

;===========Vincent monthly===================================
name4="GLDAS_FI_Feng.nc"
diri4="/Users/chang/Desktop/Small_work/permafrost/Feng/data/"

;===========Frauenfeld method ==================================
name5 = "GLDAS_FI_Oliver.nc"
diri5="/Users/chang/Desktop/Small_work/permafrost/Oliver/data/"

;=================read data ==========================

  f1  = addfile(diri1+name1,"r")
  lat=f1->lat
  lon=f1->lon
  fi_1= f1->fi
  ddt_1= f1->ddt
  ddf_1= f1->ddf


  f1  = addfile(diri2+name2,"r")
  lat=f1->lat
  lon=f1->lon
  fi_2= f1->fi
  ddt_2= f1->ddt
  ddf_2= f1->ddf


  f1  = addfile(diri3+name3,"r")
  lat=f1->lat
  lon=f1->lon
  fi_3= f1->fi
  ddt_3= f1->ddt
  ddf_3= f1->ddf


  f1  = addfile(diri4+name4,"r")
  lat=f1->lat
  lon=f1->lon
  fi_4= f1->fi
  ddt_4= f1->ddt
  ddf_4= f1->ddf


  f1  = addfile(diri5+name5,"r")
  lat=f1->lat
  lon=f1->lon
  fi_5= f1->fi
  ddt_5= f1->ddt
  ddf_5= f1->ddf
;============combined to one array=============
fi=new((/5,dimsizes(Year),dimsizes(lat),dimsizes(lon)/),"float")
ddt=new((/5,dimsizes(Year),dimsizes(lat),dimsizes(lon)/),"float")
ddf=new((/5,dimsizes(Year),dimsizes(lat),dimsizes(lon)/),"float")
fi(0,:,:,:)=fi_1(:110,:,:)
fi(1,:,:,:)=fi_2(:110,:,:)
fi(2,:,:,:)=fi_3(:110,:,:)
fi(3,:,:,:)=fi_4(:110,:,:)
fi(4,:,:,:)=fi_5(:110,:,:)

ddt(0,:,:,:)=ddt_1(:110,:,:)
ddt(1,:,:,:)=ddt_2(:110,:,:)
ddt(2,:,:,:)=ddt_3(:110,:,:)
ddt(3,:,:,:)=ddt_4(:110,:,:)
ddt(4,:,:,:)=ddt_5(:110,:,:)

ddf(0,:,:,:)=ddf_1(:110,:,:)
ddf(1,:,:,:)=ddf_2(:110,:,:)
ddf(2,:,:,:)=ddf_3(:110,:,:)
ddf(3,:,:,:)=ddf_4(:110,:,:)
ddf(4,:,:,:)=ddf_5(:110,:,:)

;==============mask. data=================

  a = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
  lsdata = a->LSMASK
  lsm = landsea_mask(lsdata,lat,lon)

  fi  = mask(fi ,lsm.eq.0,False)  
  fi  = mask(fi ,lsm.eq.2,False)  
  fi  = mask(fi ,lsm.eq.4,False)  


ilat=ind(lat.le.25)
fi(:,:,ilat,:)= fi@_FillValue
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

 nave = 9
    opt  = 0
   fi_ave_area_s     = runave_n (fi_ave_area, nave, opt, 1)
  ddf_ave_area_s   = runave_n (ddf_ave_area, nave, opt, 1)
  ddt_ave_area_s   = runave_n (ddt_ave_area, nave, opt, 1)

  wks = gsn_open_wks("pdf" ,"Time_series")      ; send graphics to PNG file
  gsn_define_colormap(wks,"amwg_blueyellowred")
  plot=new(4,graphic)


res                  = True                     ; plot mods desired
res@gsnDraw=False
res@gsnFrame = False
res@vpHeightF          = 0.15
res@vpWidthF=  0.4
res@tmYROn            = False
res@tmXTOn            = False
res@xyMarkerSizeF    = 0.006
res@xyMarkLineModes = "Lines"
res@xyMarkers      =   16 
res@trXMaxF   =  2015
res@trXMinF   =  1900
res@tmXBMode     =  "Manual"  
res@tmXBLabelDeltaF =-1
res@tmYLLabelDeltaF =-1
res@tmXBMinorOn         = False
res@xyLineThicknesses  =  (/2,2,2,2,2,2/)
res@xyDashPattern     =  "solid"
res@tmYLMajorOutwardLengthF = -0.0001
res@tmXBMajorOutwardLengthF =-0.0001
res@tmYLMinorOn         = False
res@tmBorderThicknessF  = 0.005
res@tmYLMajorThicknessF = 0.2
res@tmXBMajorThicknessF = 0.2
res@tmYLMajorLengthF  = 0.005
res@tmXBMajorLengthF  = 0.005
res@xyLineColors    = (/"grey","black","red","blue","green","orange"/) 
res@xyDashPatterns  = (/0,0,0,0,0,0/)
res@xyLineOpacities = (/1.0,1.0,0.7,1.0,1.0,1.0/)
res@tiXAxisString = "Year"; xaxis string
res@tiYAxisString = "DDF"               ; yaxis string
res@gsnLeftString = ""
 res@vpXF             = 0.4                  ; position and sizes
 res@vpYF             = 0.6  
 res@gsnLeftString  = "(b)"
 plot_ddf_0 = gsn_csm_xy (wks,Year,(/ddf_ave_area(0,:),ddf_ave_area_s(0,:),ddf_ave_area_s(1,:),ddf_ave_area_s(2,:),ddf_ave_area_s(3,:),ddf_ave_area_s(4,:)/),res)
res@tiYAxisString = "FI"               ; yaxis string
 res@vpXF             = 0.4                  ; position and sizes
 res@vpYF             = 0.3  
  res@gsnLeftString  = "(c)"
 plot_fi_0 = gsn_csm_xy (wks,Year,(/fi_ave_area(0,:),fi_ave_area_s(0,:),fi_ave_area_s(1,:),fi_ave_area_s(2,:),fi_ave_area_s(3,:),fi_ave_area_s(4,:)/),res)

 res@pmLegendDisplayMode    = "Always"            ; turn on legend
  res@pmLegendSide           = "Top"               ; Change location of 
  res@pmLegendParallelPosF   = .14                  ; move units right
  res@pmLegendOrthogonalPosF = -0.7                ; move units down
  res@pmLegendWidthF         = 0.03                ; Change width and
  res@pmLegendHeightF        = 0.08                ; height of legend.
  res@lgPerimOn              = False               ; turn off box around
  res@lgLabelFontHeightF     = .005                 ; label font height
  res@xyExplicitLegendLabels = (/"Daily dataset","Daily dataset","Monte Carlo (monthly)","Monte Carlo (annual)","Vincent","Frauenfeld"/)         ; create explicit labels


  res@tiYAxisString = "DDT"               ; yaxis string
 res@vpXF             = 0.4                  ; position and sizes
 res@vpYF             = 0.9 
  res@gsnLeftString  = "(a)"
 plot_ddt_0 = gsn_csm_xy (wks,Year,(/ddt_ave_area(0,:),ddt_ave_area_s(0,:),ddt_ave_area_s(1,:),ddt_ave_area_s(2,:),ddt_ave_area_s(3,:),ddt_ave_area_s(4,:)/),res)



  resP                     = True         ; modify the panel plot

  gsn_panel(wks,(/plot_ddt_0,plot_ddf_0,plot_fi_0 /),(/3,1/),resP) 
  
  
exit


