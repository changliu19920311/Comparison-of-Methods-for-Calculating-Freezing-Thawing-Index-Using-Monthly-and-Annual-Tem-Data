 f  = addfile("Elev_data.nc","r")
 data= f->output
 fi_elev=data(:,2,:) 
 ddt_elev=data(:,0,:) 
 ddf_elev=data(:,1,:) 
x=ispan(1,110,1)
ddt=new((/3,dimsizes(x)/),"float")
ddf=new((/3,dimsizes(x)/),"float")
fi=new((/3,dimsizes(x)/),"float")

ddt(0,0)=data(0,0,0)
wks   = gsn_open_wks ("pdf","xy")               ; send graphics to PNG file

  res                  = True                     ; plot mods desired
res@gsnDraw = False
res@gsnFrame= False
res@tmYLMinorOn         = False
res@tmYRMinorOn         = False

  res@tmXBMode          = "Explicit"              ; explicit labels
  res@tmXBValues        = ispan(1,9,1)                       ; location of labels 
  res@tmXBLabels        = (/"1","5","10","20","30","50","70","90","110"/)               ; labels themselves
  res@vpWidthF = 0.7
  res@vpHeightF = 0.3
  res@tmLabelAutoStride = True
  res@tmXTOn            = False                    ; nice stride on labels
  res@tmYLMode     =  "Manual"
  res1                        = res 
  res2                        = res  
  res3                        = res                      
  res2@xyLineThicknesses      = 2    
  res1@xyLineThicknesses      = 2.  
  res1@gsnLeftString="(a) RE (%)"

  res1@trYMinF     = 2.8		; set minimum Y-axis value
  res1@trYMaxF     =  4.4
  res2@trYMinF     = 5.0		; set minimum Y-axis value
  res2@trYMaxF     =  6.9
  res3@trYMinF     = 0.8		; set minimum Y-axis value
  res3@trYMaxF     = 1.8

   res1@tiYAxisString   ="FI"
   res2@tiYAxisString   ="DDF"
   res3@tiYAxisString   ="DDT"

  res1@xyLineColors      = "black"                   ; line color
  res3@xyLineColors      = "blue"                    ; line color
  res2@xyLineColors      = "red"                    ; line color

res1@tmYLTickStartF=res1@trYMinF
res2@tmYRTickStartF=res2@trYMinF
res3@tmYRTickStartF=res3@trYMinF
res1@tmYLPrecision  = 2
res2@tmYRPrecision  = 2
res3@tmYRPrecision  = 2
res1@tmYLTickSpacingF  = 0.4
res2@tmYRTickSpacingF  = 0.5
res3@tmYRTickSpacingF  = 0.2
 plot0  = gsn_csm_xy3 (wks,ispan(1,9,1),data(0,2,:),data(0,1,:),data(0,0,:),res1,res2,res3) ; create plot
  res1@gsnLeftString="(b) MBE"


 res1@trYMinF     = 0.0005		; set minimum Y-axis value
  res1@trYMaxF     = 0.003
  res2@trYMinF     = 0		; set minimum Y-axis value
  res2@trYMaxF     =  10
  res3@trYMinF     = 0		; set minimum Y-axis value
  res3@trYMaxF     = 10

res1@tmYLTickStartF=res1@trYMinF

res2@tmYRTickStartF=res2@trYMinF
res3@tmYRTickStartF=res3@trYMinF
res1@tmYLPrecision  = 2
res2@tmYRPrecision  = 2
res3@tmYRPrecision  = 2
res1@tmYLTickSpacingF  = 0.0005
res2@tmYRTickSpacingF  = 0.5
res3@tmYRTickSpacingF  = 0.2
 plot1  = gsn_csm_xy3 (wks,ispan(1,9,1),data(2,2,:),data(2,1,:),data(2,0,:),res1,res2,res3) ; create plot

  res1@gsnLeftString="(c) Correlation"
  res1@trYMinF     = 0.943		; set minimum Y-axis value
  res1@trYMaxF     = 0.955
  res2@trYMinF     = 0.952	; set minimum Y-axis value
  res2@trYMaxF     = 0.962
  res3@trYMinF     = 0.9914	; set minimum Y-axis value
  res3@trYMaxF     = 0.9928

res1@tmYLTickStartF=res1@trYMinF
res2@tmYRTickStartF=res2@trYMinF
res3@tmYRTickStartF=res3@trYMinF
res1@tmYLPrecision  = 3
res2@tmYRPrecision  = 3
res3@tmYRPrecision  = 4
res1@tmYLTickSpacingF  = 0.002
res1@tiXAxisString ="Years"
 plot2  = gsn_csm_xy3 (wks,ispan(1,9,1),data(1,2,:),data(1,1,:),data(1,0,:),res1,res2,res3) ; create plot


  lgres                    = True
  lgres@lgLineColors       = (/"blue","red","black"/)
  lgres@lgItemType         = "Lines"        ; show lines only (default)
  lgres@lgLabelFontHeightF = .08            ; legend label font thickness
  lgres@vpWidthF           = 0.13           ; width of legend (NDC)
  lgres@vpHeightF          = 0.10           ; height of legend (NDC)
  lgres@lgPerimOn = False
  lgres@lgPerimThicknessF  = 2.0            ; thicken the box perimeter
  lgres@lgMonoDashIndex    = True
  lgres@lgDashIndex        = 0

  labels = (/"DDT","DDF","FI"/)
  nitems = dimsizes(labels)                 ; number of legend items

; Create legend
  lbid0 = gsn_create_legend(wks,nitems,labels,lgres)
 ; lbid1 = gsn_create_legend(wks,nitems,labels,lgres)
 ; lbid2 = gsn_create_legend(wks,nitems,labels,lgres)
;
; Use gsn_add_annotation to attach this legend to our existing plot.
; This way, if we resize the plot, the legend will stay with the
; plot and be resized automatically.
;
  amres                  = True
  amres@amJust           = "BottomRight"    ; Use bottom right corner of box
                                            ; for determining its location.
  amres@amParallelPosF   = 0.5              ; Move legend to right
  amres@amOrthogonalPosF = -0.15              ; Move legend down.

  annoid0 = gsn_add_annotation(plot0,lbid0,amres)  ; add legend to plot
  annoid1 = gsn_add_annotation(plot1,lbid1,amres)  ; add legend to plot
  annoid2 = gsn_add_annotation(plot2,lbid2,amres)  ; add legend to plot
  resP                  = True                   ; modify the panel plot
resP@gsnPanelBottom = 0.1
  gsn_panel(wks,(/plot0,plot1,plot2/),(/3,1/),resP) 




