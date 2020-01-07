
var_name=(/"dif_ddf_New","dif_ddt_New","dif_fi_New","RE_ddf_New","RE_ddt_New","RE_fi_New","rr_ddf_New","rr_ddt_New","rr_fi_New"/)
out_data=new((/3,3/),"float")




;=================read the input data =================================================
do nvar=0,dimsizes(var_name)-1
  f3  = addfile("merit_New.nc","r")
  data :=  f3->$var_name(nvar)$
  data :=lonFlip(data)
  lat=data&lat
  lon=data&lon

  dims = dimsizes(data)
  lat1d := ndtooned(conform_dims(dims,lat,0))
  lon1d := ndtooned(conform_dims(dims,lon,1))
  nlatlon1d = dimsizes(lat1d)
  data_mask_1d = new(nlatlon1d,integer)

  rad    = 4.0*atan(1.0)/180.0
  re     = 6371220.0
  rr     = re*rad
  clat   = cos(lat*rad)

delete(lat)
delete(lon)


;====================mask the value using the shapefile ======================

;=================read the shapefile for permafrost=====================================
 f=addfile("permaice.shp","r")
 print(f)
  segments      =  f->segments
  geometry      =  f->geometry
  segsDims = dimsizes(segments)
  geomDims      =  dimsizes(geometry)
  lon           =  f->x                                   ;-- longitudes array of counties
  lat           =  f->y             
  geom_segIndex =  f@geom_segIndex
  geom_numSegs  =  f@geom_numSegs
  segs_xyzIndex =  f@segs_xyzIndex
  segs_numPnts  =  f@segs_numPnts
  numFeatures   =  geomDims(0)

  lines       = new(segsDims(0),graphic)   ; Array to hold polygons
  numFeatures = geomDims(0)

  unemp=f->EXTENT
  segNum = 0

    do nf=0,numFeatures-1
      if (unemp(nf).ne."") then
        startSegment = geometry(nf, geom_segIndex)
        numSegments  = geometry(nf, geom_numSegs)
        do seg=startSegment, startSegment+numSegments-1
          startPT = segments(seg, segs_xyzIndex)
          endPT   = startPT + segments(seg, segs_numPnts) - 1
          lat_sub := tofloat(lat(startPT:endPT))
          lon_sub := tofloat(lon(startPT:endPT))
          ii_latlon := ind(min(lat_sub).le.lat1d.and.lat1d.le.max(lat_sub).and.\
                           min(lon_sub).le.lon1d.and.lon1d.le.max(lon_sub))
          if(any(ismissing(ii_latlon))) then
            continue
          end if          
          nii = dimsizes(ii_latlon)
          do n=0,nii-1
            iltln = ii_latlon(n)
            if(gc_inout(lat1d(iltln),lon1d(iltln),lat_sub,lon_sub)) then
              data_mask_1d(iltln) =1
            end if
          end do
        end do
      end if
    end do

data_mask = where(onedtond(data_mask_1d,dims).eq.1,data,data@_FillValue)
printVarSummary(data_mask)
copy_VarCoords(data,data_mask)



;===============output data=======================

 data_avg = wgt_areaave(data_mask, clat, 1.0, 0)
 ny=mod(nvar,3)
 nx =toint(nvar/3)
print(nx)
print(ny)
out_data(nx,ny)=data_avg
delete(lat)
delete(lon)
delete(data)
end do 
 opt = True
 opt@title  = "        DDF        DDT       FI    "
 opt@fout   = "Annual_permafrost_New.txt"
 write_matrix (out_data, "3f12.3", opt)
exit




out_data(0,0)=dif_ddt_New_avg
out_data(0,1)=dif_ddf_New_avg
out_data(0,2)=dif_fi_New_avg

out_data(1,0)=rr_ddt_New_avg
out_data(1,1)=rr_ddf_New_avg
out_data(1,2)=rr_fi_New_avg

out_data(2,0)=re_ddt_New_avg
out_data(2,1)=re_ddf_New_avg
out_data(2,2)=re_fi_New_avg

 wks  = gsn_open_wks("pdf","test") ; send graphics to PNG file
  gsn_define_colormap(wks,"amwg_blueyellowred")
  res                     = True
  res@gsnPolar   = "NH"  
 res@cnFillOn              = True     ; turn on color fill
  res@cnLinesOn             = False    ; turn of contour lines
  res@mpMinLatF            = 25 
res@mpFillOn             = False   
res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnFillColors        =(/2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17/)
res@cnLevels := (/0.1,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.9/)
  res@pmTickMarkDisplayMode = "Always"  ; turn on tickmarks
  res@tiMainString          = "Percentage unemployment, by state"
  plot = gsn_csm_contour_map_polar(wks,data_mask,res)   ; Create map, but don't draw it yet.
delete(wks)

exit
  

