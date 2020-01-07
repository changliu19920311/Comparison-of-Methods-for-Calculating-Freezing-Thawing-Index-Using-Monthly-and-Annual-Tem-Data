
var_name=(/"dif_ddf_Feng","dif_ddt_Feng","dif_fi_Feng","RE_ddf_Feng","RE_ddt_Feng","RE_fi_Feng","rr_ddf_Feng","rr_ddt_Feng","rr_fi_Feng"/)
out_data=new((/3,3/),"float")




;=================read the input data =================================================
do nvar=0,dimsizes(var_name)-1
  f3  = addfile("merit_Feng.nc","r")
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
 opt@fout   = "Annual_permafrost_Feng.txt"
 write_matrix (out_data, "3f12.3", opt)
exit




out_data(0,0)=dif_ddt_Feng_avg
out_data(0,1)=dif_ddf_Feng_avg
out_data(0,2)=dif_fi_Feng_avg

out_data(1,0)=rr_ddt_Feng_avg
out_data(1,1)=rr_ddf_Feng_avg
out_data(1,2)=rr_fi_Feng_avg

out_data(2,0)=re_ddt_Feng_avg
out_data(2,1)=re_ddf_Feng_avg
out_data(2,2)=re_fi_Feng_avg

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
  


;-- shapefile mask resources
opt             =  True
opt@return_mask =  True    ;-- this forces the return of a 0s and 1s mask array

;-- create the mask based on the given shapefile
mask_array             =  shapefile_mask_data(var, shpname, opt)
mask_array!0           = "y" 
mask_array!1           = "x" 
mask_array@coordinates = "XLAT XLONG" 




exit


;-- create array to hold all polylines
  npoly   =  sum(geometry(:,geom_numSegs))               ;-- sum of all counties polygons
  poly    =  new(npoly,graphic)                          ;-- array of all counties polygons
  npl     =  0                                           ;-- polyline counter
  j       =  0     

  
  unemp = f->UNEMPLOY / f->PERSONS
  lon   = f->x
  lat   = f->y

  plres             = True       ; resources for polylines
  plres@gsEdgesOn   = True       ; draw border around polygons
  plres@gsEdgeColor = "black"    

  colors = (/"blue","green","yellow","red"/)

  segNum = 0
  do i=0, numFeatures-1  

     ; color assignment (probably a better way to do this?)
     if (unemp(i).ge.0.01 .and. unemp(i).lt.0.02) then
         plres@gsFillColor = colors(0)
     end if
     if (unemp(i).ge.0.02 .and. unemp(i).lt.0.03) then
         plres@gsFillColor = colors(1)
     end if
     if (unemp(i).ge.0.03 .and. unemp(i).lt.0.04) then
         plres@gsFillColor = colors(2)
     end if
     if (unemp(i).ge.0.04) then 
         plres@gsFillColor = colors(3)
     end if

     startSegment = geometry(i, geom_segIndex)
     numSegments  = geometry(i, geom_numSegs)
     do seg=startSegment, startSegment+numSegments-1
        startPT = segments(seg, segs_xyzIndex)
        endPT = startPT + segments(seg, segs_numPnts) - 1
        lines(segNum) = gsn_add_polygon(wks, plot, lon(startPT:endPT),  \
                                                   lat(startPT:endPT), plres)
        segNum = segNum + 1
     end do
  end do




CONTENT=f->CONTENT
print(CONTENT)
NUM_CODE=f->NUM_CODE
print(NUM_CODE)
EXTENT=f->EXTENT