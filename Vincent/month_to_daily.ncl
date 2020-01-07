;=== Program for computing the reconstructed daily temperature data

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"


;================set some parameters =============================================== 
out_diri ="/Volumes/Seagate Backup Plus Drive/GLDAS/GLDAS_1901_2012/GLDAS_Daily_Feng/"
out_diri1 ="/Volumes/Seagate\ Backup\ Plus\ Drive/GLDAS/GLDAS_1901_2012/GLDAS_Daily_Feng/"
int_diri="/Users/chang/Desktop/Small_work/permafrost/Feng/data/"
  ;=================read the latitude and longitude
   filename="/Volumes/Seagate Backup Plus Drive/GLDAS/GLDAS_1901_2012/monthly/tmin_GLDAS_1940.nc"
   f=addfile(filename,"r")
   lat=f->lat
   lon=f->lon
   Year=ispan(1901,2012,1)
   nlat=360
   nlon=720
   nyear=dimsizes(Year)
   Month=(/"01","02","03","04","05","06","07","08","09","10","11","12"/)
   days=(/1,31,59,90,120,151,181,212,243,273,304,334,365,396/)
   day=ispan(1,396,1)

  ;====================read the data 
  filename=int_diri+"GLDAS_M_mid.nc"
  f=addfile(filename,"r")
  Tanom=f->Tanom_mid
  Tanom!0="year"
  Tanom!1="month"
  Tanom!2="lat"
  Tanom!3="lon"

Tanom_new=new((/dimsizes(Year),14,nlat,nlon/),"float")
Tanom_new(0,0,:,:)=Tanom(0,11,:,:)
Tanom_new(1:111,0,:,:)=Tanom(0:110,11,:,:)
Tanom_new(0:110,13,:,:)=Tanom(1:111,0,:,:)
Tanom_new(111,13,:,:)=Tanom(111,0,:,:)
Tanom_new(:,1:12,:,:)=Tanom

  do nf=0,dimsizes(Year)-1
       print(nf)
       name="GLDAS_Feng_"+tostring(Year(nf))+".nc"
       system("/bin/rm "+out_diri1+name)
       setfileoption("nc","Format","LargeFile")
       fo  = addfile(out_diri+name,"c")
       data = linint1_n_Wrap (days, Tanom_new(nf,:,:,:), False, day, 0, 0)
       printVarSummary(data)
       Data = data(16:380,:,:)
       Data!0="time"
       Data&time=ispan(1,365,1)
      fo->Tanom=Data
  end do




