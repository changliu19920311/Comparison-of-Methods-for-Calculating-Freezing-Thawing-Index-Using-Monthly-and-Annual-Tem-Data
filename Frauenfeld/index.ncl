; Program for computing the frost index
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

diri1 = "/Volumes/Seagate\ Backup\ Plus\ Drive/GLDAS/GLDAS_1901_2012/monthly/"
diri3 = "/Users/chang/Desktop/Small_work/permafrost/Oliver/data/"


nlat=360
nlon=720
Year=ispan(1901,2012,1)
nyear=dimsizes(Year)
days=new((/24,360,720/),"float")
      days(0,:,:)=31
      days(1,:,:)=28
      days(2,:,:)=31
      days(3,:,:)=30
      days(4,:,:)=31      
      days(5,:,:)=30
      days(6,:,:)=31      
      days(7,:,:)=31
      days(8,:,:)=30
      days(9,:,:)=31
      days(10,:,:)=30
      days(11,:,:)=31
      days(12,:,:)=31
      days(13,:,:)=28
      days(14,:,:)=31
      days(15,:,:)=30
      days(16,:,:)=31      
      days(17,:,:)=30
      days(18,:,:)=31      
      days(19,:,:)=31
      days(20,:,:)=30
      days(21,:,:)=31
      days(22,:,:)=30
      days(23,:,:)=31
ddt=new((/nyear,nlat,nlon/),"float")
ddf=new((/nyear,nlat,nlon/),"float")

FI_rec=new((/nyear,nlat,nlon/),"float")
;===========calculate the FI for the reconstructed data======================

fils = systemfunc ("ls "+diri1+"Tmp_GLDAS_*.nc") ; file paths
  f    = addfiles (fils, "r")   
  ListSetType (f, "cat")
  t=f[:]->T
  time=f[:]->time
  utc_date = cd_calendar(time, 0)
   year   = tointeger(utc_date(:,0))  
do nf=0,dimsizes(Year)-2
print(nf)
      iyyyy=ind(year.eq.Year(nf).or.year.eq.Year(nf+1))
      tmp = t(iyyyy,:,:)
      tmp_ddt=where(tmp.le.0,0,tmp)
      Tmp_ddt=tmp_ddt*days
      ddt(nf,:,:)= dim_sum_n(Tmp_ddt(:11,:,:),0) 
      tmp_ddf=where(tmp.gt.0,0,tmp)
      Tmp_ddf=tmp_ddf*days
      ddf(nf,:,:)= -1.0*dim_sum_n(Tmp_ddf(6:17,:,:),0) 
      FI_rec(nf,:,:)=sqrt(ddf(nf,:,:))/(sqrt(ddf(nf,:,:))+sqrt(ddt(nf,:,:)))

      copy_VarCoords(tmp(0,:,:),FI_rec(0,:,:))  
      copy_VarCoords(tmp(0,:,:),ddt(0,:,:)) 
      copy_VarCoords(tmp(0,:,:),ddf(0,:,:)) 
      delete(tmp)
end do

;===========================================================================

  FI_rec!0="year"
  FI_rec&year=Year
  printVarSummary(FI_rec)
  FI_rec@long_name = "Frost index for reconstructed daily GLDAS data "


  ddt!0="year"
  ddt&year=Year
  printVarSummary(ddt)
  ddt@long_name =  "Thawing index for reconstructed daily GLDAS data "

  ddf!0="year"
  ddf&year=Year
  printVarSummary(ddf)
  ddf@long_name =  "Freezing index for reconstructed daily GLDAS data "



   system("/bin/rm "+diri3+"GLDAS_FI_Oliver.nc")
   fo      = addfile(diri3+"GLDAS_FI_Oliver.nc","c")          ; create output file
   filedimdef(fo,"time",-1,True)
   fo->fi  = FI_rec
   fo->ddt  = ddt
   fo->ddf  = ddf


exit  





