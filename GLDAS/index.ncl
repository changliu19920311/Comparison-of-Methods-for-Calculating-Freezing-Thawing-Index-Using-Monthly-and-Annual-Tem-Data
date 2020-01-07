; Program for computing the frost index
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

diri1 = "/Users/chang/Desktop/Small_work/permafrost/GLDAS/data/"
int_diri="/Volumes/Seagate Backup Plus Drive/GLDAS/GLDAS_1901_2012/Daily/"
Int_diri="/Volumes/Seagate\ Backup\ Plus\ Drive/GLDAS/GLDAS_1901_2012/Daily/"


nlat=360
nlon=720
Year=ispan(1901,2012,1)
nyear=dimsizes(Year)
FI_ori=new((/nyear,nlat,nlon/),"float")
ddt=new((/nyear,nlat,nlon/),"float")
ddf=new((/nyear,nlat,nlon/),"float")


  fils1 = systemfunc ("ls "+Int_diri+"tmin_daily_*.nc") ; file paths
  f1    = addfiles (fils1, "r")   
  ListSetType (f1, "cat")
  fils2 = systemfunc ("ls "+Int_diri+"tmax_daily_*.nc") ; file paths
  f2    = addfiles (fils2, "r")   
  ListSetType (f2, "cat")

do nf=0,dimsizes(Year)-2
print(nf)
;===========calculate the FI for the original data======================


  tmin_1 = f1[nf]->tmin
  tmin_2 = f1[nf+1]->tmin
  tmax_1 = f2[nf]->tmax
  tmax_2 = f2[nf+1]->tmax
  tmin = array_append_record (tmin_1, tmin_2, 0)
  tmax = array_append_record (tmax_1, tmax_2, 0)
  tmp=(tmin+tmax)/2
  tmp =tmp - 273.15
  copy_VarCoords(tmin,tmp)

      tmp_ddt=where(tmp.le.0,0,tmp)
      ddt(nf,:,:)= dim_sum_n(tmp_ddt(:364,:,:),0) 
      tmp_ddf=where(tmp.gt.0,0,tmp)
      ddf(nf,:,:)= -1.0*dim_sum_n(tmp_ddf(181:545,:,:),0) ;

      FI_ori(nf,:,:)=sqrt(ddf(nf,:,:))/(sqrt(ddf(nf,:,:))+sqrt(ddt(nf,:,:)))

      copy_VarCoords(tmp(0,:,:),FI_ori(0,:,:))
      copy_VarCoords(tmp(0,:,:),ddt(0,:,:)) 
      copy_VarCoords(tmp(0,:,:),ddf(0,:,:)) 

delete(tmin)
delete(tmax)
delete(tmin_1)
delete(tmin_2)
delete(tmax_1)
delete(tmax_2)
delete(tmp)
delete(tmp_ddt)
delete(tmp_ddf)
end do
 

  FI_ori!0="year"
  FI_ori&year=Year
  printVarSummary(FI_ori)
  FI_ori@long_name = "Frost index for original daily GLDAS data "


  ddt!0="year"
  ddt&year=Year
  printVarSummary(ddt)
  ddt@long_name =  "Thawing index for original daily GLDAS data "

  ddf!0="year"
  ddf&year=Year
  printVarSummary(ddf)
  ddf@long_name =  "Freezing index for original daily GLDAS data "


    
   system("/bin/rm "+diri1+"GLDAS_FI_ori.nc")
   fo      = addfile(diri1+"GLDAS_FI_ori.nc","c")          ; create output file
   filedimdef(fo,"time",-1,True)
   fo->fi  = FI_ori
   fo->ddt  = ddt
   fo->ddf  = ddf
  





