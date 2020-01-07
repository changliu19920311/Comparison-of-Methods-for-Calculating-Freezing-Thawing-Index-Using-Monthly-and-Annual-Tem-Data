; Program for computing the frost index
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

diri2 = "/Volumes/Seagate Backup Plus Drive/GLDAS/GLDAS_1901_2012/GLDAS_Daily_Feng/"
diri3 = "/Users/chang/Desktop/Small_work/permafrost/Feng/data/"
nlat=360
nlon=720
Year=ispan(1901,2012,1)
nyear=dimsizes(Year)
Month=(/"01","02","03","04","05","06","07","08","09","10","11","12"/)

FI_rec=new((/nyear,nlat,nlon/),"float")
ddt=new((/nyear,nlat,nlon/),"float")
ddf=new((/nyear,nlat,nlon/),"float")
;===========ca;culate the FI for the reconstructed data====================== 
do nf=0,dimsizes(Year)-2
print(nf)
      name="GLDAS_Feng_"+tostring(Year(nf))+".nc"
      f  = addfile(diri2+name,"r")
      tmp_1 = f->Tanom
      delete(f)

      name="GLDAS_Feng_"+tostring(Year(nf+1))+".nc"
      f  = addfile(diri2+name,"r")
      tmp_2 = f->Tanom

      tmp = array_append_record (tmp_1, tmp_2, 0)
 
         tmp_ddt=where(tmp.le.0,0,tmp)
      ddt(nf,:,:)= dim_sum_n(tmp_ddt(:364,:,:),0) 
      tmp_ddf=where(tmp.gt.0,0,tmp)
      ddf(nf,:,:)= -1.0*dim_sum_n(tmp_ddf(181:545,:,:),0) ;

     
              FI_rec(nf,:,:)=sqrt(ddf(nf,:,:))/(sqrt(ddf(nf,:,:))+sqrt(ddt(nf,:,:)))
      copy_VarCoords(tmp(0,:,:),FI_rec(0,:,:))  
      copy_VarCoords(tmp(0,:,:),ddt(0,:,:)) 
      copy_VarCoords(tmp(0,:,:),ddf(0,:,:)) 
      delete(tmp)
 delete(tmp_1)
 delete(tmp_2)


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



   system("/bin/rm "+diri3+"GLDAS_FI_Feng.nc")
   fo      = addfile(diri3+"GLDAS_FI_Feng.nc","c")          ; create output file
   filedimdef(fo,"time",-1,True)
   fo->fi  = FI_rec
   fo->ddt  = ddt
   fo->ddf  = ddf


  





