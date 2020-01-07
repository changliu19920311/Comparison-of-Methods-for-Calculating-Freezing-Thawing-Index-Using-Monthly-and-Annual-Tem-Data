;=== Program for computing the reconstructed daily temperature data

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;==============create the output file=============================

;================the method to convert month data to daily data 

out_diri="/Users/chang/Desktop/Small_work/permafrost/test/Feng/data/"
out_diri="/Users/chang/Desktop/Small_work/permafrost/test/Feng/data/"
int_diri="/Volumes/ChangLiu/GLDAS_monthly_data/"

name="GLDAS_M_mid.nc"
system("/bin/rm "+out_diri+name)
setfileoption("nc","Format","LargeFile")
fo  = addfile(out_diri+name,"c")

A=(/(/7./8.,1./8.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./),\
(/1./8.,6./8.,1./8.,0.,0.,0.,0.,0.,0.,0.,0.,0./),\
(/0.,1./8.,6./8.,1./8.,0.,0.,0.,0.,0.,0.,0.,0./),\
(/0.,0.,1./8.,6./8.,1./8.,0.,0.,0.,0.,0.,0.,0./),\
(/0.,0.,0.,1./8.,6./8.,1./8.,0.,0.,0.,0.,0.,0./),\
(/0.,0.,0.,0.,1./8.,6./8.,1./8.,0.,0.,0.,0.,0./),\
(/0.,0.,0.,0.,0.,1./8.,6./8.,1./8.,0.,0.,0.,0./),\
(/0.,0.,0.,0.,0.,0.,1./8.,6./8.,1./8.,0.,0.,0./),\
(/0.,0.,0.,0.,0.,0.,0.,1./8.,6./8.,1./8.,0.,0./),\
(/0.,0.,0.,0.,0.,0.,0.,0.,1./8.,6./8.,1./8.,0./),\
(/0.,0.,0.,0.,0.,0.,0.,0.,0.,1./8.,6./8.,1./8./),\
(/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1./8.,7./8./)\
/)

AI=inverse_matrix(A)
Year=ispan(1901,2012,1)
ntime=12
nlat=360
nlon=720
nyear=dimsizes(Year)
Tm=new((/nyear,ntime,nlat,nlon/),"float")
Tm!0="year"
Tm&year=Year

fils = systemfunc ("ls "+int_diri+"Tmp_GLDAS_*.nc") ; file paths
  f    = addfiles (fils, "r")   
  ListSetType (f, "cat")
  t=f[:]->T
  time=f[:]->time
  utc_date = cd_calendar(time, 0)
   year   = tointeger(utc_date(:,0))  


do nf=0,dimsizes(Year)-1
     iyyyy=ind(year.eq.Year(nf))
     Tm(nf,:,:,:)=t(iyyyy,:,:) 
end do
     lat=f[0]->lat
     lon=f[0]->lon

printVarSummary(Tm)
Tm_new=reshape(Tm,(/12*nyear,nlat,nlon/))
copy_VarCoords(Tm(0,0,:,:),Tm_new(0,:,:))
yyyymm = yyyymm_time(1901, 2012, "integer")
Tm_new!0="time"
Tm_new&time=yyyymm
printVarSummary(Tm_new)


 Tanom_D=Tm

do nx=0,nlat-1
  do ny=0,nlon-1
     do nt=0,nyear-1
        do nm=0,11
             T=AI(nm,0)*Tm(nt,0,nx,ny)+AI(nm,1)*Tm(nt,1,nx,ny)+AI(nm,2)*Tm(nt,2,nx,ny)+AI(nm,3)*Tm(nt,3,nx,ny)+AI(nm,4)*Tm(nt,4,nx,ny)+AI(nm,5)*Tm(nt,5,nx,ny)+AI(nm,6)*Tm(nt,6,nx,ny)+AI(nm,7)*Tm(nt,7,nx,ny)+AI(nm,8)*Tm(nt,8,nx,ny)+AI(nm,9)*Tm(nt,9,nx,ny)+AI(nm,10)*Tm(nt,10,nx,ny)+AI(nm,11)*Tm(nt,11,nx,ny)
               Tanom_D(nt,nm,nx,ny)=T
        end do
    end do
 end do
end do


fo->Tanom_mid=Tanom_D



