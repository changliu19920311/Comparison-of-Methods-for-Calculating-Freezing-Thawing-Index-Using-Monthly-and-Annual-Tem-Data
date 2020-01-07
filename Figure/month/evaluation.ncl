; evaluating the reconstucted monthly data=================
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

;================set the parameters=====================

diri1="/Users/chang/Desktop/Small_work/permafrost/Feng/data/"
name1="GLDAS_FI_Feng.nc"
diri2="/Users/chang/Desktop/Small_work/permafrost/Oliver/data/"
name2 = "GLDAS_FI_Oliver.nc"
diri3="/Users/chang/Desktop/Small_work/permafrost/New_Method/monthly/data/Frost_index_data/"
name3 = "GLDAS_FI_New_ensemble.nc"
diri4="/Users/chang/Desktop/Small_work/permafrost/GLDAS/data/"
name4="GLDAS_FI_ori.nc"

Year1=1901
Year2=2011
Year=ispan(Year1,Year2,1)


out_data_ddt=new((/3,5/),"float")
out_data_ddf=new((/3,5/),"float")
out_data_fi=new((/3,5/),"float")

;=================read the data, because the value of 2012 is missing, we only use data from 1901-2011==========================

  f1  = addfile(diri1+name1,"r")
  lat=f1->lat
  lon=f1->lon
  year_1=Year(:110)
  fi_Feng = f1->fi(:110,:,:)
  ddt_Feng = f1->ddt(:110,:,:)
  ddf_Feng = f1->ddf(:110,:,:)
  printVarSummary(fi_Feng)

  f2  = addfile(diri2+name2,"r")
  fi_Oliver= f2->fi(:110,:,:)
  ddt_Oliver= f2->ddt(:110,:,:)
  ddf_Oliver= f2->ddf(:110,:,:)
  printVarSummary(fi_Oliver)

  f3  = addfile(diri3+name3,"r")
  fi_New= f3->fi(:110,:,:)
  ddt_New= f3->ddt(:110,:,:)
  ddf_New= f3->ddf(:110,:,:)
  printVarSummary(fi_New)

  f4  = addfile(diri4+name4,"r")
  fi_GLDAS= f4->fi(:110,:,:)
  ddt_GLDAS= f4->ddt(:110,:,:)
  ddf_GLDAS= f4->ddf(:110,:,:)
  printVarSummary(fi_GLDAS)


 a = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
  lsdata = a->LSMASK
  lsm = landsea_mask(lsdata,lat,lon)

  fi_GLDAS  = mask(fi_GLDAS ,lsm.eq.0,False)  
  fi_GLDAS  = mask(fi_GLDAS ,lsm.eq.2,False)  
  fi_GLDAS  = mask(fi_GLDAS ,lsm.eq.4,False)  
 
ddt_GLDAS = where(ismissing(fi_GLDAS),fi_GLDAS@_FillValue,ddt_GLDAS )
ddf_GLDAS = where(ismissing(fi_GLDAS),fi_GLDAS@_FillValue,ddf_GLDAS )

fi_New = where(ismissing(fi_GLDAS),fi_GLDAS@_FillValue,fi_New )
ddt_New = where(ismissing(fi_GLDAS),fi_GLDAS@_FillValue,ddt_New )
ddf_New = where(ismissing(fi_GLDAS),fi_GLDAS@_FillValue,ddf_New )

fi_Oliver = where(ismissing(fi_GLDAS),fi_GLDAS@_FillValue,fi_Oliver )
ddt_Oliver = where(ismissing(fi_GLDAS),fi_GLDAS@_FillValue,ddt_Oliver )
ddf_Oliver = where(ismissing(fi_GLDAS),fi_GLDAS@_FillValue,ddf_Oliver )

fi_Feng = where(ismissing(fi_GLDAS),fi_GLDAS@_FillValue,fi_Feng )
ddt_Feng = where(ismissing(fi_GLDAS),fi_GLDAS@_FillValue,ddt_Feng )
ddf_Feng = where(ismissing(fi_GLDAS),fi_GLDAS@_FillValue,ddf_Feng )


;================= using the distribution of permafrost to mask fi ddt and ddf==================

fi_GLDAS_ave=dim_avg_n_Wrap(fi_GLDAS,0)
fi_GLDAS_ave=where(fi_GLDAS_ave.le.0.495,fi_GLDAS_ave@_FillValue,fi_GLDAS_ave)
fi_GLDAS_ave_test=conform_dims(dimsizes(fi_GLDAS),fi_GLDAS_ave,(/1,2/))

;ddt_GLDAS=where(ismissing(fi_GLDAS_ave_test),ddt_GLDAS@_FillValue,ddt_GLDAS)
;ddf_GLDAS=where(ismissing(fi_GLDAS_ave_test),ddf_GLDAS@_FillValue,ddf_GLDAS)

;fi_New=where(ismissing(fi_GLDAS_ave_test),fi_New@_FillValue,fi_New)
;ddt_New=where(ismissing(fi_GLDAS_ave_test),ddt_New@_FillValue,ddt_New)
;ddf_New=where(ismissing(fi_GLDAS_ave_test),ddf_New@_FillValue,ddf_New)

;fi_Oliver=where(ismissing(fi_GLDAS_ave_test),fi_Oliver@_FillValue,fi_Oliver)
;ddt_Oliver=where(ismissing(fi_GLDAS_ave_test),ddt_Oliver@_FillValue,ddt_Oliver)
;ddf_Oliver=where(ismissing(fi_GLDAS_ave_test),ddf_Oliver@_FillValue,ddf_Oliver)

;fi_Feng=where(ismissing(fi_GLDAS_ave_test),fi_Feng@_FillValue,fi_Feng)
;ddt_Feng=where(ismissing(fi_GLDAS_ave_test),ddt_Feng@_FillValue,ddt_Feng)
;ddf_Feng=where(ismissing(fi_GLDAS_ave_test),ddf_Feng@_FillValue,ddf_Feng)


;================calculate the relative error========================================================


re_Feng=dim_sum_n_Wrap(abs(fi_Feng-fi_GLDAS),0)
RE_fi_Feng=re_Feng/(dimsizes(Year)*fi_GLDAS_ave)*100
copy_VarCoords(fi_GLDAS_ave,RE_fi_Feng) 
re_Oliver=dim_sum_n_Wrap(abs(fi_Oliver-fi_GLDAS),0)
RE_fi_Oliver=re_Oliver/(dimsizes(Year)*fi_GLDAS_ave)*100
copy_VarCoords(fi_GLDAS_ave,RE_fi_Oliver) 
re_New=dim_sum_n_Wrap(abs(fi_New-fi_GLDAS),0)
RE_fi_New=re_New/(dimsizes(Year)*fi_GLDAS_ave)*100
copy_VarCoords(fi_GLDAS_ave,RE_fi_New) 


ddt_GLDAS_ave=dim_avg_n_Wrap(ddt_GLDAS,0)
ddt_GLDAS_ave=where(ddt_GLDAS_ave.eq.0,0.00001,ddt_GLDAS_ave)
re_Feng=dim_sum_n_Wrap(abs(ddt_Feng-ddt_GLDAS),0)
RE_ddt_Feng=re_Feng/(dimsizes(Year)*ddt_GLDAS_ave)*100
copy_VarCoords(ddt_GLDAS_ave,RE_ddt_Feng) 
re_Oliver=dim_sum_n_Wrap(abs(ddt_Oliver-ddt_GLDAS),0)

RE_ddt_Oliver=re_Oliver/(dimsizes(Year)*ddt_GLDAS_ave)*100
copy_VarCoords(ddt_GLDAS_ave,RE_ddt_Oliver) 
re_New=dim_sum_n_Wrap(abs(ddt_New-ddt_GLDAS),0)
RE_ddt_New=re_New/(dimsizes(Year)*ddt_GLDAS_ave)*100
copy_VarCoords(ddt_GLDAS_ave,RE_ddt_New) 


ddf_GLDAS_ave=dim_avg_n_Wrap(ddf_GLDAS,0)
ddf_GLDAS_ave=where(ddf_GLDAS_ave.eq.0,0.00001,ddf_GLDAS_ave)
re_Feng=dim_sum_n_Wrap(abs(ddf_Feng-ddf_GLDAS),0)
RE_ddf_Feng=re_Feng/(dimsizes(Year)*ddf_GLDAS_ave)*100
copy_VarCoords(ddf_GLDAS_ave,RE_ddf_Feng) 
re_Oliver=dim_sum_n_Wrap(abs(ddf_Oliver-ddf_GLDAS),0)
RE_ddf_Oliver=re_Oliver/(dimsizes(Year)*ddf_GLDAS_ave)*100
copy_VarCoords(ddf_GLDAS_ave,RE_ddf_Oliver) 
re_New=dim_sum_n_Wrap(abs(ddf_New-ddf_GLDAS),0)
RE_ddf_New=re_New/(dimsizes(Year)*ddf_GLDAS_ave)*100
copy_VarCoords(ddf_GLDAS_ave,RE_ddf_New) 




RE_ddf_New=where(RE_ddf_New.ge.100,RE_ddf_New@_FillValue,RE_ddf_New)
RE_ddt_New=where(RE_ddt_New.ge.100,RE_ddt_New@_FillValue,RE_ddt_New)
RE_fi_New=where(RE_fi_New.ge.100,RE_fi_New@_FillValue,RE_fi_New)

RE_ddf_Feng=where(RE_ddf_Feng.ge.100,RE_ddf_Feng@_FillValue,RE_ddf_Feng)
RE_ddt_Feng=where(RE_ddt_Feng.ge.100,RE_ddt_Feng@_FillValue,RE_ddt_Feng)
RE_fi_Feng=where(RE_fi_Feng.ge.100,RE_fi_Feng@_FillValue,RE_fi_Feng)

RE_ddf_Oliver=where(RE_ddf_Oliver.ge.100,RE_ddf_Oliver@_FillValue,RE_ddf_Oliver)
RE_ddt_Oliver=where(RE_ddt_Oliver.ge.100,RE_ddt_Oliver@_FillValue,RE_ddt_Oliver)
RE_fi_Oliver=where(RE_fi_Oliver.ge.100,RE_fi_Oliver@_FillValue,RE_fi_Oliver)


  rad    = 4.0*atan(1.0)/180.0
  re     = 6371220.0
  rr     = re*rad
  clat   = cos(lat*rad)

 RE_fi_Feng_avg = wgt_areaave(RE_fi_Feng, clat, 1.0, 0)
 RE_ddt_Feng_avg = wgt_areaave(RE_ddt_Feng, clat, 1.0, 0)
 RE_ddf_Feng_avg = wgt_areaave(RE_ddf_Feng, clat, 1.0, 0)

 RE_fi_Oliver_avg = wgt_areaave(RE_fi_Oliver, clat, 1.0, 0)
 RE_ddt_Oliver_avg = wgt_areaave(RE_ddt_Oliver, clat, 1.0, 0)
 RE_ddf_Oliver_avg = wgt_areaave(RE_ddf_Oliver, clat, 1.0, 0)

 RE_fi_New_avg = wgt_areaave(RE_fi_New, clat, 1.0, 0)
 RE_ddt_New_avg = wgt_areaave(RE_ddt_New, clat, 1.0, 0)
 RE_ddf_New_avg = wgt_areaave(RE_ddf_New, clat, 1.0, 0)

out_data_ddt(0,0)=RE_ddt_Feng_avg
out_data_ddt(1,0)=RE_ddt_Oliver_avg
out_data_ddt(2,0)=RE_ddt_New_avg

out_data_ddf(0,0)=RE_ddf_Feng_avg
out_data_ddf(1,0)=RE_ddf_Oliver_avg
out_data_ddf(2,0)=RE_ddf_New_avg

out_data_fi(0,0)=RE_fi_Feng_avg
out_data_fi(1,0)=RE_fi_Oliver_avg
out_data_fi(2,0)=RE_fi_New_avg




;====================calculate Root mean square error===============================================

rmse_fi_Feng = dim_rmsd_n(fi_Feng,fi_GLDAS,0)*100
copy_VarCoords(fi_GLDAS_ave,rmse_fi_Feng) 
rmse_fi_Oliver = dim_rmsd_n(fi_Oliver,fi_GLDAS,0)*100
copy_VarCoords(fi_GLDAS_ave,rmse_fi_Oliver) 
rmse_fi_New = dim_rmsd_n(fi_New,fi_GLDAS,0)*100
copy_VarCoords(fi_GLDAS_ave,rmse_fi_New) 

rmse_ddt_Feng = dim_rmsd_n(ddt_Feng,ddt_GLDAS,0)*100
copy_VarCoords(ddt_GLDAS_ave,rmse_ddt_Feng) 
rmse_ddt_Oliver = dim_rmsd_n(ddt_Oliver,ddt_GLDAS,0)*100
copy_VarCoords(ddt_GLDAS_ave,rmse_ddt_Oliver) 
rmse_ddt_New = dim_rmsd_n(ddt_New,ddt_GLDAS,0)*100
copy_VarCoords(ddt_GLDAS_ave,rmse_ddt_New) 

rmse_ddf_Feng = dim_rmsd_n(ddf_Feng,ddf_GLDAS,0)*100
copy_VarCoords(ddf_GLDAS_ave,rmse_ddf_Feng) 
rmse_ddf_Oliver = dim_rmsd_n(ddf_Oliver,ddf_GLDAS,0)*100
copy_VarCoords(ddf_GLDAS_ave,rmse_ddf_Oliver) 
rmse_ddf_New = dim_rmsd_n(ddf_New,ddf_GLDAS,0)*100
copy_VarCoords(ddf_GLDAS_ave,rmse_ddf_New) 


 RMSE_fi_Feng_avg = wgt_areaave(rmse_fi_Feng, clat, 1.0, 0)
 RMSE_ddt_Feng_avg = wgt_areaave(rmse_ddt_Feng, clat, 1.0, 0)
 RMSE_ddf_Feng_avg = wgt_areaave(rmse_ddf_Feng, clat, 1.0, 0)

 RMSE_fi_Oliver_avg = wgt_areaave(rmse_fi_Oliver, clat, 1.0, 0)
 RMSE_ddt_Oliver_avg = wgt_areaave(rmse_ddt_Oliver, clat, 1.0, 0)
 RMSE_ddf_Oliver_avg = wgt_areaave(rmse_ddf_Oliver, clat, 1.0, 0)

 RMSE_fi_New_avg = wgt_areaave(rmse_fi_New, clat, 1.0, 0)
 RMSE_ddt_New_avg = wgt_areaave(rmse_ddt_New, clat, 1.0, 0)
 RMSE_ddf_New_avg = wgt_areaave(rmse_ddf_New, clat, 1.0, 0)


out_data_ddt(0,1)=RMSE_ddt_Feng_avg
out_data_ddt(1,1)=RMSE_ddt_Oliver_avg
out_data_ddt(2,1)=RMSE_ddt_New_avg

out_data_ddf(0,1)=RMSE_ddf_Feng_avg
out_data_ddf(1,1)=RMSE_ddf_Oliver_avg
out_data_ddf(2,1)=RMSE_ddf_New_avg

out_data_fi(0,1)=RMSE_fi_Feng_avg
out_data_fi(1,1)=RMSE_fi_Oliver_avg
out_data_fi(2,1)=RMSE_fi_New_avg
;==============================correlation=====================================================

    n=dimsizes(Year)
    df=n-2
    rr_fi_Feng   = escorc_n(fi_Feng,fi_GLDAS,0,0) 
    rr_fi_Feng = where(rr_fi_Feng.eq.1,0.999,rr_fi_Feng)
    t    = rr_fi_Feng*sqrt((n-2)/(1-rr_fi_Feng^2))      
    p_rr_fi_Feng = student_t(t, df)
    delete(t)
    copy_VarCoords(fi_GLDAS(0,:,:),rr_fi_Feng) 
    copy_VarCoords(fi_GLDAS(0,:,:),p_rr_fi_Feng) 

    rr_fi_Oliver   = escorc_n(fi_Oliver,fi_GLDAS,0,0) 
    rr_fi_Oliver = where(rr_fi_Oliver.eq.1,0.999,rr_fi_Oliver)
    t    = rr_fi_Oliver*sqrt((n-2)/(1-rr_fi_Oliver^2))      
    p_rr_fi_Oliver = student_t(t, df)
    delete(t)
    copy_VarCoords(fi_GLDAS(0,:,:),rr_fi_Oliver) 
    copy_VarCoords(fi_GLDAS(0,:,:),p_rr_fi_Oliver) 

    rr_fi_New   = escorc_n(fi_New,fi_GLDAS,0,0) 
    rr_fi_New = where(rr_fi_New.eq.1,0.999,rr_fi_New)
    t    = rr_fi_New*sqrt((n-2)/(1-rr_fi_New^2))      
    p_rr_fi_New = student_t(t, df)
    delete(t)
    copy_VarCoords(fi_GLDAS(0,:,:),rr_fi_New) 
    copy_VarCoords(fi_GLDAS(0,:,:),p_rr_fi_New) 


    rr_ddt_Feng   = escorc_n(ddt_Feng,ddt_GLDAS,0,0) 
    rr_ddt_Feng = where(rr_ddt_Feng.eq.1,0.999,rr_ddt_Feng)
    t    = rr_ddt_Feng*sqrt((n-2)/(1-rr_ddt_Feng^2))      
    p_rr_ddt_Feng = student_t(t, df)
    delete(t)
    copy_VarCoords(ddt_GLDAS(0,:,:),rr_ddt_Feng) 
    copy_VarCoords(ddt_GLDAS(0,:,:),p_rr_ddt_Feng) 

    rr_ddt_Oliver   = escorc_n(ddt_Oliver,ddt_GLDAS,0,0) 
    rr_ddt_Oliver = where(rr_ddt_Oliver.eq.1,0.999,rr_ddt_Oliver)
    t    = rr_ddt_Oliver*sqrt((n-2)/(1-rr_ddt_Oliver^2))      
    p_rr_ddt_Oliver = student_t(t, df)
    delete(t)
    copy_VarCoords(ddt_GLDAS(0,:,:),rr_ddt_Oliver) 
    copy_VarCoords(ddt_GLDAS(0,:,:),p_rr_ddt_Oliver) 

    rr_ddt_New   = escorc_n(ddt_New,ddt_GLDAS,0,0) 
    rr_ddt_New = where(rr_ddt_New.eq.1,0.999,rr_ddt_New)
    t    = rr_ddt_New*sqrt((n-2)/(1-rr_ddt_New^2))      
    p_rr_ddt_New = student_t(t, df)
    delete(t)
    copy_VarCoords(ddt_GLDAS(0,:,:),rr_ddt_New) 
    copy_VarCoords(ddt_GLDAS(0,:,:),p_rr_ddt_New) 

    rr_ddf_Feng   = escorc_n(ddf_Feng,ddf_GLDAS,0,0) 
    rr_ddf_Feng = where(rr_ddf_Feng.eq.1,0.999,rr_ddf_Feng)
    t    = rr_ddf_Feng*sqrt((n-2)/(1-rr_ddf_Feng^2))      
    p_rr_ddf_Feng = student_t(t, df)
    delete(t)
    copy_VarCoords(ddf_GLDAS(0,:,:),rr_ddf_Feng) 
    copy_VarCoords(ddf_GLDAS(0,:,:),p_rr_ddf_Feng) 

    rr_ddf_Oliver   = escorc_n(ddf_Oliver,ddf_GLDAS,0,0) 
    rr_ddf_Oliver = where(rr_ddf_Oliver.eq.1,0.999,rr_ddf_Oliver)
    t    = rr_ddf_Oliver*sqrt((n-2)/(1-rr_ddf_Oliver^2))      
    p_rr_ddf_Oliver = student_t(t, df)
    delete(t)
    copy_VarCoords(ddf_GLDAS(0,:,:),rr_ddf_Oliver) 
    copy_VarCoords(ddf_GLDAS(0,:,:),p_rr_ddf_Oliver) 

    rr_ddf_New   = escorc_n(ddf_New,ddf_GLDAS,0,0) 
    rr_ddf_New = where(rr_ddf_New.eq.1,0.999,rr_ddf_New)
    t    = rr_ddf_New*sqrt((n-2)/(1-rr_ddf_New^2))      
    p_rr_ddf_New = student_t(t, df)
    delete(t)
    copy_VarCoords(ddf_GLDAS(0,:,:),rr_ddf_New) 
    copy_VarCoords(ddf_GLDAS(0,:,:),p_rr_ddf_New) 

    rr_fi_New_avg  = wgt_areaave(rr_fi_New, clat, 1.0, 0)
    rr_ddt_New_avg = wgt_areaave(rr_ddt_New, clat, 1.0, 0)
    rr_ddf_New_avg = wgt_areaave(rr_ddf_New, clat, 1.0, 0)
    p_rr_fi_New_avg  = wgt_areaave(p_rr_fi_New, clat, 1.0, 0)
    p_rr_ddt_New_avg = wgt_areaave(p_rr_ddt_New, clat, 1.0, 0)
    p_rr_ddf_New_avg = wgt_areaave(p_rr_ddf_New, clat, 1.0, 0)

    rr_fi_Feng_avg  = wgt_areaave(rr_fi_Feng, clat, 1.0, 0)
    rr_ddt_Feng_avg = wgt_areaave(rr_ddt_Feng, clat, 1.0, 0)
    rr_ddf_Feng_avg = wgt_areaave(rr_ddf_Feng, clat, 1.0, 0)
    p_rr_fi_Feng_avg  = wgt_areaave(p_rr_fi_Feng, clat, 1.0, 0)
    p_rr_ddt_Feng_avg = wgt_areaave(p_rr_ddt_Feng, clat, 1.0, 0)
    p_rr_ddf_Feng_avg = wgt_areaave(p_rr_ddf_Feng, clat, 1.0, 0)

    rr_fi_Oliver_avg  = wgt_areaave(rr_fi_Oliver, clat, 1.0, 0)
    rr_ddt_Oliver_avg = wgt_areaave(rr_ddt_Oliver, clat, 1.0, 0)
    rr_ddf_Oliver_avg = wgt_areaave(rr_ddf_Oliver, clat, 1.0, 0)
    p_rr_fi_Oliver_avg  = wgt_areaave(p_rr_fi_Oliver, clat, 1.0, 0)
    p_rr_ddt_Oliver_avg = wgt_areaave(p_rr_ddt_Oliver, clat, 1.0, 0)
    p_rr_ddf_Oliver_avg = wgt_areaave(p_rr_ddf_Oliver, clat, 1.0, 0)

out_data_ddt(0,2)=rr_ddt_Feng_avg
out_data_ddt(1,2)=rr_ddt_Oliver_avg
out_data_ddt(2,2)=rr_ddt_New_avg

out_data_ddf(0,2)=rr_ddf_Feng_avg
out_data_ddf(1,2)=rr_ddf_Oliver_avg
out_data_ddf(2,2)=rr_ddf_New_avg

out_data_fi(0,2)=rr_fi_Feng_avg
out_data_fi(1,2)=rr_fi_Oliver_avg
out_data_fi(2,2)=rr_fi_New_avg

out_data_ddt(0,3)=p_rr_ddt_Feng_avg
out_data_ddt(1,3)=p_rr_ddt_Oliver_avg
out_data_ddt(2,3)=p_rr_ddt_New_avg

out_data_ddf(0,3)=p_rr_ddf_Feng_avg
out_data_ddf(1,3)=p_rr_ddf_Oliver_avg
out_data_ddf(2,3)=p_rr_ddf_New_avg

out_data_fi(0,3)=p_rr_fi_Feng_avg
out_data_fi(1,3)=p_rr_fi_Oliver_avg
out_data_fi(2,3)=p_rr_fi_New_avg


;==============================Difference=====================================================

ddf_GLDAS_ave=dim_avg_n_Wrap(ddf_GLDAS,0)
ddt_GLDAS_ave=dim_avg_n_Wrap(ddt_GLDAS,0)
fi_GLDAS_ave=dim_avg_n_Wrap(fi_GLDAS,0)

ddf_Oliver_ave=dim_avg_n_Wrap(ddf_Oliver,0)
ddt_Oliver_ave=dim_avg_n_Wrap(ddt_Oliver,0)
fi_Oliver_ave=dim_avg_n_Wrap(fi_Oliver,0)

ddf_Feng_ave=dim_avg_n_Wrap(ddf_Feng,0)
ddt_Feng_ave=dim_avg_n_Wrap(ddt_Feng,0)
fi_Feng_ave=dim_avg_n_Wrap(fi_Feng,0)

ddf_New_ave=dim_avg_n_Wrap(ddf_New,0)
ddt_New_ave=dim_avg_n_Wrap(ddt_New,0)
fi_New_ave=dim_avg_n_Wrap(fi_New,0)


dif_fi_New = abs(fi_New_ave-fi_GLDAS_ave)
dif_ddf_New = abs(ddf_New_ave-ddf_GLDAS_ave)
dif_ddt_New = abs(ddt_New_ave-ddt_GLDAS_ave)

dif_fi_Feng = abs(fi_Feng_ave-fi_GLDAS_ave)
dif_ddf_Feng = abs(ddf_Feng_ave-ddf_GLDAS_ave)
dif_ddt_Feng = abs(ddt_Feng_ave-ddt_GLDAS_ave)

dif_fi_Oliver = abs(fi_Oliver_ave-fi_GLDAS_ave)
dif_ddf_Oliver = abs(ddf_Oliver_ave-ddf_GLDAS_ave)
dif_ddt_Oliver = abs(ddt_Oliver_ave-ddt_GLDAS_ave)

copy_VarCoords(ddf_GLDAS_ave,dif_fi_Oliver) 
copy_VarCoords(ddf_GLDAS_ave,dif_ddf_Oliver) 
copy_VarCoords(ddf_GLDAS_ave,dif_ddt_Oliver) 
copy_VarCoords(ddf_GLDAS_ave,dif_fi_Feng) 
copy_VarCoords(ddf_GLDAS_ave,dif_ddf_Feng) 
copy_VarCoords(ddf_GLDAS_ave,dif_ddt_Feng) 
copy_VarCoords(ddf_GLDAS_ave,dif_fi_New) 
copy_VarCoords(ddf_GLDAS_ave,dif_ddf_New) 
copy_VarCoords(ddf_GLDAS_ave,dif_ddt_New) 


 dif_fi_Feng_avg = wgt_areaave(dif_fi_Feng, clat, 1.0, 0)
 dif_ddt_Feng_avg = wgt_areaave(dif_ddt_Feng, clat, 1.0, 0)
 dif_ddf_Feng_avg = wgt_areaave(dif_ddf_Feng, clat, 1.0, 0)


 dif_fi_Oliver_avg = wgt_areaave(dif_fi_Oliver, clat, 1.0, 0)
 dif_ddt_Oliver_avg = wgt_areaave(dif_ddt_Oliver, clat, 1.0, 0)
 dif_ddf_Oliver_avg = wgt_areaave(dif_ddf_Oliver, clat, 1.0, 0)

 dif_fi_New_avg = wgt_areaave(dif_fi_New, clat, 1.0, 0)
 dif_ddt_New_avg = wgt_areaave(dif_ddt_New, clat, 1.0, 0)
 dif_ddf_New_avg = wgt_areaave(dif_ddf_New, clat, 1.0, 0)

out_data_ddt(0,4)=dif_ddt_Feng_avg
out_data_ddt(1,4)=dif_ddt_Oliver_avg
out_data_ddt(2,4)=dif_ddt_New_avg

out_data_ddf(0,4)=dif_ddf_Feng_avg
out_data_ddf(1,4)=dif_ddf_Oliver_avg
out_data_ddf(2,4)=dif_ddf_New_avg

out_data_fi(0,4)=dif_fi_Feng_avg
out_data_fi(1,4)=dif_fi_Oliver_avg
out_data_fi(2,4)=dif_fi_New_avg

print((/dif_ddt_Oliver_avg,dif_ddf_Oliver_avg,dif_fi_Oliver_avg,dif_ddt_Feng_avg,dif_ddf_Feng_avg,dif_fi_Feng_avg,dif_ddt_New_avg,dif_ddf_New_avg,dif_fi_New_avg/))

a1=(abs(dif_ddt_Oliver_avg)-abs(dif_ddt_New_avg))/abs(dif_ddt_Oliver_avg)
a2=(abs(dif_ddf_Oliver_avg)-abs(dif_ddf_New_avg))/abs(dif_ddf_Oliver_avg)
a3=(abs(dif_fi_Oliver_avg)-abs(dif_fi_New_avg))/abs(dif_fi_Oliver_avg)

b1=(abs(dif_ddt_Feng_avg)-abs(dif_ddt_New_avg))/abs(dif_ddt_Feng_avg)
b2=(abs(dif_ddf_Feng_avg)-abs(dif_ddf_New_avg))/abs(dif_ddf_Feng_avg)
b3=(abs(dif_fi_Feng_avg)-abs(dif_fi_New_avg))/abs(dif_fi_Feng_avg)

print((/a1,a2,a3/))
print((/b1,b2,b3/))

;================output data================
 opt = True
 opt@title  = "RE      RMSE      Cor   P_level    MB     "
 opt@fout   = "FI_evaluation.txt"
 write_matrix (out_data_fi, "5f12.3", opt)
 opt@fout   = "DDT_evaluation.txt"
 write_matrix (out_data_ddt, "5f12.3", opt)
 opt@fout   = "DDF_evaluation.txt"
 write_matrix (out_data_ddf, "5f12.3", opt)
;==============================plot======================================================================

fo=addfile("merit_New.nc","c")
fo->RE_fi_New=RE_fi_New
fo->RE_ddf_New=RE_ddf_New
fo->RE_ddt_New=RE_ddt_New
fo->dif_fi_New=dif_fi_New
fo->dif_ddf_New=dif_ddf_New
fo->dif_ddt_New=dif_ddt_New
fo->rr_ddf_New=rr_ddf_New
fo->rr_ddt_New=rr_ddt_New
fo->rr_fi_New=rr_fi_New


fo=addfile("merit_Oliver.nc","c")
fo->RE_fi_Oliver=RE_fi_Oliver
fo->RE_ddf_Oliver=RE_ddf_Oliver
fo->RE_ddt_Oliver=RE_ddt_Oliver
fo->dif_fi_Oliver=dif_fi_Oliver
fo->dif_ddf_Oliver=dif_ddf_Oliver
fo->dif_ddt_Oliver=dif_ddt_Oliver
fo->rr_ddf_Oliver=rr_ddf_Oliver
fo->rr_ddt_Oliver=rr_ddt_Oliver
fo->rr_fi_Oliver=rr_fi_Oliver

fo=addfile("merit_Feng.nc","c")
fo->RE_fi_Feng=RE_fi_Feng
fo->RE_ddf_Feng=RE_ddf_Feng
fo->RE_ddt_Feng=RE_ddt_Feng
fo->dif_fi_Feng=dif_fi_Feng
fo->dif_ddf_Feng=dif_ddf_Feng
fo->dif_ddt_Feng=dif_ddt_Feng
fo->rr_ddf_Feng=rr_ddf_Feng
fo->rr_ddt_Feng=rr_ddt_Feng
fo->rr_fi_Feng=rr_fi_Feng




plot_re=new(9,"graphic")
plot_rmse=new(9,"graphic")
plot_rr=new(9,"graphic")
plot_dif=new(9,"graphic")
wks = gsn_open_wks("pdf" ,"Measures")      ; send graphics to PNG file
  gsn_define_colormap(wks,"amwg_blueyellowred")

  res                       = True     ; plot mods desired
  res@cnFillOn              = True     ; turn on color fill
  res@cnLinesOn             = False    ; turn of contour lines
  plot=new(10,graphic)  
  res@gsnDraw = False 
  res@gsnFrame = False 
  res@mpCenterLonF=-45
  res@gsnPolar   = "NH"  
  res@mpMinLatF            = 20    ; range to zoom in on
  res@mpFillOn             = False              
  res@cnFillOn              = True     ; turn on color fill
  res@cnLinesOn             = False    ; turn of contour lines
  ;res@cnFillMode            = "RasterFill"    
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnFillColors        =(/2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17/)
  res@lbLabelStride         = 1 
  res@tiMainString        = " "
  res@cnLevels = (/0.1,0.25,0.5,0.75,1.,5.,10.,20.,30.,40.,45.,50.,60.,70.,80./)
  res@gsnLeftString ="RE_FI_Oliver"
  plot_re(0)=gsn_csm_contour_map_polar(wks,RE_fi_Oliver,res)
  res@gsnLeftString ="RE_FI_Feng"
  plot_re(1)=gsn_csm_contour_map_polar(wks,RE_fi_Feng,res)
  res@gsnLeftString ="RE_FI_New"
  plot_re(2)=gsn_csm_contour_map_polar(wks,RE_fi_New,res)
  res@gsnLeftString ="RE_DDT_Oliver"
  plot_re(3)=gsn_csm_contour_map_polar(wks,RE_ddt_Oliver,res)
  res@gsnLeftString ="RE_DDT_Feng"
  plot_re(4)=gsn_csm_contour_map_polar(wks,RE_ddt_Feng,res)
  res@gsnLeftString ="RE_DDT_New"
  plot_re(5)=gsn_csm_contour_map_polar(wks,RE_ddt_New,res)
  res@gsnLeftString ="RE_DDF_Oliver"
  plot_re(6)=gsn_csm_contour_map_polar(wks,RE_ddf_Oliver,res)
  res@gsnLeftString ="RE_DDF_Feng"
  plot_re(7)=gsn_csm_contour_map_polar(wks,RE_ddf_Feng,res)
  res@gsnLeftString ="RE_DDF_New"
  plot_re(8)=gsn_csm_contour_map_polar(wks,RE_ddf_New,res)



  res@cnLevels = (/0.1,0.2,0.5,1.,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0/)
  plot_rmse(0)=gsn_csm_contour_map_polar(wks,rmse_fi_Oliver,res)
  plot_rmse(1)=gsn_csm_contour_map_polar(wks,rmse_fi_Feng,res)
  plot_rmse(2)=gsn_csm_contour_map_polar(wks,rmse_fi_New,res)
  res@cnLevels = (/500.,1000.,1500.,2000.,2500.,3000.,3500.,4000.,4500.,5000.,5500.,6000.,6500.,7000.,7500./)
  plot_rmse(3)=gsn_csm_contour_map_polar(wks,rmse_ddt_Oliver,res)
  plot_rmse(4)=gsn_csm_contour_map_polar(wks,rmse_ddt_Feng,res)
  plot_rmse(5)=gsn_csm_contour_map_polar(wks,rmse_ddt_New,res)
  res@cnLevels = (/500.,1000.,1500.,2000.,2500.,3000.,3500.,4000.,4500.,5000.,5500.,6000.,6500.,7000.,7500./)
  plot_rmse(6)=gsn_csm_contour_map_polar(wks,rmse_ddf_Oliver,res)
  plot_rmse(7)=gsn_csm_contour_map_polar(wks,rmse_ddf_Feng,res)
  plot_rmse(8)=gsn_csm_contour_map_polar(wks,rmse_ddf_New,res)



  res@cnLevels = (/0.91,0.92,0.93,0.94,0.945,0.95,0.955,0.96,0.965,0.97,0.975,0.98,0.985,0.99,0.995/)

  res@gsnLeftString ="rr_FI_Oliver"
  plot_rr(0)=gsn_csm_contour_map_polar(wks,rr_fi_Oliver,res)
  res@gsnLeftString ="rr_FI_Feng"
  plot_rr(1)=gsn_csm_contour_map_polar(wks,rr_fi_Feng,res)
  res@gsnLeftString ="rr_FI_New"
  plot_rr(2)=gsn_csm_contour_map_polar(wks,rr_fi_New,res)
  res@gsnLeftString ="rr_DDT_Oliver"
  plot_rr(3)=gsn_csm_contour_map_polar(wks,rr_ddt_Oliver,res)
  res@gsnLeftString ="rr_DDT_Feng"
  plot_rr(4)=gsn_csm_contour_map_polar(wks,rr_ddt_Feng,res)
  res@gsnLeftString ="rr_DDT_New"
  plot_rr(5)=gsn_csm_contour_map_polar(wks,rr_ddt_New,res)
  res@gsnLeftString ="rr_DDF_Oliver"
  plot_rr(6)=gsn_csm_contour_map_polar(wks,rr_ddf_Oliver,res)
  res@gsnLeftString ="rr_DDF_Feng"
  plot_rr(7)=gsn_csm_contour_map_polar(wks,rr_ddf_Feng,res)
  res@gsnLeftString ="rr_DDF_New"
  plot_rr(8)=gsn_csm_contour_map_polar(wks,rr_ddf_New,res)



  res@cnLevels :=(/-0.07,-0.06,-0.05,-0.04,-0.03,-0.02,-0.01,0,0.01,0.02,0.03,0.04,0.05,0.06,0.07/)
  
  plot_dif(0)=gsn_csm_contour_map_polar(wks,dif_fi_Oliver,res)
  plot_dif(1)=gsn_csm_contour_map_polar(wks,dif_fi_Feng,res)
  plot_dif(2)=gsn_csm_contour_map_polar(wks,dif_fi_New,res)
  res@cnLevels := (/-90,-70,-50,-40,-30,-20,-10,0,10,20,30,40,50,70,90/)
  plot_dif(3)=gsn_csm_contour_map_polar(wks,dif_ddt_Oliver,res)
  plot_dif(4)=gsn_csm_contour_map_polar(wks,dif_ddt_Feng,res)
  plot_dif(5)=gsn_csm_contour_map_polar(wks,dif_ddt_New,res)
  plot_dif(6)=gsn_csm_contour_map_polar(wks,dif_ddf_Oliver,res)
  plot_dif(7)=gsn_csm_contour_map_polar(wks,dif_ddf_Feng,res)
  plot_dif(8)=gsn_csm_contour_map_polar(wks,dif_ddf_New,res)


  resP                     = True         ; modify the panel plot
  resP@gsnPanelTop   = 0.85   
  resP@gsnPanelBottom   = 0.15
  resP@gsnPanelYWhiteSpacePercent = 5
  resP@gsnPanelXWhiteSpacePercent = 5
  resP@lbLabelFontHeightF  = 0.007  
  resP@txString        = "Relative Error (RE)"
  gsn_panel(wks,plot_re,(/3,3/),resP) 
  resP@txString        = "Root Mean Square Error (RMSE)"
  gsn_panel(wks,plot_rmse,(/3,3/),resP) 
  resP@txString        = "Correlation (Cor)"
  gsn_panel(wks,plot_rr,(/3,3/),resP)
  resP@txString        = "Difference (Dif)"
  gsn_panel(wks,plot_dif,(/3,3/),resP)


exit


