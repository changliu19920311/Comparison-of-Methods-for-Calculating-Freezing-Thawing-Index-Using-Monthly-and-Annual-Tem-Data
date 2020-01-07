; Program for computing the frost index
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"



f    = addfile ("GLDAS_MC_index_ensemble.nc", "r")   
ddf=f->ddf
ddt=f->ddt
fi=f->fi
year=f->year

do i=0,10
year1=1901+i*10
year2=1901+i*10+9
if(i.eq.10)
year2=1901+i*10+11
end if
Year=ispan(year1,year2,1)
iyyyy=ind(year.ge.year1.and.year.le.year2)
   system("/bin/rm GLDAS_MC_index_"+tostring(year1)+"_"+tostring(year2)+".nc")
   fo  = addfile("GLDAS_MC_index_"+tostring(year1)+"_"+tostring(year2)+".nc","c")         
  setfileoption("nc","Format","LargeFile") 
   fo->fi  = fi(iyyyy,:,:)
   fo->ddt  = ddt(iyyyy,:,:)
   fo->ddf  = ddf(iyyyy,:,:)
print(year1)
print(year2)
delete(Year)
delete(iyyyy)
end do




  end


  





