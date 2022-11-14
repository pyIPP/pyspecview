pro mm_rd,nshot,arr,tree,na,nch,chn,sig,cnam,tor,pol,ntor,nfp,npol,oe1,oe2
; INPUTS:
;  nshot	shot number needed to pick configuration file
;  arr		for now either hf or hn
; OUTPUTS:
;  na		na as per integrate.pro
;  nch		# of active channels
;  chn		list of channel #'s
;  sig		sign of each channel
;  cnam		list of signal names
;  tor		list of toroidal angles
;  pol		list of poloidal angles
;  ntor		# of coils in toroidal array
;  nfp		# of coils oriented toroidally
;  npol		# of coils in poloidal array
 if arr eq 'hf' then begin
; returns list of channel #, names, etc
  openr,1,'/u/eric/nstx/mm/config_hf.mm'
; na=0.00457		;(m^2), see integrate.pro)
  readf,1,nfile
  shots=intarr(nfile)
  filenm=strarr(nfile)
  file=''
  for i=0,nfile-1 do begin
   readf,1,ishot,file,format='(i6,2x,a11)'; INPUTS:
;  nshot	shot number needed to pick configuration file
;  arr		for now either hf or hn
; OUTPUTS:
;  na		na as per integrate.pro
;  nch		# of active channels
;  chn		list of channel #'s
;  sig		sign of each channel
;  cnam		list of signal names
;  tor		list of toroidal angles
;  pol		list of poloidal angles
;  ntor		# of coils in toroidal array
;  nfp		# of coils oriented toroidally
;  npol		# of coils in poloidal array

;   print,long(ishot),' ',file
   if nshot ge ishot then ifi=i
   shots(i)=long(ishot)
   filenm(i)=file
  end
  close,1
; print out shot year
;  if nshot ge 100621 and nshot le 101616 then print,'year 1999, ',filenm(ifi)
;  if nshot ge 101796 and nshot le 104511 then print,'year 2000, ',filenm(ifi)
;  if nshot ge 104760 and nshot le 106490 then print,'year 2001, ',filenm(ifi)
;  if nshot ge 106766 and nshot le 109079 then print,'year 2002, ',filenm(ifi)
;  if nshot ge 109397 and nshot le 110201 then print,'year 2003, ',filenm(ifi)
;  if nshot ge 110810 and nshot le 114478 then print,'year 2004, ',filenm(ifi)
;  if nshot ge 115151 and nshot le 118162 then print,'year 2005, ',filenm(ifi)
;  if nshot ge 118929 and nshot le 121563 then print,'year 2006, ',filenm(ifi)
;  if nshot ge 122270 and nshot le 125345 then print,'year 2007, ',filenm(ifi)
;  if nshot ge 126454 and nshot le 130737 then print,'year 2008, ',filenm(ifi)
;  if nshot ge 131634 and nshot le 136160 then print,'year 2009, ',filenm(ifi)
;  if nshot ge 137110 and nshot le 142524 then print,'year 2010, ',filenm(ifi)
;  if nshot ge 143000 and nshot le 200000 then print,'year 2014, ',filenm(ifi)
  tree=''
  openr,2,'/u/eric/nstx/mm/'+filenm(ifi)
  readf,2,ntor,nfp,npol,format='(i3,i3,i3)'
;  print,ntor,nfp,npol
  readf,2,na
;  print,'na ',na
  readf,2,tree,format='(a)'
;  print,tree
  nch=ntor+nfp+npol
  chn=intarr(nch)
  cnam=strarr(nch)
  sig=fltarr(nch)
  tor=fltarr(nch)
  pol=fltarr(nch)
  point=''
  for i=0,nch-1 do begin
   if nshot lt 200000 then $
   readf,2,to,po,sign,chan,point,format='(1x,f6.1,f7.2,1x,f4.1,1x,i2,2x,a)'
   if nshot ge 200000 then $
   readf,2,to,po,sign,chan,point,format='(1x,f6.1,f7.2,1x,f5.2,1x,i2,2x,a)'
   chn(i)=chan
   cnam(i)=point
   sig(i)=sign
   tor(i)=to
   pol(i)=po
;   print,to,po,sign,chan,point
  end
  readf,2,oe1,oe2,format='(2x,i2,1x,i2)'
  close,2
 end
;no PC1 data 137591-137621
;no poloidal array data 137640-139121
;all channels 139211 - 139640
;no data 139641 - 139891
;only two cards 139892 - 139948
;all data 139949 - 139967 (no data 139951 - 139966)
;lost last card (not poloidal data?) 139973 - 140530
;
;for hn array:
;up through 2002 (101756-109079) hn data was on 907's
;109397 on it was 4MHz pci data
;in 2001 and 2002, there was also 6810 data
;_______________________________________________________________________________
;
 if arr eq 'hn' then begin
; returns list of channel #, names, etc
  nfp=0
  npol=0
  openr,1,'/u/eric/nstx/mm/config_hn.mm'
; na=0.00457		;(m^2), see integrate.pro)
  readf,1,nfile
  shots=intarr(nfile)
  filenm=strarr(nfile)
  file=''
  for i=0,nfile-1 do begin
   readf,1,ishot,file,format='(i6,2x,a11)'
;   print,long(ishot),' ',file
   if nshot ge ishot then ifi=i
   shots(i)=long(ishot)
   filenm(i)=file
  end
  close,1
; print out shot year
  tree=''
  openr,2,'/u/eric/nstx/mm/'+filenm(ifi)
  readf,2,ntor,format='(i3)'
  print,'ntor ',ntor
  readf,2,na
;  print,'na ',na
  readf,2,tree,format='(a)'
  print,tree
  nch=ntor+nfp+npol
  chn=intarr(nch)
  cnam=strarr(nch)
  sig=fltarr(nch)
  tor=fltarr(nch)
  pol=fltarr(nch)
  pol(*)=-33.47
  point=''
  for i=0,nch-1 do begin
   readf,2,to,sign,chan,point,format='(1x,f6.1,1x,f4.1,2x,i2,2x,a21)'
   chn(i)=chan
   cnam(i)=point
   sig(i)=sign
   tor(i)=to
;   print,to,sign,chan,point
  end
  readf,2,oe1,oe2,format='(2x,i2,1x,i2)'
  print,oe1,oe2
  close,2
 end
end
