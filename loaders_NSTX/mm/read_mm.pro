 pro read_mm
  read,'enter shot # ',ishot
  array=''
  dum1=''
  dummy=''
  file=''
  read,'enter choice of arrays (hf): ',array
; read,'enter filename: ',filename
  filename='config_'+array+'.mm'
  print,'opening ',filename
  openr,5,filename
  readf,5,nfile
  for i=1,nfile do begin
   readf,5,nshot,dum1,dummy,format='(i6,a2,a11)'
   if ishot ge nshot then begin
    file=dummy
   end
  end
  print,'configuration file is ',file
  openr,10,file
  readf,10,nch,ncha,format='(i2,1x,i3)'
  print,nch,ncha
  n=nch+ncha
  dev=strarr(n)
  for i=0,n-1 do begin
   readf,10,pos,sig,chn,device
   dev(i)=device
   print,pos,sig,chn,dev(i)
  end
 end
