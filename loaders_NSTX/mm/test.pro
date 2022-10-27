pro test
; shell to test mm_rd
 uread
 nshot=139048l
 rdi,'enter a shot #: ',nshot
 arr='t'
 rdc,'enter array (hf,hn): ',arr
 mm_rd,nshot,arr,na,nch,chn,sig,cnam,tor,pol,ntor,nfp,npol
 print,ntor,nfp,npol
 for i=0,nch-1 do begin
  print,i+1,chn(i),cnam(i),tor(i),pol(i),sig(i)
 end
end
