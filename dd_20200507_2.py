""" Module to read ASDEX Upgrade shotfiles

Example:

sf = dd.shotfile()
if sf.Open(diagnostic,nshot):
   timebase = sf.GetTimebase(signal)
   signal = sf.GetSignal(signal)
   sf.Close() 
"""
print('my dd')
import os
import ctypes as ct
import numpy as np
import dics

# rdm 13.08.12: lbuf is ct.c_uint32, no longer ct.c_int32


ddlib = '/afs/ipp/aug/ads/lib64/@sys/libddww8.so.8.1'
if not os.path.isfile(ddlib):
    ddlib = '/afs/ipp/aug/ads/lib64/amd64_sles11/libddww8.so.8.1'

libddww = ct.cdll.LoadLibrary(ddlib)


def to_byt(str_or_byt):
    str_out = str_or_byt.encode('utf8') if not isinstance(str_or_byt, bytes) else str_or_byt
    return str_out

    
def GetError(error):
    try:
        err = ct.c_int32(error)
    except TypeError:
        err = ct.c_int32(error.value)
    isError = libddww.xxsev_(ct.byref(err))==1
    isWarning = libddww.xxwarn_(ct.byref(err))==1
    if isError or isWarning:
        id   = ct.c_char_p(b'')
        text = ct.c_char_p(b' '*255)
        unit = ct.byref(ct.c_int32(-1))
        ctrl = ct.byref(ct.c_uint32(3))
        lid   = ct.c_uint64(0)
        ltext = ct.c_uint64(255)
        libddww.xxerrprt_(unit, text, ct.byref(err), ctrl, id, ltext, lid);
        if isError:
            print(text.value.strip())


def LastShotNr():
    """ nshot = dd.LastShotNr """
    error  = ct.c_int(0)
    _error = ct.byref(error)
    cshot  = ct.c_uint32(0)
    _cshot = ct.byref(cshot)

    result = libddww.ddlastshotnr_(_error, _cshot)

    return cshot.value


def PreviousShot(diagname, nshot, experiment='AUGD'):
    """ nshot = dd.PreviousShot(diagname, nshot) """
    diagname = to_byt(diagname)
    experiment = to_byt(experiment)
    exp    = ct.c_char_p(experiment)
    diag   = ct.c_char_p(diagname)
    cshot  = ct.c_uint32(0)
    _cshot = ct.byref(cshot)
    shot   = ct.c_uint32(nshot)
    _shot  = ct.byref(shot)
    lexp  = ct.c_uint64(len(experiment))
    ldiag = ct.c_uint64(len(diagname))

    result = libddww.ddcshotnr_(exp, diag, _shot, _cshot, lexp, ldiag)

    GetError(result)
    return cshot.value


class ddarr(np.ndarray):

    def __new__(cls, input_array, info=None):

        obj = np.asarray(input_array).view(cls)
        for attr in ('units', 'tname', 'aname', 'status', 'fmt', \
            'objtyp', 'ind2', 'descr', 'rel', 'reltyp'):
            if hasattr(info, attr):
                obj.__dict__[attr] = info.__dict__[attr]
            else:
                obj.__dict__[attr] = None

        return obj

    def __array_finalize__(self, obj):

        if obj is None: return
        for attr in ('units', 'tsame', 'aname', 'status', 'fmt', 'objtyp'):
            self.__dict__[attr] = getattr(obj, attr, None)
          # We do not need to return anything


class dd_info:
    pass


class shotfile:
    """ Shotfile class to read ASDEX Upgrade data """

    def Open(self, diagname, nshot, experiment='AUGD', edition=0):
        """ Open a shotfile \n\n """
        if hasattr(self, '_diaref'):
            self.Close()
        diagname = to_byt(diagname)
        experiment = to_byt(experiment)
        error     = ct.c_int32(0)
        _error    = ct.byref(error)
        edit      = ct.c_int32(edition)
        _edition  = ct.byref(edit)
        self._diaref = ct.byref(ct.c_int32(0))
        cshot  = ct.c_uint32(0)
        shot   = ct.c_uint32(nshot)
        diag   = ct.c_char_p(diagname)
        exp    = ct.c_char_p(experiment)
        dat    = 18*b'd'
        date   = ct.c_char_p(dat)
        _cshot = ct.byref(cshot)
        _shot  = ct.byref(shot)
        lexp   = ct.c_uint64(len(experiment))
        ldiag  = ct.c_uint64(len(diagname))
        ldate  = ct.c_uint64(len(dat))

        result = libddww.ddopen_(_error, exp, diag, _shot, _edition, self._diaref, \
                                 date, lexp, ldiag, ldate)
        GetError(error)
        if result == 0:
            self.shot = nshot
            self.diag = diagname
            self.exp  = experiment 
            self.ed   = edit.value
            self.date = date.value
        else:
            del self._diaref
            GetError(error)
        return hasattr(self, '_diaref')


    def Close(self):
        """ Close a shotfile """
        if hasattr(self, '_diaref'):
            error   = ct.c_int32(0)
            _error  = ct.byref(error)

            result = libddww.ddclose_(_error, self._diaref)

            GetError(error)
            print('DDclose: edition %d closed ' %self.ed)
            del self._diaref


    def GetInfo(self, name):
        """ Returns information about the specified signal."""
# Assuming: not more than one TB, not more than one AB 
        if hasattr(self, '_diaref'):

            name = to_byt(name)
            output = dd_info()
            head = self.GetObjectHeader(name)
            output.error = head.error
            if head.error == 0:
                output.descr   = b''.join(self.GetObjectValue(name, b'text'))
                output.buf     = head.buffer
                output.objtyp  = output.buf[0]
                output.level   = output.buf[1]
                output.status  = output.buf[2]
                output.error   = output.buf[3]
                output.address = output.buf[12]
                output.bytlen  = output.buf[13]
                output.ndim    = output.buf[16]
# Relations
                output.rel = []
                output.reltyp = []
                for objid in output.buf[4:12]:
                    if objid != 65535:
                        relname = self.GetObjectName(objid)
                        output.rel.append( relname )
                        typ = self.GetObjectValue(relname, b'objtype')
                        output.reltyp.append(typ)
                        if typ == 8:
                            output.tname = relname
                        elif typ == 13:
                            output.aname = relname
                if output.objtyp in (6, 7, 8, 13):
                    output.units   = dics.unit_d[output.buf[15]]
                    output.estatus = output.buf[17]
                    output.fmt     = output.buf[14]
                    if output.objtyp == 8:
                        output.ind_len = output.buf[21]
                    else:
                        dims = np.array(output.buf[18:22][::-1], dtype=np.int32)
                        output.ind = dims[dims > 0]
                        output.ind2 = output.ind[:output.ndim]
                        output.ind_len = np.prod(output.ind)

            return output
        return None


    def GetParameterSetInfo( self , name ):
        """ Returns information about the specified parameter set."""
        output = dd_info()
        if hasattr(self, '_diaref'):
            name = to_byt(name)
            info  = self.GetObjectValue( name , b'items' )
            
            error    = ct.c_int32(0)
            _error   = ct.byref( error )
            par_name = ct.c_char_p(name)
            nrec     = ct.c_int32(info)
            _nrec    = ct.byref( nrec )
            rname    = ct.c_char_p(b' '*8*info)
            items    = (ct.c_uint32*info)()
            _items   = ct.byref(items)
            format   = (ct.c_uint32*info)()
            _format  = ct.byref( format )
            devsig   = (ct.c_int32*info)()
            _devsig  = ct.byref(devsig)
            lname  = ct.c_uint64( len(name) )
            lrname = ct.c_uint64( 8*info )

            result = libddww.ddprinfo_(_error, self._diaref, par_name, _nrec, rname, \
                                       _items, _format, _devsig, lname, lrname)
            output.error = error.value
            if error.value != 0:
                GetError( error )
            output.N_items = nrec.value
            output.names  = []
            output.items  = []
            output.format = []
            output.devsig = []
            for i in range(info):
                tmp = rname.value[8*i:8*(i+1)]
                if tmp.strip() != b'':
                    output.names.append(tmp)
            n_pars = len(output.names)
            for j in range(n_pars):
                output.items.append(items[j])
                output.format.append(format[j])
                output.devsig.append(devsig[j])
            
        return output


    def GetParameterInfo(self, set_name, par_name):
        """ Returns information about the parameter 'par_name' of the parameter set 'set_name'."""

        if hasattr(self, '_diaref'):
            par_name = to_byt(par_name)
            set_name = to_byt(set_name)
            output = dd_info()
            error   = ct.c_int32(0)
            pset    = ct.c_char_p(set_name)
            par     = ct.c_char_p(par_name)
            item    = ct.c_uint32(0)
            format  = ct.c_uint16(0)
            _error  = ct.byref(error)
            _item   = ct.byref(item)
            _format = ct.byref(format)
            lpar = ct.c_uint64(len(par_name))
            lset = ct.c_uint64(len(set_name))

            result = libddww.dd_prinfo_(_error, self._diaref, pset, par, _item, _format, lset, lpar)

            GetError(error)
            output.error = error.value
            if error.value == 0:
                output.item = item.value
                output.fmt  = format.value
            return output
        return None


    def GetParameter(self, set_name, par_name):
        """ Returns the value of the parameter 'par_name' of the parameter set 'set_name'. """
        if hasattr(self, '_diaref'):
            par_name = to_byt(par_name)
            set_name = to_byt(set_name)
            info = self.GetParameterInfo(set_name, par_name)
            if info.error == 0:
                error     = ct.c_int32(0)
                _error    = ct.byref(error)
                setn      = ct.c_char_p(set_name)
                lset      = ct.c_uint64(len(set_name))
                par       = ct.c_char_p(par_name)
                lpar      = ct.c_uint64(len(par_name))
                physunit  = ct.c_int32(0)
                _physunit = ct.byref(physunit)
# Characters
                if info.fmt in dics.fmt2len.keys():
                    ndim = dics.fmt2len[info.fmt]
                    nlen = ndim*info.item
                    typin  = ct.c_int32(6)
                    lbuf   = ct.c_uint32(nlen)
                    buffer = ct.c_char_p(b'd'*nlen)
                    _typin = ct.byref(typin)
                    _lbuf  = ct.byref(lbuf)
                    lndim  = ct.c_uint64(ndim)

                    dummy = libddww.ddparm_(_error, self._diaref, setn, par, _typin, \
                                             _lbuf, buffer, _physunit, lset, lpar, lndim)

                    GetError(error)
                    a=[]
                    for j in range(info.item):
                        a.append(buffer.value[j*ndim:(j+1)*ndim].decode('utf-8'))
                    result = np.array(a)
                else:
                    typin = ct.c_int32(dics.fmt2typ[info.fmt])
                    lbuf = ct.c_uint32(info.item)
                    buffer = (dics.fmt2ct[info.fmt]*info.item)()
                    _typin = ct.byref(typin)
                    _lbuf = ct.byref(lbuf)
                    _buffer = ct.byref(buffer)
                    dummy = libddww.ddparm_(_error, self._diaref, setn, par, _typin, \
                                             _lbuf, _buffer, _physunit, lset, lpar)
                    result = np.frombuffer(buffer, dtype=np.dtype(buffer))[0]

                return ddarr(result, info=info)

        return None


    def GetSignal(self, signame, nbeg=0, nend=np.infty, cal=False):
        """Returns the specified signal group and if specified performes a conversion to the specified type (e.g. ct.c_float)."""
        if hasattr(self, '_diaref'):
            signame = to_byt(signame)
            info = self.GetInfo(signame)
            if info.error != 0:
                print('Error getting SignalGroup %s' %signame)
                return None
            if info.status == -1:
                print('Signal status = -1')
                return None

            if cal:
                artyp = 2
            else:
                artyp = dics.fmt2typ[info.fmt]
           # if (not cal and info.fmt in (3, 4, 5)): # byte, integer, float
             #   return self.GetObjectData(signame,nbeg=nbeg,nend=nend)

            if artyp == 6:
                char_len = dics.fmt2len[info.fmt]
            else:
                char_len = 1

            k1       = ct.c_uint32(max(1, nbeg))
            _k1      = ct.byref(k1)
            k2       = ct.c_uint32(min(info.ind[0], nend))
            _k2      = ct.byref(k2)
            
            info.ind[0] = np.subtract(k2,k1)+1

            leng   = info.ind[0]*char_len
            buffer = (dics.typ2ct[artyp]*char_len*np.prod(info.ind))()
            _buffer  = ct.byref(buffer)
            error    = ct.c_int32(0)
            _error   = ct.byref(error)
            length   = ct.c_uint32(0)
            _length  = ct.byref(length)

            ctyp     = ct.c_uint32(artyp)
            _type    = ct.byref(ctyp)
            lbuf     = ct.c_uint32(leng)
            _lbuf    = ct.byref(lbuf)
            signam   = ct.c_char_p(signame)
            physdim  = 8*b'p'
            _physdim = ct.c_char_p(physdim)
            ncal     = ct.c_int32(0)
            _ncal    = ct.byref(ncal)
            lsig     = ct.c_uint64(len(signame))
            lphysdim = ct.c_uint64(len(physdim))

            if info.objtyp == 7:
# Signal
# Calibrated Signal
                if cal:
 
                    dummy = libddww.ddccsgnl_(_error, self._diaref, signam, \
                            _k1, _k2, _type, _lbuf, _buffer, _length, _ncal, \
                            _physdim, lsig, lphysdim)
                    GetError(error)
                    if error.value != 0:
                        if error.value == 555352323:
                            print('No calibrated data, returning uncalibrated signal')
                            return self.GetObjectData(signame,nbeg=nbeg,nend=nend)
                else: # char signal
                    dummy = libddww.ddsignal_(_error, self._diaref, signam, \
                            _k1, _k2, _type, _lbuf, _buffer, _length, lsig)

                result = np.frombuffer(buffer, dtype=np.dtype(buffer)).ravel()

            elif info.objtyp == 6:
# SignalGroup

                if cal:
                    dummy = libddww.ddccsgrp_(_error, self._diaref, signam, \
                            _k1, _k2, _type, _lbuf, _buffer, _length, \
                            _ncal, _physdim, lsig, lphysdim)
                    GetError(error)
                    if error.value != 0:
                        if error.value == 556204039:
                            print('No calibrated data, returning uncalibrated signal')
                            return self.GetObjectData(signame,nbeg=nbeg,nend=nend)
                else: # char array
                    dummy = libddww.ddsgroup_(_error, self._diaref, signam, \
                            _k1, _k2, _type, _lbuf, _buffer, _length, lsig)
                    GetError(error)

                result = np.frombuffer(buffer, dtype=np.dtype(buffer))
                result = result.reshape(info.ind2, order='F')

            return ddarr(result, info=info)


    def GetAreabase(self, signame, nbeg=0, nend=np.infty):
        """ Returns the areabase to the specified signal, time is always the first independent variable."""
        if hasattr(self, '_diaref'):
            info = self.GetInfo(signame)
            if info.error == 0:
                if info.objtyp == 13:
                    return self.GetObjectData(signame,nbeg=nbeg,nend=nend)
                elif hasattr(info, 'aname'):
                    return self.GetObjectData(info.aname,nbeg=nbeg,nend=nend)
        return None


    def GetTimebase(self, signame, nbeg=0, nend=np.infty, cal=False):
        """Returns the timebase corresponding to the specified signal."""
        if hasattr(self, '_diaref'):
            signame = to_byt(signame)
            info = self.GetInfo(signame)
            if info.error != 0:
                return None
            if info.objtyp != 8:
                if hasattr(info, 'tname'):
                    return self.GetTimebase(info.tname, cal=cal)
                else:
                    return None
            else:
                if info.bytlen > 0 and info.fmt == 1005:
                    return self.GetObjectData(signame,nbeg=nbeg,nend=nend)
                else:
                    if cal:
                        artyp = 2
                    else:
                        artyp = dics.fmt2typ[info.fmt]
                    error   = ct.c_int32(0)
                    _error  = ct.byref(error)
                    lbuf  = info.ind_len
                    cfmt = dics.typ2ct[artyp]
                    signam  = ct.c_char_p(signame)
                    k1      = ct.c_uint32(max(1, nbeg))
                    _k1     = ct.byref(k1)
                    k2      = ct.c_uint32(min(info.ind_len, nend))
                    _k2     = ct.byref(k2)
                    ctyp    = ct.c_uint32(artyp)
                    _type   = ct.byref(ctyp)
                    clbuf   = ct.c_uint32(lbuf)
                    _lbuf   = ct.byref(clbuf)
                    cleng   = ct.c_uint32(0)
                    _leng   = ct.byref(cleng)
                    buffer = (cfmt*info.ind_len)()
                    _buffer = ct.byref(buffer)
                    lname   = ct.c_uint64(len(signame))
                    dummy = libddww.ddtbase_(_error, self._diaref, signam, \
                            _k1, _k2, _type, _lbuf, _buffer, _leng, lname)
                    tmp = np.frombuffer(buffer, dtype=np.dtype(buffer))
                    result = np.atleast_1d(np.squeeze(tmp))
                    return ddarr(result, info=info)


    def GetObjectData(self, obj_name, nbeg=0, nend=np.infty):
        """dd.shotfile.GetObjectData(ObjectName, fmt=ct.c_float)\n\nReturns the data part of the specified object and interpretes it as an array of the datatype specified in fmt"""

        if hasattr(self, '_diaref'):
            obj_name = to_byt(obj_name)
            info = self.GetInfo(obj_name)
            if info.error != 0:
                print('Error getting Object %s' %obj_name)
                return None
            if info.status == -1:
                return None

            cfmt = dics.fmt2ct[info.fmt]
            lbuf = info.bytlen

            error   = ct.c_int32(0)
            _error  = ct.byref(error)
            cname   = ct.c_char_p(obj_name)
            clbuf   = ct.c_uint32(lbuf)
            _lbuf   = ct.byref(clbuf)
            buffer  = (ct.c_byte*lbuf)()
            _buffer = ct.byref(buffer)
            cleng   = ct.c_uint32(0)
            _leng   = ct.byref(cleng)
            lname   = ct.c_uint64(len(obj_name) )

            dummy = libddww.ddobjdata_(_error, self._diaref, cname, \
                    _lbuf, _buffer, _leng, lname)
            result = np.frombuffer((cfmt*info.ind_len).from_buffer(buffer),dtype=np.dtype(cfmt))

            GetError(error)

            if info.objtyp == 6:
                result = result.reshape(info.ind2, order='F')
            elif info.objtyp == 13:
                result = result.reshape(info.ind)
            return ddarr(result, info=info)
        return None


    def GetListInfo(self):
        if hasattr(self, '_diaref'):
            error   = ct.c_int32(0)
            _error  = ct.byref(error)
            lbuf    = ct.c_uint32(255)
            _lbuf   = ct.byref(lbuf)
            buf     = (ct.c_int32*255)()
            _buf    = ct.byref(buf)
            length  = ct.c_int32(0)
            _length = ct.byref(length) 
            result  = libddww.ddflist_(_error, self._diaref, _lbuf, _buf, _length)
            GetError(error)
            if length.value != 0:
                return np.int(buf[0:np.int(length)])
        return None


    def GetObjectValue(self, name, field):
        """dd.shotfile.GetObjectValue(Name, Field)\n\nReturns the value specified in Field of the object Name."""
        if hasattr(self, '_diaref'):
            name = to_byt(name)
            field = to_byt(field)
            ovalue = ct.c_int32(0)
            if field in (b'relations', b'special'):
                ovalue = (ct.c_int32*8)()
            if field in (b'format', b'indices'):
                ovalue = (ct.c_int32*3)()
            if field == b'dataformat':
                ovalue = ct.c_uint16(0)
            if field == b'text':
                ovalue = (ct.c_char*128)()
            _value  = ct.byref(ovalue)
            error   = ct.c_int32(0)
            _error  = ct.byref(error)
            _name   = ct.c_char_p(name)
            _field  = ct.c_char_p(field)
            lname   = ct.c_uint64(len(name))
            lfield  = ct.c_uint64(len(field))

            result = libddww.ddobjval_(_error, self._diaref, _name, _field, _value, \
                                       lname, lfield)

            GetError(error)
            if np.size(ovalue) == 1:
                return ovalue.value
            else:
                return np.frombuffer(ovalue, dtype=np.dtype(ovalue))[0]
        return None


    def GetObjectHeader(self, name):
        if hasattr(self, '_diaref'):
            name = to_byt(name)
            output = dd_info()
            text = 64*b't'
            error   = ct.c_int32(0)
            _error  = ct.byref(error)
            _name   = ct.c_char_p(name)
            buffer  = (ct.c_int32*26)()
            _buffer = ct.byref(buffer)
            _text   = ct.c_char_p(text)
            lname = ct.c_uint64(len(name))
            ltext = ct.c_uint64(len(text))

            result = libddww.ddobjhdr_(_error, self._diaref, _name, _buffer, _text, lname, ltext)

            GetError(error)
            output.error = error.value
            if error.value == 0:
                output.buffer = buffer[0:26]
                output.text   = _text.value
            return output
        return None


    def GetObjectName(self, obj):
        if hasattr(self, '_diaref'):
            name = 8*b'1'
            error   = ct.c_int32(0)
            _error  = ct.byref(error)
            _name   = ct.c_char_p(name)
            obje    = ct.c_int32(obj)
            _obje   = ct.byref(obje)
            lname = ct.c_uint64(len(name))

            result = libddww.ddobjname_(_error, self._diaref, _obje, _name, lname)

            if error.value == 0:
                return _name.value.replace(b' ',b'')
        return -1


    def GetNames(self):
        if hasattr(self, '_diaref'):
            i = 1
            ok = True
            result = []
            while self.GetObjectName(i) != -1:
                name = self.GetObjectName(i)
                name = name.decode('utf8') if isinstance(name, bytes) else name
                result.append(name)
                i += 1
            return result
        else:
            return None
