# enable local tools
import os
from math import pi,floor,atan2,sin,cos,sqrt
import pylab

try:
    import casatools
    from casatasks import casalog
    
    mytb=casatools.table()
    myme=casatools.measures()
    mymd=casatools.msmetadata()
except ImportError:
    from taskinit import *

    mytb=tbtool()
    myme=metool()
    mymd=msmdtool()

def polfromgain(vis,tablein,caltable,paoffset):

    casalog.origin('polfromgain')

    casalog.post("Deriving calibrator linear polarization from gain ratios.")

    try:

        if ((type(vis)==str) & (os.path.exists(vis))):
            mymd.open(vis)
        else:
            raise Exception('Visibility data set not found - please verify the name')

        nant=mymd.nantennas()
        nfld=mymd.nfields()
        fldnames=mymd.fieldnames()
        nspw=mymd.nspw()

        rempol=False
        if ((type(tablein)==str) & (os.path.exists(tablein))):
            if type(caltable)==str and len(caltable)>0:

                if os.path.exists(caltable):
                    raise Exception('Output caltable='+caltable+' exists.  Choose another name or delete it.')

                casalog.post("New caltable, "+caltable+", corrected for linear polarization, will be generated.")
                mytb.open(tablein)
                myout=mytb.copy(newtablename=caltable,deep=True)
                mytb.close()
                myout.close()
                rempol=True
            else:
                casalog.post("No new caltable will be generated")
                caltable=tablein
        else:
            raise Exception('input calibration table not found - please verify the name')
            

        if paoffset!=0.0:
            casalog.post("NB: default band position angle will be offset by "+str(paoffset)+"deg.")

        # Field coords
        mytb.open(caltable+'/FIELD')
        dirs=mytb.getcol('DELAY_DIR')[:,0,:]
        mytb.close()
        
        # Must retrieve nominal feed angles from MS.FEED!
        mytb.open(vis+'/FEED')
        nfeed=mytb.nrows()
        fang=mytb.getcol('RECEPTOR_ANGLE')
        fspw=mytb.getcol('SPECTRAL_WINDOW_ID')
        fant=mytb.getcol('ANTENNA_ID')
        rang=pylab.zeros((nant,nspw));
        for ifeed in range(nfeed):
            rang[fant[ifeed],fspw[ifeed]]=fang[0,ifeed]
        mytb.close()

        R=pylab.zeros((nspw,nfld))
        Q=pylab.zeros((nspw,nfld))
        U=pylab.zeros((nspw,nfld))
        mask=pylab.zeros((nspw,nfld),dtype=bool)

        IQUV={}
        nomod=not rempol
        mytb.open(caltable,nomodify=nomod)
        uflds=pylab.unique(mytb.getcol('FIELD_ID'))
        uspws=pylab.unique(mytb.getcol('SPECTRAL_WINDOW_ID'))
        for ifld in uflds:
            rah=dirs[0,ifld]*12.0/pi
            decr=dirs[1,ifld]
            IQUV[fldnames[ifld]]={}
            for ispw in uspws:

                r=pylab.zeros(nant)
                q=pylab.zeros(nant)
                u=pylab.zeros(nant)
                antok=pylab.zeros(nant,dtype=bool)

                for iant in range(nant):
                    qstring='FIELD_ID=='+str(ifld)+' && SPECTRAL_WINDOW_ID=='+str(ispw)+' && ANTENNA1=='+str(iant)
                    st=mytb.query(query=qstring)
                    nrows=st.nrows()
                    if nrows > 0:

                        times=st.getcol('TIME')
                        gains=st.getcol('CPARAM')
                        flags=st.getcol('FLAG')
                        flags=pylab.logical_or(flags[0,0,:],flags[1,0,:])  # 1D

                        # Escape if insufficient data
                        if (nrows-pylab.sum(flags))<3:
                            antok[iant]=False
                            st.close()
                            continue


                        # parang
                        parang=pylab.zeros(len(times))
                
                        apos=mymd.antennaposition(iant)
                        latr=myme.measure(apos,'WGS84')['m1']['value']
                        myme.doframe(apos)
                        har=pylab.zeros(nrows)
                        for itim in range(len(times)):
                            tm=myme.epoch('UTC',str(times[itim])+'s')
                            last=myme.measure(tm,'LAST')['m0']['value']
                            last-=floor(last)  # days
                            last*=24.0  # hours
                            ha=last-rah  # hours
                            har[itim]=ha*2.0*pi/24.0  # radians
                    
                        parang=pylab.arctan2( (cos(latr)*pylab.sin(har)),
                                             (sin(latr)*cos(decr)-cos(latr)*sin(decr)*pylab.cos(har)) )

                        parang+=rang[iant,ispw]
                        parang+=(paoffset*pi/180.)       # manual feed pa offset
                    
                        # indep var matrix
                        A=pylab.ones((nrows,3))
                        A[:,1]=pylab.cos(2*parang)
                        A[:,2]=pylab.sin(2*parang)
                        A[flags,:]=0.0  # zero flagged rows
                    
                        # squared gain amplitude ratio
                        amps=pylab.absolute(gains)
                        amps[amps==0.0]=1.0
                        gratio2=pylab.square(amps[0,0,:]/amps[1,0,:])
                        gratio2[flags]=0.0  # zero flagged samples
                
                        fit=pylab.lstsq(A,gratio2)

                        r[iant]=fit[0][0]
                        q[iant]=fit[0][1]/r[iant]/2.0
                        u[iant]=fit[0][2]/r[iant]/2.0
                        r[iant]=sqrt(r[iant])
                        p=sqrt(q[iant]**2+u[iant]**2)
                        x=0.5*atan2(u[iant],q[iant])*180/pi
                    
                        antok[iant]=True;

                        #print 'Fld='+fldnames[ifld],'Spw='+str(ispw),'Ant='+str(iant), '(PA offset='+str(rang[iant,ispw]*180/pi+paoffset)+'deg)','Gx/Gy='+str(r[iant]),'q='+str(q[iant]),'u='+str(u[iant]),'p='+str(p),'x='+str(x)
                        casalog.post('Fld='+fldnames[ifld]+' Spw='+str(ispw)+' Ant='+str(iant)+' (PA offset='+str(rang[iant,ispw]*180/pi+paoffset)+'deg)'+' q='+str(q[iant])+' u='+str(u[iant])+' p='+str(p)+' x='+str(x)+' Gx/Gy='+str(sqrt(r[iant])))

                        if rempol:
                            if p<1.0:
                                Qpsi=q[iant]*pylab.cos(2*parang) + u[iant]*pylab.sin(2*parang)
                                gains[0,0,:]/=pylab.sqrt(1.0+Qpsi)
                                gains[1,0,:]/=pylab.sqrt(1.0-Qpsi)
                                st.putcol('CPARAM',gains)
                            else:
                                st.close()
                                raise Exception('Spurious fractional polarization!')

                    st.close()

                nantok=pylab.sum(antok)
                if nantok>0:
                    Q[ispw,ifld]=pylab.sum(q)/nantok
                    U[ispw,ifld]=pylab.sum(u)/nantok
                    R[ispw,ifld]=pylab.sum(r)/nantok
                    mask[ispw,ifld]=True
                
                    P=sqrt(Q[ispw,ifld]**2+U[ispw,ifld]**2)
                    X=0.5*atan2(U[ispw,ifld],Q[ispw,ifld])*180/pi

                #print 'Fld='+fldnames[ifld],'Spw='+str(ispw),'Ant=*', '(PA offset='+str(rang[iant,ispw]*180/pi+paoffset)+'deg)','Gx/Gy='+str(R[ispw,ifld]),'Q='+str(Q[ispw,ifld]),'U='+str(U[ispw,ifld]),'P='+str(P),'X='+str(X)

                    casalog.post('Fld='+fldnames[ifld]+' Spw='+str(ispw)+' Ant=*'+' (PA offset='+str(rang[iant,ispw]*180/pi+paoffset)+'deg)'+' Q='+str(Q[ispw,ifld])+' U='+str(U[ispw,ifld])+' P='+str(P)+' X='+str(X))

                    IQUV[fldnames[ifld]]['Spw'+str(ispw)]=[1.0,Q[ispw,ifld],U[ispw,ifld],0.0]

                else:
                    mask[ispw,ifld]=False

            if sum(mask[:,ifld])>0:
                casalog.post('For field='+fldnames[ifld]+' there are '+str(sum(mask[:,ifld]))+' good spws.')
                Qm=pylab.mean(Q[mask[:,ifld],ifld])
                Um=pylab.mean(U[mask[:,ifld],ifld])
                IQUV[fldnames[ifld]]['SpwAve']=[1.0,Qm,Um,0.0]
                Qe=pylab.std(Q[mask[:,ifld],ifld])
                Ue=pylab.std(U[mask[:,ifld],ifld])
                Pm=sqrt(Qm**2+Um**2)
                Xm=0.5*atan2(Um,Qm)*180/pi
                casalog.post('Spw mean: Fld='+fldnames[ifld]+' Q='+str(Qm)+' U='+str(Um)+' P='+str(Pm)+' X='+str(Xm))
                
        mytb.close()
        mymd.close()

        casalog.post("NB: Returning dictionary containing fractional Stokes results.")
        return IQUV

    except Exception as instance:
        print('*** Error ***',instance)
        mytb.close()
        mymd.close()
        raise

