import numpy as np

def ps_scale_for_lat(lat):
    '''
    This function calculates the length scaling factor for a polar stereographic
       projection (ie. SSM/I grid) which can be used to correct area calculations. 
       The scaling factor is defined (Snyder, 1982) as:
    
        k = (mc/m)*(t/tc), where
            
        m = cos(lat)/sqrt(1 - e2*sin(lat)^2)
        t = tan(Pi/4 - lat/2)/((1 - e*sin(lat))/(1 + e*sin(lat)))^(e/2)
        e2 = 0.006693883 is the earth eccentricity (Hughes ellipsoid)
        e = sqrt(e2)
        mc = m at the reference latitude (70 degrees)
        tc = t at the reference latitude (70 degrees)
            
        The ratio mc/tc is precalculated and stored in the variable mc_rc.
        See https://pubs/usgs.gov/pp/1395/report.pdf

        Areas calculated in polar stereographic need to be multiplied by
        k**2
        
        ref:Snyder, 1982, Map Projections used by the U.S. Geological Survey 
            Page 160
            https://pubs/usgs.gov/pp/1395/report.pdf  
    '''

    if lat > 0:
        hemisphere=1
    else:
        hemisphere=-1


    if ( hemisphere==1):
        #mc=0.3430355367529768
        #tc=0.17744189713382205
        mc_tc = 1.9332273960882438
    else:
        #mc/tc calculated at -71 degrees
        #mc=0.32654678137878085
        #tc=0.16840732452916296
        mc_tc=1.939029565915543 ; # for 71 deg, southern PS

    # for WGS84, a=6378137, 1/f = 298.257223563 -> 1-sqrt(1-e^2) = f
    # -> 1-(1-f)^2 = e2 =    0.006694379990141
    e2=0.006694379990141;  #from WGS84 parameters
    e = np.sqrt(e2);

    latr = np.abs(lat*np.pi/180.)
    if  isinstance(latr, (np.ndarray)):
        latr[lat==-90]=(-89.9999*np.pi/180.)
    slat = np.sin(latr)
    clat = np.cos(latr)

    m = clat/np.sqrt(1. - e2*slat**2)
    #print(m)
    t = np.tan(np.pi/4. - latr/2.)/((1.-e*slat) / (1.+e*slat))**(e/2.);
    #print(t)
    k = t/m *mc_tc

    scale=(1./k);
    return scale
    #
    # check: at S pole (this comes out right!)
    # if False:
    #     a=6378137.
    #     slat=np.sin(-90*np.pi/180)
    #     # radius of curvature
    #     Rc=a*(1-e2)/(1-e2*slat**2)**(3/2)
    #     L=Rc*.01*np.pi/180;
    #     A_e=L**2
    #     D=pc.data().from_dict({'latitude':np.array([-90, -89.99, -89.99]), 'longitude':np.array([0., 0., 90.])}).get_xy(EPSG=3031)
    #     A_ps = (D.x[2]-D.x[0])*(D.y[1]-D.y[0])
    #     print(A_e/A_ps)
    #     print(pc.ps_scale_for_lat(-89.9999)**2)
