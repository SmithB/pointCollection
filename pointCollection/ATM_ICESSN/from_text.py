import re
import datetime
import csv
import numpy as np
import pointCollection as pc


def get_ATM_time(self):
    m=re.search(r"ATM.*_(\d{4})(\d{2})(\d{2})_(\d{2})(\d{2})(\d{2}).*.csv", self.filename)
    year_offset=0
    if m is None:
        m=re.search(r"ATM.*_(\d{2})(\d{2})(\d{2})_(\d{2})(\d{2})(\d{2}).*", self.filename)
        year_offset=2000
    m=[int(m.group(ind+1)) for ind in range(3)]
    m[0] += year_offset
    t0=datetime.datetime(*m)
    # The time is relative to the day of the granule -> discard h, m, s
    #datetime.timedelta(hours=this_time[3], minutes=this_time[4], seconds=this_time[5])
    t0 -= datetime.datetime(2000, 1, 1, 0, 0, 0)
    t0=np.float64(t0.days)+np.float64(t0.seconds)/24./3600
    self.assign(days_J2k = t0 + self.utc_seconds_of_day.astype(np.float64)/24./3600.)

def from_text(self, thefile, EPSG=None, format='csv'):
    ITRF_re=re.compile('(ITRF\d\d)')
    lines=[]
    ITRF=None
    try:
        with open(thefile,'r') as fh:
            format='csv'
            csvreader = csv.reader(fh)
            for line in csvreader:
                if len(line)==0:
                    continue
                if line[0][0]=='#':
                    last_comment=line
                    if 'ITRF' in line[-1]:
                        ITRF=ITRF_re.search(line[-1]).group(1)
                else:
                   lines += [[*map(float, line)]]
    except Exception:
        format='ssv'
        with open(thefile,'r') as fh:
            last_comment = '# UTC_Seconds_Of_Day, Latitude(deg), Longitude(deg), WGS84_Ellipsoid_Height(m), South-to-North_Slope, West-to-East_Slope, RMS_Fit(cm), Number_Of_ATM_Measurments_Used, Number_Of_ATM_Measurements_Removed, Distance_Of_Block_To_The_Right_Of_Aircraft(m), Track_Identifier'
            last_comment = last_comment.replace(',','').replace('# ','').split()
            for line in fh:
                try:
                    lines += [[*map(float, line.rstrip().split())]]
                except Exception:
                    print(f'problem with {thefile}:\n{line}')
                    pass
    lines=np.c_[lines]
    fields={}
    data={}

    for col, field  in enumerate(last_comment):
        field=field.replace('# ','')
        short_name=field.split('(')[0].replace('-','_').lower().lstrip()
        fields[short_name]={}
        if '(cm)' in field:
            fields[short_name]['scale']=0.01
        else:
            fields[short_name]['scale']=1
        fields[short_name]['long_name']=field
        data[short_name]=lines[:, col]
    self=pc.data().from_dict(data)
    self.filename = thefile

    if EPSG is not None:
        self.get_xy(EPSG=EPSG)
        self.EPSG=EPSG
    self.attrs['ITRF']=ITRF

    get_ATM_time(self)
    # convert to vanilla format
    self.assign(z=self.wgs84_ellipsoid_height,
                time=2000+self.days_J2k/365.25,
                sigma=self.rms_fit)
    return self
