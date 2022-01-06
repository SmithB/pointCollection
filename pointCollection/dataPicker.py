#! /usr/bin/env python

from scipy.interpolate import interpn
import numpy as np
import matplotlib.pyplot as plt
import pointCollection as pc
import time

class dataPicker(object):
    def __init__(self, pt_data, scatter_args=None, 
                 img_data=None, img_args=None, 
                 field='z', 
                 handles=None,  
                 correction=None,
                 plot_fn=None,
                 file_args=None, W=1.e4):
        if handles is None:
            handles = self.__init_ui__(pt_data, scatter_args, img_data, img_args)
        self.last_click_time=0.0
        self.max_click_time = 0.2
        self.D=pt_data
        self.handles=handles
        self.field=field
        self.messages=[[]]
        self.last_pt=[[]]
        self.last_data=None
        self.correction=correction
        if plot_fn is None:
            self.plot_fn=default_data_picker_plot_function
        self.x_atc_ctr=np.NaN
        if file_args is None:
            self.file_args={}
        else:
            self.file_args=file_args
        self.W=W
        self.cid=self.handles['fig'].canvas.mpl_connect('button_press_event', self.buttondown)
        self.cid=self.handles['fig'].canvas.mpl_connect('button_release_event', self.buttonup)
     
    def __init_ui__(pt_data,scatter_args, img_data, img_args, field='z'):
        fig=plt.figure()
        hax=fig.subplots(1,2)
        handles={'fig':fig, 'map_ax':hax[0], 'plot_ax':hax[1]}
        
        if img_data is not None:
            img_data.show(ax=handles['map_ax'], **img_args)
        
        if pt_data is not None:
            hax[0].scatter(pt_data.x, pt_data.y, 2, \
                           c=getattr(pt_data,field), **scatter_args)
        
        return handles
    
    def buttondown(self, event):
        # We'll make sure there's a short time between button down and up, so that we don't 
        # catch mouse zooms
        if not event.inaxes in [self.handles['map_ax']]:
            return
        self.last_click_time=time.time()
    
    def buttonup(self, event):
        # buttonup indicates picked points
        try:
            if not event.inaxes in [self.handles['map_ax']]:
                self.messages += ['tile_picker: last point not in map axis']
                return
            dt_click = time.time()-self.last_click_time
            if time.time()-self.last_click_time > self.max_click_time:
                self.messages += [f'too much time has elapsed : {dt_click}']
                return
            xy0=(event.xdata, event.ydata)
            tx = 'xy =[%f,%f]' % xy0
            self.handles['plot_ax'].set_title(tx)


            dist=np.abs((self.D.x-xy0[0])+1j*(self.D.y-xy0[1]))
            Dsub=self.D[dist<self.W/2]
            self.last_data=Dsub            
            self.plot_fn(xy0, Dsub, self.handles['plot_ax'])

            
        except Exception as e:
            self.messages += [e]
            plt.gca().set_title('ERROR')
        self.handles['plot_ax'].figure.canvas.draw()
    
    def clear_lines(self):
        lines=list(self.handles['plot_ax'].collections)
        for num in range(len(list(self.handles['plot_ax'].collections))):
            self.handles['plot_ax'].collections.pop(0)
        for num in range(len(list(self.handles['plot_ax'].lines))):
            self.handles['plot_ax'].lines.pop(0)
        self.handles['plot_ax'].figure.canvas.draw()

def default_data_picker_plot_function(xy0, D, hax):

    #hax.cla()
    hax.plot(D.time, D.z-D.dem_h, '.')

    plt.gca().set_xlim([np.nanmin(D.time), np.nanmax(D.time)])
    plt.gca().set_ylim([np.nanmin(D.z-D.dem_h), np.nanmax(D.z-D.dem_h)])
    

    
