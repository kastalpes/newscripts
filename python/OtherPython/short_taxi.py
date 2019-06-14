from GL import *
import yt
yt.enable_parallelism()
class short_taxi():
    """This is a package of information that various tools need
    about a simulation"""
    def __init__(self,directory=".",name="NAME",outname=None,frames='all',fields=['density']):
        self.directory=directory
        self.name=name
        self.frames=frames
        self.fields=fields
        if outname is None:
            self.outname = name
        else:
            self.outname = outname
    def return_frames(self):
        return self.frames
    def load(self,frame):
        self.ds = yt.load("%s/DD%04d/data%04d"%(self.directory,frame,frame))
        return self.ds
    def get_region(self,frame):
        self.load(frame)
        self.reg = self.ds.all_data()
        return self.reg
