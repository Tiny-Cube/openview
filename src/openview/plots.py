import numpy as np
from PySide6.QtCore import Qt, QPointF
from PySide6.QtGui import QPainter, QColor, QPen, QPainterPath
from PySide6.QtWidgets import QWidget

class Plot(QWidget):
    def __init__(self,title,parent=None):
        super().__init__(parent); self.title=title; self.x=np.array([]); self.y=np.array([])
        self.xlabel=''; self.ylabel=''; self.setMinimumHeight(180)

    def set_data(self,x,y,xlabel='',ylabel=''):
        self.x=np.asarray(x); self.y=np.asarray(y); self.xlabel=xlabel; self.ylabel=ylabel; self.update()

    def paintEvent(self,event):
        p=QPainter(self); p.setRenderHint(QPainter.RenderHint.Antialiasing)
        p.fillRect(self.rect(),QColor('#ffffff')); p.setPen(QColor('#203047')); p.drawText(15,25,self.title)
        if not len(self.x) or len(self.x)!=len(self.y):
            p.setPen(QColor('#52647a')); p.drawText(self.rect(),Qt.AlignmentFlag.AlignCenter,'No data in this file'); return
        left,top,w,h=65,45,max(1,self.width()-90),max(1,self.height()-100)
        xmin,xmax=float(self.x.min()),float(self.x.max()); ymin,ymax=float(self.y.min()),float(self.y.max())
        if xmax==xmin: xmax=xmin+1
        if ymax==ymin: ymax=ymin+1e-6
        pad=(ymax-ymin)*.08; ymin-=pad; ymax+=pad
        p.setPen(QPen(QColor('#dce4ed'),1))
        for i in range(5):
            y=top+h*i/4; p.drawLine(QPointF(left,y),QPointF(left+w,y))
            p.setPen(QColor('#52647a')); p.drawText(3,int(y+4),f'{ymax-(ymax-ymin)*i/4:.5g}'); p.setPen(QColor('#dce4ed'))
        x=left+(self.x-xmin)/(xmax-xmin)*w; y=top+h-(self.y-ymin)/(ymax-ymin)*h
        path=QPainterPath(QPointF(float(x[0]),float(y[0])))
        for px,py in zip(x[1:],y[1:]): path.lineTo(float(px),float(py))
        p.setPen(QPen(QColor('#087f78'),2)); p.drawPath(path)
        p.setPen(QColor('#52647a')); p.drawText(left,top+h+18,f'{xmin:g}'); p.drawText(left+w-45,top+h+18,f'{xmax:g}')
        p.drawText(left,top+h+40,self.xlabel+'   '+self.ylabel)

def spectrum(frequencies,intensities,width=12):
    if not frequencies: return [],[]
    f=np.asarray(frequencies); intensities=np.asarray(intensities if len(intensities)==len(f) else np.ones(len(f)))
    x=np.linspace(min(-100.,float(f.min())-100),max(4000.,float(f.max())+100),2400)
    y=np.zeros_like(x)
    for freq,intensity in zip(f,intensities): y+=intensity*(width**2/((x-freq)**2+width**2))
    return x,y
