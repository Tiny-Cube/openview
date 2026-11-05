"""Vector-drawn element coordination palette; independent of original artwork."""
import math
from PySide6.QtCore import Qt, QSize, QPointF, Signal, QTimer
from PySide6.QtGui import QIcon, QPixmap, QPainter, QPen, QColor, QPolygonF
from PySide6.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QLabel, QToolButton, QButtonGroup, QScrollArea, QFrame
from .fragments import FRAGMENTS, available_fragments, element_options, fragment_name, CATALOG
from .model import SYMBOLS


def fragment_icon(key,number=6):
    pix=QPixmap(92,58);pix.fill(Qt.GlobalColor.transparent)
    p=QPainter(pix);p.setRenderHint(QPainter.RenderHint.Antialiasing)
    font=p.font();font.setPixelSize(15);p.setFont(font)
    ink=QColor('#273747');p.setPen(QPen(ink,1.6))
    if key=='atom':
        p.drawText(pix.rect(),Qt.AlignmentFlag.AlignCenter,SYMBOLS[number]+'\nAtom');p.end();return QIcon(pix)
    center=QPointF(46,29)
    rays={'sing':[(math.pi,1)],'dit':[(math.pi,3)],'trid':[(math.pi,2)],
          'trisd':[(math.pi,1),(math.pi/3,2)],'triaa':[(math.pi,1.5),(math.pi/3,1.5)],
          'bent':[(math.pi*.19,1),(math.pi*.81,1)],
          'tetsss':[(math.pi,1),(-math.pi/3,1),(math.pi/3,1)],
          'st':[(math.pi,1),(0,3)],'dd':[(math.pi,2),(0,2)],
          'ssd':[(math.pi,2),(-math.pi/3,1),(math.pi/3,1)],
          'saa':[(math.pi,1),(-math.pi/3,1.5),(math.pi/3,1.5)],
          'tet':[(-math.pi/2,1),(math.pi,1),(.85,1),(.25,1)]}.get(key)
    if rays is None:
        rays=[]
        for (x,y,z),order in zip(FRAGMENTS[key].directions,FRAGMENTS[key].orders):
            px=.94*x+.34*z;py=.15*x-.91*y-.38*z
            rays.append((math.atan2(py,px),order))
    for i,(angle,order) in enumerate(rays):
        d=QPointF(math.cos(angle),math.sin(angle));normal=QPointF(-d.y(),d.x())
        start=center+d*10;end=center+d*28
        if key in ('tet','tetsss') and i==2:
            p.setBrush(ink);p.drawPolygon(QPolygonF([start,end+normal*4,end-normal*4]));p.setBrush(Qt.BrushStyle.NoBrush)
        elif key=='tet' and i==3:
            for step in range(5):
                f=(step+1)/5;at=start+(end-start)*f;p.drawLine(at-normal*4*f,at+normal*4*f)
        else:
            offsets=[0] if order==1 else [-2,2] if order in (1.5,2) else [-3,0,3]
            for j,off in enumerate(offsets):
                p.setPen(QPen(ink,1.4,Qt.PenStyle.DashLine if order==1.5 and j==1 else Qt.PenStyle.SolidLine))
                p.drawLine(start+normal*off,end+normal*off)
    p.setPen(ink)
    pairs={'trid':2,'trisd':1,'triaa':1,'bent':2,'tetsss':1}.get(key,0)
    p.setBrush(ink)
    for pair in range(pairs):
        for offset in (-2,2):p.drawEllipse(QPointF(43+pair*12+offset,16),1,1)
    symbol=SYMBOLS[number];p.drawText(round(46-p.fontMetrics().horizontalAdvance(symbol)/2),35,symbol);p.end();return QIcon(pix)


class ElementPalette(QWidget):
    chosen=Signal(str)
    def __init__(self,parent=None):
        super().__init__(parent)
        layout=QVBoxLayout(self);layout.setContentsMargins(0,0,0,0);layout.setSpacing(3)
        self.caption=QLabel();self.caption.setWordWrap(True);layout.addWidget(self.caption)
        self.scroll=QScrollArea();self.scroll.setWidgetResizable(True);self.scroll.setFrameShape(QFrame.Shape.NoFrame)
        self.scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded);self.scroll.setFixedHeight(72);layout.addWidget(self.scroll)
        self.number=None;self.buttons={};self.group=None;self.set_element(6,'tet')

    def set_element(self,number,key=None):
        key=key or element_options(number)['default']
        if self.number==number:
            self.set_fragment(key);return
        self.number=number
        old=self.scroll.takeWidget()
        if old:old.deleteLater()
        if self.group:self.group.deleteLater()
        panel=QWidget();row=QHBoxLayout(panel);row.setContentsMargins(0,0,0,0);row.setSpacing(3)
        self.group=QButtonGroup(self);self.group.setExclusive(True);self.buttons={}
        for fragment in available_fragments(number):
            b=QToolButton();b.setIcon(fragment_icon(fragment.key,number));b.setIconSize(QSize(66,46));b.setFixedSize(72,52)
            name=fragment_name(number,fragment.key)
            b.setCheckable(True);b.setToolTip(name);b.setAccessibleName(name)
            b.setStyleSheet('QToolButton:checked {background:#c5e9e2;border:2px solid #187d70;border-radius:4px;}')
            b.clicked.connect(lambda checked=False,key=fragment.key:self.chosen.emit(key))
            self.group.addButton(b);self.buttons[fragment.key]=b;row.addWidget(b)
        row.addStretch();self.scroll.setWidget(panel);self.set_fragment(key)

    def set_fragment(self,key):
        self.buttons[key].setChecked(True)
        name=fragment_name(self.number,key)
        self.caption.setText('Select '+element_options(self.number)['name']+' Fragment:  '+name.split(' · ',1)[1])
        if SYMBOLS[self.number] not in CATALOG['elements']:self.caption.setText(self.caption.text()+' · no original coordination preset')
        QTimer.singleShot(0,lambda:self.scroll.ensureWidgetVisible(self.buttons[key],4,0) if key in self.buttons else None)


CarbonPalette=ElementPalette
