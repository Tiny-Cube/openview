from PySide6.QtCore import QAbstractTableModel, Qt, Signal
from .model import SYMBOLS, element

class AtomTable(QAbstractTableModel):
    changed=Signal(int,int,object)
    def __init__(self,parent=None):
        super().__init__(parent); self.molecule=None; self.frame=0
    def set_molecule(self,mol):
        self.beginResetModel(); self.molecule=mol; self.frame=0; self.endResetModel()
    def rowCount(self,parent=None):
        if parent is not None and parent.isValid():return 0
        return len(self.molecule.numbers) if self.molecule is not None else 0
    def columnCount(self,parent=None):return 0 if parent is not None and parent.isValid() else 5
    def headerData(self,section,orientation,role=Qt.ItemDataRole.DisplayRole):
        if role==Qt.ItemDataRole.DisplayRole:
            return ['#','Element','X / Å','Y / Å','Z / Å'][section] if orientation==Qt.Orientation.Horizontal else str(section+1)
    def data(self,index,role=Qt.ItemDataRole.DisplayRole):
        if not index.isValid() or self.molecule is None: return None
        if role in (Qt.ItemDataRole.DisplayRole,Qt.ItemDataRole.EditRole):
            r,c=index.row(),index.column()
            if c==0:return r+1
            if c==1:return 'Open' if self.molecule.numbers[r]==0 else SYMBOLS[self.molecule.numbers[r]]
            return f'{self.molecule.frames[self.frame,r,c-2]:.6f}'
    def flags(self,index):
        flags=super().flags(index)
        if index.column()>0: flags|=Qt.ItemFlag.ItemIsEditable
        return flags
    def setData(self,index,value,role=Qt.ItemDataRole.EditRole):
        if not index.isValid() or self.molecule is None or role!=Qt.ItemDataRole.EditRole or index.column()<1:return False
        try:
            parsed=element(value) if index.column()==1 else float(value)
            import math
            if not math.isfinite(parsed) or (index.column()>1 and abs(parsed)>100000): return False
        except (ValueError,TypeError):return False
        self.changed.emit(index.row(),index.column(),parsed)
        self.dataChanged.emit(index,index); return True
    def refresh(self,frame):
        self.frame=frame
        if self.rowCount():self.dataChanged.emit(self.index(0,0),self.index(self.rowCount()-1,4))
