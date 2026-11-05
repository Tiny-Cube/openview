from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (QDialog, QVBoxLayout, QLabel, QDoubleSpinBox,
    QSlider, QCheckBox, QDialogButtonBox)
from .geometry import geometry_value, edited_geometry


class GeometryDialog(QDialog):
    previewChanged = Signal(object)

    def __init__(self, positions, bonds, indices, parent=None):
        super().__init__(parent)
        self.positions = positions.copy(); self.bonds = bonds; self.indices = list(indices)
        self.candidate = self.positions.copy()
        count = len(indices)
        kind = {2: 'Bond length', 3: 'Bond angle', 4: 'Dihedral angle'}[count]
        self.setWindowTitle(kind); self.resize(440, 230)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel('Atoms: ' + ' → '.join(str(i + 1) for i in indices)))
        hint = QLabel('The first side stays fixed. Adjust the last atom or its attached group; OK commits one undoable edit, Cancel restores the structure.')
        hint.setWordWrap(True); layout.addWidget(hint)
        current = geometry_value(positions, indices)
        self.value = QDoubleSpinBox(); self.value.setDecimals(3)
        self.value.setRange(*({2: (.01, max(20., current * 2)), 3: (.1, 179.9), 4: (-180., 180.)}[count]))
        self.value.setSuffix(' Å' if count == 2 else '°'); self.value.setValue(current)
        layout.addWidget(self.value)
        self.slider = QSlider(Qt.Orientation.Horizontal)
        self.slider.setRange(round(self.value.minimum()*100), round(self.value.maximum()*100))
        self.slider.setValue(round(current*100)); layout.addWidget(self.slider)
        self.group = QCheckBox('Move attached group'); self.group.setChecked(True); layout.addWidget(self.group)
        self.error = QLabel(); self.error.setWordWrap(True); layout.addWidget(self.error)
        self.buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        self.buttons.accepted.connect(self.accept); self.buttons.rejected.connect(self.reject); layout.addWidget(self.buttons)
        self.value.valueChanged.connect(self.update_preview)
        self.group.toggled.connect(self.update_preview)
        self.slider.valueChanged.connect(lambda value: self.value.setValue(value/100))

    def update_preview(self):
        self.slider.blockSignals(True); self.slider.setValue(round(self.value.value()*100)); self.slider.blockSignals(False)
        try:
            self.candidate = edited_geometry(self.positions, self.bonds, self.indices, self.value.value(), self.group.isChecked())
        except ValueError as error:
            self.error.setText(str(error)); self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(False)
            self.previewChanged.emit(self.positions)
        else:
            self.error.clear(); self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(True)
            self.previewChanged.emit(self.candidate)
