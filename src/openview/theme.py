"""Explicit light colors, including native popup and editor palette roles."""
from PySide6.QtGui import QColor, QPalette

STYLE='''
QWidget { color:#203047; font-size:13px; }
QMainWindow,QDialog { background:#f1f5f9; color:#203047; }
QToolBar { background:#ffffff; border:0; border-bottom:1px solid #ced8e3; padding:8px; spacing:9px; }
QToolButton,QPushButton { color:#203047; background:#ffffff; border:1px solid #bccbda; border-radius:5px; padding:7px 12px; }
QPushButton:hover,QToolButton:hover { background:#e7f0fa; border-color:#6b9dcc; }
QPushButton:checked,QToolButton:checked { background:#dcecfb; color:#124b80; border:1px solid #4686bf; }
QPushButton:disabled,QToolButton:disabled { color:#788594; background:#edf1f5; }
QLabel#brand { font-size:23px; font-weight:700; color:#087f78; }
QLabel#muted { color:#52647a; }
QComboBox,QSpinBox,QDoubleSpinBox,QLineEdit,QPlainTextEdit { color:#203047; background:#ffffff; border:1px solid #bccbda; border-radius:4px; padding:6px; selection-background-color:#d7e8fa; selection-color:#163a61; }
QComboBox { combobox-popup:0; }
QCheckBox::indicator { width:14px; height:14px; border:1px solid #8295a8; border-radius:2px; background:#ffffff; }
QCheckBox::indicator:checked { background:#16877f; border-color:#086b67; }
QComboBox:disabled,QSpinBox:disabled,QDoubleSpinBox:disabled { color:#788594; background:#edf1f5; }
QComboBox QAbstractItemView { color:#203047; background:#ffffff; alternate-background-color:#f3f7fb; selection-background-color:#d7e8fa; selection-color:#163a61; border:1px solid #9eb3c8; outline:0; }
QComboBox QAbstractItemView::item { min-height:24px; padding:3px 6px; }
QTabWidget::pane { border:1px solid #ced8e3; background:#ffffff; }
QTabBar::tab { background:#e8eef5; padding:9px 13px; color:#52647a; }
QTabBar::tab:selected { color:#086b67; background:#ffffff; border-top:2px solid #168c82; }
QTableView { color:#203047; background:#ffffff; alternate-background-color:#f2f6fb; gridline-color:#e0e7ef; selection-background-color:#d7e8fa; selection-color:#163a61; border:0; }
QTableView QLineEdit { background:#ffffff; color:#203047; }
QHeaderView::section { background:#e6edf5; color:#334e68; padding:6px; border:0; }
QGroupBox { border:1px solid #ced8e3; border-radius:5px; margin-top:12px; padding:12px; }
QGroupBox::title { subcontrol-origin:margin; left:10px; color:#087f78; }
QStatusBar { background:#e9eff6; color:#41566d; }
QSlider::groove:horizontal { height:5px; background:#c8d5e3; }
QSlider::handle:horizontal { width:14px; background:#12877e; margin:-5px 0; border-radius:6px; }
QSplitter::handle { background:#d5dfe9; width:3px; }
QMenu { color:#203047; background:#ffffff; border:1px solid #bccbda; }
QMenu::item:selected { color:#163a61; background:#d7e8fa; }
QToolTip { color:#203047; background:#ffffff; border:1px solid #bccbda; padding:4px; }
'''

def apply_theme(app):
    app.setStyle('Fusion')
    palette=QPalette()
    colors={'Window':'#f1f5f9','WindowText':'#203047','Base':'#ffffff',
            'AlternateBase':'#f2f6fb','Text':'#203047','Button':'#ffffff',
            'ButtonText':'#203047','Highlight':'#d7e8fa','HighlightedText':'#163a61',
            'ToolTipBase':'#ffffff','ToolTipText':'#203047','PlaceholderText':'#64748b',
            'Light':'#ffffff','Midlight':'#e6edf5','Mid':'#bccbda','Dark':'#8295a8',
            'Shadow':'#52647a','Link':'#165da0','LinkVisited':'#6b4696'}
    for role,color in colors.items():palette.setColor(getattr(QPalette.ColorRole,role),QColor(color))
    for role in ('WindowText','Text','ButtonText'):
        palette.setColor(QPalette.ColorGroup.Disabled,getattr(QPalette.ColorRole,role),QColor('#788594'))
    app.setPalette(palette);app.setStyleSheet(STYLE)
