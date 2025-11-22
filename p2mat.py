#!/usr/bin/env python
######################################################
# Created: August, 2024
# By: Md Kamruzzaman (@mkamruz)
# Objective:
#          1. Create a global variable named PKGDIR
######################################################
from __future__ import print_function

# Import packages
import numpy as np
import pandas as pd
import time
import sys
from utils.mspp import MaterialStructuralPropertyPrediction
from utils import helper

import os

# UI related libraries
from PyQt5.Qt import Qt
from PyQt5.QtGui import QPixmap, QIcon
from PyQt5.QtWidgets import (
    QApplication,
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QTextEdit,
    QPushButton,
    QComboBox,
    QFormLayout,
    QDialog,
    QPlainTextEdit,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
    QFileDialog,
    QProgressBar,
    QAbstractItemView
)

class Popup(QDialog):
    def __init__(self, *args, **kwargs):
        super(Popup, self).__init__(*args, **kwargs)
        self.setLayout(QVBoxLayout())
        self.label = QLabel(self)
        self.layout().addWidget(self.label)
        

class MatPropPred(QWidget):
    def __init__(self, width=800, height=600, batch_size=128, verbose=0):
        super(MatPropPred, self).__init__()
        
        self.setWindowTitle("Material property prediction")
        self.resize(width, height)
        # Set the application icon
        # Replace 'icon.png' with the actual name/path of your icon file
        self.setWindowIcon(QIcon('icon.png'))
        self.__create_UI__()
        self.popup = None
        self.batch = batch_size
        self.verbose = verbose
        self.prop_data = None
        self.width = width
        self.height = height

    def __show_popup__(self):
        if self.popup is None:
            self.popup = Popup()
        self.popup.exec_()

    def __close_popup__(self):
        if isinstance(self.popup, QDialog):
            self.popup.close()

    def __get_smiles__(self):
        return self.smile_txt.toPlainText().split("\n")

    def __create_prop_table__(self):
        self.result_table.setHidden(False)
        
        # Row count
        self.result_table.setRowCount(self.prop_data.shape[0]+1)

        # Column count
        self.result_table.setColumnCount(self.prop_data.shape[1])
    
        # Set header
        for idx, col_name in enumerate(self.prop_data.columns):
            if 'MP' in col_name:
                col_name = "Melting point (K)"
            elif 'BP' in col_name:
                col_name = "Boiling point (K)"
            elif 'HC' in col_name:
                col_name = "Heat capacity"
            elif 'dH' in col_name:
                col_name = "Heat H2"
            elif 'Saturated_smiles' == col_name:
                col_name = "Saturated smiles"
            elif 'H2_uptake' == col_name:
                col_name = "H2 uptake"
                
            self.result_table.setHorizontalHeaderItem(idx, QTableWidgetItem(col_name))

        for row in range(self.prop_data.shape[0]):
            for col in range(self.prop_data.shape[1]):
                if col in range(1, 5):
                    self.result_table.setItem(row,col,QTableWidgetItem(str(np.round(self.prop_data.iloc[row, col], 3))))
                else:
                    self.result_table.setItem(row,col,QTableWidgetItem(str(self.prop_data.iloc[row, col])))

        # Fit table horizontally
        self.result_table.horizontalHeader().setStretchLastSection(True) 
        self.result_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch) 
    
    def __prop_prod_btn__(self, s):
        print(f"Button click={s}")
        self.result_details_lbl.setPlainText("")
        self.result_details_lbl.setHidden(True)
        self.result_table.setHidden(True)
        self.result_lbl.setHidden(True)
        self.button2.setHidden(True)
        self.pbar.setHidden(True)
        
        #self.__show_popup__()

        txt = self.__get_smiles__()
        #print(txt)

        if len(''.join(txt))>0:
            self.pbar.setMaximum(len(txt))
            self.resize(self.width, self.height+100)
            self.pbar.setHidden(False)
            time.sleep(0.05)
            print('Starting loop')
            wrong_smiles = []
            self.prop_data = None
            
            for i, smile in enumerate(txt):
                time.sleep(0.05)
                QApplication.processEvents() 

                try:
                    mspp = MaterialStructuralPropertyPrediction(verbose=self.verbose, batch=self.batch)
                    wrong_smile = mspp.prediction_from_smile([smile])
                    if len(wrong_smile)>0:
                        wrong_smiles.extend(wrong_smile)

                    tmp_data = mspp.get_data()
                    if tmp_data is None or len(tmp_data)==0:
                        print(f"No data for smile: {smile}")
                        continue

                    if self.prop_data is None:
                        self.prop_data = mspp.get_data()
                    else:
                        self.prop_data = pd.concat([self.prop_data, mspp.get_data()], axis=0, ignore_index=True)

                except Exception as e:
                    print(f'Error on prediction: {e}')
                finally:
                    del[mspp]

                self.pbar.setValue(i+1)
        
            if len(wrong_smiles)>0:
                self.result_lbl.setHidden(False)
                self.result_details_lbl.setHidden(False)
                self.result_details_lbl.appendPlainText("Following SMILEs are in wrong format:")
                for i,s in enumerate(wrong_smiles):
                    self.result_details_lbl.appendPlainText(f"{i+1}) {s}")
    
            if self.prop_data is not None and len(self.prop_data)>0:
                print(self.prop_data.shape)
                print(self.prop_data.columns)
                self.button2.setHidden(False)
                self.result_lbl.setHidden(False)
                self.resize(self.width, self.height+300)
                self.__create_prop_table__()
      
                
    
        #self.__close_popup__()

    def __save_file__(self, s):
        print(f"Button click={s}")
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getSaveFileName(self,"Save File","","csv Files(*.csv)",options = options)
        
        if file_name and len(file_name)>0 and file_name[-4:] != '.csv':
            file_name = file_name + ".csv"
            if self.prop_data is not None and len(self.prop_data)>0:
                self.prop_data.to_csv(file_name, index=False)
            
        print(f"File saved to {file_name}")
    
        
    def __create_UI__(self):
        # Create all app object
        logo_lbl = QLabel("")
        logo_path = helper.resource_path(os.path.join("design", "snl-logo-inline.svg"))
        print(f"Logo path: {logo_path}")
        pixmap = QPixmap(logo_path)
        logo_lbl.setPixmap(pixmap)
        logo_lbl.resize(pixmap.width(), pixmap.height())
        
        smile_lbl = QLabel("Enter a SMILE string\n(Multiple SMILEs in new line)")
        self.smile_txt = QTextEdit()
        
        #prop_lbl = QLabel("Select a property")
        self.prop_choice = QComboBox()
        self.prop_choice.addItems(["Melting point"])

        self.pbar = QProgressBar()
        
        #button_style = "background-color: #fff; color: #06c !important; border: 1px solid #06c; 
        #text-align: center; border-radius: 8px; font-size: 1em; padding: 1.25rem 3rem 1.25rem 3rem; width: 120px; height: 40px"
        
        self.button1 = QPushButton("Predict properties")
        self.button1.clicked.connect(self.__prop_prod_btn__)
        
        self.button2 = QPushButton("Save data")
        self.button2.clicked.connect(self.__save_file__)
        
        self.result_lbl = QLabel("Result")
        self.result_details_lbl = QPlainTextEdit("")
        self.result_details_lbl.setReadOnly(True)
        self.result_table = QTableWidget()

        # Disable editing
        self.result_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        
        # Disable selection
        self.result_table.setSelectionMode(QAbstractItemView.NoSelection)
        
        # Disable focus
        self.result_table.setFocusPolicy(Qt.NoFocus)

        self.result_lbl.setHidden(True)
        self.result_details_lbl.setHidden(True)
        self.result_table.setHidden(True)
        self.button2.setHidden(True)
        self.pbar.setHidden(True)
        
        # All design is here
        row_logo = QHBoxLayout()
        row_logo.addWidget(logo_lbl)
        
        row_smile = QHBoxLayout()
        row_smile.addWidget(smile_lbl) #, alignment=Qt.AlignLeft)
        row_smile.addWidget(self.smile_txt) #, alignment=Qt.AlignCenter)
        
        row_button_1 = QHBoxLayout()
        #row_button_1.addWidget(self.pbar, alignment=Qt.AlignLeft)
        row_button_1.addWidget(self.button1, alignment=Qt.AlignRight)

        row_pb = QHBoxLayout()
        row_pb.addWidget(self.pbar)

        result_form = QFormLayout()
        result_form.addRow(self.result_lbl)
        result_form.addRow(self.result_details_lbl)
        result_form.addRow(self.result_table)

        row_button_2 = QHBoxLayout()
        row_button_2.addWidget(self.button2, alignment=Qt.AlignRight)

        
        # Events
        master_layout = QVBoxLayout()
        
        #master_layout.addLayout(form_layout)
        master_layout.addLayout(row_logo)
        master_layout.addLayout(row_smile)
        master_layout.addLayout(row_button_1)
        master_layout.addLayout(row_pb)
        master_layout.addLayout(result_form)
        master_layout.addLayout(row_button_2)
        
        self.setLayout(master_layout)

if __name__=="__main__":
    app = QApplication(sys.argv)
    app.setWindowIcon(QIcon('icon.png'))

    css_file = helper.resource_path(os.path.join("design", "myapp.css"))
    print(f"CSS file path: {css_file}")

    stylesheet = helper.read_and_replace_css(css_file)
    app.setStyleSheet(stylesheet)
    
    window = MatPropPred()
    window.show()
    sys.exit(app.exec_())
    
