######################################################
# Created: August, 2024
# By: Md Kamruzzaman (@mkamruz)
# Objective:
#          1. Drawing SMILEs string
######################################################
import base64
from cairosvg import svg2png
import os
import hashlib
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image, ImageDraw, ImageFont

class DepictSmile:
    def __init__(self, save_file_path='./', size=(300, 300), dpi=200, verbose=0):
        self.molSize = size
        self.file_path = save_file_path
        self.dpi = dpi
        self.verbose = verbose

    def __moltosvg__(self, mol, kekulize = True):
        mc = Chem.Mol(mol.ToBinary())
        if kekulize:
            try:
                Chem.Kekulize(mc)
            except:
                mc = Chem.Mol(mol.ToBinary())
        if not mc.GetNumConformers():
            rdDepictor.Compute2DCoords(mc)
        drawer = rdMolDraw2D.MolDraw2DSVG(self.molSize[0], self.molSize[1])
        drawer.DrawMolecule(mc)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        return svg.replace('svg:','')
    
    #svg_img = SVG(moltosvg(m))
    
    def __SHA1__(self, msg: str) -> str:
        return hashlib.sha1(msg.encode()).hexdigest()

    def __get_file_nme__(self, smiles):
        return self.__SHA1__(smiles)

    def __add_text_to_image__(self, img_src, title):
        img = Image.open(img_src, 'r')
        draw = ImageDraw.Draw(img)
        w, h = img.size
        title_len = len(title)
        term = title_len//40
        si=0
        ei=0
        for i in range(term+1):
            si=ei
            ei+=41
            ei = min(ei, title_len)
            
            #font = ImageFont.truetype(size=16)
            text_w, text_h = (50*(term-i+1))-self.molSize[0], self.molSize[1]-20  #draw.textsize(title, font)
    
            X_POSN = h - text_h
            Y_POSN = w//2 - text_w//2  # Or (w - text_w) // 2
        
            draw.text((X_POSN, Y_POSN), title[si:ei], (0,0,0))
            
    
        img.save(img_src)
    
    def save_svg_to_png(self, smiles):
        m = Chem.MolFromSmiles(smiles)
        svg_img = self.__moltosvg__(m)
        fn = self.__get_file_nme__(smiles)
        fp = os.path.join(self.file_path, f'{fn}.png')
        
        svg2png(bytestring=svg_img, dpi=self.dpi, write_to=fp)
        self.__add_text_to_image__(fp, smiles, )
        
        if self.verbose>0:
            print(f"File saved to {fp}")
    
    def get_decoded_image(self, smiles):
        fn = self.__get_file_nme__(smiles)
        img = os.path.join(self.file_path, f'{fn}.png')
        enc = None
        
        if not os.path.exists(img):
            self.save_svg_to_png(smiles)
        
        with open(img, "rb") as image_file:
            encoded_string= base64.b64encode(image_file.read())
            enc=encoded_string.decode('utf-8')
        return enc