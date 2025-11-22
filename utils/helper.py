######################################################
# Created: August, 2024
# By: Md Kamruzzaman (@mkamruz)
# Objective:
#          1. Helper methods
######################################################
import os 
import sys

def read_and_replace_css(css_path):
    with open(css_path, "r") as file:
        css_content = file.read()
        # Replace placeholder with actual package directory
        #css_content = css_content.replace("__PACKAGE_DIR__", pkg_path)
    return css_content

def get_pkg_dir():
    """ Get absolute path to resource, works for dev and for PyInstaller """
    if getattr(sys, 'frozen', False):
        # The application is frozen (bundled)
        base_path = sys._MEIPASS
    else:
        # The application is not frozen
        base_path = os.path.abspath(".")

    return base_path

def resource_path(relative_path):
    return os.path.join(get_pkg_dir(), relative_path)