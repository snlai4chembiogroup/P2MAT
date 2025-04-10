######################################################
# Created: August, 2024
# By: Md Kamruzzaman (@mkamruz)
# Objective:
#          1. Helper methods
######################################################

def read_and_replace_css(css_path, pkg_path):
    with open(css_path, "r") as file:
        css_content = file.read()
        # Replace placeholder with actual package directory
        css_content = css_content.replace("__PACKAGE_DIR__", pkg_path)
    return css_content