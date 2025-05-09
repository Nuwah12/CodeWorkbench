import os

folder = "/mnt/data2/noah/fish_image/images_maxProj"

for filename in os.listdir(folder):
    if ".tif" in filename:
        newname = filename.replace("_MaxProj", "")
        #newname = newname.replace("Location", "Mask")
        os.rename(os.path.join(folder, filename), os.path.join(folder, newname))
