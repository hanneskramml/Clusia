import owncloud
from config import RPi


class Cloud:
    def __init__(self, url=RPi.CLOUD_URL, password=None):
        self.oc = owncloud.Client.from_public_link(url, folder_password=password)

    def upload(self, file):
        try:
            self.oc.drop_file(file)
            print("File uploaded to ucloud: {}".format(file))

        except Exception as e:
            print("ERROR uploading file ({}): {}".format(file, e))
