'''Get test data for tests and/or examples'''

# based on https://github.com/GlacioHack/xdem/blob/d91bf1cc9b3f36d77f3729649bc8e9edc6b42f9f/xdem/examples.py#L33

import os
import tarfile
import tempfile
import urllib.request

def download_example_data() -> None:
    """
    Fetch the GOES ABI example files.
    """

    # Static commit hash (to be updated as needed)
    commit = "bc4e02ee93a9a8ca10738df3bcec2b829c838d69"
    # The URL from which to download the tarball
    url = f"https://github.com/spestana/goes-ortho-data/tarball/main#commit={commit}"

    # Path and filename for tarball
    tar_path = "./resources/data.tar.gz"

    response = urllib.request.urlopen(url)
    # If the response was right, download the tarball
    if response.getcode() == 200:
        with open(tar_path, "wb") as outfile:
            outfile.write(response.read())
    else:
        raise ValueError(f"Example GOES ABI data fetch gave non-200 response: {response.status_code}")

    # Extract the tarball
    with tarfile.open(tar_path) as tar:
        tar.extractall("./resources/")